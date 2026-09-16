//! XTC trajectory format.
//!
//! This is a direct port of `gromacs/fileio/xtcio.cpp` and the compressed
//! coordinate code in `gromacs/fileio/libxdrf.cpp`.
//!
//! Frame layout (all big endian):
//!
//! ```text
//! int    magic      1995 (or 2023 for >298261617 atoms)
//! int    natoms
//! int    step
//! float  time
//! float  boxm[3][3]
//! <compressed coordinates: int natoms, float prec, ...>
//! ```

use crate::frame::{is_triclinic, Frame, Matrix, Rvec};
use crate::xdr::{Reader, Result, Writer, XdrError};

pub const XTC_MAGIC: i32 = 1995;
pub const XTC_NEW_MAGIC: i32 = 2023;
pub const XTC_1995_MAX_NATOMS: usize = 298261617;

const FIRSTIDX: usize = 9;
const MAGICINTS: [i32; 73] = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 8, // 0-9
    10, 12, 16, 20, 25, 32, 40, 50, 64, 80, // 10-19
    101, 128, 161, 203, 256, 322, 406, 512, 645, 812, // 20-29
    1024, 1290, 1625, 2048, 2580, 3250, 4096, 5060, 6501, 8192, // 30-39
    10321, 13003, 16384, 20642, 26007, 32768, 41285, 52015, 65536, 82570, // 40-49
    104031, 131072, 165140, 208063, 262144, 330280, 416127, 524287, 660561, 832255, // 50-59
    1048576, 1321122, 1664510, 2097152, 2642245, 3329021, 4194304, 5284491, 6658042,
    8388607, // 60-69
    10568983, 13316085, 16777216, // 70-72
];
const LASTIDX: usize = MAGICINTS.len();

/// Largest float that can still be converted to an in-range 32 bit integer.
fn max_absolute_int() -> f32 {
    (i32::MAX as f32).next_down_positive()
}

trait NextDown {
    fn next_down_positive(self) -> f32;
}

impl NextDown for f32 {
    fn next_down_positive(self) -> f32 {
        // nextafterf(x, 0) for a positive value.
        f32::from_bits(self.to_bits() - 1)
    }
}

struct WriteBuffer {
    index: usize,
    lastbits: i32,
    lastbyte: u32,
    data: Vec<u8>,
}

struct ReadBuffer<'a> {
    index: usize,
    lastbits: i32,
    lastbyte: u32,
    data: &'a [u8],
}

/// `sendbits`: appends `num` to the bit stream using `num_of_bits` bits.
fn sendbits(buffer: &mut WriteBuffer, mut num_of_bits: i32, num: i32) {
    let mut lastbits = buffer.lastbits;
    let mut lastbyte = buffer.lastbyte;
    while num_of_bits >= 8 {
        lastbyte = (lastbyte << 8) | ((num >> (num_of_bits - 8)) as u32 & 0xff);
        if buffer.index >= buffer.data.len() {
            buffer.data.resize(buffer.index + 1, 0);
        }
        buffer.data[buffer.index] = (lastbyte >> lastbits) as u8;
        buffer.index += 1;
        num_of_bits -= 8;
    }
    if num_of_bits > 0 {
        lastbyte = (lastbyte << num_of_bits) | ((num as u32) & ((1u32 << num_of_bits) - 1));
        lastbits += num_of_bits;
        if lastbits >= 8 {
            lastbits -= 8;
            if buffer.index >= buffer.data.len() {
                buffer.data.resize(buffer.index + 1, 0);
            }
            buffer.data[buffer.index] = (lastbyte >> lastbits) as u8;
            buffer.index += 1;
        }
    }
    buffer.lastbits = lastbits;
    buffer.lastbyte = lastbyte;
    if lastbits > 0 {
        if buffer.index >= buffer.data.len() {
            buffer.data.resize(buffer.index + 1, 0);
        }
        buffer.data[buffer.index] = (lastbyte << (8 - lastbits)) as u8;
    }
}

/// `receivebits`: extracts `num_of_bits` bits from the stream.
fn receivebits(buffer: &mut ReadBuffer, mut num_of_bits: i32) -> i32 {
    let mask: u32 = if num_of_bits >= 32 {
        u32::MAX
    } else {
        (1u32 << num_of_bits) - 1
    };
    let mut lastbits = buffer.lastbits;
    let mut lastbyte = buffer.lastbyte;
    let mut num: u32 = 0;

    while num_of_bits >= 8 {
        lastbyte = (lastbyte << 8) | buffer.data[buffer.index] as u32;
        buffer.index += 1;
        num |= (lastbyte >> lastbits) << (num_of_bits - 8);
        num_of_bits -= 8;
    }
    if num_of_bits > 0 {
        if lastbits < num_of_bits {
            lastbits += 8;
            lastbyte = (lastbyte << 8) | buffer.data[buffer.index] as u32;
            buffer.index += 1;
        }
        lastbits -= num_of_bits;
        num |= (lastbyte >> lastbits) & ((1u32 << num_of_bits) - 1);
    }
    buffer.lastbits = lastbits;
    buffer.lastbyte = lastbyte;
    (num & mask) as i32
}

/// `sizeofint`: number of bits needed to store an integer with the given max size.
fn sizeofint(size: i32) -> i32 {
    let mut num: i32 = 1;
    let mut num_of_bits: i32 = 0;
    while size >= num && num_of_bits < 32 {
        num_of_bits += 1;
        num <<= 1;
    }
    num_of_bits
}

/// `sizeofints`: bit size of a set of compressed small integers.
fn sizeofints(num_of_ints: usize, sizes: &[u32]) -> i32 {
    let mut bytes = [0u32; 32];
    let mut num_of_bytes: usize = 1;
    bytes[0] = 1;
    let mut num_of_bits: i32 = 0;

    for i in 0..num_of_ints {
        let mut tmp: u32 = 0;
        let mut bytecnt: usize = 0;
        while bytecnt < num_of_bytes {
            tmp = bytes[bytecnt]
                .wrapping_mul(sizes[i])
                .wrapping_add(tmp);
            bytes[bytecnt] = tmp & 0xff;
            tmp >>= 8;
            bytecnt += 1;
        }
        while tmp != 0 {
            bytes[bytecnt] = tmp & 0xff;
            tmp >>= 8;
            bytecnt += 1;
        }
        num_of_bytes = bytecnt;
    }
    num_of_bytes -= 1;
    let mut num: u32 = 1;
    while bytes[num_of_bytes] >= num {
        num_of_bits += 1;
        num *= 2;
    }
    num_of_bits + (num_of_bytes as i32) * 8
}

fn sendints(buffer: &mut WriteBuffer, num_of_ints: usize, num_of_bits: i32, sizes: &[u32], nums: &[u32]) {
    let mut bytes = [0u32; 32];
    let mut tmp = nums[0];
    let mut num_of_bytes: usize = 0;
    loop {
        bytes[num_of_bytes] = tmp & 0xff;
        num_of_bytes += 1;
        tmp >>= 8;
        if tmp == 0 {
            break;
        }
    }
    for i in 1..num_of_ints {
        if nums[i] >= sizes[i] {
            panic!(
                "major breakdown in sendints num {} doesn't match size {}",
                nums[i], sizes[i]
            );
        }
        tmp = nums[i];
        let mut bytecnt = 0;
        while bytecnt < num_of_bytes {
            tmp = bytes[bytecnt]
                .wrapping_mul(sizes[i])
                .wrapping_add(tmp);
            bytes[bytecnt] = tmp & 0xff;
            tmp >>= 8;
            bytecnt += 1;
        }
        while tmp != 0 {
            bytes[bytecnt] = tmp & 0xff;
            tmp >>= 8;
            bytecnt += 1;
        }
        num_of_bytes = bytecnt;
    }

    if num_of_bits >= (num_of_bytes as i32) * 8 {
        for i in 0..num_of_bytes {
            sendbits(buffer, 8, bytes[i] as i32);
        }
        sendbits(buffer, num_of_bits - (num_of_bytes as i32) * 8, 0);
    } else {
        for i in 0..num_of_bytes - 1 {
            sendbits(buffer, 8, bytes[i] as i32);
        }
        let remaining = num_of_bits - ((num_of_bytes - 1) as i32) * 8;
        sendbits(buffer, remaining, bytes[num_of_bytes - 1] as i32);
    }
}

fn receiveints(buffer: &mut ReadBuffer, num_of_ints: usize, mut num_of_bits: i32, sizes: &[u32], nums: &mut [i32]) {
    let mut bytes = [0u32; 32];
    let mut num_of_bytes: usize = 0;
    while num_of_bits > 8 {
        bytes[num_of_bytes] = receivebits(buffer, 8) as u32;
        num_of_bytes += 1;
        num_of_bits -= 8;
    }
    if num_of_bits > 0 {
        bytes[num_of_bytes] = receivebits(buffer, num_of_bits) as u32;
        num_of_bytes += 1;
    }
    for i in (1..num_of_ints).rev() {
        if sizes[i] == 0 {
            panic!("Cannot read trajectory, file possibly corrupted.");
        }
        let mut num: u32 = 0;
        for j in (0..num_of_bytes).rev() {
            num = (num << 8) | bytes[j];
            let p = num / sizes[i];
            bytes[j] = p;
            num -= p * sizes[i];
        }
        nums[i] = num as i32;
    }
    nums[0] = (bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24)) as i32;
}

/// Writes compressed coordinates (`xdr3dfcoord`), e.g. everything after the boxm.
#[allow(unused_assignments)]
fn write_3dfcoord(w: &mut Writer, x: &[Rvec], prec: f32, magic_number: i32) {
    let size = x.len();
    w.int(size as i32);
    if size <= 9 {
        for coord in x {
            for c in coord {
                w.float(*c);
            }
        }
        return;
    }
    w.float(prec);

    let size3 = size * 3;
    let flat: Vec<f32> = x.iter().flat_map(|c| c.iter().copied()).collect();
    let mut ip = vec![0i32; size3];
    let mut minint = [i32::MAX; 3];
    let mut maxint = [i32::MIN; 3];
    let mut mindiff = i32::MAX;
    let mut oldlint = [0i32; 3];
    let mut errval = 1;

    for i in 0..size {
        for d in 0..3 {
            let v = flat[i * 3 + d];
            let mut lf = if v >= 0.0 { v * prec + 0.5 } else { v * prec - 0.5 };
            if lf.abs() > max_absolute_int() {
                errval = 0;
            }
            if lf > i32::MAX as f32 {
                lf = i32::MAX as f32;
            }
            let lint = lf as i32;
            if lint < minint[d] {
                minint[d] = lint;
            }
            if lint > maxint[d] {
                maxint[d] = lint;
            }
            ip[i * 3 + d] = lint;
        }
        if i > 0 {
            let diff = (oldlint[0] - ip[i * 3]).abs()
                + (oldlint[1] - ip[i * 3 + 1]).abs()
                + (oldlint[2] - ip[i * 3 + 2]).abs();
            if diff < mindiff {
                mindiff = diff;
            }
        }
        oldlint = [ip[i * 3], ip[i * 3 + 1], ip[i * 3 + 2]];
    }

    for d in 0..3 {
        w.int(minint[d]);
    }
    for d in 0..3 {
        w.int(maxint[d]);
    }

    if (maxint[0] as f32) - (minint[0] as f32) >= max_absolute_int()
        || (maxint[1] as f32) - (minint[1] as f32) >= max_absolute_int()
        || (maxint[2] as f32) - (minint[2] as f32) >= max_absolute_int()
    {
        errval = 0;
    }

    let sizeint: Vec<u32> = (0..3)
        .map(|d| (maxint[d] - minint[d] + 1) as u32)
        .collect();
    let mut bitsizeint = [0i32; 3];
    let bitsize: i32;
    if (sizeint[0] | sizeint[1] | sizeint[2]) > 0xffffff {
        bitsizeint[0] = sizeofint(sizeint[0] as i32);
        bitsizeint[1] = sizeofint(sizeint[1] as i32);
        bitsizeint[2] = sizeofint(sizeint[2] as i32);
        bitsize = 0;
    } else {
        bitsize = sizeofints(3, &sizeint);
    }

    let mut smallidx = FIRSTIDX;
    while smallidx < LASTIDX && MAGICINTS[smallidx] < mindiff {
        smallidx += 1;
    }
    w.int(smallidx as i32);

    let maxidx = std::cmp::min(LASTIDX, smallidx + 8);
    let minidx = maxidx - 8;
    let mut smaller = MAGICINTS[std::cmp::max(FIRSTIDX, smallidx.saturating_sub(1))] / 2;
    let mut smallnum = MAGICINTS[smallidx] / 2;
    let mut sizesmall = [MAGICINTS[smallidx] as u32; 3];
    let larger = MAGICINTS[maxidx] / 2;

    let mut buffer = WriteBuffer {
        index: 0,
        lastbits: 0,
        lastbyte: 0,
        data: vec![0u8; (size3 as f64 * 1.2) as usize + 32],
    };

    let mut prevcoord;
    prevcoord = [0i32; 3];
    let mut prevrun: i32 = -1;
    let mut i: usize = 0;
    let mut tmpcoord = [0u32; 30];

    while i < size {
        let mut is_small = 0;
        // thiscoord is a window into ip starting at 3*i
        let mut ci = i * 3;
        let is_smaller: i32;
        if smallidx < maxidx
            && i >= 1
            && (ip[ci] - prevcoord[0]).abs() < larger
            && (ip[ci + 1] - prevcoord[1]).abs() < larger
            && (ip[ci + 2] - prevcoord[2]).abs() < larger
        {
            is_smaller = 1;
        } else if smallidx > minidx {
            is_smaller = -1;
        } else {
            is_smaller = 0;
        }

        if i + 1 < size {
            if (ip[ci] - ip[ci + 3]).abs() < smallnum
                && (ip[ci + 1] - ip[ci + 4]).abs() < smallnum
                && (ip[ci + 2] - ip[ci + 5]).abs() < smallnum
            {
                for d in 0..3 {
                    let t = ip[ci + d];
                    ip[ci + d] = ip[ci + 3 + d];
                    ip[ci + 3 + d] = t;
                }
                is_small = 1;
            }
        }

        tmpcoord[0] = (ip[ci] - minint[0]) as u32;
        tmpcoord[1] = (ip[ci + 1] - minint[1]) as u32;
        tmpcoord[2] = (ip[ci + 2] - minint[2]) as u32;
        if bitsize == 0 {
            sendbits(&mut buffer, bitsizeint[0], tmpcoord[0] as i32);
            sendbits(&mut buffer, bitsizeint[1], tmpcoord[1] as i32);
            sendbits(&mut buffer, bitsizeint[2], tmpcoord[2] as i32);
        } else {
            sendints(&mut buffer, 3, bitsize, &sizeint, &tmpcoord[0..3]);
        }
        prevcoord = [ip[ci], ip[ci + 1], ip[ci + 2]];
        ci += 3;
        i += 1;

        let mut run = 0usize;
        let mut is_smaller = is_smaller;
        if is_small == 0 && is_smaller == -1 {
            is_smaller = 0;
        }
        while is_small != 0 && run < 8 * 3 {
            if is_smaller == -1 {
                let sq = (ip[ci] - prevcoord[0]).pow(2)
                    + (ip[ci + 1] - prevcoord[1]).pow(2)
                    + (ip[ci + 2] - prevcoord[2]).pow(2);
                if sq >= smaller * smaller {
                    is_smaller = 0;
                }
            }
            tmpcoord[run] = (ip[ci] - prevcoord[0] + smallnum) as u32;
            run += 1;
            tmpcoord[run] = (ip[ci + 1] - prevcoord[1] + smallnum) as u32;
            run += 1;
            tmpcoord[run] = (ip[ci + 2] - prevcoord[2] + smallnum) as u32;
            run += 1;

            prevcoord = [ip[ci], ip[ci + 1], ip[ci + 2]];
            i += 1;
            ci += 3;
            is_small = 0;
            if i < size
                && (ip[ci] - prevcoord[0]).abs() < smallnum
                && (ip[ci + 1] - prevcoord[1]).abs() < smallnum
                && (ip[ci + 2] - prevcoord[2]).abs() < smallnum
            {
                is_small = 1;
            }
        }

        if run as i32 != prevrun || is_smaller != 0 {
            prevrun = run as i32;
            sendbits(&mut buffer, 1, 1);
            sendbits(&mut buffer, 5, run as i32 + is_smaller + 1);
        } else {
            sendbits(&mut buffer, 1, 0);
        }
        let mut k = 0;
        while k < run {
            sendints(
                &mut buffer,
                3,
                smallidx as i32,
                &sizesmall,
                &tmpcoord[k..k + 3],
            );
            k += 3;
        }
        if is_smaller != 0 {
            smallidx = (smallidx as i32 + is_smaller) as usize;
            if is_smaller < 0 {
                smallnum = smaller;
                smaller = MAGICINTS[smallidx - 1] / 2;
            } else {
                smaller = smallnum;
                smallnum = MAGICINTS[smallidx] / 2;
            }
            sizesmall = [MAGICINTS[smallidx] as u32; 3];
        }
    }

    if buffer.lastbits != 0 {
        buffer.index += 1;
    }
    if magic_number == XTC_NEW_MAGIC {
        w.int64(buffer.index as i64);
    } else {
        w.int(buffer.index as i32);
    }
    w.opaque(&buffer.data[..buffer.index]);

    if errval == 0 {
        // GROMACS returns 0 here; the frame is still written, so we keep going
        // but note that the data was out of the representable range.
    }
}

/// Reads compressed coordinates (`xdr3dfcoord`) and returns `(coords, precision)`.
fn read_3dfcoord(r: &mut Reader, magic_number: i32) -> Result<(Vec<Rvec>, f32)> {
    let lsize = r.int()?;
    if lsize < 0 {
        return Err(XdrError::Invalid(format!("negative atom count {lsize}")));
    }
    let size = lsize as usize;
    if size <= 9 {
        let mut x = Vec::with_capacity(size);
        for _ in 0..size {
            x.push([r.float()?, r.float()?, r.float()?]);
        }
        return Ok((x, -1.0));
    }
    let precision = r.float()?;
    let minint = [r.int()?, r.int()?, r.int()?];
    let maxint = [r.int()?, r.int()?, r.int()?];
    let sizeint: Vec<u32> = (0..3)
        .map(|d| (maxint[d] - minint[d] + 1) as u32)
        .collect();
    let mut bitsizeint = [0i32; 3];
    let bitsize: i32;
    if (sizeint[0] | sizeint[1] | sizeint[2]) > 0xffffff {
        bitsizeint[0] = sizeofint(sizeint[0] as i32);
        bitsizeint[1] = sizeofint(sizeint[1] as i32);
        bitsizeint[2] = sizeofint(sizeint[2] as i32);
        bitsize = 0;
    } else {
        bitsize = sizeofints(3, &sizeint);
    }
    let mut smallidx = r.int()? as usize;
    let mut smaller = MAGICINTS[std::cmp::max(FIRSTIDX, smallidx.saturating_sub(1))] / 2;
    let mut smallnum = MAGICINTS[smallidx] / 2;
    let mut sizesmall = [MAGICINTS[smallidx] as u32; 3];

    let buffer_size = if magic_number == XTC_NEW_MAGIC {
        r.int64()? as usize
    } else {
        r.int()? as usize
    };
    // `xdr_opaque` pads the encoded block to a multiple of four bytes, so the
    // padding has to be skipped here as well.
    let raw = r.opaque(buffer_size)?;

    let mut buffer = ReadBuffer {
        index: 0,
        lastbits: 0,
        lastbyte: 0,
        data: raw,
    };

    let inv_precision = 1.0f32 / precision;
    let mut x: Vec<Rvec> = vec![[0.0; 3]; size];
    let mut thiscoord = [0i32; 3];
    let mut prevcoord: [i32; 3];
    let mut run: i32 = 0;
    let mut i: usize = 0;
    // `out` mirrors the monotonic `lfp` pointer of the C implementation: the
    // decoded coordinates are not written out in decoding order.
    let mut out: usize = 0;
    while i < size {
        if bitsize == 0 {
            thiscoord[0] = receivebits(&mut buffer, bitsizeint[0]);
            thiscoord[1] = receivebits(&mut buffer, bitsizeint[1]);
            thiscoord[2] = receivebits(&mut buffer, bitsizeint[2]);
        } else {
            let mut tmp = [0i32; 3];
            receiveints(&mut buffer, 3, bitsize, &sizeint, &mut tmp);
            thiscoord = tmp;
        }
        i += 1;
        for d in 0..3 {
            thiscoord[d] += minint[d];
        }
        prevcoord = thiscoord;

        let flag = receivebits(&mut buffer, 1);
        let mut is_smaller = 0;
        if flag == 1 {
            run = receivebits(&mut buffer, 5);
            is_smaller = run % 3;
            run -= is_smaller;
            is_smaller -= 1;
        }
        if run > 0 {
            let mut k = 0;
            while k < run {
                let mut tmp = [0i32; 3];
                receiveints(&mut buffer, 3, smallidx as i32, &sizesmall, &mut tmp);
                thiscoord = tmp;
                i += 1;
                for d in 0..3 {
                    thiscoord[d] += prevcoord[d] - smallnum;
                }
                if k == 0 {
                    // interchange first with second atom for better compression
                    std::mem::swap(&mut thiscoord[0], &mut prevcoord[0]);
                    std::mem::swap(&mut thiscoord[1], &mut prevcoord[1]);
                    std::mem::swap(&mut thiscoord[2], &mut prevcoord[2]);
                    x[out] = [
                        prevcoord[0] as f32 * inv_precision,
                        prevcoord[1] as f32 * inv_precision,
                        prevcoord[2] as f32 * inv_precision,
                    ];
                    out += 1;
                } else {
                    prevcoord = thiscoord;
                }
                x[out] = [
                    thiscoord[0] as f32 * inv_precision,
                    thiscoord[1] as f32 * inv_precision,
                    thiscoord[2] as f32 * inv_precision,
                ];
                out += 1;
                k += 3;
            }
        } else {
            x[out] = [
                thiscoord[0] as f32 * inv_precision,
                thiscoord[1] as f32 * inv_precision,
                thiscoord[2] as f32 * inv_precision,
            ];
            out += 1;
        }

        smallidx = (smallidx as i32 + is_smaller) as usize;
        if is_smaller < 0 {
            smallnum = smaller;
            if smallidx > FIRSTIDX {
                smaller = MAGICINTS[smallidx - 1] / 2;
            } else {
                smaller = 0;
            }
        } else if is_smaller > 0 {
            smaller = smallnum;
            smallnum = MAGICINTS[smallidx] / 2;
        }
        sizesmall = [MAGICINTS[smallidx] as u32; 3];
    }

    Ok((x, precision))
}

/// Writes a single XTC frame to `w`.
pub fn write_frame(w: &mut Writer, frame: &Frame, prec: f32) {
    let natoms = frame.x.as_ref().map(|x| x.len()).unwrap_or(frame.natoms);
    let magic = if natoms > XTC_1995_MAX_NATOMS {
        XTC_NEW_MAGIC
    } else {
        XTC_MAGIC
    };
    w.int(magic);
    w.int(natoms as i32);
    w.int(frame.step.unwrap_or(0) as i32);
    w.float(frame.time.unwrap_or(0.0) as f32);

    let boxm = frame.boxm.unwrap_or([[0.0; 3]; 3]);
    for row in boxm.iter() {
        for c in row.iter() {
            w.float(*c);
        }
    }
    let empty = Vec::new();
    let x = frame.x.as_ref().unwrap_or(&empty);
    write_3dfcoord(w, x, prec, magic);
}

/// Reads one XTC frame; returns `None` at clean end of file.
pub fn read_frame(r: &mut Reader) -> Result<Option<Frame>> {
    if r.remaining() < 4 {
        return Ok(None);
    }
    let magic = r.int()?;
    if magic != XTC_MAGIC && magic != XTC_NEW_MAGIC {
        return Err(XdrError::Invalid(format!(
            "Magic Number Error in XTC file (read {magic}, should be {XTC_MAGIC} or {XTC_NEW_MAGIC})"
        )));
    }
    let natoms = r.int()? as usize;
    let step = r.int()? as i64;
    let time = r.float()? as f64;

    let mut boxm = [[0.0f32; 3]; 3];
    for row in boxm.iter_mut() {
        for c in row.iter_mut() {
            *c = r.float()?;
        }
    }
    let (x, prec) = read_3dfcoord(r, magic)?;

    let mut frame = Frame::new(natoms);
    frame.step = Some(step);
    frame.time = Some(time);
    frame.boxm = Some(boxm);
    frame.x = Some(x);
    frame.prec = Some(prec);
    Ok(Some(frame))
}

/// Reads the raw bytes of one XTC frame from a stream.
///
/// The frame is assembled in a small buffer and then decoded with
/// [`read_frame`], which lets the tools process trajectories that are far
/// larger than memory.  Returns `Ok(None)` at a clean end of file.
/// Reads the fixed part of one frame: the frame header and the header of the
/// compressed coordinate block.
///
/// Returns the bytes read together with the number of payload bytes that
/// follow them, or `None` at a clean end of file.  Splitting a frame this way
/// lets [`read_frame_bytes`] decode it and [`frame_count`] skip the payload of
/// a trajectory that is never decoded.
fn read_frame_prefix<R: std::io::Read>(r: &mut R) -> Result<Option<(Vec<u8>, usize)>> {
    use crate::xdr::{read_exact, read_or_eof, xdr_pad};

    let mut buf: Vec<u8> = Vec::with_capacity(4096);

    let mut magic_bytes = [0u8; 4];
    if !read_or_eof(r, &mut magic_bytes)? {
        return Ok(None);
    }
    buf.extend_from_slice(&magic_bytes);
    let magic = i32::from_be_bytes(magic_bytes);
    if magic != XTC_MAGIC && magic != XTC_NEW_MAGIC {
        return Err(XdrError::Invalid(format!(
            "Magic Number Error in XTC file (read {magic}, should be {XTC_MAGIC} or {XTC_NEW_MAGIC})"
        )));
    }

    // natoms, step, time and the 3x3 box.
    let mut fixed = [0u8; 4 + 4 + 4 + 36];
    read_exact(r, &mut fixed)?;
    buf.extend_from_slice(&fixed);

    // Size field of the compressed block.
    let mut size_bytes = [0u8; 4];
    read_exact(r, &mut size_bytes)?;
    buf.extend_from_slice(&size_bytes);
    let lsize = i32::from_be_bytes(size_bytes);
    if lsize < 0 {
        return Err(XdrError::Invalid(format!("negative atom count {lsize}")));
    }
    let size = lsize as usize;
    if size <= 9 {
        // Stored as plain single precision coordinates.
        return Ok(Some((buf, size * 3 * 4)));
    }

    // precision, minint[3], maxint[3] and smallidx.
    let mut mid = [0u8; 4 + 12 + 12 + 4];
    read_exact(r, &mut mid)?;
    buf.extend_from_slice(&mid);

    let buffer_size = if magic == XTC_NEW_MAGIC {
        let mut b = [0u8; 8];
        read_exact(r, &mut b)?;
        buf.extend_from_slice(&b);
        i64::from_be_bytes(b) as usize
    } else {
        let mut b = [0u8; 4];
        read_exact(r, &mut b)?;
        buf.extend_from_slice(&b);
        i32::from_be_bytes(b) as usize
    };

    // The payload is opaque XDR data and therefore padded to four bytes.
    Ok(Some((buf, buffer_size + xdr_pad(buffer_size))))
}

pub fn read_frame_bytes<R: std::io::Read>(r: &mut R) -> Result<Option<Vec<u8>>> {
    use crate::xdr::read_exact;

    let Some((mut buf, payload)) = read_frame_prefix(r)? else {
        return Ok(None);
    };
    // The padding of the opaque payload is kept so that the buffer can be
    // decoded as a whole.
    let mut data = vec![0u8; payload];
    read_exact(r, &mut data)?;
    buf.extend_from_slice(&data);
    Ok(Some(buf))
}

/// Counts the frames of an XTC file.
///
/// Only the fixed part of every frame is read and the compressed coordinates
/// are skipped with a seek, so a trajectory much larger than memory is counted
/// without decoding a single frame.
pub fn frame_count(path: &str) -> Result<usize> {
    use std::io::{BufReader, Seek, SeekFrom};

    let file = std::fs::File::open(path)
        .map_err(|e| XdrError::Invalid(format!("cannot read {path}: {e}")))?;
    let mut r = BufReader::with_capacity(1 << 16, file);
    let mut frames = 0usize;
    while let Some((_, payload)) = read_frame_prefix(&mut r)? {
        r.seek(SeekFrom::Current(payload as i64))
            .map_err(|e| XdrError::Invalid(format!("cannot read {path}: {e}")))?;
        frames += 1;
    }
    Ok(frames)
}

/// True when `boxm` needs all nine values on disk (used by the GRO writer too).
pub fn box_is_triclinic(boxm: &Matrix) -> bool {
    is_triclinic(boxm)
}
