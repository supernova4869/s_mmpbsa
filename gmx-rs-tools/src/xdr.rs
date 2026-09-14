//! XDR (external data representation) helpers.
//!
//! Mirrors `gromacs/fileio/xdr_serializer.cpp` and the XDR primitives that
//! GROMACS relies on: big endian integers/floats and 4-byte aligned opaque
//! data.  GROMACS writes strings as `int(len + 1)` (an historical quirk to
//! stay compatible with the terminating null byte) followed by the standard
//! XDR string (`int(len)`, the characters and padding).

use std::fmt;

/// Errors produced while decoding a binary file.
#[derive(Debug)]
pub enum XdrError {
    /// The underlying data ended prematurely.
    Truncated(&'static str),
    /// The data was structurally invalid.
    Invalid(String),
}

impl fmt::Display for XdrError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            XdrError::Truncated(what) => write!(f, "unexpected end of data while reading {what}"),
            XdrError::Invalid(msg) => write!(f, "{msg}"),
        }
    }
}

impl std::error::Error for XdrError {}

pub type Result<T> = std::result::Result<T, XdrError>;

/// Reads exactly `buf.len()` bytes.
pub fn read_exact<R: std::io::Read>(r: &mut R, buf: &mut [u8]) -> Result<()> {
    r.read_exact(buf)
        .map_err(|e| XdrError::Invalid(format!("read error: {e}")))
}

/// Reads exactly `buf.len()` bytes, reporting a clean end of file.
///
/// Returns `Ok(false)` when the stream is already at its end.
pub fn read_or_eof<R: std::io::Read>(r: &mut R, buf: &mut [u8]) -> Result<bool> {
    let mut filled = 0;
    while filled < buf.len() {
        match r.read(&mut buf[filled..]) {
            Ok(0) => {
                if filled == 0 {
                    return Ok(false);
                }
                return Err(XdrError::Truncated("binary data"));
            }
            Ok(n) => filled += n,
            Err(e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
            Err(e) => return Err(XdrError::Invalid(format!("read error: {e}"))),
        }
    }
    Ok(true)
}

/// Number of padding bytes needed to reach the next multiple of four.
pub fn xdr_pad(len: usize) -> usize {
    (4 - (len % 4)) % 4
}

/// Sequential big-endian reader over a byte slice.
pub struct Reader<'a> {
    data: &'a [u8],
    pos: usize,
}

impl<'a> Reader<'a> {
    pub fn new(data: &'a [u8]) -> Self {
        Reader { data, pos: 0 }
    }

    pub fn position(&self) -> usize {
        self.pos
    }

    pub fn seek(&mut self, pos: usize) {
        self.pos = pos;
    }

    pub fn remaining(&self) -> usize {
        self.data.len() - self.pos
    }

    pub fn bytes(&mut self, n: usize) -> Result<&'a [u8]> {
        if self.pos + n > self.data.len() {
            return Err(XdrError::Truncated("opaque data"));
        }
        let slice = &self.data[self.pos..self.pos + n];
        self.pos += n;
        Ok(slice)
    }

    /// Reads `n` bytes and skips the XDR padding.
    pub fn opaque(&mut self, n: usize) -> Result<&'a [u8]> {
        let slice = self.bytes(n)?;
        self.pos += xdr_pad(n);
        Ok(slice)
    }

    pub fn u32(&mut self) -> Result<u32> {
        let b = self.bytes(4)?;
        Ok(u32::from_be_bytes([b[0], b[1], b[2], b[3]]))
    }

    pub fn int(&mut self) -> Result<i32> {
        Ok(self.u32()? as i32)
    }

    pub fn int64(&mut self) -> Result<i64> {
        let b = self.bytes(8)?;
        Ok(i64::from_be_bytes([
            b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7],
        ]))
    }

    pub fn uchar(&mut self) -> Result<u8> {
        Ok(self.u32()? as u8)
    }

    pub fn ushort(&mut self) -> Result<u16> {
        Ok(self.u32()? as u16)
    }

    pub fn bool(&mut self) -> Result<bool> {
        Ok(self.int()? != 0)
    }

    pub fn float(&mut self) -> Result<f32> {
        Ok(f32::from_bits(self.u32()?))
    }

    pub fn double(&mut self) -> Result<f64> {
        let b = self.bytes(8)?;
        Ok(f64::from_bits(u64::from_be_bytes([
            b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7],
        ])))
    }

    /// Reads a `real`, i.e. a float or double depending on the file precision.
    pub fn real(&mut self, double_precision: bool) -> Result<f64> {
        if double_precision {
            self.double()
        } else {
            Ok(self.float()? as f64)
        }
    }

    /// Reads a GROMACS serialized string (length+1, length, chars, padding).
    pub fn gmx_string(&mut self) -> Result<String> {
        let _len_plus_one = self.int()?;
        self.xdr_string()
    }

    /// Reads a plain XDR string (`int(len)`, chars, padding).
    pub fn xdr_string(&mut self) -> Result<String> {
        let len = self.int()?;
        if len < 0 {
            return Err(XdrError::Invalid(format!("negative string length {len}")));
        }
        let len = len as usize;
        let raw = self.opaque(len)?;
        // The serialized buffer contains a terminating null we do not want.
        let end = raw.iter().position(|&c| c == 0).unwrap_or(raw.len());
        Ok(String::from_utf8_lossy(&raw[..end]).into_owned())
    }

    pub fn int_array(&mut self, n: usize) -> Result<Vec<i32>> {
        let mut v = Vec::with_capacity(n);
        for _ in 0..n {
            v.push(self.int()?);
        }
        Ok(v)
    }

    pub fn real_array(&mut self, n: usize, double_precision: bool) -> Result<Vec<f64>> {
        let mut v = Vec::with_capacity(n);
        for _ in 0..n {
            v.push(self.real(double_precision)?);
        }
        Ok(v)
    }
}

/// Growable big-endian writer producing XDR encoded bytes.
#[derive(Default)]
pub struct Writer {
    pub data: Vec<u8>,
}

impl Writer {
    pub fn new() -> Self {
        Writer { data: Vec::new() }
    }

    pub fn len(&self) -> usize {
        self.data.len()
    }

    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    pub fn into_vec(self) -> Vec<u8> {
        self.data
    }

    pub fn bytes(&mut self, b: &[u8]) {
        self.data.extend_from_slice(b);
    }

    pub fn opaque(&mut self, b: &[u8]) {
        self.data.extend_from_slice(b);
        self.data.extend(std::iter::repeat(0u8).take(xdr_pad(b.len())));
    }

    pub fn int(&mut self, v: i32) {
        self.data.extend_from_slice(&v.to_be_bytes());
    }

    pub fn u32(&mut self, v: u32) {
        self.data.extend_from_slice(&v.to_be_bytes());
    }

    pub fn int64(&mut self, v: i64) {
        self.data.extend_from_slice(&v.to_be_bytes());
    }

    pub fn byte(&mut self, v: u8) {
        self.int(v as i32);
    }

    pub fn ushort(&mut self, v: u16) {
        self.int(v as i32);
    }

    pub fn bool(&mut self, v: bool) {
        self.int(if v { 1 } else { 0 });
    }

    pub fn float(&mut self, v: f32) {
        self.u32(v.to_bits());
    }

    pub fn double(&mut self, v: f64) {
        self.data.extend_from_slice(&v.to_bits().to_be_bytes());
    }

    pub fn real(&mut self, v: f64, double_precision: bool) {
        if double_precision {
            self.double(v);
        } else {
            self.float(v as f32);
        }
    }

    /// Writes a GROMACS serialized string (length+1, length, chars, padding).
    pub fn gmx_string(&mut self, s: &str) {
        self.int(s.len() as i32 + 1);
        self.int(s.len() as i32);
        self.data.extend_from_slice(s.as_bytes());
        self.data.extend(std::iter::repeat(0u8).take(xdr_pad(s.len())));
    }

    /// Patches a previously written 64 bit integer at `offset`.
    pub fn patch_int64(&mut self, offset: usize, v: i64) {
        self.data[offset..offset + 8].copy_from_slice(&v.to_be_bytes());
    }
}
