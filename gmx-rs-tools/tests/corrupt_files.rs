//! Regression tests for corrupt and truncated trajectory / run input files.
//!
//! Every parser here decodes untrusted input, so a malformed file must come
//! back as an `Err` — never as a panic, an integer overflow or an oversized
//! allocation.

use gmx_rs_tools::frame::Frame;
use gmx_rs_tools::gro;
use gmx_rs_tools::tpr::{self, TpxHeader};
use gmx_rs_tools::trr;
use gmx_rs_tools::xdr::{self, Writer, XdrError};
use gmx_rs_tools::xtc;
use gmx_rs_tools::trx::frame_count;

fn sample_frame(n: usize) -> Frame {
    let mut f = Frame::new(n);
    f.step = Some(7);
    f.time = Some(1.25);
    f.boxm = Some([
        [2.5, 0.0, 0.0],
        [0.0, 2.5, 0.0],
        [0.0, 0.0, 2.5],
    ]);
    f.x = Some((0..n).map(|i| [i as f32 * 0.01; 3]).collect());
    f
}

fn patch_i32(bytes: &mut [u8], offset: usize, value: i32) {
    bytes[offset..offset + 4].copy_from_slice(&value.to_be_bytes());
}

// --- TRR ---------------------------------------------------------------------

/// Byte offsets inside a TRR frame written by `trr::write_frame`:
/// magic, version string, seven section sizes, x/v/f sizes, natoms, ...
const TRR_BOX_SIZE: usize = 32;
const TRR_NATOMS: usize = 64;

#[test]
fn trr_zero_natoms_with_coordinates_errors() {
    let mut w = Writer::new();
    trr::write_frame(&mut w, &sample_frame(4));
    let mut bytes = w.into_vec();
    // box_size = 0, natoms = 0, x_size != 0: the float-size probe used to
    // divide by natoms * 3 == 0.
    patch_i32(&mut bytes, TRR_BOX_SIZE, 0);
    patch_i32(&mut bytes, TRR_NATOMS, 0);
    let mut r = xdr::Reader::new(&bytes);
    assert!(trr::read_frame(&mut r).is_err());
}

#[test]
fn trr_negative_natoms_errors() {
    let mut w = Writer::new();
    trr::write_frame(&mut w, &sample_frame(4));
    let mut bytes = w.into_vec();
    patch_i32(&mut bytes, TRR_NATOMS, -3);
    let mut r = xdr::Reader::new(&bytes);
    assert!(trr::read_frame(&mut r).is_err());
}

#[test]
fn trr_negative_section_size_errors() {
    let mut w = Writer::new();
    trr::write_frame(&mut w, &sample_frame(4));
    let mut bytes = w.into_vec();
    patch_i32(&mut bytes, TRR_BOX_SIZE, -(36));
    let mut r = xdr::Reader::new(&bytes);
    assert!(trr::read_frame(&mut r).is_err());
}

#[test]
fn trr_truncated_frame_errors() {
    let mut w = Writer::new();
    trr::write_frame(&mut w, &sample_frame(64));
    let bytes = w.into_vec();
    let mut r = xdr::Reader::new(&bytes[..bytes.len() / 2]);
    assert!(matches!(
        trr::read_frame(&mut r),
        Err(XdrError::Truncated(_))
    ));
}

#[test]
fn trr_frame_count_on_truncated_file_errors() {
    let mut w = Writer::new();
    trr::write_frame(&mut w, &sample_frame(40));
    trr::write_frame(&mut w, &sample_frame(40));
    let bytes = w.into_vec();
    let path = std::env::temp_dir().join("gmx-rs-tools-truncated.trr");
    // Cut into the payload of the second frame: the header parses, but the
    // seek past the declared payload lands beyond the end of the file.
    std::fs::write(&path, &bytes[..bytes.len() - 8]).unwrap();
    let result = frame_count(path.to_str().unwrap());
    let _ = std::fs::remove_file(&path);
    assert!(matches!(result, Err(XdrError::Truncated(_))));
}

// --- XTC ---------------------------------------------------------------------

/// Byte offsets inside an XTC frame written by `xtc::write_frame` (1995 magic).
const XTC_SMALLIDX: usize = 84;

#[test]
fn xtc_out_of_range_smallidx_errors() {
    for bad in [0x7FFF7FFF, -1] {
        let mut w = Writer::new();
        xtc::write_frame(&mut w, &sample_frame(32), 1000.0);
        let mut bytes = w.into_vec();
        patch_i32(&mut bytes, XTC_SMALLIDX, bad);
        let mut r = xdr::Reader::new(&bytes);
        assert!(
            xtc::read_frame(&mut r).is_err(),
            "smallidx {bad} must be rejected"
        );
    }
}

#[test]
fn xtc_empty_bitstream_errors() {
    // A frame header claiming ten atoms with an empty compressed block: the
    // bit reader used to index past the payload and panic.
    let mut w = Writer::new();
    w.int(xtc::XTC_MAGIC);
    w.int(10);
    w.int(0);
    w.float(0.0);
    for c in 0..9 {
        w.float(if c % 4 == 0 { 2.5 } else { 0.0 });
    }
    w.int(10); // lsize
    w.float(1000.0); // precision
    for _ in 0..3 {
        w.int(0);
    } // minint
    for _ in 0..3 {
        w.int(0);
    } // maxint
    w.int(9); // smallidx (FIRSTIDX)
    w.int(0); // empty compressed block
    let bytes = w.into_vec();
    let mut r = xdr::Reader::new(&bytes);
    assert!(matches!(
        xtc::read_frame(&mut r),
        Err(XdrError::Truncated(_))
    ));
}

#[test]
fn xtc_truncated_frame_errors() {
    let mut w = Writer::new();
    xtc::write_frame(&mut w, &sample_frame(64), 1000.0);
    let bytes = w.into_vec();
    let mut r = xdr::Reader::new(&bytes[..bytes.len() / 2]);
    assert!(xtc::read_frame(&mut r).is_err());
}

#[test]
fn xtc_negative_block_size_errors() {
    let mut w = Writer::new();
    xtc::write_frame(&mut w, &sample_frame(32), 1000.0);
    let mut bytes = w.into_vec();
    patch_i32(&mut bytes, 88, -5);
    let mut r = xdr::Reader::new(&bytes);
    assert!(xtc::read_frame(&mut r).is_err());
}

#[test]
fn xtc_frame_count_on_truncated_file_errors() {
    let mut w = Writer::new();
    xtc::write_frame(&mut w, &sample_frame(40), 1000.0);
    xtc::write_frame(&mut w, &sample_frame(40), 1000.0);
    let bytes = w.into_vec();
    let path = std::env::temp_dir().join("gmx-rs-tools-truncated.xtc");
    std::fs::write(&path, &bytes[..bytes.len() - 8]).unwrap();
    let result = frame_count(path.to_str().unwrap());
    let _ = std::fs::remove_file(&path);
    assert!(matches!(result, Err(XdrError::Truncated(_))));
}

// --- TPR ---------------------------------------------------------------------

fn header_with(natoms: i32, b_x: bool) -> TpxHeader {
    TpxHeader {
        version_string: "VERSION 2026.3".into(),
        precision: 4,
        is_double: false,
        file_version: 138,
        file_generation: 29,
        tag: String::new(),
        natoms,
        ngtc: 0,
        fep_state: 0,
        lambda: 0.0,
        b_ir: false,
        b_top: false,
        b_x,
        b_v: false,
        b_f: false,
        b_box: false,
        size_of_tpr_body: 0,
    }
}

#[test]
fn tpr_huge_natoms_with_small_body_errors() {
    // One billion atoms claimed, eight bytes of body: the coordinate section
    // used to `Vec::with_capacity` the full count first.
    let header = header_with(1_000_000_000, true);
    let result = tpr::parse_body(&header, &[0u8; 8]);
    assert!(matches!(result, Err(XdrError::Truncated(_))));
}

#[test]
fn tpr_negative_natoms_errors() {
    let header = header_with(-1, true);
    let result = tpr::parse_body(&header, &[0u8; 8]);
    assert!(result.is_err());
}

// --- GRO ---------------------------------------------------------------------

#[test]
fn gro_negative_natoms_errors() {
    let path = std::env::temp_dir().join("gmx-rs-tools-negative.gro");
    std::fs::write(&path, "test\n-5\n").unwrap();
    let result = gro::read_all(path.to_str().unwrap());
    let _ = std::fs::remove_file(&path);
    assert!(result.is_err());
}

/// A gro line whose residue/atom name columns contain a multi-byte UTF-8
/// character must not panic on the (byte-indexed) fixed-column slices.
#[test]
fn gro_multibyte_characters_do_not_panic() {
    // Column layout: resid(5) resname(5) atomname(5) resnr(5) then 3x8.3f.
    let line = format!("{:>5}{}{}{:>5}{:8.3}{:8.3}{:8.3}", 1, "ütfé", "åtom", 1, 0.0, 0.0, 0.0);
    let content = format!("multibyte\n1\n{}\n   2.00000   2.00000   2.00000\n", line);
    assert!(line.len() >= 39);
    let path = std::env::temp_dir().join("gmx-rs-tools-multibyte.gro");
    std::fs::write(&path, content).unwrap();
    let result = gro::read_all(path.to_str().unwrap());
    let _ = std::fs::remove_file(&path);
    assert!(result.is_ok(), "a multi-byte gro line must parse, not panic");
}

// --- XDR helper --------------------------------------------------------------

#[test]
fn read_exact_chunked_reports_truncation() {
    let mut data = std::io::Cursor::new(vec![1u8; 10]);
    assert!(xdr::read_exact_chunked(&mut data, 10).is_ok());
    let mut data = std::io::Cursor::new(vec![1u8; 10]);
    assert!(matches!(
        xdr::read_exact_chunked(&mut data, 11),
        Err(XdrError::Truncated(_))
    ));
}
