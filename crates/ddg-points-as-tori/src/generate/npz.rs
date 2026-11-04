//! A writer of NumPy `.npz` archives (uncompressed, like `np.savez`) for the training and
//! validation files of SPEC.md §10.3. `numpy.load` and our readers ([`crate::npy::parse_npz_stored`],
//! `train::npz`, which also reads stored members) open the result. Members are little-endian `<f8`,
//! `<i4` or `|b1`, C order, `.npy` format version 1.0. No ZIP64: a member or the archive of 4 GiB or
//! more is an error (split the dataset into files of 500 shapes as the reference does).

use crate::PatError;
use std::io::Write;
use std::path::Path;

/// The values of a member.
#[derive(Debug, Clone, Copy)]
pub enum NpzOut<'a> {
    /// `<f8`.
    F64(&'a [f64]),
    /// `<i4`.
    I32(&'a [i32]),
    /// `|b1`.
    Bool(&'a [bool]),
}

fn crc_table() -> [u32; 256] {
    let mut t = [0u32; 256];
    for (i, e) in t.iter_mut().enumerate() {
        let mut c = i as u32;
        for _ in 0..8 {
            c = if c & 1 != 0 { 0xEDB8_8320 ^ (c >> 1) } else { c >> 1 };
        }
        *e = c;
    }
    t
}

fn crc32(table: &[u32; 256], data: &[u8]) -> u32 {
    let mut c = 0xFFFF_FFFFu32;
    for &b in data {
        c = table[((c ^ b as u32) & 0xFF) as usize] ^ (c >> 8);
    }
    !c
}

/// The bytes of an `.npy` file for `data` with the given shape.
pub fn npy_bytes(shape: &[usize], data: NpzOut<'_>) -> Result<Vec<u8>, PatError> {
    let (descr, n) = match data {
        NpzOut::F64(d) => ("<f8", d.len()),
        NpzOut::I32(d) => ("<i4", d.len()),
        NpzOut::Bool(d) => ("|b1", d.len()),
    };
    if shape.iter().product::<usize>() != n {
        return Err(PatError::SizeMismatch { what: format!("shape {shape:?} does not match {n} values") });
    }
    let shape_s = match shape.len() {
        1 => format!("({},)", shape[0]),
        _ => format!("({})", shape.iter().map(|s| s.to_string()).collect::<Vec<_>>().join(", ")),
    };
    let mut header = format!("{{'descr': '{descr}', 'fortran_order': False, 'shape': {shape_s}, }}");
    // total length 10 + header (+ newline) divisible by 64
    while (10 + header.len() + 1) % 64 != 0 {
        header.push(' ');
    }
    header.push('\n');
    let mut out = Vec::with_capacity(10 + header.len() + n * 8);
    out.extend_from_slice(b"\x93NUMPY\x01\x00");
    out.extend_from_slice(&(header.len() as u16).to_le_bytes());
    out.extend_from_slice(header.as_bytes());
    match data {
        NpzOut::F64(d) => d.iter().for_each(|x| out.extend_from_slice(&x.to_le_bytes())),
        NpzOut::I32(d) => d.iter().for_each(|x| out.extend_from_slice(&x.to_le_bytes())),
        NpzOut::Bool(d) => out.extend(d.iter().map(|&b| b as u8)),
    }
    Ok(out)
}

/// Write the archive `path` with the members `(name, shape, data)` (names get `.npy`).
pub fn write_npz(path: &Path, members: &[(&str, Vec<usize>, NpzOut<'_>)]) -> Result<(), PatError> {
    let err = |e: std::io::Error| PatError::Format { reason: format!("{}: {e}", path.display()) };
    let table = crc_table();
    let mut file = std::io::BufWriter::new(std::fs::File::create(path).map_err(err)?);
    let mut central: Vec<u8> = Vec::new();
    let mut offset: u64 = 0;
    for (name, shape, data) in members {
        let bytes = npy_bytes(shape, *data)?;
        let fname = format!("{name}.npy");
        if bytes.len() as u64 >= u32::MAX as u64 || offset >= u32::MAX as u64 {
            return Err(PatError::Format { reason: format!("{}: member {name} or the archive exceeds 4 GiB (no ZIP64 writer)", path.display()) });
        }
        let crc = crc32(&table, &bytes);
        let mut h = Vec::new();
        h.extend_from_slice(&0x0403_4b50u32.to_le_bytes());
        h.extend_from_slice(&20u16.to_le_bytes()); // version needed
        h.extend_from_slice(&0u16.to_le_bytes()); // flags
        h.extend_from_slice(&0u16.to_le_bytes()); // stored
        h.extend_from_slice(&0u16.to_le_bytes()); // time
        h.extend_from_slice(&33u16.to_le_bytes()); // date 1980-01-01
        h.extend_from_slice(&crc.to_le_bytes());
        h.extend_from_slice(&(bytes.len() as u32).to_le_bytes());
        h.extend_from_slice(&(bytes.len() as u32).to_le_bytes());
        h.extend_from_slice(&(fname.len() as u16).to_le_bytes());
        h.extend_from_slice(&0u16.to_le_bytes());
        h.extend_from_slice(fname.as_bytes());
        file.write_all(&h).map_err(err)?;
        file.write_all(&bytes).map_err(err)?;
        central.extend_from_slice(&0x0201_4b50u32.to_le_bytes());
        central.extend_from_slice(&20u16.to_le_bytes()); // made by
        central.extend_from_slice(&20u16.to_le_bytes()); // needed
        central.extend_from_slice(&0u16.to_le_bytes());
        central.extend_from_slice(&0u16.to_le_bytes());
        central.extend_from_slice(&0u16.to_le_bytes());
        central.extend_from_slice(&33u16.to_le_bytes());
        central.extend_from_slice(&crc.to_le_bytes());
        central.extend_from_slice(&(bytes.len() as u32).to_le_bytes());
        central.extend_from_slice(&(bytes.len() as u32).to_le_bytes());
        central.extend_from_slice(&(fname.len() as u16).to_le_bytes());
        central.extend_from_slice(&0u16.to_le_bytes()); // extra
        central.extend_from_slice(&0u16.to_le_bytes()); // comment
        central.extend_from_slice(&0u16.to_le_bytes()); // disk
        central.extend_from_slice(&0u16.to_le_bytes()); // internal attrs
        central.extend_from_slice(&0u32.to_le_bytes()); // external attrs
        central.extend_from_slice(&(offset as u32).to_le_bytes());
        central.extend_from_slice(fname.as_bytes());
        offset += (h.len() + bytes.len()) as u64;
    }
    if offset + central.len() as u64 >= u32::MAX as u64 {
        return Err(PatError::Format { reason: format!("{}: archive exceeds 4 GiB (no ZIP64 writer)", path.display()) });
    }
    file.write_all(&central).map_err(err)?;
    let mut e = Vec::new();
    e.extend_from_slice(&0x0605_4b50u32.to_le_bytes());
    e.extend_from_slice(&0u16.to_le_bytes());
    e.extend_from_slice(&0u16.to_le_bytes());
    e.extend_from_slice(&(members.len() as u16).to_le_bytes());
    e.extend_from_slice(&(members.len() as u16).to_le_bytes());
    e.extend_from_slice(&(central.len() as u32).to_le_bytes());
    e.extend_from_slice(&(offset as u32).to_le_bytes());
    e.extend_from_slice(&0u16.to_le_bytes());
    file.write_all(&e).map_err(err)?;
    file.flush().map_err(err)
}
