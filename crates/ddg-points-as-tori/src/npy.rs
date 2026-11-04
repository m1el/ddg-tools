//! A minimal reader of NumPy `.npy` files and of uncompressed (`np.savez`) `.npz` archives:
//! little-endian `f8`, `f4`, `i4`, `i8`, `b1`, C order. Used by the parity tests and tools.
//! (Compressed archives are read by the training data reader, WP7.)

use crate::PatError;

/// A loaded array, converted to f64 (`i4`/`i8`/`b1` exactly, `f4` exactly).
#[derive(Debug, Clone, PartialEq)]
pub struct NpyArray {
    /// Shape.
    pub shape: Vec<usize>,
    /// Values, C order.
    pub data: Vec<f64>,
    /// The `descr` string (e.g. `<f8`).
    pub descr: String,
    /// For `f4` arrays, the f32 values bit-exact.
    pub f32_data: Option<Vec<f32>>,
}

fn ferr(reason: String) -> PatError {
    PatError::Format { reason }
}

fn field<'a>(header: &'a str, key: &str) -> Option<&'a str> {
    let pat = format!("'{key}':");
    let at = header.find(&pat)? + pat.len();
    Some(header[at..].trim_start())
}

/// Parse `.npy` bytes.
pub fn parse_npy(bytes: &[u8]) -> Result<NpyArray, PatError> {
    if bytes.len() < 10 || &bytes[..6] != b"\x93NUMPY" {
        return Err(ferr("not an .npy file (bad magic)".into()));
    }
    let major = bytes[6];
    let (hlen, off) = match major {
        1 => (u16::from_le_bytes([bytes[8], bytes[9]]) as usize, 10),
        2 | 3 => {
            if bytes.len() < 12 {
                return Err(ferr("truncated .npy header".into()));
            }
            (u32::from_le_bytes([bytes[8], bytes[9], bytes[10], bytes[11]]) as usize, 12)
        }
        v => return Err(ferr(format!(".npy version {v} not supported"))),
    };
    let header = std::str::from_utf8(bytes.get(off..off + hlen).ok_or_else(|| ferr("truncated .npy header".into()))?)
        .map_err(|_| ferr("non-UTF-8 .npy header".into()))?;
    let descr = field(header, "descr").ok_or_else(|| ferr(format!("no descr in {header}")))?;
    let descr: String = descr.trim_start_matches('\'').chars().take_while(|&c| c != '\'').collect();
    let fo = field(header, "fortran_order").ok_or_else(|| ferr(format!("no fortran_order in {header}")))?;
    if fo.starts_with("True") {
        return Err(ferr("Fortran-order arrays are not supported".into()));
    }
    let sh = field(header, "shape").ok_or_else(|| ferr(format!("no shape in {header}")))?;
    let sh = sh.trim_start_matches('(');
    let sh = &sh[..sh.find(')').ok_or_else(|| ferr(format!("bad shape in {header}")))?];
    let shape: Vec<usize> = sh
        .split(',')
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .map(|s| s.parse::<usize>().map_err(|_| ferr(format!("bad shape entry {s:?}"))))
        .collect::<Result<_, _>>()?;
    let count: usize = shape.iter().product();
    let body = &bytes[off + hlen..];
    let width = match descr.as_str() {
        "<f8" | "<i8" => 8,
        "<f4" | "<i4" => 4,
        "|b1" => 1,
        d => return Err(ferr(format!("dtype {d} not supported"))),
    };
    if body.len() < count * width {
        return Err(ferr(format!("{} data bytes for {count} values of {descr}", body.len())));
    }
    let mut f32_data = None;
    let data: Vec<f64> = match descr.as_str() {
        "<f8" => body.chunks_exact(8).take(count).map(|c| f64::from_le_bytes(c.try_into().unwrap_or([0; 8]))).collect(),
        "<i8" => body.chunks_exact(8).take(count).map(|c| i64::from_le_bytes(c.try_into().unwrap_or([0; 8])) as f64).collect(),
        "<f4" => {
            let v: Vec<f32> = body.chunks_exact(4).take(count).map(|c| f32::from_le_bytes(c.try_into().unwrap_or([0; 4]))).collect();
            let d = v.iter().map(|&x| x as f64).collect();
            f32_data = Some(v);
            d
        }
        "<i4" => body.chunks_exact(4).take(count).map(|c| i32::from_le_bytes(c.try_into().unwrap_or([0; 4])) as f64).collect(),
        _ => body.iter().take(count).map(|&b| (b != 0) as u8 as f64).collect(),
    };
    Ok(NpyArray { shape, data, descr, f32_data })
}

/// Read an `.npy` file.
pub fn read_npy(path: &std::path::Path) -> Result<NpyArray, PatError> {
    let b = std::fs::read(path).map_err(|e| ferr(format!("{}: {e}", path.display())))?;
    parse_npy(&b).map_err(|e| ferr(format!("{}: {e}", path.display())))
}

/// The members of an uncompressed `.npz` archive (`np.savez`), by name without `.npy`.
pub fn parse_npz_stored(bytes: &[u8]) -> Result<Vec<(String, NpyArray)>, PatError> {
    let mut out = Vec::new();
    let mut p = 0usize;
    let u16at = |b: &[u8], i: usize| u16::from_le_bytes([b[i], b[i + 1]]) as usize;
    let u32at = |b: &[u8], i: usize| u32::from_le_bytes([b[i], b[i + 1], b[i + 2], b[i + 3]]) as usize;
    while p + 30 <= bytes.len() && bytes[p..p + 4] == [0x50, 0x4b, 0x03, 0x04] {
        let method = u16at(bytes, p + 8);
        let mut csize = u32at(bytes, p + 18);
        let nlen = u16at(bytes, p + 26);
        let xlen = u16at(bytes, p + 28);
        let name = String::from_utf8_lossy(&bytes[p + 30..p + 30 + nlen]).to_string();
        let extra = &bytes[p + 30 + nlen..p + 30 + nlen + xlen];
        if csize == 0xffff_ffff {
            // zip64 extra field: id 1, then uncompressed and compressed sizes
            let mut e = 0;
            while e + 4 <= extra.len() {
                let id = u16at(extra, e);
                let sz = u16at(extra, e + 2);
                if id == 1 && sz >= 16 {
                    csize = u64::from_le_bytes(extra[e + 12..e + 20].try_into().unwrap_or([0; 8])) as usize;
                }
                e += 4 + sz;
            }
        }
        if method != 0 {
            return Err(ferr(format!("npz member {name} is compressed (method {method}); not supported here")));
        }
        let start = p + 30 + nlen + xlen;
        let data = bytes.get(start..start + csize).ok_or_else(|| ferr(format!("npz member {name} truncated")))?;
        let arr = parse_npy(data).map_err(|e| ferr(format!("npz member {name}: {e}")))?;
        out.push((name.trim_end_matches(".npy").to_string(), arr));
        p = start + csize;
    }
    if out.is_empty() {
        return Err(ferr("not an uncompressed npz archive".into()));
    }
    Ok(out)
}
