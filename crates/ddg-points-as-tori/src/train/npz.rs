//! A typed reader of NumPy `.npz` archives, compressed (`np.savez_compressed`, deflate) or not
//! (SPEC.md §10.3). Members are read one at a time straight from the archive into typed vectors
//! (`f64`, `i32`, `bool`), little-endian and C order; nothing is converted to a wider type, so the
//! 262 MB `neighbors` member of a 500-shape file stays 262 MB.
//!
//! Every error names the file and the member.

use crate::PatError;
use std::io::Read;
use std::path::{Path, PathBuf};

/// The values of one member.
#[derive(Debug, Clone, PartialEq)]
pub enum NpzData {
    /// `<f8`.
    F64(Vec<f64>),
    /// `<f4`.
    F32(Vec<f32>),
    /// `<i4`.
    I32(Vec<i32>),
    /// `<i8`.
    I64(Vec<i64>),
    /// `|b1`.
    Bool(Vec<bool>),
}

impl NpzData {
    fn len(&self) -> usize {
        match self {
            NpzData::F64(v) => v.len(),
            NpzData::F32(v) => v.len(),
            NpzData::I32(v) => v.len(),
            NpzData::I64(v) => v.len(),
            NpzData::Bool(v) => v.len(),
        }
    }

    fn dtype(&self) -> &'static str {
        match self {
            NpzData::F64(_) => "<f8",
            NpzData::F32(_) => "<f4",
            NpzData::I32(_) => "<i4",
            NpzData::I64(_) => "<i8",
            NpzData::Bool(_) => "|b1",
        }
    }
}

/// One member: its shape and values.
#[derive(Debug, Clone, PartialEq)]
pub struct NpzArray {
    /// Shape (C order).
    pub shape: Vec<usize>,
    /// Values.
    pub data: NpzData,
}

/// An open archive.
pub struct Npz {
    path: PathBuf,
    archive: zip::ZipArchive<std::io::BufReader<std::fs::File>>,
}

fn ferr(path: &Path, member: &str, reason: impl std::fmt::Display) -> PatError {
    PatError::Format { reason: format!("{} member {member}: {reason}", path.display()) }
}

fn header_field<'a>(header: &'a str, key: &str) -> Option<&'a str> {
    let pat = format!("'{key}':");
    let at = header.find(&pat)? + pat.len();
    Some(header[at..].trim_start())
}

/// Parse an `.npy` header dictionary: `(descr, shape)`.
fn parse_header(header: &str) -> Result<(String, Vec<usize>), String> {
    let descr = header_field(header, "descr").ok_or_else(|| format!("no descr in {header}"))?;
    let descr: String = descr.trim_start_matches('\'').chars().take_while(|&c| c != '\'').collect();
    let fo = header_field(header, "fortran_order").ok_or_else(|| format!("no fortran_order in {header}"))?;
    if fo.starts_with("True") {
        return Err("Fortran-order arrays are not supported".into());
    }
    let sh = header_field(header, "shape").ok_or_else(|| format!("no shape in {header}"))?;
    let sh = sh.trim_start_matches('(');
    let end = sh.find(')').ok_or_else(|| format!("bad shape in {header}"))?;
    let shape = sh[..end]
        .split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(|s| s.parse::<usize>().map_err(|_| format!("bad shape entry {s:?}")))
        .collect::<Result<Vec<_>, _>>()?;
    Ok((descr, shape))
}

fn read_exact_vec<R: Read>(r: &mut R, n: usize) -> std::io::Result<Vec<u8>> {
    let mut v = vec![0u8; n];
    r.read_exact(&mut v)?;
    Ok(v)
}

/// Read one `.npy` stream.
fn read_npy<R: Read>(r: &mut R) -> Result<NpzArray, String> {
    let mut pre = [0u8; 8];
    r.read_exact(&mut pre).map_err(|e| format!("truncated header: {e}"))?;
    if &pre[..6] != b"\x93NUMPY" {
        return Err("not an .npy stream (bad magic)".into());
    }
    let hlen = match pre[6] {
        1 => {
            let mut b = [0u8; 2];
            r.read_exact(&mut b).map_err(|e| format!("truncated header: {e}"))?;
            u16::from_le_bytes(b) as usize
        }
        2 | 3 => {
            let mut b = [0u8; 4];
            r.read_exact(&mut b).map_err(|e| format!("truncated header: {e}"))?;
            u32::from_le_bytes(b) as usize
        }
        v => return Err(format!(".npy version {v} not supported")),
    };
    let hb = read_exact_vec(r, hlen).map_err(|e| format!("truncated header: {e}"))?;
    let header = std::str::from_utf8(&hb).map_err(|_| "non-UTF-8 header".to_string())?;
    let (descr, shape) = parse_header(header)?;
    let count = shape.iter().try_fold(1usize, |a, &b| a.checked_mul(b)).ok_or("shape overflows")?;
    let width = match descr.as_str() {
        "<f8" | "<i8" => 8,
        "<f4" | "<i4" => 4,
        "|b1" => 1,
        d => return Err(format!("dtype {d} not supported")),
    };
    let bytes = count.checked_mul(width).ok_or("size overflows")?;
    let body = read_exact_vec(r, bytes).map_err(|e| format!("{count} values of {descr} expected: {e}"))?;
    let data = match descr.as_str() {
        "<f8" => NpzData::F64(body.chunks_exact(8).map(|c| f64::from_le_bytes(c.try_into().unwrap_or([0; 8]))).collect()),
        "<i8" => NpzData::I64(body.chunks_exact(8).map(|c| i64::from_le_bytes(c.try_into().unwrap_or([0; 8]))).collect()),
        "<f4" => NpzData::F32(body.chunks_exact(4).map(|c| f32::from_le_bytes(c.try_into().unwrap_or([0; 4]))).collect()),
        "<i4" => NpzData::I32(body.chunks_exact(4).map(|c| i32::from_le_bytes(c.try_into().unwrap_or([0; 4]))).collect()),
        _ => {
            if let Some(&b) = body.iter().find(|&&b| b > 1) {
                return Err(format!("bool value byte {b}"));
            }
            NpzData::Bool(body.iter().map(|&b| b == 1).collect())
        }
    };
    Ok(NpzArray { shape, data })
}

impl Npz {
    /// Open an archive.
    pub fn open(path: &Path) -> Result<Npz, PatError> {
        let f = std::fs::File::open(path).map_err(|e| PatError::Format { reason: format!("{}: {e}", path.display()) })?;
        let archive = zip::ZipArchive::new(std::io::BufReader::with_capacity(1 << 20, f))
            .map_err(|e| PatError::Format { reason: format!("{}: not a zip/npz archive: {e}", path.display()) })?;
        Ok(Npz { path: path.to_path_buf(), archive })
    }

    /// Member names without `.npy`, in archive order.
    pub fn names(&self) -> Vec<String> {
        self.archive.file_names().map(|n| n.trim_end_matches(".npy").to_string()).collect()
    }

    /// The file this archive was opened from.
    pub fn path(&self) -> &Path {
        &self.path
    }

    /// Read member `name` (without `.npy`).
    pub fn read(&mut self, name: &str) -> Result<NpzArray, PatError> {
        let path = self.path.clone();
        let mut f = self.archive.by_name(&format!("{name}.npy")).map_err(|e| ferr(&path, name, format!("missing: {e}")))?;
        let arr = read_npy(&mut f).map_err(|e| ferr(&path, name, e))?;
        debug_assert_eq!(arr.data.len(), arr.shape.iter().product::<usize>());
        Ok(arr)
    }

    fn typed<T>(&mut self, name: &str, want: &'static str, f: impl FnOnce(NpzData) -> Option<Vec<T>>) -> Result<(Vec<usize>, Vec<T>), PatError> {
        let a = self.read(name)?;
        let got = a.data.dtype();
        let shape = a.shape;
        match f(a.data) {
            Some(v) => Ok((shape, v)),
            None => Err(ferr(&self.path, name, format!("dtype {got}, expected {want}"))),
        }
    }

    /// Read an `<f8` member.
    pub fn read_f64(&mut self, name: &str) -> Result<(Vec<usize>, Vec<f64>), PatError> {
        self.typed(name, "<f8", |d| if let NpzData::F64(v) = d { Some(v) } else { None })
    }

    /// Read an `<i4` member.
    pub fn read_i32(&mut self, name: &str) -> Result<(Vec<usize>, Vec<i32>), PatError> {
        self.typed(name, "<i4", |d| if let NpzData::I32(v) = d { Some(v) } else { None })
    }

    /// Read a `|b1` member.
    pub fn read_bool(&mut self, name: &str) -> Result<(Vec<usize>, Vec<bool>), PatError> {
        self.typed(name, "|b1", |d| if let NpzData::Bool(v) = d { Some(v) } else { None })
    }
}

/// Write an uncompressed `.npz` (`np.savez` layout) of `f64`/`i32`/`bool` members: test fixtures
/// and the arrays handed to the JAX comparison script.
pub fn write_npz(path: &Path, members: &[(&str, &NpzArray)]) -> Result<(), PatError> {
    let werr = |e: std::io::Error| PatError::Format { reason: format!("{}: {e}", path.display()) };
    let f = std::fs::File::create(path).map_err(werr)?;
    let mut z = zip::ZipWriter::new(f);
    let opts = zip::write::FileOptions::default().compression_method(zip::CompressionMethod::Stored).large_file(true);
    for (name, a) in members {
        let n: usize = a.shape.iter().product();
        if n != a.data.len() {
            return Err(PatError::SizeMismatch { what: format!("npz member {name}: shape {:?} with {} values", a.shape, a.data.len()) });
        }
        let shape = match a.shape.len() {
            1 => format!("({},)", a.shape[0]),
            _ => format!("({})", a.shape.iter().map(|s| s.to_string()).collect::<Vec<_>>().join(", ")),
        };
        let mut header = format!("{{'descr': '{}', 'fortran_order': False, 'shape': {shape}, }}", a.data.dtype());
        while (10 + header.len() + 1) % 64 != 0 {
            header.push(' ');
        }
        header.push('\n');
        let mut b = Vec::with_capacity(10 + header.len() + 8 * n);
        b.extend_from_slice(b"\x93NUMPY\x01\x00");
        b.extend_from_slice(&(header.len() as u16).to_le_bytes());
        b.extend_from_slice(header.as_bytes());
        match &a.data {
            NpzData::F64(v) => v.iter().for_each(|x| b.extend_from_slice(&x.to_le_bytes())),
            NpzData::F32(v) => v.iter().for_each(|x| b.extend_from_slice(&x.to_le_bytes())),
            NpzData::I32(v) => v.iter().for_each(|x| b.extend_from_slice(&x.to_le_bytes())),
            NpzData::I64(v) => v.iter().for_each(|x| b.extend_from_slice(&x.to_le_bytes())),
            NpzData::Bool(v) => v.iter().for_each(|&x| b.push(x as u8)),
        }
        z.start_file(format!("{name}.npy"), opts).map_err(|e| PatError::Format { reason: format!("{}: {e}", path.display()) })?;
        std::io::Write::write_all(&mut z, &b).map_err(werr)?;
    }
    z.finish().map_err(|e| PatError::Format { reason: format!("{}: {e}", path.display()) })?;
    Ok(())
}
