//! Wavefront OBJ I/O for curve networks: `v` lines are vertices, `l` lines with two or more
//! indices are polylines (`SPEC.md` §9.1, §16). Faces and everything else are ignored, so that a
//! mesh file reads as its vertices without edges (and then fails with `IsolatedVertex`).

use std::path::Path;

use crate::error::RcError;
use crate::network::CurveNetwork;

/// Parses the text of a curve OBJ file and validates it (topology and geometry, `SPEC.md` §2.3).
///
/// `l a b c` contributes the edges `(a, b)` and `(b, c)`; a closed loop repeats its first index
/// at the end or is written as several two-index `l` lines. Indices are 1-based, negative
/// indices count from the end of the vertex list read so far, and `a/b/c` forms use `a`.
///
/// # Errors
/// [`RcError::Parse`] with the 1-based line number for malformed `v` or `l` lines or indices out
/// of range; the validation errors of [`CurveNetwork::new`] and
/// [`CurveNetwork::validate_positions`] otherwise.
pub fn parse_curve_obj(text: &str) -> Result<(CurveNetwork, Vec<[f64; 3]>), RcError> {
    let mut positions: Vec<[f64; 3]> = Vec::new();
    let mut edges: Vec<[u32; 2]> = Vec::new();
    for (ln, raw) in text.lines().enumerate() {
        let line_no = ln + 1;
        let line = raw.trim();
        let mut it = line.split_whitespace();
        let Some(tag) = it.next() else { continue };
        let parse_err = |message: String| RcError::Parse {
            line: line_no,
            message,
        };
        match tag {
            "v" => {
                let mut c = [0.0f64; 3];
                for (k, slot) in c.iter_mut().enumerate() {
                    let tok = it
                        .next()
                        .ok_or_else(|| parse_err(format!("vertex has only {k} coordinates")))?;
                    *slot = tok
                        .parse::<f64>()
                        .map_err(|_| parse_err(format!("cannot read coordinate {tok:?}")))?;
                }
                positions.push(c);
            }
            "l" => {
                let mut idx: Vec<u32> = Vec::new();
                for tok in it {
                    let head = tok.split('/').next().unwrap_or("");
                    let i: i64 = head
                        .parse()
                        .map_err(|_| parse_err(format!("cannot read vertex index {tok:?}")))?;
                    let n = positions.len() as i64;
                    let zero_based = if i > 0 {
                        i - 1
                    } else if i < 0 {
                        n + i
                    } else {
                        return Err(parse_err("vertex index 0 (OBJ indices start at 1)".into()));
                    };
                    if zero_based < 0 || zero_based >= n {
                        return Err(parse_err(format!(
                            "vertex index {i} is out of range ({n} vertices read so far)"
                        )));
                    }
                    idx.push(zero_based as u32);
                }
                if idx.len() < 2 {
                    return Err(parse_err(format!(
                        "polyline has {} vertices, at least 2 are needed",
                        idx.len()
                    )));
                }
                for w in idx.windows(2) {
                    edges.push([w[0], w[1]]);
                }
            }
            _ => {}
        }
    }
    let network = CurveNetwork::new(positions.len(), edges)?;
    network.validate_positions(&positions)?;
    Ok((network, positions))
}

/// Reads and validates a curve OBJ file (see [`parse_curve_obj`]).
///
/// # Errors
/// [`RcError::Io`] if the file cannot be read, otherwise the errors of [`parse_curve_obj`].
pub fn read_curve_obj(path: impl AsRef<Path>) -> Result<(CurveNetwork, Vec<[f64; 3]>), RcError> {
    let path = path.as_ref();
    let text = std::fs::read_to_string(path).map_err(|e| RcError::Io {
        path: path.display().to_string(),
        message: e.to_string(),
    })?;
    parse_curve_obj(&text)
}

/// Writes the network as `v` lines and one two-index `l` line per edge. Coordinates are printed
/// with the shortest representation that reads back exactly.
///
/// # Errors
/// [`RcError::PositionCountMismatch`] for a wrong number of positions, [`RcError::Io`] on write
/// failure.
pub fn write_curve_obj(
    path: impl AsRef<Path>,
    network: &CurveNetwork,
    positions: &[[f64; 3]],
) -> Result<(), RcError> {
    use std::fmt::Write;
    if positions.len() != network.n_vertices() {
        return Err(RcError::PositionCountMismatch {
            expected: network.n_vertices(),
            found: positions.len(),
        });
    }
    let mut s = String::new();
    for p in positions {
        let _ = writeln!(s, "v {} {} {}", p[0], p[1], p[2]);
    }
    for &[a, b] in network.edges() {
        let _ = writeln!(s, "l {} {}", a + 1, b + 1);
    }
    let path = path.as_ref();
    std::fs::write(path, s).map_err(|e| RcError::Io {
        path: path.display().to_string(),
        message: e.to_string(),
    })
}
