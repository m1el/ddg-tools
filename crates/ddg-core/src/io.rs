//! Reading and writing Wavefront OBJ meshes, and reading the cone files of the optimal-cones
//! dataset.
//!
//! [`ObjMesh`] keeps faces as polygons in the same layout Blender uses (and
//! [`crate::triangulate`] consumes): one vertex index per polygon corner plus the start offset
//! of every polygon.

use std::fmt;
use std::fmt::Write as _;
use std::path::Path;

use crate::topology::NONE;

/// A file could not be read or parsed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FileError {
    /// The operating system reported an error.
    Io {
        /// The file.
        path: String,
        /// The system's message.
        message: String,
    },
    /// The content is malformed.
    Parse {
        /// The file, if the text came from one.
        path: Option<String>,
        /// 1-based line number.
        line: usize,
        /// What is wrong with that line.
        message: String,
    },
    /// The arguments of a writer are inconsistent.
    InvalidInput {
        /// What is wrong.
        message: String,
    },
}

impl fmt::Display for FileError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FileError::Io { path, message } => write!(f, "{path}: {message}"),
            FileError::Parse { path: Some(path), line, message } => write!(f, "{path}:{line}: {message}"),
            FileError::Parse { path: None, line, message } => write!(f, "line {line}: {message}"),
            FileError::InvalidInput { message } => f.write_str(message),
        }
    }
}

impl std::error::Error for FileError {}

fn parse_error(line: usize, message: String) -> FileError {
    FileError::Parse { path: None, line, message }
}

fn with_path(error: FileError, path: &Path) -> FileError {
    match error {
        FileError::Parse { line, message, .. } => {
            FileError::Parse { path: Some(path.display().to_string()), line, message }
        }
        other => other,
    }
}

fn read_file(path: &Path) -> Result<String, FileError> {
    std::fs::read_to_string(path)
        .map_err(|e| FileError::Io { path: path.display().to_string(), message: e.to_string() })
}

/// The geometry and polygons of an OBJ file. Normals, groups, materials and other statements
/// are ignored.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ObjMesh {
    /// Vertex positions (`v` lines; a fourth coordinate or vertex colours are ignored).
    pub positions: Vec<[f64; 3]>,
    /// Texture coordinates (`vt` lines; a missing `v` coordinate is 0, a third one is ignored).
    pub texcoords: Vec<[f64; 2]>,
    /// For every polygon corner, its 0-based vertex index.
    pub corner_verts: Vec<u32>,
    /// For every polygon corner, its 0-based index into [`ObjMesh::texcoords`], or
    /// [`NONE`](crate::topology::NONE) if the corner has none.
    pub corner_texcoords: Vec<u32>,
    /// For every polygon, the index of its first corner. Polygon `f` owns the corners
    /// `face_starts[f] .. face_starts[f + 1]` (up to `corner_verts.len()` for the last one).
    pub face_starts: Vec<u32>,
}

impl ObjMesh {
    /// Number of polygons.
    pub fn n_faces(&self) -> usize {
        self.face_starts.len()
    }

    /// The corner index range of polygon `f`.
    pub fn face_corners(&self, f: usize) -> std::ops::Range<usize> {
        let start = self.face_starts[f] as usize;
        let end = self.face_starts.get(f + 1).map_or(self.corner_verts.len(), |&e| e as usize);
        start..end
    }

    /// The faces as triangles, or `None` if any polygon is not a triangle.
    pub fn triangles(&self) -> Option<Vec<[u32; 3]>> {
        if self.corner_verts.len() != 3 * self.n_faces() {
            return None;
        }
        let mut out = Vec::with_capacity(self.n_faces());
        for f in 0..self.n_faces() {
            let r = self.face_corners(f);
            if r.len() != 3 {
                return None;
            }
            out.push([self.corner_verts[r.start], self.corner_verts[r.start + 1], self.corner_verts[r.start + 2]]);
        }
        Some(out)
    }

    /// The texture coordinate of every polygon corner, or `None` if some corner has none.
    pub fn corner_uvs(&self) -> Option<Vec<[f64; 2]>> {
        self.corner_texcoords
            .iter()
            .map(|&t| (t != NONE).then(|| self.texcoords[t as usize]))
            .collect()
    }
}

/// Reads an OBJ file. See [`parse_obj`].
pub fn read_obj(path: impl AsRef<Path>) -> Result<ObjMesh, FileError> {
    let path = path.as_ref();
    parse_obj(&read_file(path)?).map_err(|e| with_path(e, path))
}

/// Parses OBJ text.
///
/// Supported: `v`, `vt`, and polygonal `f` statements with corners of the form `v`, `v/vt`,
/// `v//vn` and `v/vt/vn`; 1-based and negative (relative) indices; comments; line continuation
/// with a trailing backslash. Normals (`vn` lines and normal indices) and every other statement
/// are skipped.
///
/// # Errors
/// [`FileError::Parse`] with the line number for: a malformed or non-finite number, a `v` with
/// fewer than three coordinates, a face with fewer than three corners, an index that is zero or
/// refers to a vertex or texture coordinate that does not exist.
pub fn parse_obj(text: &str) -> Result<ObjMesh, FileError> {
    let mut mesh = ObjMesh::default();
    // Positive indices may legally refer to elements defined later in the file, so they are
    // range-checked at the end; remember where each face came from for the message.
    let mut face_lines: Vec<usize> = Vec::new();
    let mut joined = String::new();
    let mut pending_start: Option<usize> = None;

    for (i, raw) in text.lines().enumerate() {
        let number = i + 1;
        let raw = raw.trim_end();
        let (line, line_number): (&str, usize) = if let Some(stripped) = raw.strip_suffix('\\') {
            joined.push_str(stripped);
            joined.push(' ');
            pending_start.get_or_insert(number);
            continue;
        } else if let Some(start) = pending_start.take() {
            joined.push_str(raw);
            (joined.as_str(), start)
        } else {
            (raw, number)
        };
        let result = parse_obj_line(line, line_number, &mut mesh, &mut face_lines);
        joined.clear();
        result?;
    }
    if let Some(start) = pending_start {
        let line = std::mem::take(&mut joined);
        parse_obj_line(&line, start, &mut mesh, &mut face_lines)?;
    }

    for (f, &line) in face_lines.iter().enumerate() {
        for c in mesh.face_corners(f) {
            let v = mesh.corner_verts[c];
            if v as usize >= mesh.positions.len() {
                return Err(parse_error(
                    line,
                    format!("face refers to vertex {}, but the file defines {} vertices", v + 1, mesh.positions.len()),
                ));
            }
            let t = mesh.corner_texcoords[c];
            if t != NONE && t as usize >= mesh.texcoords.len() {
                return Err(parse_error(
                    line,
                    format!(
                        "face refers to texture coordinate {}, but the file defines {}",
                        t + 1,
                        mesh.texcoords.len()
                    ),
                ));
            }
        }
    }
    Ok(mesh)
}

fn parse_obj_line(
    line: &str,
    line_number: usize,
    mesh: &mut ObjMesh,
    face_lines: &mut Vec<usize>,
) -> Result<(), FileError> {
    let line = line.split('#').next().unwrap_or("");
    let mut tokens = line.split_ascii_whitespace();
    let Some(keyword) = tokens.next() else {
        return Ok(());
    };
    let number = |token: &str| -> Result<f64, FileError> {
        match token.parse::<f64>() {
            Ok(x) if x.is_finite() => Ok(x),
            Ok(_) => Err(parse_error(line_number, format!("non-finite number '{token}'"))),
            Err(_) => Err(parse_error(line_number, format!("cannot parse '{token}' as a number"))),
        }
    };
    match keyword {
        "v" => {
            let mut p = [0.0; 3];
            for (k, slot) in p.iter_mut().enumerate() {
                let token = tokens.next().ok_or_else(|| {
                    parse_error(line_number, format!("vertex needs 3 coordinates, found {k}"))
                })?;
                *slot = number(token)?;
            }
            mesh.positions.push(p);
        }
        "vt" => {
            let u = tokens
                .next()
                .ok_or_else(|| parse_error(line_number, "texture coordinate needs at least 1 value".to_string()))?;
            let u = number(u)?;
            let v = tokens.next().map(number).transpose()?.unwrap_or(0.0);
            mesh.texcoords.push([u, v]);
        }
        "f" => {
            let start = mesh.corner_verts.len();
            for token in tokens {
                let mut parts = token.split('/');
                let v = parts.next().unwrap_or("");
                let v = resolve_index(v, mesh.positions.len(), "vertex", line_number)?;
                let t = match parts.next() {
                    None | Some("") => NONE,
                    Some(t) => resolve_index(t, mesh.texcoords.len(), "texture coordinate", line_number)?,
                };
                mesh.corner_verts.push(v);
                mesh.corner_texcoords.push(t);
            }
            let size = mesh.corner_verts.len() - start;
            if size < 3 {
                return Err(parse_error(line_number, format!("face has {size} corners, at least 3 are required")));
            }
            if mesh.corner_verts.len() >= NONE as usize {
                return Err(parse_error(line_number, "more than 2^32 - 2 face corners".to_string()));
            }
            mesh.face_starts.push(start as u32);
            face_lines.push(line_number);
        }
        _ => {}
    }
    Ok(())
}

/// Converts a 1-based or negative OBJ index to a 0-based one. `defined` is the number of
/// elements defined so far, which is what negative indices count back from.
fn resolve_index(token: &str, defined: usize, what: &str, line: usize) -> Result<u32, FileError> {
    let index: i64 = token
        .parse()
        .map_err(|_| parse_error(line, format!("cannot parse '{token}' as a {what} index")))?;
    if index > 0 {
        if index >= i64::from(NONE) {
            return Err(parse_error(line, format!("{what} index {index} is too large")));
        }
        Ok((index - 1) as u32)
    } else if index < 0 {
        let back = index.unsigned_abs();
        if back > defined as u64 {
            return Err(parse_error(
                line,
                format!("relative {what} index {index} reaches before the first of the {defined} defined so far"),
            ));
        }
        Ok((defined as u64 - back) as u32)
    } else {
        Err(parse_error(line, format!("{what} index 0 is invalid: OBJ indices start at 1")))
    }
}

/// Formats a polygon mesh as OBJ text.
///
/// `corner_verts` and `face_starts` use the layout of [`ObjMesh`]. If `corner_uvs` is given it
/// must hold one UV per polygon corner; one `vt` line is written per corner and referenced as
/// `v/vt`. Numbers are written with the shortest representation that reads back exactly.
///
/// # Errors
/// [`FileError::InvalidInput`] if the arrays are inconsistent (bad face offsets, a face with
/// fewer than three corners, an out-of-range vertex index, wrong number of UVs) or a number is
/// not finite.
pub fn obj_to_string(
    positions: &[[f64; 3]],
    corner_verts: &[u32],
    face_starts: &[u32],
    corner_uvs: Option<&[[f64; 2]]>,
) -> Result<String, FileError> {
    let invalid = |message: String| Err(FileError::InvalidInput { message });
    if let Some(uvs) = corner_uvs
        && uvs.len() != corner_verts.len()
    {
        return invalid(format!("{} corner UVs given for {} corners", uvs.len(), corner_verts.len()));
    }
    if let Some(v) = positions.iter().position(|p| !p.iter().all(|x| x.is_finite())) {
        return invalid(format!("position of vertex {v} is not finite"));
    }
    if let Some(c) = corner_uvs.and_then(|uvs| uvs.iter().position(|p| !p.iter().all(|x| x.is_finite()))) {
        return invalid(format!("UV of corner {c} is not finite"));
    }
    if face_starts.is_empty() != corner_verts.is_empty() {
        return invalid(format!("{} faces but {} corners", face_starts.len(), corner_verts.len()));
    }
    let mut out = String::new();
    out.push_str("# written by ddg-core\n");
    for p in positions {
        let _ = writeln!(out, "v {} {} {}", p[0], p[1], p[2]);
    }
    if let Some(uvs) = corner_uvs {
        for uv in uvs {
            let _ = writeln!(out, "vt {} {}", uv[0], uv[1]);
        }
    }
    for (f, &start) in face_starts.iter().enumerate() {
        let start = start as usize;
        let end = face_starts.get(f + 1).map_or(corner_verts.len(), |&e| e as usize);
        if (f == 0 && start != 0) || end > corner_verts.len() || end < start + 3 {
            return invalid(format!(
                "face {f} has corner range {start}..{end}: faces must start at 0, be contiguous and have at \
                 least 3 corners ({} corners in total)",
                corner_verts.len()
            ));
        }
        out.push('f');
        for c in start..end {
            let v = corner_verts[c];
            if v as usize >= positions.len() {
                return invalid(format!(
                    "corner {c} of face {f} refers to vertex {v}, but there are {} positions",
                    positions.len()
                ));
            }
            if corner_uvs.is_some() {
                let _ = write!(out, " {}/{}", v + 1, c + 1);
            } else {
                let _ = write!(out, " {}", v + 1);
            }
        }
        out.push('\n');
    }
    Ok(out)
}

/// Writes a polygon mesh to an OBJ file. See [`obj_to_string`].
pub fn write_obj(
    path: impl AsRef<Path>,
    positions: &[[f64; 3]],
    corner_verts: &[u32],
    face_starts: &[u32],
    corner_uvs: Option<&[[f64; 2]]>,
) -> Result<(), FileError> {
    let text = obj_to_string(positions, corner_verts, face_starts, corner_uvs)?;
    let path = path.as_ref();
    std::fs::write(path, text).map_err(|e| FileError::Io { path: path.display().to_string(), message: e.to_string() })
}

/// A cone singularity: a vertex and its cone angle.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cone {
    /// 0-based vertex index.
    pub vertex: u32,
    /// Cone angle `µ = 2π − (angle sum around the vertex in the flattened metric)`, in radians.
    /// Positive for a cone tip.
    pub angle: f64,
}

/// Reads a cone file of the optimal-cones dataset. See [`parse_cones`].
pub fn read_cones(path: impl AsRef<Path>, n_vertices: usize) -> Result<Vec<Cone>, FileError> {
    let path = path.as_ref();
    parse_cones(&read_file(path)?, n_vertices).map_err(|e| with_path(e, path))
}

/// Parses a cone file: one `index angle` pair per line, where `index` is a **1-based** vertex
/// index and `angle` the cone angle in radians. Blank lines and `#` comments are skipped. The
/// returned cones use 0-based indices and keep the file's order.
///
/// # Errors
/// [`FileError::Parse`] with the line number for a malformed line, an index that is 0 or larger
/// than `n_vertices`, a non-finite angle, or a vertex listed twice.
pub fn parse_cones(text: &str, n_vertices: usize) -> Result<Vec<Cone>, FileError> {
    let mut cones = Vec::new();
    let mut first_line: std::collections::HashMap<u32, usize> = std::collections::HashMap::new();
    for (i, raw) in text.lines().enumerate() {
        let line_number = i + 1;
        let line = raw.split('#').next().unwrap_or("");
        let mut tokens = line.split_ascii_whitespace();
        let Some(index) = tokens.next() else {
            continue;
        };
        let (Some(angle), None) = (tokens.next(), tokens.next()) else {
            return Err(parse_error(line_number, format!("expected 'index angle', found '{}'", raw.trim())));
        };
        let index: u64 = index
            .parse()
            .map_err(|_| parse_error(line_number, format!("cannot parse '{index}' as a vertex index")))?;
        if index == 0 {
            return Err(parse_error(line_number, "vertex index 0 is invalid: cone files are 1-based".to_string()));
        }
        if index > n_vertices as u64 {
            return Err(parse_error(
                line_number,
                format!("vertex index {index} exceeds the number of vertices ({n_vertices})"),
            ));
        }
        let angle = match angle.parse::<f64>() {
            Ok(a) if a.is_finite() => a,
            _ => return Err(parse_error(line_number, format!("cannot parse '{angle}' as a finite angle"))),
        };
        let vertex = (index - 1) as u32;
        if let Some(previous) = first_line.insert(vertex, line_number) {
            return Err(parse_error(
                line_number,
                format!("vertex {index} already has a cone (line {previous})"),
            ));
        }
        cones.push(Cone { vertex, angle });
    }
    Ok(cones)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_all_corner_forms_and_negative_indices() {
        let text = "\
# a comment
o thing
v 0 0 0
v 1 0 0 1.0
v 1 1 0   # trailing comment
v 0 1 0
vt 0 0
vt 1 0 0
vt 1 1
vt 0.5
vn 0 0 1
g group
f 1 2 3
f 1/1 2/2 3/3 4/4
f 1//1 3//1 4//1
f -4/-4/-1 -3/-3/-1 -1/1/1
s off
f 1 2 \\
  3 4
";
        let m = parse_obj(text).unwrap();
        assert_eq!(m.positions.len(), 4);
        assert_eq!(m.positions[2], [1.0, 1.0, 0.0]);
        assert_eq!(m.texcoords, vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.5, 0.0]]);
        assert_eq!(m.n_faces(), 5);
        assert_eq!(m.face_starts, vec![0, 3, 7, 10, 13]);
        assert_eq!(m.corner_verts, vec![0, 1, 2, 0, 1, 2, 3, 0, 2, 3, 0, 1, 3, 0, 1, 2, 3]);
        assert_eq!(
            m.corner_texcoords,
            vec![NONE, NONE, NONE, 0, 1, 2, 3, NONE, NONE, NONE, 0, 1, 0, NONE, NONE, NONE, NONE]
        );
        assert_eq!(m.face_corners(1), 3..7);
        assert_eq!(m.face_corners(4), 13..17);
        assert_eq!(m.triangles(), None);
        assert_eq!(m.corner_uvs(), None);
    }

    #[test]
    fn triangles_and_corner_uvs() {
        let m = parse_obj("v 0 0 0\nv 1 0 0\nv 0 1 0\nvt 0 0\nvt 1 0\nvt 0 1\nf 1/1 2/2 3/3\n").unwrap();
        assert_eq!(m.triangles(), Some(vec![[0, 1, 2]]));
        assert_eq!(m.corner_uvs(), Some(vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]));
    }

    #[test]
    fn parse_errors_name_the_line() {
        let err = |text: &str| parse_obj(text).unwrap_err().to_string();
        assert_eq!(err("v 0 0\n"), "line 1: vertex needs 3 coordinates, found 2");
        assert_eq!(err("v 0 0 0\nv 0 x 0\n"), "line 2: cannot parse 'x' as a number");
        assert_eq!(err("v 0 0 nan\n"), "line 1: non-finite number 'nan'");
        assert_eq!(err("v 0 0 0\nv 1 0 0\nf 1 2\n"), "line 3: face has 2 corners, at least 3 are required");
        assert_eq!(
            err("v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 0\n"),
            "line 4: vertex index 0 is invalid: OBJ indices start at 1"
        );
        assert_eq!(
            err("v 0 0 0\nv 1 0 0\nv 0 1 0\n\nf 1 2 4\n"),
            "line 5: face refers to vertex 4, but the file defines 3 vertices"
        );
        assert_eq!(
            err("v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1/1 2/1 3/1\n"),
            "line 4: face refers to texture coordinate 1, but the file defines 0"
        );
        assert_eq!(
            err("v 0 0 0\nv 1 0 0\nv 0 1 0\nf -1 -2 -4\n"),
            "line 4: relative vertex index -4 reaches before the first of the 3 defined so far"
        );
        assert_eq!(err("v 0 0 0\nf 1 a 1\n"), "line 2: cannot parse 'a' as a vertex index");
        let e = read_obj("/nonexistent/dir/mesh.obj").unwrap_err();
        assert!(matches!(e, FileError::Io { .. }));
        assert!(e.to_string().starts_with("/nonexistent/dir/mesh.obj: "));
    }

    #[test]
    fn positive_indices_may_refer_forward() {
        let m = parse_obj("f 1 2 3\nv 0 0 0\nv 1 0 0\nv 0 1 0\n").unwrap();
        assert_eq!(m.triangles(), Some(vec![[0, 1, 2]]));
    }

    #[test]
    fn obj_round_trip_with_corner_uvs() {
        let positions = vec![[0.1, -2.5e-7, 3.0], [1.0 / 3.0, 2.0, 1e20], [0.0, 1.0, -0.0], [5.0, 5.5, 5.25]];
        let corner_verts = vec![0, 1, 2, 0, 2, 3, 1];
        let face_starts = vec![0, 3];
        let uvs: Vec<[f64; 2]> = (0..7).map(|c| [c as f64 / 7.0, 1.0 - (c as f64).sqrt()]).collect();
        let text = obj_to_string(&positions, &corner_verts, &face_starts, Some(&uvs)).unwrap();
        let m = parse_obj(&text).unwrap();
        assert_eq!(m.positions, positions);
        assert_eq!(m.corner_verts, corner_verts);
        assert_eq!(m.face_starts, face_starts);
        assert_eq!(m.corner_uvs().unwrap(), uvs);

        let text = obj_to_string(&positions, &corner_verts, &face_starts, None).unwrap();
        let m = parse_obj(&text).unwrap();
        assert_eq!(m.positions, positions);
        assert_eq!(m.corner_verts, corner_verts);
        assert!(m.texcoords.is_empty());
        assert!(m.corner_texcoords.iter().all(|&t| t == NONE));

        let dir = std::env::temp_dir().join(format!("ddg-core-io-test-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("round_trip.obj");
        write_obj(&path, &positions, &corner_verts, &face_starts, Some(&uvs)).unwrap();
        let m = read_obj(&path).unwrap();
        assert_eq!(m.corner_uvs().unwrap(), uvs);
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn writer_rejects_inconsistent_input() {
        let p = vec![[0.0; 3]; 3];
        let err = |r: Result<String, FileError>| r.unwrap_err().to_string();
        assert_eq!(err(obj_to_string(&p, &[0, 1, 2], &[0], Some(&[[0.0; 2]; 2]))), "2 corner UVs given for 3 corners");
        assert_eq!(
            err(obj_to_string(&p, &[0, 1, 3], &[0], None)),
            "corner 2 of face 0 refers to vertex 3, but there are 3 positions"
        );
        assert!(err(obj_to_string(&p, &[0, 1, 2, 0], &[0, 2], None)).starts_with("face 0 has corner range 0..2"));
        assert!(err(obj_to_string(&p, &[0, 1, 2], &[1], None)).starts_with("face 0 has corner range 1..3"));
        assert_eq!(
            err(obj_to_string(&[[0.0, f64::NAN, 0.0]], &[], &[], None)),
            "position of vertex 0 is not finite"
        );
        assert_eq!(obj_to_string(&[], &[], &[], None).unwrap(), "# written by ddg-core\n");
    }

    #[test]
    fn cones_are_one_based() {
        let cones = parse_cones("3 0.5\n\n# comment\n1 -1.25e-1  # trailing\n", 3).unwrap();
        assert_eq!(cones, vec![Cone { vertex: 2, angle: 0.5 }, Cone { vertex: 0, angle: -0.125 }]);
        let err = |text: &str| parse_cones(text, 3).unwrap_err().to_string();
        assert_eq!(err("0 1.0\n"), "line 1: vertex index 0 is invalid: cone files are 1-based");
        assert_eq!(err("1 1.0\n4 1.0\n"), "line 2: vertex index 4 exceeds the number of vertices (3)");
        assert_eq!(err("1 1.0 2.0\n"), "line 1: expected 'index angle', found '1 1.0 2.0'");
        assert_eq!(err("1\n"), "line 1: expected 'index angle', found '1'");
        assert_eq!(err("x 1.0\n"), "line 1: cannot parse 'x' as a vertex index");
        assert_eq!(err("1 inf\n"), "line 1: cannot parse 'inf' as a finite angle");
        assert_eq!(err("2 1.0\n1 0.5\n2 0.25\n"), "line 3: vertex 2 already has a cone (line 1)");
    }
}
