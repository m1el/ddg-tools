//! Training data (SPEC.md §10.3; NOTES-REFERENCE.md §10): the reference's npz files, validated.
//!
//! Every file has the members `positions (S,P,3) f8`, `normals (S,P,3) f8`, `neighbors (S,P,k) i4`,
//! `masks (S,P) b1`, `mesh_vertices (S,V,3) f8`, `mesh_faces (S,F,3) i4`, `face_masks (S,F) b1`.
//! Checks (an error names the file, member, shape and row):
//! - the shapes agree;
//! - `masks` and `face_masks` are prefixes;
//! - for valid rows, neighbor column 0 is the row itself, the other entries are valid point
//!   indices or `−1`, and the `−1` entries are a suffix;
//! - valid faces index stored vertices;
//! - valid positions and normals are finite, normals non-zero.
//!
//! Normals of length 2 (the pretraining blob and the validation files, N§10) are normalized at load
//! and counted (D23) unless [`DataOptions::raw_normals`] (parity with the reference, whose fitting
//! uses the raw normal).

use super::npz::Npz;
use crate::PatError;
use std::path::Path;

/// Options of [`load_file`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DataOptions {
    /// Keep the stored normals unnormalized (the reference's behavior; D23). Default false.
    pub raw_normals: bool,
    /// A valid normal with `|‖n‖ − 1| > normal_tolerance` is counted as non-unit (and normalized
    /// unless `raw_normals`). Default 1e-3.
    pub normal_tolerance: f64,
}

impl Default for DataOptions {
    fn default() -> Self {
        DataOptions { raw_normals: false, normal_tolerance: 1e-3 }
    }
}

/// One training shape: a padded point cloud with neighborhoods and its watertight mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct TrainShape {
    /// Source file name (without directories) and index of the shape in it.
    pub source: (String, usize),
    /// Stored rows `P` (2048).
    pub rows: usize,
    /// Valid points `n ≤ P` (the mask prefix).
    pub n: usize,
    /// Neighbors per row `k` (64).
    pub k: usize,
    /// `P` positions (padding rows as stored, zeros).
    pub positions: Vec<[f64; 3]>,
    /// `P` normals: unit for valid rows unless `raw_normals`; padding as stored.
    pub normals: Vec<[f64; 3]>,
    /// `P·k` neighbor indices, row major, `−1` for none.
    pub neighbors: Vec<i32>,
    /// Mesh vertices (all stored rows; faces index them).
    pub vertices: Vec<[f64; 3]>,
    /// Valid faces only.
    pub faces: Vec<[u32; 3]>,
}

impl TrainShape {
    /// The neighbor row of point `i` (`k` entries, `−1` padded).
    pub fn neighbor_row(&self, i: usize) -> &[i32] {
        &self.neighbors[i * self.k..(i + 1) * self.k]
    }

    /// Number of valid neighbors in row `i` (the prefix before the first `−1`).
    pub fn valid_neighbors(&self, i: usize) -> usize {
        self.neighbor_row(i).iter().take_while(|&&j| j >= 0).count()
    }
}

/// What loading found.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct DataReport {
    /// Shapes read.
    pub shapes: usize,
    /// Shapes with `n < P` (padded clouds).
    pub padded_shapes: usize,
    /// Valid rows whose neighborhood has fewer than `k` valid entries.
    pub short_neighborhoods: usize,
    /// Valid normals with `|‖n‖ − 1| > tolerance` (normalized unless `raw_normals`).
    pub nonunit_normals: usize,
    /// Smallest and largest valid normal length seen.
    pub normal_length_range: (f64, f64),
    /// Shapes with no valid face of positive area (queries cannot be sampled on them).
    pub meshless_shapes: usize,
}

fn shape_err(path: &Path, what: String) -> PatError {
    PatError::Format { reason: format!("{}: {what}", path.display()) }
}

fn expect_shape(path: &Path, member: &str, got: &[usize], want: &[usize]) -> Result<(), PatError> {
    if got != want {
        return Err(shape_err(path, format!("member {member} has shape {got:?}, expected {want:?}")));
    }
    Ok(())
}

fn prefix_len(path: &Path, member: &str, s: usize, m: &[bool]) -> Result<usize, PatError> {
    let n = m.iter().take_while(|&&b| b).count();
    if let Some(off) = m[n..].iter().position(|&b| b) {
        return Err(shape_err(path, format!("member {member}, shape {s}: mask is not a prefix (false at {n}, true at {})", n + off)));
    }
    Ok(n)
}

/// Load and validate one npz file.
pub fn load_file(path: &Path, opts: &DataOptions) -> Result<(Vec<TrainShape>, DataReport), PatError> {
    if !(opts.normal_tolerance.is_finite() && opts.normal_tolerance >= 0.0) {
        return Err(PatError::InvalidOption { what: format!("DataOptions::normal_tolerance = {}", opts.normal_tolerance) });
    }
    let mut z = Npz::open(path)?;
    let (ps, pos) = z.read_f64("positions")?;
    if ps.len() != 3 || ps[2] != 3 {
        return Err(shape_err(path, format!("member positions has shape {ps:?}, expected (S, P, 3)")));
    }
    let (s_count, p) = (ps[0], ps[1]);
    let (ns, nrm) = z.read_f64("normals")?;
    expect_shape(path, "normals", &ns, &[s_count, p, 3])?;
    let (ms, masks) = z.read_bool("masks")?;
    expect_shape(path, "masks", &ms, &[s_count, p])?;
    let (nbs, nbr) = z.read_i32("neighbors")?;
    if nbs.len() != 3 || nbs[0] != s_count || nbs[1] != p || nbs[2] == 0 {
        return Err(shape_err(path, format!("member neighbors has shape {nbs:?}, expected ({s_count}, {p}, k)")));
    }
    let k = nbs[2];
    let (vs, verts) = z.read_f64("mesh_vertices")?;
    if vs.len() != 3 || vs[0] != s_count || vs[2] != 3 {
        return Err(shape_err(path, format!("member mesh_vertices has shape {vs:?}, expected ({s_count}, V, 3)")));
    }
    let v_rows = vs[1];
    let (fs, faces) = z.read_i32("mesh_faces")?;
    if fs.len() != 3 || fs[0] != s_count || fs[2] != 3 {
        return Err(shape_err(path, format!("member mesh_faces has shape {fs:?}, expected ({s_count}, F, 3)")));
    }
    let f_rows = fs[1];
    let (fms, fmask) = z.read_bool("face_masks")?;
    expect_shape(path, "face_masks", &fms, &[s_count, f_rows])?;

    let fname = path.file_name().map(|f| f.to_string_lossy().to_string()).unwrap_or_default();
    let mut report = DataReport { normal_length_range: (f64::INFINITY, 0.0), ..Default::default() };
    let mut shapes = Vec::with_capacity(s_count);
    let v3 = |v: &[f64], i: usize| [v[3 * i], v[3 * i + 1], v[3 * i + 2]];
    for s in 0..s_count {
        let n = prefix_len(path, "masks", s, &masks[s * p..(s + 1) * p])?;
        let nf = prefix_len(path, "face_masks", s, &fmask[s * f_rows..(s + 1) * f_rows])?;
        let mut positions: Vec<[f64; 3]> = (0..p).map(|i| v3(&pos, s * p + i)).collect();
        let mut normals: Vec<[f64; 3]> = (0..p).map(|i| v3(&nrm, s * p + i)).collect();
        for i in 0..n {
            if !positions[i].iter().all(|x| x.is_finite()) {
                return Err(shape_err(path, format!("shape {s}: position {i} is not finite")));
            }
            let nn = normals[i];
            let l = (nn[0] * nn[0] + nn[1] * nn[1] + nn[2] * nn[2]).sqrt();
            if !l.is_finite() || l == 0.0 {
                return Err(shape_err(path, format!("shape {s}: normal {i} has length {l}")));
            }
            report.normal_length_range.0 = report.normal_length_range.0.min(l);
            report.normal_length_range.1 = report.normal_length_range.1.max(l);
            if (l - 1.0).abs() > opts.normal_tolerance {
                report.nonunit_normals += 1;
            }
            if !opts.raw_normals {
                normals[i] = [nn[0] / l, nn[1] / l, nn[2] / l];
            }
        }
        // padding rows stay as stored (zeros); they never reach the network unless
        // `reference_padding`
        positions.truncate(p);
        let neighbors = nbr[s * p * k..(s + 1) * p * k].to_vec();
        for i in 0..n {
            let row = &neighbors[i * k..(i + 1) * k];
            if row[0] != i as i32 {
                return Err(shape_err(path, format!("shape {s}: neighbor row {i} starts with {}, not with itself", row[0])));
            }
            let valid = row.iter().take_while(|&&j| j >= 0).count();
            if let Some(&j) = row[valid..].iter().find(|&&j| j != -1) {
                return Err(shape_err(path, format!("shape {s}: neighbor row {i} has entry {j} after a −1")));
            }
            if let Some(&j) = row[..valid].iter().find(|&&j| j as usize >= n) {
                return Err(shape_err(path, format!("shape {s}: neighbor row {i} has index {j} ≥ n = {n}")));
            }
            if valid < k {
                report.short_neighborhoods += 1;
            }
        }
        let vertices: Vec<[f64; 3]> = (0..v_rows).map(|i| v3(&verts, s * v_rows + i)).collect();
        let mut fv = Vec::with_capacity(nf);
        let mut area_pos = false;
        for f in 0..nf {
            let t = &faces[(s * f_rows + f) * 3..(s * f_rows + f) * 3 + 3];
            if t.iter().any(|&j| j < 0 || j as usize >= v_rows) {
                return Err(shape_err(path, format!("shape {s}: face {f} = {t:?} indexes outside 0..{v_rows}")));
            }
            let tri = [t[0] as u32, t[1] as u32, t[2] as u32];
            let [a, b, c] = tri.map(|j| vertices[j as usize]);
            if !a.iter().chain(&b).chain(&c).all(|x| x.is_finite()) {
                return Err(shape_err(path, format!("shape {s}: face {f} has a non-finite vertex")));
            }
            let e1 = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
            let e2 = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
            let cr = [e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]];
            area_pos |= cr.iter().any(|&x| x != 0.0);
            fv.push(tri);
        }
        report.meshless_shapes += (!area_pos) as usize;
        report.padded_shapes += (n < p) as usize;
        shapes.push(TrainShape { source: (fname.clone(), s), rows: p, n, k, positions, normals, neighbors, vertices, faces: fv });
    }
    report.shapes = s_count;
    if report.normal_length_range.0 > report.normal_length_range.1 {
        report.normal_length_range = (0.0, 0.0);
    }
    Ok((shapes, report))
}

impl DataReport {
    /// Merge the report of another file.
    pub fn merge(&mut self, o: &DataReport) {
        self.shapes += o.shapes;
        self.padded_shapes += o.padded_shapes;
        self.short_neighborhoods += o.short_neighborhoods;
        self.nonunit_normals += o.nonunit_normals;
        self.meshless_shapes += o.meshless_shapes;
        if o.shapes > 0 && o.normal_length_range.1 > 0.0 {
            if self.normal_length_range.1 == 0.0 {
                self.normal_length_range = o.normal_length_range;
            } else {
                self.normal_length_range.0 = self.normal_length_range.0.min(o.normal_length_range.0);
                self.normal_length_range.1 = self.normal_length_range.1.max(o.normal_length_range.1);
            }
        }
    }
}
