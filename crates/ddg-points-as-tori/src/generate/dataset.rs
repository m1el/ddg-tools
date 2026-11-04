//! Training and validation files in the authors' layout (SPEC.md §10.3, NOTES-REFERENCE §10):
//!
//! | member | dtype | shape |
//! |---|---|---|
//! | `positions`, `normals` | f8 | `(S, P, 3)` |
//! | `neighbors` | i4 | `(S, P, k)`: self first, then the `k − 1` nearest (our kd-tree, `(d², index)` order), `−1` padding |
//! | `masks` | b1 | `(S, P)`: true prefix of the `n ≤ P` valid points |
//! | `mesh_vertices` | f8 | `(S, V, 3)`: zero padding |
//! | `mesh_faces` | i4 | `(S, F, 3)`: `−1` padding |
//! | `face_masks` | b1 | `(S, F)` |
//!
//! with `P = V = F = 2048`, `k = 64`. Unlike the reference's blob files, normals are unit length
//! (its length-2 normals are an `+=` bug; D23). Shapes are generated in parallel, each from its
//! own seeded stream, so a file is bitwise independent of the thread count.

use super::blob::{BlobOptions, MeshData, generate_blob};
use super::bvh::TriBvh;
use super::npz::{NpzOut, write_npz};
use super::rng::Rng;
use super::sampling::{farthest_point_sampling, sample_mesh_uniform};
use super::scanner::{ScanOptions, simulate_scan};
use crate::PatError;
use crate::spatial::{KdTree, KdTreeOptions, fitting_neighborhoods};
use rayon::prelude::*;
use std::path::Path;

/// How the cloud of a shape is sampled.
#[derive(Debug, Clone, PartialEq)]
pub enum CloudMode {
    /// Mode 0: 50/50 between area-uniform samples and farthest points of `4 n` uniform samples.
    Uniform,
    /// Mode −1: the simulated scanner.
    Scan(ScanOptions),
}

/// One shape: a cloud and its mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct ShapeRecord {
    /// Cloud positions (`n ≤ P`).
    pub positions: Vec<[f64; 3]>,
    /// Unit normals.
    pub normals: Vec<[f64; 3]>,
    /// Mesh vertices.
    pub vertices: Vec<[f64; 3]>,
    /// Mesh faces.
    pub faces: Vec<[u32; 3]>,
}

/// Sizes of a file.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Layout {
    /// Rows of points `P` (2048).
    pub points: usize,
    /// Neighbors per row `k` (64).
    pub k: usize,
    /// Rows of mesh vertices `V` (2048).
    pub vertices: usize,
    /// Rows of mesh faces `F` (2048).
    pub faces: usize,
}

impl Default for Layout {
    fn default() -> Self {
        Layout { points: 2048, k: 64, vertices: 2048, faces: 2048 }
    }
}

/// Sample a cloud of at most `n` points on a closed mesh (`seed`: the stream).
pub fn sample_cloud(mesh: &MeshData, mode: &CloudMode, n: usize, seed: u64) -> Result<(Vec<[f64; 3]>, Vec<[f64; 3]>), PatError> {
    let mut rng = Rng::new(seed);
    match mode {
        CloudMode::Uniform => {
            if rng.f64() < 0.5 {
                sample_mesh_uniform(&mesh.vertices, &mesh.faces, n, &mut rng)
            } else {
                let (p, nr) = sample_mesh_uniform(&mesh.vertices, &mesh.faces, 4 * n, &mut rng)?;
                let idx = farthest_point_sampling(&p, n, 0);
                Ok((idx.iter().map(|&i| p[i]).collect(), idx.iter().map(|&i| nr[i]).collect()))
            }
        }
        CloudMode::Scan(opts) => {
            let bvh = TriBvh::new(&mesh.vertices, &mesh.faces)?;
            let mut o = opts.clone();
            o.max_points = o.max_points.min(n);
            let s = simulate_scan(&bvh, &o, rng.next_u64())?;
            Ok((s.positions, s.normals))
        }
    }
}

/// Write the shapes to `path`. Shapes with more points, vertices or faces than the layout allows
/// are an error naming the shape.
pub fn write_training_npz(path: &Path, shapes: &[ShapeRecord], layout: &Layout) -> Result<(), PatError> {
    let Layout { points: p, k, vertices: v, faces: f } = *layout;
    let s = shapes.len();
    let mut positions = vec![0.0f64; s * p * 3];
    let mut normals = vec![0.0f64; s * p * 3];
    let mut neighbors = vec![-1i32; s * p * k];
    let mut masks = vec![false; s * p];
    let mut mverts = vec![0.0f64; s * v * 3];
    let mut mfaces = vec![-1i32; s * f * 3];
    let mut fmasks = vec![false; s * f];
    for (i, sh) in shapes.iter().enumerate() {
        let n = sh.positions.len();
        if n == 0 || n > p || sh.normals.len() != n || sh.vertices.len() > v || sh.faces.len() > f || sh.faces.is_empty() {
            return Err(PatError::SizeMismatch {
                what: format!("shape {i}: {n} points / {} normals / {} vertices / {} faces for layout {layout:?}", sh.normals.len(), sh.vertices.len(), sh.faces.len()),
            });
        }
        for (j, (a, b)) in sh.positions.iter().zip(&sh.normals).enumerate() {
            positions[(i * p + j) * 3..(i * p + j) * 3 + 3].copy_from_slice(a);
            normals[(i * p + j) * 3..(i * p + j) * 3 + 3].copy_from_slice(b);
            masks[i * p + j] = true;
        }
        let tree = KdTree::new(&sh.positions, KdTreeOptions::default())?;
        let nb = fitting_neighborhoods(&tree, &sh.positions, k);
        for j in 0..n {
            for (c, &x) in nb.row(j).iter().enumerate() {
                neighbors[(i * p + j) * k + c] = x as i32;
            }
        }
        for (j, a) in sh.vertices.iter().enumerate() {
            mverts[(i * v + j) * 3..(i * v + j) * 3 + 3].copy_from_slice(a);
        }
        for (j, t) in sh.faces.iter().enumerate() {
            for c in 0..3 {
                mfaces[(i * f + j) * 3 + c] = t[c] as i32;
            }
            fmasks[i * f + j] = true;
        }
    }
    write_npz(
        path,
        &[
            ("positions", vec![s, p, 3], NpzOut::F64(&positions)),
            ("normals", vec![s, p, 3], NpzOut::F64(&normals)),
            ("neighbors", vec![s, p, k], NpzOut::I32(&neighbors)),
            ("masks", vec![s, p], NpzOut::Bool(&masks)),
            ("mesh_vertices", vec![s, v, 3], NpzOut::F64(&mverts)),
            ("mesh_faces", vec![s, f, 3], NpzOut::I32(&mfaces)),
            ("face_masks", vec![s, f], NpzOut::Bool(&fmasks)),
        ],
    )
}

/// Generate `count` blob shapes (numbers `first..first + count` of the run `seed`) with clouds in
/// `mode`, in parallel on the current rayon pool. Deterministic and independent of the thread
/// count. Fails with the first error by shape index.
pub fn generate_blob_records(seed: u64, first: u64, count: usize, mode: &CloudMode, layout: &Layout, opts: &BlobOptions) -> Result<Vec<ShapeRecord>, PatError> {
    (0..count as u64)
        .into_par_iter()
        .map(|i| {
            let index = first + i;
            let (mesh, _) = generate_blob(seed, index, opts)?;
            let (positions, normals) = sample_cloud(&mesh, mode, layout.points, Rng::for_item(seed ^ 0xB10B, index).next_u64())?;
            Ok(ShapeRecord { positions, normals, vertices: mesh.vertices, faces: mesh.faces })
        })
        .collect()
}
