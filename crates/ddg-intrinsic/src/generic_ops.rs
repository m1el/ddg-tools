//! Operators written once over [`Triangulation`] (Integer Coordinates contract, WP11).
//!
//! The cotan Laplacian and the lumped masses depend on the edge lengths and the connectivity
//! alone, so they are the same for every representation of an intrinsic triangulation. These
//! free functions compute them from the trait; on signposts they are **bitwise identical** to
//! [`IntrinsicTriangulation::cotan_laplacian`] and [`IntrinsicTriangulation::lumped_mass`]
//! (tested), because they perform the same operations in the same order.

use ddg_core::SparseMatrix;
use ddg_core::geometry::corner_cotan;

use crate::interface::Triangulation;

/// The cotan Laplacian of any triangulation: `n_vertices × n_vertices`, symmetric positive
/// semidefinite with zero row sums (see [`IntrinsicTriangulation::cotan_laplacian`] for the
/// conventions, which this function shares).
///
/// [formula: cotan-laplacian]
pub fn cotan_laplacian<T: Triangulation + ?Sized>(t: &T) -> SparseMatrix {
    let complex = t.complex();
    let n = complex.n_vertices();
    let n_faces = complex.n_faces();
    let lengths = t.edge_lengths();
    let mut triplets: Vec<(u32, u32, f64)> = Vec::with_capacity(6 * n_faces + n);
    for v in 0..n as u32 {
        triplets.push((v, v, 0.0));
    }
    for h in 0..complex.n_halfedges() {
        let (i, j) = (complex.from_vertex(h) as u32, complex.to_vertex(h) as u32);
        if i == j {
            continue;
        }
        let l = crate::interface::face_lengths(complex, lengths, h / 3);
        let k = h % 3;
        let cot = corner_cotan(l[k], l[(k + 1) % 3], l[(k + 2) % 3])
            .expect("every face keeps the triangle inequalities");
        let w = 0.5 * cot;
        triplets.push((i, j, -w));
        triplets.push((j, i, -w));
    }
    let mut l = SparseMatrix::from_triplets(n, n, &triplets)
        .expect("vertex indices of a valid triangulation are in range");
    let col_ptr = l.col_ptr().to_vec();
    for j in 0..n {
        let range = col_ptr[j]..col_ptr[j + 1];
        let mut sum = 0.0;
        let mut diagonal_slot = range.start;
        for k in range {
            if l.row_idx()[k] as usize == j {
                diagonal_slot = k;
            } else {
                sum += l.values()[k];
            }
        }
        l.values_mut()[diagonal_slot] = -sum;
    }
    l
}

/// The lumped vertex areas `m_i = ⅓ Σ A_f` of any triangulation (a face that uses a vertex
/// twice counts twice for it); they sum to the total area.
///
/// [formula: lumped-mass]
pub fn lumped_mass<T: Triangulation + ?Sized>(t: &T) -> Vec<f64> {
    let complex = t.complex();
    let mut mass = vec![0.0; complex.n_vertices()];
    for (f, face) in complex.faces.iter().enumerate() {
        let third = t.face_area(f) / 3.0;
        for &v in face {
            mass[v as usize] += third;
        }
    }
    mass
}

/// The lumped masses divided by the total area (they sum to one).
///
/// [formula: lumped-mass]
pub fn lumped_mass_normalized<T: Triangulation + ?Sized>(t: &T) -> Vec<f64> {
    let mut mass = lumped_mass(t);
    let total = t.total_area();
    for m in &mut mass {
        *m /= total;
    }
    mass
}
