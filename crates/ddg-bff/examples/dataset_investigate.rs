//! What do the UV maps stored in the optimal-cones dataset contain, beyond what the cone
//! files say? Diagnostics behind the notes in the README.
//!
//! ```text
//! cargo run --release -p ddg-bff --example dataset_investigate [mesh ...]
//! ```
//!
//! Per mesh:
//!
//! * `seam p50/max` — relative difference of the lengths of the two sides of the cut edges in
//!   the stored map (median and maximum): the size of the six-digit rounding of the file.
//! * `cone diff` — the largest difference, modulo 2π, between a published cone angle and the
//!   angle the stored map realizes at that cone, measured on the boundary polygon of the cut
//!   map (flipped triangles do not affect it), and the vertex where it occurs.
//! * `unlisted` — the sum of the cone angles that the stored map realizes at the vertices of
//!   its cut that are not listed as cones. On a closed mesh this must make up for the
//!   Gauss–Bonnet discrepancy of the published angles (`GB discr`).
//! * `rms` — distance of our layout (same cut, published angles) from the stored one, and
//!   `rms*` the same with the cone angles of the stored map wherever they differ from the
//!   published ones by more than 1e-4.

#[path = "../tests/common/mod.rs"]
mod common;

use common::{compare_layouts, dataset_names, load, stored_cone_angles};
use ddg_bff::{
    Cone, ConeFlattenInput, ConeFlattenOptions, CutSource, GaussBonnetPolicy, ScaleFactorOptions,
    cut_mesh, flatten_with_cones,
};

fn main() {
    let Some(mut names) = dataset_names() else {
        eprintln!("dataset not found at {}", common::dataset_dir().display());
        std::process::exit(1);
    };
    let wanted: Vec<String> = std::env::args().skip(1).collect();
    if !wanted.is_empty() {
        names.retain(|n| wanted.contains(n));
    }
    println!(
        "{:<20} {:>8} {:>8} | {:>10} {:>7} | {:>9} {:>9} | {:>8} {:>8}",
        "mesh",
        "seam p50",
        "seam max",
        "cone diff",
        "vertex",
        "unlisted",
        "GB discr",
        "rms",
        "rms*"
    );
    let tau = std::f64::consts::TAU;
    for name in &names {
        let data = load(name);
        let (mesh, uv) = (&data.mesh, &data.reference_uvs);
        let length = |h: usize| {
            let (a, b) = (uv[h], uv[mesh.next(h)]);
            (b[0] - a[0]).hypot(b[1] - a[1])
        };
        let mut seams: Vec<f64> = data
            .reference_cut
            .iter()
            .map(|&e| {
                let h = mesh.edge_halfedge(e as usize);
                let t = mesh.twin(h).expect("cut edges are interior");
                (length(h) - length(t)).abs() / length(h).max(length(t))
            })
            .collect();
        seams.sort_by(|a, b| a.total_cmp(b));

        // Cone angles of the stored map at the listed cones.
        let stored = stored_cone_angles(&data);
        let (mut worst, mut worst_vertex) = (0.0f64, 0u32);
        let mut cones: Vec<Cone> = data.cones.clone();
        for (cone, &angle) in cones.iter_mut().zip(&stored) {
            let difference = angle - cone.angle;
            if difference.abs() > worst.abs() {
                (worst, worst_vertex) = (difference, cone.vertex);
            }
            if difference.abs() > 1e-4 {
                cone.angle = angle;
            }
        }

        // Cone angles of the stored map at the other vertices of its cut.
        let cut = cut_mesh(mesh, &data.reference_cut).expect("stored cut");
        let mut vertex_uv = vec![[0.0; 2]; cut.mesh.n_vertices()];
        for (f, face) in cut.mesh.faces().iter().enumerate() {
            for k in 0..3 {
                vertex_uv[face[k] as usize] = uv[3 * f + k];
            }
        }
        let mut interior = vec![0.0; mesh.n_vertices()];
        let mut on_cut = vec![false; mesh.n_vertices()];
        for boundary in cut.mesh.boundary_loops() {
            let n = boundary.len();
            for i in 0..n {
                let (a, b, c) = (
                    vertex_uv[boundary[(i + n - 1) % n] as usize],
                    vertex_uv[boundary[i] as usize],
                    vertex_uv[boundary[(i + 1) % n] as usize],
                );
                let (ux, uy, vx, vy) = (b[0] - a[0], b[1] - a[1], c[0] - b[0], c[1] - b[1]);
                let original = cut.vertex_to_original[boundary[i] as usize] as usize;
                interior[original] += tau / 2.0 - (ux * vy - uy * vx).atan2(ux * vx + uy * vy);
                on_cut[original] = true;
            }
        }
        let mut listed = vec![false; mesh.n_vertices()];
        for cone in &data.cones {
            listed[cone.vertex as usize] = true;
        }
        let unlisted: f64 = (0..mesh.n_vertices())
            .filter(|&v| on_cut[v] && !listed[v] && !mesh.is_boundary_vertex(v))
            .map(|v| {
                let angle = (tau - interior[v]).rem_euclid(tau);
                if angle > tau / 2.0 {
                    angle - tau
                } else {
                    angle
                }
            })
            .sum();
        let discrepancy = if mesh.boundary_loops().is_empty() {
            let sum: f64 = data.cones.iter().map(|c| c.angle).sum();
            format!("{:9.1e}", sum - tau * mesh.euler_characteristic() as f64)
        } else {
            "        -".to_string()
        };

        // Our layout with the published and with the stored cone angles.
        let (completed, _) = data.completed_reference_cut();
        let residual = |cones: &[Cone]| {
            flatten_with_cones(
                &ConeFlattenInput {
                    mesh,
                    geometry: &data.geometry,
                    cones,
                    scale_factors: None,
                    cut: CutSource::Edges(&completed),
                    anchors: &[],
                },
                &ConeFlattenOptions {
                    scale_factors: ScaleFactorOptions {
                        gauss_bonnet: GaussBonnetPolicy::SpreadOverVertices,
                        ..ScaleFactorOptions::default()
                    },
                    ..ConeFlattenOptions::default()
                },
            )
            .map(|flat| {
                compare_layouts(&flat.face_chart, flat.n_charts(), &flat.corner_uvs, uv).rms
            })
        };
        let published = residual(&data.cones).unwrap_or(f64::NAN);
        let corrected = if cones == data.cones {
            "       =".to_string()
        } else {
            format!("{:8.1e}", residual(&cones).unwrap_or(f64::NAN))
        };
        println!(
            "{:<20} {:8.1e} {:8.1e} | {:10.2e} {:>7} | {:9.1e} {} | {:8.1e} {}",
            name,
            seams[seams.len() / 2],
            seams[seams.len() - 1],
            worst,
            worst_vertex,
            unlisted,
            discrepancy,
            published,
            corrected
        );
    }
}
