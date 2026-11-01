//! Anisotropic coarsening study (paper Sec. 8.6.2, Fig. 27): for a mesh and a principal
//! direction, the largest valid `τ`, and for a list of strengths the elongation of the coarse
//! faces and their shape. Coarse meshes come from `coarsen_by_area` (a stand-in order).
//!
//! ```text
//! cargo run --release -p ddg-intrinsic-simplification --example anisotropy_study -- MESH.obj \
//!     max|min TARGET_VERTICES TAU...
//! ```

use ddg_core::io::read_obj;
use ddg_core::{Geometry, TriMesh};
use ddg_intrinsic_simplification::CoarseMeshOptions;
use ddg_intrinsic_simplification::anisotropy::{
    AnisotropyOptions, PrincipalDirection, anisotropic_coarse_mesh, coarse_face_alignment,
    coarsen_by_area, principal_direction_field,
};
use ddg_intrinsic_simplification::curvature::CurvatureOptions;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 5 {
        eprintln!("usage: anisotropy_study MESH.obj max|min TARGET TAU...");
        std::process::exit(2);
    }
    let which = if args[2] == "min" {
        PrincipalDirection::Min
    } else {
        PrincipalDirection::Max
    };
    let target: usize = args[3].parse().unwrap_or(500);
    let taus: Vec<f64> = args[4..].iter().filter_map(|s| s.parse().ok()).collect();
    let obj = read_obj(&args[1]).unwrap_or_else(|e| {
        eprintln!("{e}");
        std::process::exit(1)
    });
    let mut faces = Vec::new();
    for f in 0..obj.n_faces() {
        let c = &obj.corner_verts[obj.face_corners(f)];
        for k in 1..c.len().saturating_sub(1) {
            faces.push([c[0], c[k], c[k + 1]]);
        }
    }
    let positions = obj.positions;
    let mesh = TriMesh::new(positions.len(), faces.clone()).unwrap_or_else(|e| {
        eprintln!("{e}");
        std::process::exit(1)
    });
    let lengths = Geometry::from_positions(&mesh, &positions)
        .unwrap_or_else(|e| {
            eprintln!("{e}");
            std::process::exit(1)
        })
        .edge_lengths()
        .to_vec();
    let field = principal_direction_field(&positions, &faces, which, &CurvatureOptions::default())
        .unwrap_or_else(|e| {
            eprintln!("{e}");
            std::process::exit(1)
        });
    let other = match which {
        PrincipalDirection::Max => &field.curvatures.d_min,
        PrincipalDirection::Min => &field.curvatures.d_max,
    };
    let a =
        ddg_intrinsic_simplification::anisotropy::edge_alignment(&mesh, &lengths, &field.angles)
            .unwrap();
    let tau_max = ddg_intrinsic_simplification::anisotropy::max_valid_tau(&mesh, &lengths, &a);
    println!(
        "{} {:?}: {} vertices, {} failed curvature fits, tau_max = {tau_max:.4}",
        args[1],
        which,
        positions.len(),
        field.curvatures.n_failed()
    );
    println!(
        "  tau   used  clamped  faces  angle(u)  angle(other)  aspect  3D-minangle  scaled-minangle"
    );
    for &tau in &taus {
        let opts = AnisotropyOptions {
            tau,
            clamp_tau: true,
            ..Default::default()
        };
        let (mut c, rep) = match anisotropic_coarse_mesh(
            &mesh,
            &lengths,
            &field.angles,
            &[],
            CoarseMeshOptions::default(),
            &opts,
        ) {
            Ok(x) => x,
            Err(e) => {
                println!("  {tau:.2}: {e}");
                continue;
            }
        };
        let left = coarsen_by_area(&mut c, target).unwrap();
        let cert = c.certify(false).is_ok();
        let out = c.compact().unwrap();
        let normals = &field.curvatures.normal;
        let u = coarse_face_alignment(&out, &positions, &field.directions, normals, 4);
        let o = coarse_face_alignment(&out, &positions, other, normals, 4);
        // shape of the chord triangles in 3D and of the coarse faces in the scaled metric
        let ff = out.faces_fine();
        let min_angle = |l: [f64; 3]| {
            let ang = |o: f64, a: f64, b: f64| {
                ((a * a + b * b - o * o) / (2.0 * a * b))
                    .clamp(-1.0, 1.0)
                    .acos()
            };
            ang(l[0], l[1], l[2])
                .min(ang(l[1], l[2], l[0]))
                .min(ang(l[2], l[0], l[1]))
        };
        let d3 = |p: u32, q: u32| {
            let (p, q) = (positions[p as usize], positions[q as usize]);
            ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt()
        };
        let m3: f64 = ff
            .iter()
            .map(|f| min_angle([d3(f[1], f[2]), d3(f[2], f[0]), d3(f[0], f[1])]).to_degrees())
            .sum::<f64>()
            / ff.len() as f64;
        let ms: f64 = out
            .corner_angles
            .iter()
            .map(|a| a.iter().cloned().fold(f64::INFINITY, f64::min).to_degrees())
            .sum::<f64>()
            / out.corner_angles.len() as f64;
        println!(
            "  {tau:.2}  {:.3}  {:5}  {:5}  {:7.1}  {:11.1}  {:6.2}  {:10.1}  {:14.1}   ({left} vertices, certify {})",
            rep.tau_used,
            rep.clamped,
            u.n_faces,
            u.mean_angle_deg,
            o.mean_angle_deg,
            u.mean_aspect,
            m3,
            ms,
            if cert { "ok" } else { "FAILED" }
        );
    }
}
