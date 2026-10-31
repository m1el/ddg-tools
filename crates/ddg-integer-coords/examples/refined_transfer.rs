//! Transfer after refining by insertions (not the paper's Fig. 14, which is `fig14.rs`): a function
//! brought to the input mesh by (a) copying the values at the shared vertices and (b) the `L²`
//! transfer of Eq. 13-14.
//!
//! The paper's experiment: an intrinsic Delaunay triangulation is refined (Sec. 4.2), a smooth
//! function is sampled at its vertices, and the function on the input mesh is compared with the
//! exact one; at the marked size the transfer was 28x more accurate than the original mesh and
//! the refined triangulation 11.5x (SPEC 10.1). We have no refinement yet (WP6), so the
//! triangulation is refined here by inserting vertices (the circumcenter, or the centroid if the
//! circumcenter is outside, of the face of the largest area) and Delaunay flipping every 25
//! insertions, to a given multiple of `|V⁰|` extra vertices. The exact function is
//! `u(x) = sin(6 pi x / D) + 0.5 cos(4 pi y / D + 0.4) + 0.25 sin(10 pi z / D)` with `D` the
//! diagonal of the bounding box, sampled at the vertices of `T¹` (the positions of inserted
//! vertices are the exact positions on the input surface). Errors are `L²` norms over the input
//! surface of `u − (piecewise linear function)`, by quadrature on the pieces of the common
//! subdivision, relative to `‖u‖`.
//!
//! Columns: `copy` is `u` sampled at the input vertices (the input mesh's own interpolant, what
//! copying the shared vertices gives), `transfer` the `L²` transfer of the samples at `T¹`'s
//! vertices, `T1` the error of the refined triangulation's own piecewise linear function.
//!
//! ```text
//! cargo run --release -p ddg-integer-coords --example refined_transfer
//! ```
//!
//! Writes PLY files with vertex colours to `runs_tmp/wp7-refined/`.

use std::path::PathBuf;
use std::time::Instant;

use ddg_core::io::read_obj;
use ddg_core::{Geometry, TriMesh};
use ddg_integer_coords::{IntegerTriangulation, Transfer};
use ddg_intrinsic::{DelaunayOptions, Triangulation};

struct Mesh {
    name: String,
    mesh: TriMesh,
    positions: Vec<[f64; 3]>,
    geometry: Geometry,
}

fn load(relative: &str) -> Option<Mesh> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../data")
        .join(relative);
    let obj = read_obj(&path).ok()?;
    let mut faces = Vec::new();
    for f in 0..obj.n_faces() {
        let c = &obj.corner_verts[obj.face_corners(f)];
        for k in 1..c.len() - 1 {
            faces.push([c[0], c[k], c[k + 1]]);
        }
    }
    let mesh = TriMesh::new(obj.positions.len(), faces).ok()?;
    let geometry = Geometry::from_positions(&mesh, &obj.positions).ok()?;
    Some(Mesh {
        name: path.file_stem()?.to_str()?.to_string(),
        mesh,
        positions: obj.positions,
        geometry,
    })
}

fn heron(l: [f64; 3]) -> f64 {
    let s = 0.5 * (l[0] + l[1] + l[2]);
    (s * (s - l[0]) * (s - l[1]) * (s - l[2])).max(0.0).sqrt()
}

/// Barycentric coordinates of the circumcenter of a face with these lengths (halfedge `k` runs
/// from corner `k` to `k + 1`, so `l[k]` is opposite corner `k + 2`), or the centroid.
fn center(l: [f64; 3]) -> [f64; 3] {
    // Squared length opposite corner c: l[(c + 1) % 3].
    let opp = |c: usize| l[(c + 1) % 3] * l[(c + 1) % 3];
    let w: [f64; 3] =
        std::array::from_fn(|c| opp(c) * (opp((c + 1) % 3) + opp((c + 2) % 3) - opp(c)));
    let s: f64 = w.iter().sum();
    let b = [w[0] / s, w[1] / s, w[2] / s];
    if s > 0.0 && b.iter().all(|&x| x > 0.05) {
        b
    } else {
        [1.0 / 3.0; 3]
    }
}

fn refine(tri: &mut IntegerTriangulation, extra: usize) -> usize {
    let mut done = 0;
    let mut failed = std::collections::HashSet::new();
    let mut since_flip = 0;
    let mut stuck = 0;
    while done < extra && stuck < 1000 {
        let nf = tri.complex().n_faces();
        let best = (0..nf)
            .filter(|f| !failed.contains(f))
            .max_by(|&a, &b| heron(tri.face_lengths(a)).total_cmp(&heron(tri.face_lengths(b))));
        let Some(f) = best else {
            failed.clear();
            stuck += 1;
            continue;
        };
        match tri.insert_vertex(f, center(tri.face_lengths(f))) {
            Ok(_) => {
                done += 1;
                since_flip += 1;
            }
            Err(_) => {
                failed.insert(f);
                stuck += 1;
            }
        }
        if since_flip == 25 || done == extra {
            tri.flip_to_delaunay(&DelaunayOptions::default())
                .expect("delaunay");
            failed.clear();
            since_flip = 0;
        }
    }
    done
}

fn ramp(t: f64) -> [u8; 3] {
    let t = t.clamp(0.0, 1.0);
    [
        (255.0 * t) as u8,
        (255.0 * (1.0 - (2.0 * t - 1.0).abs())) as u8,
        (255.0 * (1.0 - t)) as u8,
    ]
}

fn write_ply(path: &std::path::Path, m: &Mesh, values: &[f64], lo: f64, hi: f64) {
    use std::fmt::Write;
    let mut s = String::new();
    writeln!(s, "ply\nformat ascii 1.0\nelement vertex {}\nproperty float x\nproperty float y\nproperty float z\nproperty uchar red\nproperty uchar green\nproperty uchar blue\nelement face {}\nproperty list uchar int vertex_indices\nend_header", m.positions.len(), m.mesh.n_faces()).unwrap();
    for (p, v) in m.positions.iter().zip(values) {
        let c = ramp((v - lo) / (hi - lo));
        writeln!(s, "{} {} {} {} {} {}", p[0], p[1], p[2], c[0], c[1], c[2]).unwrap();
    }
    for f in m.mesh.faces() {
        writeln!(s, "3 {} {} {}", f[0], f[1], f[2]).unwrap();
    }
    let _ = std::fs::write(path, s);
}

fn run(m: &Mesh, out_dir: &std::path::Path) -> String {
    use std::fmt::Write;
    let mut out = String::new();
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for p in &m.positions {
        for c in 0..3 {
            lo[c] = lo[c].min(p[c]);
            hi[c] = hi[c].max(p[c]);
        }
    }
    let d = (0..3).map(|c| (hi[c] - lo[c]).powi(2)).sum::<f64>().sqrt();
    let pi = std::f64::consts::PI;
    let u = move |p: [f64; 3]| {
        (6.0 * pi * p[0] / d).sin()
            + 0.5 * (4.0 * pi * p[1] / d + 0.4).cos()
            + 0.25 * (10.0 * pi * p[2] / d).sin()
    };
    writeln!(
        out,
        "{}: |V0| = {}, |F0| = {}",
        m.name,
        m.positions.len(),
        m.mesh.n_faces()
    )
    .unwrap();
    writeln!(
        out,
        "{:>7} {:>8} {:>8} | {:>10} {:>10} {:>10} | {:>8} {:>8} | {:>8} {:>8}",
        "extra", "|V1|", "|S|", "copy", "transfer", "T1", "copy/tr", "T1/tr", "t_build", "t_solve"
    )
    .unwrap();
    let mut tri = IntegerTriangulation::new(&m.mesh, &m.geometry).expect("init");
    tri.flip_to_delaunay(&DelaunayOptions::default())
        .expect("delaunay");
    let mut inserted = 0;
    let mut timing;
    for frac in [0.0, 0.1, 0.25, 0.5, 1.0] {
        let target = (frac * m.positions.len() as f64) as usize;
        inserted += refine(&mut tri, target - inserted);
        let t0 = Instant::now();
        let sub = tri.common_subdivision().expect("subdivision");
        let t_sub = t0.elapsed().as_secs_f64();
        let tr = Transfer::new(&tri, &sub).expect("transfer");
        let t_build = t0.elapsed().as_secs_f64();
        let r = tr.report();
        timing = format!(
            "  timings at the last level (s): common subdivision {t_sub:.4}, M_S {:.4}, P0 P1 {:.4}, M0 M1 {:.4}, factorize both {:.4}, residual {:.1e}",
            r.t_mass_s.as_secs_f64(),
            r.t_interpolation.as_secs_f64(),
            r.t_congruence.as_secs_f64(),
            r.t_factor.as_secs_f64(),
            tr.to_input_report(&vec![1.0; tr.n_intrinsic_vertices()])
                .unwrap()
                .1
                .residual
        );
        let pos = tr.subdivision_positions(&m.positions).unwrap();
        let n1 = tr.n_intrinsic_vertices();
        let f: Vec<f64> = pos[..n1].iter().map(|&p| u(p)).collect();
        let t0 = Instant::now();
        let fhat = tr.to_input(&f).unwrap();
        let t_solve = t0.elapsed().as_secs_f64();
        let fcopy = tr.copy_to_input(&f).unwrap();
        let err = |g: &dyn Fn(&mut [f64])| {
            let mut on_s = vec![0.0; pos.len()];
            g(&mut on_s);

            tr.l2_error_to_function(&pos, &on_s, u).unwrap()
        };
        let e_copy = err(&|y| tr.p0().mul_vec(&fcopy, y).unwrap());
        let e_tr = err(&|y| tr.p0().mul_vec(&fhat, y).unwrap());
        let e_t1 = err(&|y| tr.p1().mul_vec(&f, y).unwrap());
        let norm = tr
            .l2_error_to_function(&pos, &vec![0.0; pos.len()], u)
            .unwrap();
        writeln!(
            out,
            "{:>7} {:>8} {:>8} | {:>10.3e} {:>10.3e} {:>10.3e} | {:>8.2} {:>8.2} | {:>7.3}s {:>7.4}s",
            inserted,
            n1,
            pos.len(),
            e_copy / norm,
            e_tr / norm,
            e_t1 / norm,
            e_copy / e_tr,
            e_t1 / e_tr,
            t_build,
            t_solve
        )
        .unwrap();
        if frac == 1.0 {
            writeln!(out, "{timing}").unwrap();
            let exact: Vec<f64> = m.positions.iter().map(|&p| u(p)).collect();
            let (a, b) = exact
                .iter()
                .fold((f64::MAX, f64::MIN), |(a, b), &v| (a.min(v), b.max(v)));
            let _ = std::fs::create_dir_all(out_dir);
            write_ply(
                &out_dir.join(format!("{}-exact.ply", m.name)),
                m,
                &exact,
                a,
                b,
            );
            write_ply(
                &out_dir.join(format!("{}-copy.ply", m.name)),
                m,
                &fcopy,
                a,
                b,
            );
            write_ply(
                &out_dir.join(format!("{}-transfer.ply", m.name)),
                m,
                &fhat,
                a,
                b,
            );
            let diff: Vec<f64> = fhat
                .iter()
                .zip(&exact)
                .map(|(x, y)| (x - y).abs())
                .collect();
            let dcopy: Vec<f64> = fcopy
                .iter()
                .zip(&exact)
                .map(|(x, y)| (x - y).abs())
                .collect();
            let top = diff.iter().chain(&dcopy).cloned().fold(0.0, f64::max);
            write_ply(
                &out_dir.join(format!("{}-error-transfer.ply", m.name)),
                m,
                &diff,
                0.0,
                top,
            );
            write_ply(
                &out_dir.join(format!("{}-error-copy.ply", m.name)),
                m,
                &dcopy,
                0.0,
                top,
            );
        }
    }
    out
}

fn main() {
    let out_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../runs_tmp/wp7-refined");
    let meshes: Vec<Mesh> = [
        "models/spot/spot/spot_triangulated.obj",
        "optimal-cones/optimal-cones-data/data/camelbad.obj",
    ]
    .iter()
    .filter_map(|p| load(p))
    .collect();
    if meshes.is_empty() {
        eprintln!("no meshes in data/; nothing to do");
        return;
    }
    // One after the other: the timings are single-thread numbers.
    for m in &meshes {
        println!("{}", run(m, &out_dir));
    }
}
