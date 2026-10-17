//! Differential test of `map_to_sphere` against the reference's `BFF::mapToSphere`.
//!
//! Reads the output of the driver `bff_ref sphere` (`reference/bff_ref.cpp`: the position of
//! every corner on the sphere, the positions the reference used and the vertex it removed) and
//! maps the same mesh with this crate, with the same pole. `reference/run_sphere.sh` runs
//! everything.
//!
//! ```text
//! reference_sphere_diff write DIR             # writes the generated test spheres to DIR
//! reference_sphere_diff compare out.txt ...   # one line per driver output
//! ```
//!
//! The two programs center differently (the reference weights faces by their area on the
//! sphere and stops at `|µ| < 1e-3`; this crate weights them by their area on the surface, see
//! the README), so the reference's output is first centered the way this crate centers
//! ([`mobius_center`], tolerance `1e-12`); then both are canonical up to a rotation, and the
//! residual is measured after the best rotation, as chordal distance on the unit sphere.
//!
//! Columns: `V`; the removed vertex; for the reference as it comes, measured by this crate's
//! [`sphere_distortion`]: quasi-conformal error (mean / max), area distortion (area-weighted
//! rms of the log scale factor), flipped faces and the center of mass `|µ|` of the surface area;
//! the same for this crate with the reference's diagonal shift; the residual (rms / max) of this
//! crate with the shift and with exact solves against the re-centered reference.

use ddg_bff::quality::fit_similarity;
use ddg_bff::{
    CenteringOptions, DiskBff, DiskOptions, PoleChoice, REFERENCE_DIAGONAL_SHIFT,
    REFERENCE_MOLLIFY_DELTA, SphereDistortion, SphereOptions, map_to_sphere, mobius_center,
    mollify_edge_lengths, sphere_distortion,
};
use ddg_core::generators::{icosphere, perturb};
use ddg_core::geometry::edge_length;
use ddg_core::io::write_obj;
use ddg_core::{Geometry, TriMesh};

struct Reference {
    positions: Vec<[f64; 3]>,
    faces: Vec<[u32; 3]>,
    corners: Vec<[f64; 3]>,
    /// The unit-disk layout of the punctured sphere, per corner (zero in the pole's star).
    disk: Vec<[f64; 2]>,
    pole: u32,
}

fn parse(path: &str) -> Reference {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let mut r = Reference {
        positions: Vec::new(),
        faces: Vec::new(),
        corners: Vec::new(),
        disk: Vec::new(),
        pole: u32::MAX,
    };
    for line in text.lines() {
        let mut tokens = line.split_ascii_whitespace();
        let tag = tokens.next();
        if !matches!(tag, Some("f" | "p" | "pole" | "d")) {
            continue;
        }
        let x: Vec<f64> = tokens
            .map(|t| {
                t.parse()
                    .unwrap_or_else(|_| panic!("{path}: bad number {t}"))
            })
            .collect();
        match tag {
            Some("f") => {
                r.faces.push([x[0] as u32, x[1] as u32, x[2] as u32]);
                r.corners
                    .extend([[x[3], x[4], x[5]], [x[6], x[7], x[8]], [x[9], x[10], x[11]]]);
            }
            Some("p") => r.positions.push([x[0], x[1], x[2]]),
            Some("d") => r.disk.extend([[x[0], x[1]], [x[2], x[3]], [x[4], x[5]]]),
            Some("pole") => r.pole = x[0] as u32,
            _ => {}
        }
    }
    assert!(
        r.pole != u32::MAX,
        "{path}: no pole line (not a sphere run?)"
    );
    r
}

fn norm(p: [f64; 3]) -> f64 {
    (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt()
}

/// The polar factor of `Σ q pᵀ` (the best rotation or reflection taking `p` to `q`), applied:
/// rms and max of `|R p_i − q_i|`, and the sign of `det R`.
fn orthogonal_fit(p: &[[f64; 3]], q: &[[f64; 3]]) -> (f64, f64, f64) {
    let mut x = [[0.0; 3]; 3];
    for (a, b) in p.iter().zip(q) {
        for i in 0..3 {
            for j in 0..3 {
                x[i][j] += b[i] * a[j];
            }
        }
    }
    let det = |m: &[[f64; 3]; 3]| {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    };
    for _ in 0..100 {
        let d = det(&x);
        let cof = |r0: usize, r1: usize, c0: usize, c1: usize| {
            x[r0][c0] * x[r1][c1] - x[r0][c1] * x[r1][c0]
        };
        let it = [
            [cof(1, 2, 1, 2), -cof(1, 2, 0, 2), cof(1, 2, 0, 1)],
            [-cof(0, 2, 1, 2), cof(0, 2, 0, 2), -cof(0, 2, 0, 1)],
            [cof(0, 1, 1, 2), -cof(0, 1, 0, 2), cof(0, 1, 0, 1)],
        ];
        let next: [[f64; 3]; 3] =
            std::array::from_fn(|i| std::array::from_fn(|j| 0.5 * (x[i][j] + it[i][j] / d)));
        let change = (0..9)
            .map(|k| (next[k / 3][k % 3] - x[k / 3][k % 3]).abs())
            .fold(0.0, f64::max);
        x = next;
        if change < 1e-15 {
            break;
        }
    }
    let mut sum = 0.0;
    let mut worst: f64 = 0.0;
    for (a, b) in p.iter().zip(q) {
        let r: [f64; 3] = std::array::from_fn(|i| (0..3).map(|j| x[i][j] * a[j]).sum());
        let e = norm([r[0] - b[0], r[1] - b[1], r[2] - b[2]]);
        sum += e * e;
        worst = worst.max(e);
    }
    ((sum / p.len() as f64).sqrt(), worst, det(&x).signum())
}

fn center_norm(mesh: &TriMesh, geometry: &Geometry, positions: &[[f64; 3]]) -> f64 {
    let mut copy = positions.to_vec();
    let options = CenteringOptions {
        max_iterations: 0,
        ..CenteringOptions::default()
    };
    mobius_center(mesh, geometry, &mut copy, &options)
        .unwrap()
        .final_center
}

fn summary(d: &SphereDistortion) -> String {
    format!(
        "{:.4} / {:.3} | {:.3} | {}",
        d.faces.quasi_conformal.mean,
        d.faces.quasi_conformal.max,
        d.faces.log_scale_deviation,
        d.faces.flipped_triangles + d.faces.degenerate_triangles
    )
}

fn compare(path: &str) {
    let r = parse(path);
    let name = std::path::Path::new(path)
        .file_name()
        .unwrap()
        .to_string_lossy()
        .trim_end_matches(".out.txt")
        .to_string();
    let mesh = TriMesh::new(r.positions.len(), r.faces.clone()).unwrap();
    // The geometry the reference builds: edge lengths from its positions, mollified.
    let lengths: Vec<f64> = mesh
        .edges()
        .iter()
        .map(|&[a, b]| edge_length(r.positions[a as usize], r.positions[b as usize]))
        .collect();
    let mollified = mollify_edge_lengths(&mesh, &lengths, REFERENCE_MOLLIFY_DELTA).unwrap();
    let geometry = Geometry::from_edge_lengths(&mesh, mollified.edge_lengths).unwrap();

    // One point per vertex; the corners of a vertex must agree.
    let mut reference = vec![[f64::NAN; 3]; mesh.n_vertices()];
    let mut corner_spread: f64 = 0.0;
    for (c, &p) in r.corners.iter().enumerate() {
        let v = r.faces[c / 3][c % 3] as usize;
        if reference[v][0].is_nan() {
            reference[v] = p;
        } else {
            let q = reference[v];
            corner_spread = corner_spread.max(norm([p[0] - q[0], p[1] - q[1], p[2] - q[2]]));
        }
    }
    // The reference's sphere may be mirrored relative to this crate's orientation.
    let mut mirrored = false;
    let mut ref_distortion = sphere_distortion(&mesh, &geometry, &reference).unwrap();
    if ref_distortion.covering_degree < 0.0 {
        mirrored = true;
        for p in &mut reference {
            p[0] = -p[0];
        }
        ref_distortion = sphere_distortion(&mesh, &geometry, &reference).unwrap();
    }
    let ref_mu = center_norm(&mesh, &geometry, &reference);
    // The punctured disk alone: this crate's unit-disk flattening of the mesh without the
    // pole's star, with the reference's shift, against the driver's, after the best similarity
    // (or the best similarity composed with a reflection), relative to the layout's diameter.
    let kept: Vec<u32> = (0..mesh.n_faces() as u32)
        .filter(|&f| !mesh.faces()[f as usize].contains(&r.pole))
        .collect();
    let sub = mesh.submesh(&kept).unwrap();
    let sub_lengths: Vec<f64> = sub
        .mesh
        .edges()
        .iter()
        .map(|&[a, b]| {
            let e = mesh
                .find_edge(
                    sub.vertex_to_parent[a as usize] as usize,
                    sub.vertex_to_parent[b as usize] as usize,
                )
                .unwrap();
            geometry.edge_lengths()[e]
        })
        .collect();
    let sub_geometry = Geometry::from_edge_lengths(&sub.mesh, sub_lengths).unwrap();
    let mut disk = DiskBff::new(
        &sub.mesh,
        &sub_geometry,
        DiskOptions {
            diagonal_shift: REFERENCE_DIAGONAL_SHIFT,
            ..DiskOptions::default()
        },
    )
    .unwrap();
    let ours_disk = disk.flatten_to_disk(10).unwrap();
    // The punctured mesh as a file of its own, for the reference's own disk mode.
    let punctured: Vec<[f64; 3]> = sub
        .vertex_to_parent
        .iter()
        .map(|&v| r.positions[v as usize])
        .collect();
    write_mesh(
        std::path::PathBuf::from(path.replace(".out.txt", "-punctured.obj")),
        &punctured,
        sub.mesh.faces(),
    );
    let (mut p_ours, mut p_ref) = (Vec::new(), Vec::new());
    for (k, &f) in kept.iter().enumerate() {
        for j in 0..3 {
            p_ours.push(ours_disk.uv[sub.mesh.faces()[k][j] as usize]);
            p_ref.push(r.disk[3 * f as usize + j]);
        }
    }
    let (disk_rms, disk_max) = fit_plane(&p_ours, &p_ref);
    let mut centered = reference.clone();
    let tight = CenteringOptions {
        tolerance: 1e-12,
        ..CenteringOptions::default()
    };
    let recentered = mobius_center(&mesh, &geometry, &mut centered, &tight).unwrap();
    assert!(recentered.converged(), "{name}: {recentered:?}");

    let ours = |shift: f64| {
        let options = SphereOptions {
            pole: PoleChoice::Vertex(r.pole),
            diagonal_shift: shift,
            centering: Some(tight),
            ..SphereOptions::default()
        };
        map_to_sphere(&mesh, &geometry, &options).unwrap()
    };
    let shifted = ours(REFERENCE_DIAGONAL_SHIFT);
    let exact = ours(0.0);
    let (rms_s, max_s, det_s) = orthogonal_fit(&shifted.positions, &centered);
    let (rms_e, max_e, _) = orthogonal_fit(&exact.positions, &centered);
    println!(
        "| {name} | {} | {} | {} | {:.0e} | {} | {:.0e} | {:.1e} / {:.1e} | {:.1e} / {:.1e} | {:.1e} / {:.1e} |{}{}",
        mesh.n_vertices(),
        r.pole,
        summary(&ref_distortion),
        ref_mu,
        summary(&shifted.report.distortion),
        shifted.report.centering.unwrap().final_center,
        rms_s,
        max_s,
        rms_e,
        max_e,
        disk_rms,
        disk_max,
        if mirrored { " reference mirrored;" } else { "" },
        if det_s < 0.0 || corner_spread > 1e-12 {
            format!(" det {det_s}, corner spread {corner_spread:.1e}")
        } else {
            String::new()
        },
    );
}

/// The best similarity, or similarity composed with a reflection, from `p` to `q`: rms and max
/// residual relative to the diameter of `q`.
fn fit_plane(p: &[[f64; 2]], q: &[[f64; 2]]) -> (f64, f64) {
    let (lo, hi) = q.iter().fold(
        ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]),
        |(lo, hi), x| {
            (
                [lo[0].min(x[0]), lo[1].min(x[1])],
                [hi[0].max(x[0]), hi[1].max(x[1])],
            )
        },
    );
    let diameter = (hi[0] - lo[0]).hypot(hi[1] - lo[1]);
    let direct = fit_similarity(p, q).unwrap();
    let mirrored: Vec<[f64; 2]> = p.iter().map(|x| [-x[0], x[1]]).collect();
    let flipped = fit_similarity(&mirrored, q).unwrap();
    let best = if flipped.rms < direct.rms {
        flipped
    } else {
        direct
    };
    (best.rms / diameter, best.max / diameter)
}

/// The driver's `disk` output for a punctured sphere against this crate's `flatten_to_disk`
/// with the reference's shift.
fn compare_disk(path: &str) {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let (mut positions, mut faces, mut uv) = (Vec::new(), Vec::new(), Vec::new());
    for line in text.lines() {
        let mut tokens = line.split_ascii_whitespace();
        let tag = tokens.next();
        if !matches!(tag, Some("f" | "p")) {
            continue;
        }
        let x: Vec<f64> = tokens.map(|t| t.parse().unwrap()).collect();
        if tag == Some("f") {
            faces.push([x[0] as u32, x[1] as u32, x[2] as u32]);
            uv.extend([[x[3], x[4]], [x[5], x[6]], [x[7], x[8]]]);
        } else {
            positions.push([x[0], x[1], x[2]]);
        }
    }
    let mesh = TriMesh::new(positions.len(), faces.clone()).unwrap();
    let lengths: Vec<f64> = mesh
        .edges()
        .iter()
        .map(|&[a, b]| edge_length(positions[a as usize], positions[b as usize]))
        .collect();
    let mollified = mollify_edge_lengths(&mesh, &lengths, REFERENCE_MOLLIFY_DELTA).unwrap();
    let geometry = Geometry::from_edge_lengths(&mesh, mollified.edge_lengths).unwrap();
    let mut disk = DiskBff::new(
        &mesh,
        &geometry,
        DiskOptions {
            diagonal_shift: REFERENCE_DIAGONAL_SHIFT,
            ..DiskOptions::default()
        },
    )
    .unwrap();
    let ours = disk.flatten_to_disk(10).unwrap();
    let mine: Vec<[f64; 2]> = faces
        .iter()
        .flat_map(|f| f.map(|v| ours.uv[v as usize]))
        .collect();
    let (rms, max) = fit_plane(&mine, &uv);
    let name = std::path::Path::new(path)
        .file_name()
        .unwrap()
        .to_string_lossy()
        .trim_end_matches("-punctured.out.txt")
        .to_string();
    println!("| {name} | {} | {rms:.1e} / {max:.1e} |", mesh.n_vertices());
}

fn write_mesh(path: std::path::PathBuf, positions: &[[f64; 3]], faces: &[[u32; 3]]) {
    let corners: Vec<u32> = faces.iter().flatten().copied().collect();
    let starts: Vec<u32> = (0..faces.len() as u32).map(|f| 3 * f).collect();
    write_obj(&path, positions, &corners, &starts, None)
        .unwrap_or_else(|e| panic!("{}: {e}", path.display()));
}

fn write(dir: &str) {
    let dir = std::path::Path::new(dir);
    let g = icosphere(3).unwrap();
    write_mesh(dir.join("generated-icosphere3.obj"), &g.positions, &g.faces);
    let mut g = icosphere(4).unwrap();
    for p in &mut g.positions {
        p[0] *= 2.0;
        p[2] *= 0.5;
    }
    perturb(&mut g.positions, [0.01; 3], 7);
    write_mesh(dir.join("generated-ellipsoid.obj"), &g.positions, &g.faces);
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        Some("write") if args.len() == 2 => write(&args[1]),
        Some("compare") => {
            println!(
                "| mesh | V | pole | reference: QC mean / max | area rms u | flipped | reference |µ| | this crate (shift 1e-8): QC mean / max | area rms u | flipped | |µ| | residual, shift 1e-8 (rms / max) | residual, exact (rms / max) | punctured disk, shift 1e-8 (rms / max) |"
            );
            println!("|---|---:|---:|---|---:|---:|---:|---|---:|---:|---:|---|---|---|");
            for path in &args[1..] {
                compare(path);
            }
        }
        Some("disk") => {
            println!(
                "| punctured mesh | V | this crate (shift 1e-8) against the driver's disk mode (rms / max) |"
            );
            println!("|---|---:|---|");
            for path in &args[1..] {
                compare_disk(path);
            }
        }
        _ => {
            eprintln!(
                "usage: reference_sphere_diff write DIR | compare out.txt ... | disk out.txt ..."
            );
            std::process::exit(2);
        }
    }
}
