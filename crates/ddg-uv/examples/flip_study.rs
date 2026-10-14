//! Where the flipped triangles of the end-to-end unwrap of a cut come from.
//!
//! ```text
//! cargo run --release -p ddg-uv --example flip_study -- CASE TRIANGULATION [MIN_OFFSET ...]
//! CASE:          sphere | cube | sheared | needles | camelbad | spot | atlas:W (spot, 8 clusters,
//!                ripple wavelength W)
//! TRIANGULATION: caller | flips | refined
//! ```
//!
//! Runs the chain (start, flow, `Session::unwrap`) and prints, for every flipped triangle of the
//! unwrap, its 3D area against the mean of the inserted mesh, its smallest 3D angle, its
//! shortest edge against its longest, and for each of its vertices where it comes from and how
//! close to a vertex or a side of the old triangle it is (`min(t, 1 - t)` for a point of an edge,
//! the smallest barycentric coordinate for a point inside a triangle). Then it repeats the
//! unwrap of the same flow state with other values of `min_offset` and prints the counts.
//! `DDG_MODELS` and `DDG_DATASET` locate spot and camelbad.
#![allow(clippy::needless_range_loop)]
use std::path::PathBuf;

use ddg_core::io::read_obj;
use ddg_uv::PolygonMesh;
use ddg_uv::cutting::{
    ConePlacement, CutUnwrapOptions, DEFAULT_CROSSING_OFFSET, FlowKind, IntrinsicOptions, Session,
    Start, Triangulation, VertexSource,
};
use ddg_uv::primitives::{grid_cube, uv_sphere};
use ddg_variational_cuts::{
    Energy, PAPER_FEATURE_SMOOTHING, ScaleWeight, StepSize, StopRule, length_weight_for_ripples,
};

fn sheared_grid(n: usize, sigma: f64) -> PolygonMesh {
    let h = 1.0 / n as f64;
    let mut positions = Vec::new();
    for j in 0..=n {
        for i in 0..=n {
            positions.push([(i as f64 + sigma * j as f64) * h, j as f64 * h, 0.0]);
        }
    }
    let at = |i: usize, j: usize| (j * (n + 1) + i) as u32;
    let mut faces = Vec::new();
    for j in 0..n {
        for i in 0..n {
            let (a, b, c, d) = (at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1));
            faces.push([a, b, c]);
            faces.push([a, c, d]);
        }
    }
    PolygonMesh::from_triangles(positions, &faces)
}

fn needles_and_caps(n: usize) -> PolygonMesh {
    let d = 1.0 / n as f64;
    let (ha, hb) = (0.05 * d, 3.0 * d);
    let pairs = ((1.0 / (ha + hb)).round() as usize).max(1);
    let scale = 1.0 / (pairs as f64 * (ha + hb));
    let rows = 2 * pairs + 1;
    let mut y = vec![0.0f64];
    for j in 0..2 * pairs {
        let last = y[j];
        y.push(last + scale * if j % 2 == 0 { ha } else { hb });
    }
    *y.last_mut().unwrap() = 1.0;
    let mut positions = Vec::new();
    let mut row_x: Vec<Vec<f64>> = Vec::new();
    let mut start = Vec::new();
    for j in 0..rows {
        let xs: Vec<f64> = if j % 2 == 0 {
            (0..=n).map(|i| i as f64 * d).collect()
        } else {
            std::iter::once(0.0)
                .chain((0..n).map(|i| (i as f64 + 0.5) * d))
                .chain(std::iter::once(1.0))
                .collect()
        };
        start.push(positions.len() as u32);
        for &x in &xs {
            positions.push([x, y[j], 0.0]);
        }
        row_x.push(xs);
    }
    let mut faces = Vec::new();
    for j in 0..rows - 1 {
        let (xa, xb) = (&row_x[j], &row_x[j + 1]);
        let (mut i, mut m) = (0usize, 0usize);
        while i < xa.len() - 1 || m < xb.len() - 1 {
            let bottom = if i == xa.len() - 1 {
                false
            } else if m == xb.len() - 1 {
                true
            } else {
                xa[i + 1] <= xb[m + 1]
            };
            let (a, b) = (start[j] + i as u32, start[j + 1] + m as u32);
            if bottom {
                faces.push([a, a + 1, b]);
                i += 1;
            } else {
                faces.push([a, b + 1, b]);
                m += 1;
            }
        }
    }
    PolygonMesh::from_triangles(positions, &faces)
}

fn load(path: PathBuf, polygons: bool) -> PolygonMesh {
    let obj = read_obj(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    if polygons {
        let polys: Vec<Vec<u32>> = (0..obj.n_faces())
            .map(|f| obj.corner_verts[obj.face_corners(f)].to_vec())
            .collect();
        PolygonMesh::from_polygons(obj.positions.clone(), &polys)
    } else {
        PolygonMesh::from_triangles(obj.positions.clone(), &obj.triangles().unwrap())
    }
}

fn area(mesh: &PolygonMesh) -> f64 {
    mesh.triangles.iter().map(|t| tri_area(mesh, *t)).sum()
}

fn p3(mesh: &PolygonMesh, c: u32) -> [f64; 3] {
    mesh.positions[mesh.corner_vertices[c as usize] as usize]
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn norm(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

fn tri_area(mesh: &PolygonMesh, t: [u32; 3]) -> f64 {
    let (a, b, c) = (p3(mesh, t[0]), p3(mesh, t[1]), p3(mesh, t[2]));
    let (u, v) = (sub(b, a), sub(c, a));
    let x = [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ];
    0.5 * norm(x)
}

struct Outcome {
    flipped: usize,
    n_clamped: usize,
    merged: usize,
    displacement: f64,
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let (case, tri_name) = (args[0].clone(), args[1].clone());
    let offsets: Vec<f64> = if args.len() > 2 {
        args[2..].iter().map(|s| s.parse().unwrap()).collect()
    } else {
        vec![1e-6, 1e-4, 1e-3, 1e-2]
    };
    let models = std::env::var_os("DDG_MODELS").map_or_else(
        || PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/models/spot/spot"),
        PathBuf::from,
    );
    let dataset = std::env::var_os("DDG_DATASET").map_or_else(
        || {
            PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("../../data/optimal-cones/optimal-cones-data/data")
        },
        PathBuf::from,
    );
    let (name, wavelength) = match case.split_once(':') {
        Some((n, w)) => (n.to_string(), w.parse::<f64>().unwrap()),
        None => (case.clone(), 0.5),
    };
    let (mesh, clusters, steps) = match name.as_str() {
        "sphere" => (uv_sphere(32, 16, 1.0).unwrap(), None, 40),
        "cube" => (grid_cube(6, 2.0), None, 40),
        "sheared" => (sheared_grid(32, 2.0), None, 8),
        "sheared165" => (sheared_grid(32, 1.0 / 15f64.to_radians().tan()), None, 8),
        "needles" => (needles_and_caps(16), None, 8),
        "camelbad" => (load(dataset.join("camelbad.obj"), false), None, 20),
        "spot" => (load(models.join("spot_quadrangulated.obj"), true), None, 40),
        "atlas" => (
            load(models.join("spot_quadrangulated.obj"), true),
            Some(8),
            60,
        ),
        other => panic!("unknown case {other}"),
    };
    let merge: f64 = std::env::var("MERGE").map_or(0.0, |v| v.parse().unwrap());
    let make_tri = |o: f64| match tri_name.as_str() {
        "caller" => Triangulation::Caller,
        "flips" => Triangulation::Intrinsic(IntrinsicOptions {
            min_offset: o,
            merge_tolerance: merge,
            ..IntrinsicOptions::flips_only()
        }),
        "refined" => Triangulation::Intrinsic(IntrinsicOptions {
            min_offset: o,
            merge_tolerance: merge,
            ..IntrinsicOptions::refined()
        }),
        other => panic!("unknown triangulation {other}"),
    };
    let flat = matches!(name.as_str(), "sheared" | "sheared165" | "needles");
    let center = {
        let mut c = [0.0f64; 3];
        for p in &mesh.positions {
            for d in 0..3 {
                c[d] += p[d] / mesh.positions.len() as f64;
            }
        }
        (0..mesh.positions.len())
            .min_by(|&a, &b| {
                let dist = |v: usize| {
                    (0..3)
                        .map(|d| (mesh.positions[v][d] - c[d]).powi(2))
                        .sum::<f64>()
                };
                dist(a).total_cmp(&dist(b))
            })
            .unwrap() as u32
    };
    let run = |offset: f64, verbose: bool| -> Outcome {
        let tri = make_tri(offset);
        let (n, start) = match clusters {
            Some(n) => (n, Start::NormalClusters { max_iterations: 50 }),
            None => (
                2,
                Start::GeodesicCircle {
                    vertex: center,
                    radius: 0.2 * area(&mesh).sqrt(),
                },
            ),
        };
        let mut session =
            Session::with_triangulation(&mesh.input(), n, &start, FlowKind::Auto, &tri).unwrap();
        let options = session.recommended_options(StepSize::default());
        let energy = if clusters.is_some() {
            Energy::ScaledDirichlet(ScaleWeight::FeatureSize {
                r_eps: 1.0,
                smoothing: PAPER_FEATURE_SMOOTHING,
            })
        } else {
            Energy::Hencky
        };
        let mut problem = session.new_problem(energy, 0.0);
        let omega = session.roughening_rate(&problem, &options.eval).unwrap();
        problem.length_weight = if flat {
            1.0
        } else {
            length_weight_for_ripples(omega, if clusters.is_some() { wavelength } else { 0.25 })
        };
        let report = session
            .run(
                &problem,
                &options,
                &StopRule {
                    max_steps: steps,
                    motion_tolerance: if std::env::var_os("NOSTOP").is_some() {
                        0.0
                    } else {
                        1e-5
                    },
                    ..StopRule::default()
                },
            )
            .unwrap();
        if std::env::var_os("REPORT").is_some() {
            println!(
                "RUN {} steps, stop {:?}, final energy {:.6}",
                report.steps.len(),
                report.stop,
                report.final_energy.total
            );
            for (k, r) in report.steps.iter().take(4).enumerate() {
                println!(
                    "STEP {k}: tau {:.3e}, max_speed {:.3e}, max_length_speed {:.3e}, omega {:.3e}, cut_motion {:.3e}, energy {:.6}, cut points {}",
                    r.tau,
                    r.max_speed,
                    r.max_length_speed,
                    r.roughening_rate,
                    r.cut_motion,
                    r.energy.total,
                    r.n_cut_points
                );
            }
        }
        let (cut, result) = session
            .unwrap(
                DEFAULT_CROSSING_OFFSET,
                &CutUnwrapOptions {
                    cones: ConePlacement::Optimal,
                    ..Default::default()
                },
            )
            .unwrap();
        let m = &result.inserted.mesh;
        let mean = area(m) / m.triangles.len() as f64;
        let mut flipped = 0;
        for (t, tri) in m.triangles.iter().enumerate() {
            let uv = tri.map(|c| result.map.uvs[c as usize]);
            let a = (uv[1][0] - uv[0][0]) * (uv[2][1] - uv[0][1])
                - (uv[1][1] - uv[0][1]) * (uv[2][0] - uv[0][0]);
            if a >= 0.0 {
                continue;
            }
            flipped += 1;
            if !verbose {
                continue;
            }
            let p = tri.map(|c| p3(m, c));
            let e = [
                norm(sub(p[1], p[0])),
                norm(sub(p[2], p[1])),
                norm(sub(p[0], p[2])),
            ];
            let ang = |i: usize| {
                let (a, b, c) = (e[i], e[(i + 1) % 3], e[(i + 2) % 3]);
                ((b * b + c * c - a * a) / (2.0 * b * c))
                    .clamp(-1.0, 1.0)
                    .acos()
            };
            let min_angle = (0..3).map(ang).fold(f64::INFINITY, f64::min).to_degrees();
            let verts = tri.map(|c| m.corner_vertices[c as usize]);
            let mut desc = Vec::new();
            let mut min_offset = f64::INFINITY;
            for v in verts {
                let (text, off) = match result.inserted.vertex_source[v as usize] {
                    VertexSource::Old(o) => (format!("old{o}"), f64::INFINITY),
                    VertexSource::EdgePoint(i) => {
                        let t = cut.edge_points[i as usize].t;
                        (format!("edge(t={t:.2e})"), t.min(1.0 - t))
                    }
                    VertexSource::Interior { path, index } => {
                        let t = cut.paths[path as usize].interior[index as usize].t;
                        (format!("diag(t={t:.2e})"), t.min(1.0 - t))
                    }
                    VertexSource::Junction(j) => {
                        let b = cut.junctions[j as usize].barycentric;
                        let lo = b.iter().copied().fold(f64::INFINITY, f64::min);
                        (format!("inner(bary_min={lo:.2e})"), lo)
                    }
                };
                min_offset = min_offset.min(off);
                desc.push(text);
            }
            let shortest = e.iter().copied().fold(f64::INFINITY, f64::min);
            let longest = e.iter().copied().fold(0.0, f64::max);
            let class = if min_offset < 1e-2 && desc.iter().any(|d| d.starts_with("old")) {
                "a point close to a corner/side with a vertex of the mesh"
            } else if shortest / longest < 0.05 {
                "short piece between two cut vertices"
            } else if min_offset < 1e-2 {
                "a point close to a side"
            } else {
                "other (ear clipping or the flattening)"
            };
            println!(
                "FLIPPED tri {t}: area/mean {:.3e}, min angle {min_angle:.2} deg, shortest/longest edge {:.3e}, vertices {}, offset {min_offset:.2e}: {class}",
                tri_area_m(m, tri) / mean,
                shortest / longest,
                desc.join(" "),
            );
        }
        let (_, rep) = session
            .describe_cut_with_report(DEFAULT_CROSSING_OFFSET)
            .unwrap();
        Outcome {
            flipped,
            n_clamped: cut.n_clamped_crossings,
            merged: rep.n_merged_nodes,
            displacement: rep.max_displacement_relative,
        }
    };
    println!("CASE {case} {tri_name}");
    for (i, &o) in offsets.iter().enumerate() {
        let out = run(o, i == 0);
        println!(
            "SWEEP {case} {tri_name} merge {merge} min_offset {o:.0e}: {} flipped, {} moved points, {} bends merged, displacement {:.2e} of an edge",
            out.flipped, out.n_clamped, out.merged, out.displacement
        );
    }
}

fn tri_area_m(m: &PolygonMesh, t: &[u32; 3]) -> f64 {
    tri_area(m, *t)
}
