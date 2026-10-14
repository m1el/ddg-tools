//! Research worker for `scripts/investigate_camel.py uv-opt`.
//! Reads lines of zero-based `vertex:angle` pairs on stdin, flattens with BFF,
//! and prints area distortion both including and excluding reversed faces.
//! `published` measures the dataset UVs. Does not change production behavior.
use ddg_bff::{ConeFlattenInput, ConeFlattenOptions, CutSource, flatten_with_cones};
use ddg_core::io::{Cone, read_obj, write_obj};
use ddg_core::operators::angle_defects;
use ddg_core::{Geometry, TriMesh};
use ddg_optimal_cones::{
    BoundaryTreatment, ConeProblem, SolveOptions, Status, solve_with_boundary,
};
use std::io::{self, BufRead, Write};
use std::path::PathBuf;

fn measure(mesh: &TriMesh, geo: &Geometry, tips: &[bool], uv: &[[f64; 2]]) -> String {
    let areas = geo.face_areas();
    let signed: Vec<f64> = (0..mesh.n_faces())
        .map(|f| {
            let [p, q, r] = [uv[3 * f], uv[3 * f + 1], uv[3 * f + 2]];
            0.5 * ((q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0]))
        })
        .collect();
    let orientation = if signed.iter().filter(|b| **b < 0.).count() > mesh.n_faces() / 2 {
        -1.
    } else {
        1.
    };
    let mut sums = [[0.; 3]; 2];
    let mut flips = 0;
    let mut spike_flips = 0;
    let mut degenerate = 0;
    for (f, (&a, &b)) in areas.iter().zip(&signed).enumerate() {
        if b == 0. || !b.is_finite() {
            degenerate += 1;
            continue;
        }
        let flipped = b * orientation < 0.;
        if flipped {
            flips += 1;
            if mesh.faces()[f].iter().any(|v| tips[*v as usize]) {
                spike_flips += 1;
            }
        }
        let q = 0.5 * (b.abs() / a).ln();
        let w = a / geo.total_area();
        for (k, s) in sums.iter_mut().enumerate() {
            if k == 1 && flipped {
                continue;
            }
            s[0] += w;
            s[1] += w * q;
            s[2] += w * q * q;
        }
    }
    let e: Vec<f64> = sums.iter().map(|s| s[2] - s[1] * s[1] / s[0]).collect();
    format!(
        "{:.16e} {:.16e} {} {} {}",
        e[0], e[1], flips, spike_flips, degenerate
    )
}

fn main() {
    let name = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "camelshotnoise".into());
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../data/optimal-cones/optimal-cones-data/data");
    let obj = read_obj(dir.join(format!("{name}.obj"))).unwrap();
    let mesh = TriMesh::new(obj.positions.len(), obj.triangles().unwrap()).unwrap();
    let geo = Geometry::from_positions(&mesh, &obj.positions).unwrap();
    let tips: Vec<bool> = angle_defects(&mesh, &geo)
        .unwrap()
        .iter()
        .map(|k| *k > 6.)
        .collect();
    let mut proxy_positions = obj.positions.clone();
    for v in 0..mesh.n_vertices() {
        if tips[v] {
            let ring: Vec<usize> = mesh.vertex_neighbors(v).collect();
            proxy_positions[v] = std::array::from_fn(|k| {
                ring.iter().map(|w| obj.positions[*w][k]).sum::<f64>() / ring.len() as f64
            });
        }
    }
    let proxy = Geometry::from_positions(&mesh, &proxy_positions).unwrap();
    let options = ConeFlattenOptions {
        overlap_samples_per_triangle: 0.,
        ..Default::default()
    };
    for line in io::stdin().lock().lines() {
        let line = line.unwrap();
        if line == "published" {
            println!(
                "{}",
                measure(&mesh, &geo, &tips, &obj.corner_uvs().unwrap())
            );
        } else if line.starts_with("generate ") {
            let args: Vec<&str> = line.split_whitespace().collect();
            let lambda: f64 = args[1].parse().unwrap();
            let boundary = if args[2] == "free" {
                BoundaryTreatment::Free
            } else {
                BoundaryTreatment::Fixed
            };
            let problem = ConeProblem::new(&mesh, &proxy, lambda);
            let solve_options = SolveOptions {
                merge_cone_clusters: true,
                max_merge_cluster_size: usize::MAX,
                ..Default::default()
            };
            match solve_with_boundary(&problem, boundary, &solve_options) {
                Ok(s) if s.report.status == Status::Converged => {
                    let mut costs = vec![f64::INFINITY; mesh.n_vertices()];
                    for c in &s.cones {
                        costs[c.vertex as usize] = 0.;
                    }
                    let refit = ConeProblem {
                        cone_costs: Some(costs),
                        ..ConeProblem::new(&mesh, &proxy, 1.)
                    };
                    match solve_with_boundary(
                        &refit,
                        BoundaryTreatment::Free,
                        &SolveOptions::default(),
                    ) {
                        Ok(r) if r.report.status == Status::Converged => {
                            println!(
                                "cones {}",
                                r.cones
                                    .iter()
                                    .map(|c| format!("{}:{:.16e}", c.vertex, c.angle))
                                    .collect::<Vec<_>>()
                                    .join(" ")
                            );
                        }
                        _ => println!("error refit did not converge"),
                    }
                }
                _ => println!("error placement did not converge"),
            }
        } else {
            let (output, cone_text) = if let Some(rest) = line.strip_prefix("save ") {
                let (path, pairs) = rest.split_once(' ').unwrap();
                (Some(path), pairs)
            } else {
                (None, line.as_str())
            };
            let cones: Vec<Cone> = cone_text
                .split_whitespace()
                .map(|item| {
                    let (v, a) = item.split_once(':').unwrap();
                    Cone {
                        vertex: v.parse().unwrap(),
                        angle: a.parse().unwrap(),
                    }
                })
                .collect();
            let input = ConeFlattenInput {
                mesh: &mesh,
                geometry: &geo,
                cones: &cones,
                scale_factors: None,
                cut: CutSource::Plan {
                    may_cut: None,
                    must_cut: &[],
                },
                anchors: &[],
            };
            match flatten_with_cones(&input, &options) {
                Ok(flat) => {
                    if let Some(path) = output {
                        write_obj(
                            path,
                            &obj.positions,
                            &obj.corner_verts,
                            &obj.face_starts,
                            Some(&flat.corner_uvs),
                        )
                        .unwrap();
                    }
                    println!("{}", measure(&mesh, &geo, &tips, &flat.corner_uvs));
                }
                Err(e) => println!("error {e}"),
            }
        }
        io::stdout().flush().unwrap();
    }
}
