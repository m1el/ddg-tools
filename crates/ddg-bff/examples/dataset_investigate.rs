//! Diagnostics: where does the stored UV map of a dataset mesh differ from ours?
//! (Scratch tool; prints per-mesh statistics.)

#[path = "../tests/common/mod.rs"]
mod common;

use common::{dataset_names, load};
use ddg_bff::quality::{log_scale_factor, signed_double_area};
use ddg_bff::{
    ConeFlattenInput, ConeFlattenOptions, CutSource, GaussBonnetPolicy, ScaleFactorOptions,
    flatten_with_cones,
};

fn quantiles(mut x: Vec<f64>) -> [f64; 5] {
    x.sort_by(|a, b| a.total_cmp(b));
    let n = x.len();
    let at = |q: f64| x[((n - 1) as f64 * q).round() as usize];
    [at(0.0), at(0.01), at(0.5), at(0.99), at(1.0)]
}

fn main() {
    let mut names = dataset_names().expect("dataset");
    let wanted: Vec<String> = std::env::args().skip(1).collect();
    if !wanted.is_empty() {
        names.retain(|n| wanted.contains(n));
    }
    println!(
        "{:<20} | {:>9} {:>9} | {:>9} {:>9} {:>9} | {:>9} {:>9} | {:>9} {:>9} {:>9}",
        "mesh", "bdry min", "bdry max", "seam p50", "seam p99", "seam max",
        "du rms", "du max", "cone p50", "cone max", "noncone"
    );
    println!("(then: boundary window max | per-vertex stored/ours log length ratio: min, p1, p99, max)");
    for name in &names {
        let data = load(name);
        let mesh = &data.mesh;
        let uv = &data.reference_uvs;
        let len = |h: usize| {
            let (a, b) = (uv[h], uv[mesh.next(h)]);
            (b[0] - a[0]).hypot(b[1] - a[1])
        };
        // 1. Original boundary: stored UV length / 3D length (u = 0 there means constant).
        let mut boundary: Vec<f64> = (0..mesh.n_halfedges())
            .filter(|&h| mesh.twin(h).is_none())
            .map(|h| (len(h) / data.geometry.edge_lengths()[mesh.edge(h)]).ln())
            .collect();
        let (bmin, bmax) = if boundary.is_empty() {
            (f64::NAN, f64::NAN)
        } else {
            let q = quantiles(std::mem::take(&mut boundary));
            (q[0] - q[2], q[4] - q[2])
        };
        // 2. Seam mismatch of the stored map.
        let seams: Vec<f64> = data
            .reference_cut
            .iter()
            .map(|&e| {
                let h = mesh.edge_halfedge(e as usize);
                let t = mesh.twin(h).unwrap();
                (len(h) - len(t)).abs() / len(h).max(len(t))
            })
            .collect();
        let sq = quantiles(seams);
        // 3. Per-face log scale factor: stored vs ours (same cut).
        let ours = flatten_with_cones(
            &ConeFlattenInput {
                mesh,
                geometry: &data.geometry,
                cones: &data.cones,
                scale_factors: None,
                cut: CutSource::Edges(&data.reference_cut),
                anchors: &[],
            },
            &ConeFlattenOptions {
                scale_factors: ScaleFactorOptions {
                    gauss_bonnet: GaussBonnetPolicy::SpreadOverVertices,
                    ..ScaleFactorOptions::default()
                },
                ..ConeFlattenOptions::default()
            },
        );
        let Ok(ours) = ours else {
            println!("{name:<20} | flatten failed: {}", ours.unwrap_err());
            continue;
        };
        let areas = data.geometry.face_areas();
        let face_u = |layout: &[[f64; 2]], f: usize| {
            log_scale_factor(
                areas[f],
                0.5 * signed_double_area(layout[3 * f], layout[3 * f + 1], layout[3 * f + 2]).abs(),
            )
        };
        let mut diff: Vec<f64> = (0..mesh.n_faces())
            .map(|f| face_u(uv, f) - face_u(&ours.corner_uvs, f))
            .collect();
        let total: f64 = areas.iter().sum();
        let mean: f64 = diff.iter().zip(areas).map(|(d, a)| d * a).sum::<f64>() / total;
        for d in &mut diff {
            *d -= mean;
        }
        let rms = (diff.iter().zip(areas).map(|(d, a)| d * d * a).sum::<f64>() / total).sqrt();
        let dmax = diff.iter().fold(0.0f64, |m, d| m.max(d.abs()));
        // 4. Realized angle sums of the stored map at cones and at other interior vertices.
        let mut sums = vec![0.0; mesh.n_vertices()];
        for f in 0..mesh.n_faces() {
            for k in 0..3 {
                let (a, b, c) = (uv[3 * f + k], uv[3 * f + (k + 1) % 3], uv[3 * f + (k + 2) % 3]);
                let (ux, uy, vx, vy) = (b[0] - a[0], b[1] - a[1], c[0] - a[0], c[1] - a[1]);
                sums[mesh.faces()[f][k] as usize] += (ux * vy - uy * vx).atan2(ux * vx + uy * vy);
            }
        }
        let mut is_cone = vec![false; mesh.n_vertices()];
        let mut cone_err = Vec::new();
        for c in &data.cones {
            is_cone[c.vertex as usize] = true;
            cone_err.push((2.0 * std::f64::consts::PI - sums[c.vertex as usize] - c.angle).abs());
        }
        // Cone angles from the boundary polygon alone (not affected by flipped or degenerate
        // triangles): at the copies p of a cone on the cut, the interior angles of the
        // boundary polygon sum to 2 pi - mu modulo 2 pi.
        let cut_mesh_ = ddg_bff::cut_mesh(mesh, &data.reference_cut).unwrap();
        let polygon_angle = |layout: &[[f64; 2]]| -> Vec<f64> {
            let mut vertex_uv = vec![[0.0; 2]; cut_mesh_.mesh.n_vertices()];
            for (f, face) in cut_mesh_.mesh.faces().iter().enumerate() {
                for k in 0..3 {
                    vertex_uv[face[k] as usize] = layout[3 * f + k];
                }
            }
            let mut interior = vec![0.0; mesh.n_vertices()];
            for boundary in cut_mesh_.mesh.boundary_loops() {
                let n = boundary.len();
                for i in 0..n {
                    let (a, b, c) = (
                        vertex_uv[boundary[(i + n - 1) % n] as usize],
                        vertex_uv[boundary[i] as usize],
                        vertex_uv[boundary[(i + 1) % n] as usize],
                    );
                    let (ux, uy, vx, vy) = (b[0] - a[0], b[1] - a[1], c[0] - b[0], c[1] - b[1]);
                    let turn = (ux * vy - uy * vx).atan2(ux * vx + uy * vy);
                    interior[cut_mesh_.vertex_to_original[boundary[i] as usize] as usize] +=
                        std::f64::consts::PI - turn;
                }
            }
            interior
        };
        let wrap = |x: f64| {
            let t = std::f64::consts::TAU;
            let y = x.rem_euclid(t);
            if y > t / 2.0 { y - t } else { y }
        };
        let stored_polygon = polygon_angle(uv);
        let ours_polygon = polygon_angle(&ours.corner_uvs);
        let mut worst_polygon = (0.0f64, 0u32, 0.0f64);
        for c in &data.cones {
            let v = c.vertex as usize;
            let stored_error = wrap(2.0 * std::f64::consts::PI - stored_polygon[v] - c.angle);
            let ours_error = wrap(2.0 * std::f64::consts::PI - ours_polygon[v] - c.angle);
            if stored_error.abs() > worst_polygon.0.abs() {
                worst_polygon = (stored_error, c.vertex, ours_error);
            }
        }
        println!(
            "  {name}: largest (stored cone angle - listed) mod 2pi from the boundary polygon: {:.3e} at vertex {} (ours there: {:.1e})",
            worst_polygon.0, worst_polygon.1, worst_polygon.2
        );
        // Vertices on the stored cut that are not listed as cones but have a cone angle in
        // the stored map.
        {
            let copies = cut_mesh_.copies_per_vertex();
            let mut listed = vec![false; mesh.n_vertices()];
            for c in &data.cones {
                listed[c.vertex as usize] = true;
            }
            let mut on_cut = vec![false; mesh.n_vertices()];
            for &e in &data.reference_cut {
                let [a, b] = mesh.edge_vertices(e as usize);
                on_cut[a] = true;
                on_cut[b] = true;
            }
            let mut unlisted: Vec<(usize, f64)> = (0..mesh.n_vertices())
                .filter(|&v| on_cut[v] && !listed[v] && !mesh.is_boundary_vertex(v))
                .map(|v| (v, wrap(2.0 * std::f64::consts::PI - stored_polygon[v])))
                .collect();
            let total: f64 = unlisted.iter().map(|x| x.1).sum();
            let n_cut_vertices = unlisted.len();
            unlisted.retain(|x| x.1.abs() > 1e-4);
            unlisted.sort_by(|a, b| b.1.abs().total_cmp(&a.1.abs()));
            let _ = copies;
            println!(
                "  {name}: {} unlisted vertices on the stored cut, their stored cone angles sum to {:.2e}; {} exceed 1e-4: {:?}",
                n_cut_vertices,
                total,
                unlisted.len(),
                &unlisted[..unlisted.len().min(4)]
            );
        }
        // Does the stored map agree better with a flattening that uses the cone angles
        // measured from the stored map where they differ from the listed ones?
        {
            let mut cones = data.cones.clone();
            let mut changed = 0;
            for c in &mut cones {
                let v = c.vertex as usize;
                let error = wrap(2.0 * std::f64::consts::PI - stored_polygon[v] - c.angle);
                if error.abs() > 1e-4 {
                    c.angle += error;
                    changed += 1;
                }
            }
            // Optionally also treat unlisted cut vertices with a clear cone angle in the
            // stored map as cones (threshold from the environment).
            if let Ok(t) = std::env::var("UNLISTED_THRESHOLD") {
                let t: f64 = t.parse().unwrap();
                let mut listed = vec![false; mesh.n_vertices()];
                for c in &data.cones {
                    listed[c.vertex as usize] = true;
                }
                let mut on_cut = vec![false; mesh.n_vertices()];
                for &e in &data.reference_cut {
                    let [a, b] = mesh.edge_vertices(e as usize);
                    on_cut[a] = true;
                    on_cut[b] = true;
                }
                for v in 0..mesh.n_vertices() {
                    if on_cut[v] && !listed[v] && !mesh.is_boundary_vertex(v) {
                        let angle = wrap(2.0 * std::f64::consts::PI - stored_polygon[v]);
                        if angle.abs() > t {
                            cones.push(ddg_bff::Cone {
                                vertex: v as u32,
                                angle,
                            });
                            changed += 1;
                        }
                    }
                }
            }
            let options = ConeFlattenOptions {
                scale_factors: ScaleFactorOptions {
                    gauss_bonnet: GaussBonnetPolicy::SpreadOverVertices,
                    ..ScaleFactorOptions::default()
                },
                ..ConeFlattenOptions::default()
            };
            let before = common::compare_layouts(
                &ours.face_chart,
                ours.n_charts(),
                &ours.corner_uvs,
                uv,
            );
            if changed > 0 {
                let again = flatten_with_cones(
                    &ConeFlattenInput {
                        mesh,
                        geometry: &data.geometry,
                        cones: &cones,
                        scale_factors: None,
                        cut: CutSource::Edges(&data.reference_cut),
                        anchors: &[],
                    },
                    &options,
                )
                .unwrap();
                let after = common::compare_layouts(
                    &again.face_chart,
                    again.n_charts(),
                    &again.corner_uvs,
                    uv,
                );
                println!(
                    "  {name}: with {changed} cone angle(s) taken from the stored map: residual rms {:.1e} -> {:.1e}, max {:.1e} -> {:.1e}",
                    before.rms, after.rms, before.max, after.max
                );
            }
        }
        if std::env::var("CONE_TABLE").is_ok() {
            println!("  cone: vertex, listed angle, stored realized - listed, ours realized - listed, wedges");
            for (c, r) in data.cones.iter().zip(&ours.cones) {
                println!(
                    "    {:6} {:10.6} {:12.3e} {:12.3e} {:3}",
                    c.vertex,
                    c.angle,
                    2.0 * std::f64::consts::PI - sums[c.vertex as usize] - c.angle,
                    r.realized_angle - c.angle,
                    r.n_wedges
                );
            }
        }
        let cq = quantiles(cone_err);
        let noncone = (0..mesh.n_vertices())
            .filter(|&v| !is_cone[v] && !mesh.is_boundary_vertex(v))
            .map(|v| (2.0 * std::f64::consts::PI - sums[v]).abs())
            .fold(0.0, f64::max);
        // 5. Stored length / our length along the original boundary, smoothed over windows of
        //    consecutive edges to average out the rounding of the file, and the same along
        //    all edges by distance class (near / far from the cut).
        let ours_len = |h: usize| {
            let (a, b) = (ours.corner_uvs[h], ours.corner_uvs[mesh.next(h)]);
            (b[0] - a[0]).hypot(b[1] - a[1])
        };
        let all: Vec<f64> = (0..mesh.n_halfedges()).map(|h| (len(h) / ours_len(h)).ln()).collect();
        let median = quantiles(all.iter().copied().filter(|x| x.is_finite()).collect())[2];
        let mut window_max = 0.0f64;
        for boundary in mesh.boundary_loops() {
            let ratios: Vec<f64> = boundary
                .iter()
                .map(|&v| all[mesh.outgoing_halfedge(v as usize)] - median)
                .collect();
            let w = 16.min(ratios.len());
            for start in 0..ratios.len() {
                let mean: f64 =
                    (0..w).map(|k| ratios[(start + k) % ratios.len()]).sum::<f64>() / w as f64;
                window_max = window_max.max(mean.abs());
            }
        }
        if std::env::var("BOUNDARY_PROFILE").is_ok() {
            for boundary in mesh.boundary_loops() {
                let w = 16;
                let stored: Vec<f64> = boundary
                    .iter()
                    .map(|&v| {
                        let h = mesh.outgoing_halfedge(v as usize);
                        (len(h) / data.geometry.edge_lengths()[mesh.edge(h)]).ln()
                    })
                    .collect();
                let mine: Vec<f64> = boundary
                    .iter()
                    .map(|&v| {
                        let h = mesh.outgoing_halfedge(v as usize);
                        (ours_len(h) / data.geometry.edge_lengths()[mesh.edge(h)]).ln()
                    })
                    .collect();
                let ms = quantiles(stored.clone())[2];
                let mm = quantiles(mine.clone())[2];
                println!("  boundary loop of {} edges: window means of ln(uv length / 3D length), stored | ours", boundary.len());
                for start in (0..boundary.len()).step_by(w) {
                    let end = (start + w).min(boundary.len());
                    let a: f64 = stored[start..end].iter().sum::<f64>() / (end - start) as f64 - ms;
                    let b: f64 = mine[start..end].iter().sum::<f64>() / (end - start) as f64 - mm;
                    println!("    {:4} {:9.2e} | {:9.2e}", start, a, b);
                }
            }
        }
        // The same smoothing over the faces around each vertex: a per-vertex estimate of
        // (stored u) - (our u), and its largest value.
        let mut vertex_sum = vec![0.0; mesh.n_vertices()];
        let mut vertex_count = vec![0.0; mesh.n_vertices()];
        for h in 0..mesh.n_halfedges() {
            if all[h].is_finite() {
                for v in [mesh.from_vertex(h), mesh.to_vertex(h)] {
                    vertex_sum[v] += all[h] - median;
                    vertex_count[v] += 1.0;
                }
            }
        }
        if std::env::var("DISTANCE_PROFILE").is_ok() {
            // Mean of the per-vertex estimate by edge distance from the boundary / the cones.
            let bfs = |sources: Vec<usize>| {
                let mut d = vec![usize::MAX; mesh.n_vertices()];
                let mut queue = std::collections::VecDeque::new();
                for s in sources {
                    d[s] = 0;
                    queue.push_back(s);
                }
                while let Some(v) = queue.pop_front() {
                    for w in mesh.vertex_neighbors(v) {
                        if d[w] == usize::MAX {
                            d[w] = d[v] + 1;
                            queue.push_back(w);
                        }
                    }
                }
                d
            };
            let from_boundary = bfs(mesh.boundary_vertices().iter().map(|&v| v as usize).collect());
            let from_cones = bfs(data.cones.iter().map(|c| c.vertex as usize).collect());
            for (label, d) in [("boundary", &from_boundary), ("cones", &from_cones)] {
                println!("  mean and rms of (stored - ours) log length ratio by distance from the {label}:");
                let max = d.iter().copied().filter(|&x| x != usize::MAX).max().unwrap_or(0);
                for k in (0..=max).filter(|k| *k < 6 || k % 5 == 0) {
                    let values: Vec<f64> = (0..mesh.n_vertices())
                        .filter(|&v| d[v] == k && vertex_count[v] > 0.0)
                        .map(|v| vertex_sum[v] / vertex_count[v])
                        .collect();
                    if values.is_empty() {
                        continue;
                    }
                    let mean = values.iter().sum::<f64>() / values.len() as f64;
                    let rms = (values.iter().map(|x| (x - mean) * (x - mean)).sum::<f64>()
                        / values.len() as f64)
                        .sqrt();
                    println!("    {:3} n={:5} mean {:9.2e} rms {:9.2e}", k, values.len(), mean, rms);
                }
            }
        }
        let vq = quantiles(
            (0..mesh.n_vertices())
                .filter(|&v| vertex_count[v] > 0.0)
                .map(|v| vertex_sum[v] / vertex_count[v])
                .collect(),
        );
        println!(
            "{:<20} | {:9.1e} {:9.1e} | {:9.1e} {:9.1e} {:9.1e} | {:9.1e} {:9.1e} | {:9.1e} {:9.1e} {:9.1e} | {:9.1e} | {:9.1e} {:9.1e} {:9.1e} {:9.1e}",
            name, bmin, bmax, sq[2], sq[3], sq[4], rms, dmax, cq[2], cq[4], noncone, window_max,
            vq[0], vq[1], vq[3], vq[4]
        );
    }
}
