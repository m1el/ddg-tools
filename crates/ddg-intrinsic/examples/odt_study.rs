//! The intrinsic optimal Delaunay triangulation (paper Sec. 4.3, Fig. 14, 17) against
//! Delaunay flipping and Delaunay refinement, at comparable vertex counts.
//!
//! ```text
//! cargo run --release -p ddg-intrinsic --example odt_study -- [options] MESH ...
//!   MESH            an OBJ file, or one of the generated meshes: sheared135, sheared165,
//!                   cells1x16, jitter33 (flat, 32 cells per side)
//!   --modes LIST    comma-separated subset of input,idt,idr25,idr30,odt,odt25,odt30
//!                   (default all). odtNN: iODT with the target length chosen by bisection so
//!                   that the vertex count is closest to that of iDR at NN degrees; odt: the
//!                   default options (target = mean input edge).
//!   --no-cond       skip the condition numbers
//!   --iterations N  iteration limit of iODT (default: the default options)
//! ```
//!
//! One `row |` line per mode: vertices, inserted, iterations and status, the time of the
//! operation alone (the calibration runs are not counted), angles (min, 1 %, median, 99 %, max;
//! percent of corners below 30° and above 90°), areas (max/min, coefficient of variation), the
//! condition number of the pencil (L, M) (cotan Laplacian, lumped masses: `λ_max / λ_2`), and
//! the consistency of the result (largest angle-sum error at a vertex, largest relative
//! mismatch of an edge trace over the input mesh). Timings are single-thread; run one process
//! at a time when they matter.

use std::time::Instant;

use ddg_core::generators::{SplitMix64, grid_disk};
use ddg_core::io::read_obj;
use ddg_core::sparse::{Method, Parallelism, Symbolic};
use ddg_core::{Geometry, TriMesh};
use ddg_intrinsic::{
    EdgeTraceOptions, IntrinsicTriangulation, OptimalDelaunayOptions, OptimalDelaunayReport,
    RefineOptions, TargetLength,
};

fn load(path: &str) -> Result<(TriMesh, Geometry), String> {
    let (positions, faces) = match path {
        "sheared135" | "sheared165" | "cells1x16" | "jitter33" => generated(path)?,
        _ => {
            let obj = read_obj(path).map_err(|e| e.to_string())?;
            let mut faces = Vec::new();
            for f in 0..obj.n_faces() {
                let corners = &obj.corner_verts[obj.face_corners(f)];
                for k in 1..corners.len().saturating_sub(1) {
                    faces.push([corners[0], corners[k], corners[k + 1]]);
                }
            }
            (obj.positions, faces)
        }
    };
    let mesh = TriMesh::new(positions.len(), faces).map_err(|e| e.to_string())?;
    let geometry = Geometry::from_positions(&mesh, &positions).map_err(|e| e.to_string())?;
    Ok((mesh, geometry))
}

/// Positions and triangles.
type Soup = (Vec<[f64; 3]>, Vec<[u32; 3]>);

/// Flat test meshes with the same surface (a parallelogram or rectangle) and bad triangles.
fn generated(name: &str) -> Result<Soup, String> {
    let n = 32;
    let (nx, ny) = if name == "cells1x16" { (n, 16 * n) } else { (n, n) };
    let g = grid_disk(nx, ny).map_err(|e| e.to_string())?;
    let mut positions = g.positions;
    match name {
        "sheared135" => positions.iter_mut().for_each(|p| p[0] += p[1]),
        "sheared165" => {
            let s = 1.0 / 15f64.to_radians().tan();
            positions.iter_mut().for_each(|p| p[0] += s * p[1]);
        }
        "cells1x16" => positions.iter_mut().for_each(|p| p[1] *= 16.0),
        "jitter33" => {
            let h = 1.0 / n as f64;
            let mut rng = SplitMix64::new(5);
            for (i, p) in positions.iter_mut().enumerate() {
                let (x, y) = (i / (ny + 1), i % (ny + 1));
                if x > 0 && x < nx && y > 0 && y < ny {
                    p[0] += 0.33 * h * rng.uniform(-1.0, 1.0);
                    p[1] += 0.33 * h * rng.uniform(-1.0, 1.0);
                }
            }
        }
        _ => {}
    }
    Ok((positions, g.faces))
}

/// `λ_max / λ_2` of `L x = λ M x` (cotan Laplacian, lumped mass): power iteration for the
/// largest eigenvalue, inverse iteration with a sparse Cholesky factorization of `L + εM` for
/// the smallest non-zero one (constants removed in the `M` inner product).
fn condition(tri: &IntrinsicTriangulation) -> Result<(f64, f64), String> {
    let n = tri.n_vertices();
    let l = tri.cotan_laplacian();
    let m = tri.lumped_mass();
    let mut rng = SplitMix64::new(1);
    let total: f64 = m.iter().sum();
    let deflate = |x: &mut [f64]| {
        let c: f64 = x.iter().zip(&m).map(|(a, b)| a * b).sum::<f64>() / total;
        x.iter_mut().for_each(|a| *a -= c);
    };
    let rayleigh = |x: &[f64]| -> f64 {
        let mut y = vec![0.0; n];
        l.mul_vec(x, &mut y).unwrap();
        let num: f64 = x.iter().zip(&y).map(|(a, b)| a * b).sum();
        let den: f64 = x.iter().zip(&m).map(|(a, b)| a * a * b).sum();
        num / den
    };
    // Largest: power iteration on M^{-1} L.
    let mut x: Vec<f64> = (0..n).map(|_| rng.uniform(-1.0, 1.0)).collect();
    let mut y = vec![0.0; n];
    for _ in 0..400 {
        l.mul_vec(&x, &mut y).map_err(|e| e.to_string())?;
        for i in 0..n {
            x[i] = y[i] / m[i];
        }
        let norm = x.iter().map(|a| a * a).sum::<f64>().sqrt();
        x.iter_mut().for_each(|a| *a /= norm);
    }
    let lambda_max = rayleigh(&x);
    // Smallest non-zero: inverse iteration.
    let eps = 1e-8 * lambda_max;
    let mut k = l.clone();
    k.add_to_diagonal(&m.iter().map(|w| eps * w).collect::<Vec<_>>())
        .map_err(|e| e.to_string())?;
    let factor = Symbolic::analyze(&k, Method::Llt)
        .and_then(|s| s.factor(&k, Parallelism::Sequential))
        .map_err(|e| e.to_string())?;
    let mut x: Vec<f64> = (0..n).map(|_| rng.uniform(-1.0, 1.0)).collect();
    deflate(&mut x);
    let mut lambda_2 = f64::INFINITY;
    for _ in 0..200 {
        let mut rhs: Vec<f64> = x.iter().zip(&m).map(|(a, b)| a * b).collect();
        factor
            .solve_in_place(&mut rhs, Parallelism::Sequential)
            .map_err(|e| e.to_string())?;
        deflate(&mut rhs);
        let norm = rhs.iter().map(|a| a * a).sum::<f64>().sqrt();
        x = rhs.into_iter().map(|a| a / norm).collect();
        let next = rayleigh(&x);
        if (next - lambda_2).abs() < 1e-10 * next {
            lambda_2 = next;
            break;
        }
        lambda_2 = next;
    }
    Ok((lambda_max, lambda_2))
}

struct Run {
    tri: IntrinsicTriangulation,
    seconds: f64,
    note: String,
}

fn odt(
    base: &IntrinsicTriangulation,
    target: TargetLength,
    iterations: Option<usize>,
) -> Result<(IntrinsicTriangulation, OptimalDelaunayReport, f64), String> {
    let mut tri = base.clone();
    let mut options = OptimalDelaunayOptions {
        target_length: target,
        ..OptimalDelaunayOptions::default()
    };
    if let Some(n) = iterations {
        options.max_iterations = n;
    }
    let start = Instant::now();
    let report = tri.optimal_delaunay(&options).map_err(|e| e.to_string())?;
    Ok((tri, report, start.elapsed().as_secs_f64()))
}

fn odt_note(r: &OptimalDelaunayReport) -> String {
    let halved: usize = r.iterations.iter().map(|it| it.n_halved).sum();
    let stuck: usize = r.iterations.iter().map(|it| it.n_not_moved).sum();
    let refused: usize = r.iterations.iter().map(|it| it.n_split_refused).sum();
    format!(
        "{} it {:?}, target {:.4e}, last motion {:.2e} (rel {:.2e}), halved {halved}, \
         not moved {stuck}, splits refused {refused}, not Delaunay {}",
        r.n_iterations,
        r.status,
        r.target_length.unwrap_or(0.0),
        r.last_max_motion,
        r.last_max_relative_motion,
        r.delaunay.n_not_delaunay
    )
}

/// iODT with the target length whose vertex count is closest to `goal`.
fn odt_matched(
    base: &IntrinsicTriangulation,
    goal: usize,
    iterations: Option<usize>,
) -> Result<Run, String> {
    let (mut lo, mut hi) = (0.2f64, 6.0f64); // factors of the mean input edge
    let mut best: Option<(usize, IntrinsicTriangulation, OptimalDelaunayReport, f64, f64)> = None;
    for _ in 0..12 {
        let k = (lo * hi).sqrt();
        let (tri, report, seconds) = odt(base, TargetLength::MeanInputEdge(k), iterations)?;
        let count = tri.n_vertices();
        let miss = count.abs_diff(goal);
        if best.as_ref().is_none_or(|b| miss < b.0) {
            best = Some((miss, tri, report, seconds, k));
        }
        if count > goal {
            lo = k;
        } else {
            hi = k;
        }
        if miss * 100 < goal {
            break;
        }
    }
    let (_, tri, report, seconds, k) = best.ok_or("no run")?;
    Ok(Run {
        tri,
        seconds,
        note: format!("factor {k:.3}: {}", odt_note(&report)),
    })
}

fn main() {
    let mut modes: Vec<String> = ["input", "idt", "idr25", "idr30", "odt", "odt25", "odt30"]
        .map(String::from)
        .to_vec();
    let mut cond = true;
    let mut iterations = None;
    let mut meshes = Vec::new();
    let mut args = std::env::args().skip(1);
    while let Some(a) = args.next() {
        match a.as_str() {
            "--modes" => {
                modes = args
                    .next()
                    .expect("--modes needs a value")
                    .split(',')
                    .map(String::from)
                    .collect()
            }
            "--no-cond" => cond = false,
            "--iterations" => {
                iterations = Some(args.next().and_then(|v| v.parse().ok()).expect("a number"))
            }
            _ => meshes.push(a),
        }
    }
    println!(
        "row | mesh | mode | V | inserted | seconds | min° | 1%° | median° | 99%° | max° | <30° % | \
         >90° % | area max/min | area cv | λ_max/λ_2 | angle-sum err | trace mismatch | note"
    );
    for path in &meshes {
        let name = std::path::Path::new(path)
            .file_stem()
            .map_or(path.clone(), |s| s.to_string_lossy().to_string());
        let (mesh, geometry) = match load(path) {
            Ok(x) => x,
            Err(e) => {
                println!("error | {name} | {e}");
                continue;
            }
        };
        let fresh = IntrinsicTriangulation::new(&mesh, &geometry).expect("a valid mesh");
        let mut counts = std::collections::HashMap::new();
        for mode in &modes {
            let run: Result<Run, String> = (|| match mode.as_str() {
                "input" => Ok(Run {
                    tri: fresh.clone(),
                    seconds: 0.0,
                    note: String::new(),
                }),
                "idt" => {
                    let mut tri = fresh.clone();
                    let start = Instant::now();
                    let r = tri
                        .flip_to_delaunay(&Default::default())
                        .map_err(|e| e.to_string())?;
                    Ok(Run {
                        tri,
                        seconds: start.elapsed().as_secs_f64(),
                        note: format!("{} flips", r.n_flips),
                    })
                }
                "idr25" | "idr30" => {
                    let degrees: f64 = mode[3..].parse().unwrap();
                    let mut tri = fresh.clone();
                    let start = Instant::now();
                    let r = tri
                        .refine(&RefineOptions {
                            min_angle: degrees.to_radians(),
                            ..RefineOptions::default()
                        })
                        .map_err(|e| e.to_string())?;
                    Ok(Run {
                        tri,
                        seconds: start.elapsed().as_secs_f64(),
                        note: format!("{:?}, {} below the bound", r.status, r.below_bound.len()),
                    })
                }
                "odt" => {
                    let (tri, report, seconds) = odt(
                        &fresh,
                        OptimalDelaunayOptions::default().target_length,
                        iterations,
                    )?;
                    Ok(Run {
                        tri,
                        seconds,
                        note: odt_note(&report),
                    })
                }
                "odt25" | "odt30" => {
                    let goal = *counts
                        .get(&format!("idr{}", &mode[3..]))
                        .ok_or("run the matching idr mode first")?;
                    odt_matched(&fresh, goal, iterations)
                }
                other => Err(format!("unknown mode {other}")),
            })();
            let run = match run {
                Ok(r) => r,
                Err(e) => {
                    println!("error | {name} | {mode} | {e}");
                    continue;
                }
            };
            let tri = &run.tri;
            counts.insert(mode.clone(), tri.n_vertices());
            let s = tri.shape_statistics();
            let deg = f64::to_degrees;
            let kappa = if cond {
                match condition(tri) {
                    Ok((max, min)) => format!("{:.3e}", max / min),
                    Err(e) => format!("({e})"),
                }
            } else {
                "-".to_string()
            };
            let check = tri.check().map_err(|e| e.to_string());
            let traces = tri
                .trace_edges(&EdgeTraceOptions::default())
                .map(|t| format!("{:.1e}", t.report.max_relative_mismatch))
                .unwrap_or_else(|e| format!("FAILED: {e}"));
            println!(
                "row | {name} | {mode} | {} | {} | {:.3} | {:.3} | {:.2} | {:.2} | {:.2} | {:.2} | \
                 {:.2} | {:.2} | {:.3e} | {:.3} | {kappa} | {} | {traces} | {}",
                tri.n_vertices(),
                tri.n_vertices() - tri.n_input_vertices(),
                run.seconds,
                deg(s.angles.min),
                deg(s.angles.quantiles[0]),
                deg(s.angles.quantiles[2]),
                deg(s.angles.quantiles[4]),
                deg(s.angles.max),
                100.0 * s.angles_below[2],
                100.0 * s.angles_above[0],
                s.area_ratio,
                s.area_cv,
                check.map_or_else(|e| format!("FAILED: {e}"), |c| format!("{:.1e}", c.max_angle_sum_error)),
                run.note,
            );
            // Histogram of the corner angles (10° bins), for the figure.
            let mut bins = [0usize; 18];
            for h in 0..tri.n_halfedges() {
                let b = (deg(tri.corner_angle(h)) / 10.0) as usize;
                bins[b.min(17)] += 1;
            }
            println!("hist | {name} | {mode} | {bins:?}");
        }
    }
}
