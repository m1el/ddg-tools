//! Errors and speedups of the applications on the coarse mesh against the fine solutions
//! (paper Sec. 8.5, Fig. 21-25; SPEC §10.1). Single thread for every timing except the line
//! marked "parallel"; speedups are ratios of wall times on the same machine.
//!
//! ```text
//! cargo run --release -p ddg-intrinsic-simplification --example applications -- \
//!     MESH.obj [--fractions 0.5,0.1,0.02] [--sources 5] [--rows 8] [--voronoi 20] \
//!     [--exp 100] [--threads 8]
//! ```
//!
//! The coarse mesh comes from `make_coarse` below, which uses the stand-in order
//! `greedy_coarsening` (smallest memoryless cost). **Swap that one function for `simplify` (the
//! ICE order) to rerun the numbers with the paper's method**: nothing else here depends on how the
//! coarse mesh was made (it takes a `CoarseMesh` with the fine endpoints of the exponential-map
//! test tracked as extra points).

use std::time::Instant;

use ddg_core::io::read_obj;
use ddg_core::operators;
use ddg_core::sparse::ThreadPool;
use ddg_core::{Geometry, TriMesh};
use ddg_intrinsic_apps::exact::{ExactOptions, ExactSource, exact_distance};
use ddg_intrinsic_simplification::applications::{
    CoarseApplications, CoarseHeat, DistanceError, PoissonSolver, greedy_coarsening,
};
use ddg_intrinsic_simplification::{CoarseMesh, CoarseMeshOptions, CoarsePoint, InputPoint};

/// The one place that decides how the coarse mesh is made.
fn make_coarse(
    positions: &[[f64; 3]],
    faces: &[[u32; 3]],
    points: &[InputPoint],
    target: usize,
) -> CoarseMesh {
    let mut mesh =
        CoarseMesh::from_positions(positions, faces, points, CoarseMeshOptions::default())
            .unwrap_or_else(|e| panic!("{e}"));
    greedy_coarsening(&mut mesh, target).unwrap_or_else(|e| panic!("{e}"));
    mesh
}

fn arg<T: std::str::FromStr>(args: &[String], name: &str, default: T) -> T {
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1))
        .and_then(|x| x.parse().ok())
        .unwrap_or(default)
}

fn median(v: &mut [f64]) -> f64 {
    v.sort_by(|a, b| a.total_cmp(b));
    if v.is_empty() {
        f64::NAN
    } else {
        v[v.len() / 2]
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let Some(path) = args.get(1) else {
        eprintln!(
            "usage: applications MESH.obj [--fractions a,b,c] [--sources n] [--rows n] [--voronoi n] [--exp n] [--threads n]"
        );
        std::process::exit(2);
    };
    let fractions: Vec<f64> = args
        .iter()
        .position(|a| a == "--fractions")
        .and_then(|i| args.get(i + 1))
        .map(|s| s.split(',').filter_map(|x| x.parse().ok()).collect())
        .unwrap_or_else(|| vec![0.5, 0.1, 0.02]);
    let n_sources: usize = arg(&args, "--sources", 5);
    let n_rows: usize = arg(&args, "--rows", 8);
    let n_voronoi: usize = arg(&args, "--voronoi", 20);
    let n_exp: usize = arg(&args, "--exp", 100);
    let threads: usize = arg(&args, "--threads", 8);
    let obj = read_obj(path).unwrap_or_else(|e| {
        eprintln!("{path}: {e}");
        std::process::exit(1);
    });
    let mut faces = Vec::new();
    for f in 0..obj.n_faces() {
        let c = &obj.corner_verts[obj.face_corners(f)];
        for k in 1..c.len().saturating_sub(1) {
            faces.push([c[0], c[k], c[k + 1]]);
        }
    }
    let positions = obj.positions;
    let n = positions.len();
    let mesh = TriMesh::new(n, faces.clone()).unwrap_or_else(|e| {
        eprintln!("{path}: {e}");
        std::process::exit(1);
    });
    let geo = Geometry::from_positions(&mesh, &positions).unwrap_or_else(|e| {
        eprintln!("{path}: {e}");
        std::process::exit(1);
    });
    let opts = ExactOptions::default();
    println!("# {path}: {n} vertices, {} faces", faces.len());

    // ----- fine references (computed once) -----
    let step = (n / n_sources.max(1)).max(1);
    let sources: Vec<usize> = (0..n_sources).map(|k| (k * step + 7) % n).collect();
    let mut fine_dist = Vec::new();
    let t = Instant::now();
    for &s in &sources {
        fine_dist.push(exact_distance(&mesh, &geo, &[s], &opts).unwrap().distances);
    }
    let t_fine_exact = t.elapsed().as_secs_f64() / sources.len() as f64;
    let fl = operators::cotan_laplacian(&mesh, &geo).unwrap();
    let fm = operators::lumped_mass(&mesh, &geo).unwrap();
    // Poisson data: a smooth function of the position (normalized coordinates).
    let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
    for p in &positions {
        for k in 0..3 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    let unit = |p: &[f64; 3], k: usize| (p[k] - lo[k]) / (hi[k] - lo[k]).max(1e-300);
    let f: Vec<f64> = positions
        .iter()
        .map(|p| (4.0 * unit(p, 0)).sin() + unit(p, 1) * unit(p, 2))
        .collect();
    let t = Instant::now();
    let fine_solver = PoissonSolver::new(&fl, &fm, &[]).unwrap();
    let (xf, _) = fine_solver.solve(&f, &[]).unwrap();
    let t_fine_poisson = t.elapsed().as_secs_f64();
    // Voronoi sources and all-pairs sampled rows.
    let vstep = (n / n_voronoi.max(1)).max(1);
    let vsources: Vec<usize> = (0..n_voronoi).map(|k| (k * vstep + 3) % n).collect();
    let t = Instant::now();
    let vfine: Vec<Vec<f64>> = vsources
        .iter()
        .map(|&s| exact_distance(&mesh, &geo, &[s], &opts).unwrap().distances)
        .collect();
    let t_fine_voronoi = t.elapsed().as_secs_f64();
    let rstep = (n / n_rows.max(1)).max(1);
    let rows: Vec<usize> = (0..n_rows).map(|k| (k * rstep + 11) % n).collect();
    let t = Instant::now();
    let fine_rows: Vec<Vec<f64>> = rows
        .iter()
        .map(|&s| exact_distance(&mesh, &geo, &[s], &opts).unwrap().distances)
        .collect();
    // All pairs on the fine mesh: |V| single-source runs; extrapolated from the sampled rows.
    let t_fine_allpairs = t.elapsed().as_secs_f64() / rows.len() as f64 * n as f64;
    // Exponential map on the fine mesh (the identity coarse mesh): the endpoints are tracked.
    let identity =
        CoarseMesh::from_positions(&positions, &faces, &[], CoarseMeshOptions::default()).unwrap();
    let fine_app = CoarseApplications::new(&identity).unwrap();
    let fout = fine_app.output();
    let diameter = fine_dist
        .iter()
        .flat_map(|d| d.iter().copied())
        .fold(0.0, f64::max);
    let exp_len = 0.25 * diameter;
    let exp_starts: Vec<(usize, f64)> = (0..n_exp)
        .map(|k| {
            (
                (k * 977 + 13) % n,
                ((k * 37 + 5) % 100) as f64 / 100.0 * std::f64::consts::TAU,
            )
        })
        .filter(|&(v, _)| !fout.boundary[v])
        .collect();
    let mut exp_end_points = Vec::new();
    let mut fine_exp_ok = Vec::new();
    let t = Instant::now();
    for &(v, phi) in &exp_starts {
        let angle = phi * fout.angle_sum[v] / std::f64::consts::TAU;
        let r = fine_app.exp_map_vertex(v as u32, angle, exp_len).unwrap();
        match r.point {
            CoarsePoint::Face { face, bary } if r.traced >= exp_len * 0.999 => {
                exp_end_points.push(InputPoint { face, bary });
                fine_exp_ok.push((v, phi));
            }
            _ => {}
        }
    }
    let t_fine_exp = t.elapsed().as_secs_f64() / exp_starts.len().max(1) as f64;

    println!(
        "# fine: exact distance {:.1} ms/source, Poisson {:.1} ms, Voronoi ({} sources) {:.2} s, all-pairs (extrapolated) {:.1} s, exp map {:.2} ms ({} traces)",
        1e3 * t_fine_exact,
        1e3 * t_fine_poisson,
        vsources.len(),
        t_fine_voronoi,
        t_fine_allpairs,
        1e3 * t_fine_exp,
        exp_starts.len()
    );
    println!(
        "{:>7} {:>7} | {:>8} {:>8} | {:>8} {:>8} {:>8} | {:>8} {:>8} | {:>8} {:>8} {:>8} | {:>7} {:>8} | {:>8} {:>8}",
        "frac",
        "|V̄|",
        "exact%",
        "speedup",
        "heat%",
        "speedup",
        "prefact",
        "pois%",
        "speedup",
        "apairs%",
        "speedup",
        "par.spd",
        "vor.err%",
        "speedup",
        "exp/L",
        "speedup"
    );
    let pool1 = ThreadPool::new(1).unwrap();
    let pool_n = ThreadPool::new(threads).unwrap();
    for &fraction in &fractions {
        let target = ((n as f64) * fraction).ceil() as usize;
        let coarse = make_coarse(&positions, &faces, &exp_end_points, target);
        coarse.certify(false).unwrap();
        let app = CoarseApplications::new(&coarse).unwrap();
        let nc = app.output().n_vertices();
        // Single-source geodesic, exact and heat (setup included for heat, excluded as "prefact").
        let (mut e_exact, mut e_heat) = (0.0, 0.0);
        let (mut t_exact, mut t_heat, mut t_heat_pre) = (0.0, 0.0, 0.0);
        for (k, &s) in sources.iter().enumerate() {
            let t = Instant::now();
            let (_, d) = app.geodesic_exact(s, &opts).unwrap();
            t_exact += t.elapsed().as_secs_f64();
            e_exact += DistanceError::of(&d, &fine_dist[k]).unwrap().rel_l1;
            let t = Instant::now();
            let heat = CoarseHeat::new(&app, 1.0).unwrap();
            let t2 = Instant::now();
            let sx = app.source_of(s).unwrap();
            let dbar = heat.distance(&[sx]).unwrap();
            let d = app.prolong(&dbar).unwrap();
            t_heat_pre += t2.elapsed().as_secs_f64();
            t_heat += t.elapsed().as_secs_f64();
            e_heat += DistanceError::of(&d, &fine_dist[k]).unwrap().rel_l1;
        }
        let ns = sources.len() as f64;
        // Poisson.
        let t = Instant::now();
        let fbar = app.sample(&f).unwrap();
        let (_, x, _) = app.poisson(&fbar, &[], &[]).unwrap();
        let t_poisson = t.elapsed().as_secs_f64();
        let e_poisson = DistanceError::of(&x, &xf).unwrap().rel_l2;
        // All pairs: single thread, then parallel; error on the sampled rows.
        let t = Instant::now();
        let dbar = app.all_pairs(&pool1, &opts).unwrap();
        let t_dbar1 = t.elapsed().as_secs_f64();
        let t = Instant::now();
        let dbar_n = app.all_pairs(&pool_n, &opts).unwrap();
        let t_dbar_n = t.elapsed().as_secs_f64();
        assert_eq!(dbar, dbar_n, "parallel all pairs differs from serial");
        let t = Instant::now();
        let _ = app.prolongate_rows(&dbar, &rows[..1]).unwrap();
        // The dense product costs |V| rows of this kind.
        let t_prolong_dense = t.elapsed().as_secs_f64() * n as f64;
        let part = app.prolongate_rows(&dbar, &rows).unwrap();
        let (mut es, mut rs) = (0.0, 0.0);
        for (k, r) in fine_rows.iter().enumerate() {
            for (a, b) in part.row(k).iter().zip(r) {
                es += (a - b).abs();
                rs += b;
            }
        }
        let e_pairs = es / rs;
        // Voronoi.
        let t = Instant::now();
        let vor = app.voronoi(&vsources, &pool1, &opts).unwrap();
        let t_vor = t.elapsed().as_secs_f64();
        let mut wrong = 0usize;
        for i in 0..n {
            let mut best = 0;
            for k in 1..vfine.len() {
                if vfine[k][i] < vfine[best][i] {
                    best = k;
                }
            }
            if best as u32 != vor.fine_labels[i] {
                wrong += 1;
            }
        }
        // Exponential map: coarse endpoint vs the tracked fine endpoint. Distance between the two
        // coarse points: exact distance from the traced point to the corners of the tracked
        // endpoint's face, plus the planar distance (an upper bound), over the traced length.
        let mut ratios = Vec::new();
        let t = Instant::now();
        let mut traced = Vec::new();
        for &(v, phi) in &fine_exp_ok {
            let Some(cv) = app.output().coarse_of_fine[v] else {
                continue;
            };
            let theta = app.output().angle_sum[cv as usize];
            let angle = phi * theta / std::f64::consts::TAU;
            traced.push((cv, angle, app.exp_map_vertex(cv, angle, exp_len).unwrap()));
        }
        let t_exp = t.elapsed().as_secs_f64() / traced.len().max(1) as f64;
        let first_extra = n; // the endpoints are the extra tracked points, in order
        for (this, &(v, _)) in fine_exp_ok.iter().enumerate() {
            let Some(cv) = app.output().coarse_of_fine[v] else {
                continue;
            };
            let Some(entry) = traced.iter().find(|(c, _, _)| *c == cv) else {
                continue;
            };
            let end = &entry.2;
            let target = app.output().points[first_extra + this];
            let src = match end.point {
                CoarsePoint::Vertex(u) => ExactSource::Vertex(u as usize),
                CoarsePoint::Face { face, bary } => ExactSource::Point {
                    face: face as usize,
                    bary,
                },
            };
            let dd = app.exact_distance(&[src], &opts).unwrap();
            let dist = match target {
                CoarsePoint::Vertex(u) => dd[u as usize],
                CoarsePoint::Face { face, bary } => {
                    let c = app.output().faces[face as usize];
                    let q = app.face_position(face, bary).unwrap();
                    let lay = app.face_layout(face).unwrap();
                    (0..3)
                        .map(|k| {
                            dd[c[k] as usize]
                                + ((q[0] - lay[k][0]).powi(2) + (q[1] - lay[k][1]).powi(2)).sqrt()
                        })
                        .fold(f64::INFINITY, f64::min)
                }
            };
            ratios.push(dist / exp_len);
        }
        println!(
            "{:>7.3} {:>7} | {:>8.3} {:>8.0} | {:>8.3} {:>8.0} {:>8.0} | {:>8.3} {:>8.1} | {:>8.3} {:>8.0} {:>8.0} | {:>7.2} {:>8.0} | {:>8.3} {:>8.1}",
            nc as f64 / n as f64,
            nc,
            100.0 * e_exact / ns,
            t_fine_exact / (t_exact / ns),
            100.0 * e_heat / ns,
            t_fine_exact / (t_heat / ns),
            t_fine_exact / (t_heat_pre / ns),
            100.0 * e_poisson,
            t_fine_poisson / t_poisson,
            100.0 * e_pairs,
            t_fine_allpairs / (t_dbar1 + t_prolong_dense),
            t_fine_allpairs / (t_dbar_n + t_prolong_dense),
            100.0 * wrong as f64 / n as f64,
            t_fine_voronoi / t_vor,
            median(&mut ratios),
            t_fine_exp / t_exp,
        );
    }
    println!(
        "# exact%/heat%: mean relative L1 error of the prolongated distance over the fine vertices ({} sources); heat speedup includes the factorizations, 'prefact' excludes them",
        sources.len()
    );
    println!(
        "# pois%: relative L2 error; apairs%: relative L1 error on {} sampled rows, fine time extrapolated from them, the speedup includes the dense product P D̄ Pᵀ (extrapolated from one row); par.spd with {threads} threads for D̄",
        rows.len()
    );
    println!(
        "# vor.err%: fine vertices with a different nearest source ({} sources); exp/L: median distance between the coarse endpoint and the tracked fine endpoint (an upper bound, see source) over the length L = a quarter of the largest sampled distance",
        vsources.len()
    );
}
