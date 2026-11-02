//! The pieces of one hierarchical solve (`SPEC.md` §8.3, §6.3; WP17): the block cluster tree
//! build of `B`, `B⁰`, one product with 1–4 vectors, the vertex tree of `L^V_{2−s}`, one
//! preconditioner application with 1–2 vectors. A plain `main` printing one CSV row per mesh
//! and thread count (`--header` prints the column names).
//!
//! ```text
//! cargo bench -p ddg-repulsive-surfaces --bench product -- [--header] [--mesh SPEC]
//!     [--threads 1,4,16] [--reps N]
//! ```
//!
//! `SPEC`: `torus:FACES` (default `torus:45000`) or an OBJ path. Times are the median of `reps`.

use std::sync::Arc;
use std::time::Instant;

use ddg_core::generators::{perturb, torus};
use ddg_repulsive_surfaces::constraints::{BarycenterMode, Constraint, ConstraintSet};
use ddg_repulsive_surfaces::energy::exact::make_pool;
use ddg_repulsive_surfaces::geometry::mean_edge_length;
use ddg_repulsive_surfaces::metric::hierarchical::{ClusterTreeOptions, FaceMetricTree, VertexFractionalTree};
use ddg_repulsive_surfaces::solve::LinearOperator;
use ddg_repulsive_surfaces::solve::Preconditioner;
use ddg_repulsive_surfaces::solve::iterative::{HierarchicalOptions, MetricExtras, hierarchical_saddle};
use ddg_repulsive_surfaces::{Exponents, SobolevWeights, Surface, read_obj};

fn load(spec: &str) -> (Surface, Vec<[f64; 3]>) {
    if let Some(f) = spec.strip_prefix("torus:") {
        let f: f64 = f.parse().expect("faces");
        let nm = (f / 5.0).sqrt().round() as usize;
        let nmaj = (f / (2.0 * nm as f64)).round() as usize;
        let m = torus(nmaj, nm, 1.5, 0.6).expect("torus");
        let s = Surface::new(m.positions.len(), m.faces).expect("surface");
        let mut p = m.positions;
        let l = mean_edge_length(&s, &p);
        perturb(&mut p, [0.05 * l; 3], 7);
        return (s, p);
    }
    read_obj(spec).unwrap_or_else(|e| panic!("{spec}: {e}"))
}

fn median_ms(reps: usize, mut f: impl FnMut()) -> f64 {
    let mut t: Vec<f64> = (0..reps)
        .map(|_| {
            let t0 = Instant::now();
            f();
            t0.elapsed().as_secs_f64() * 1e3
        })
        .collect();
    t.sort_by(f64::total_cmp);
    t[t.len() / 2]
}

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let mut mesh = "torus:45000".to_string();
    let mut threads = vec![1usize, 4, 16];
    let mut reps = 5;
    let mut i = 0;
    while i < argv.len() {
        match argv[i].as_str() {
            "--header" => println!(
                "mesh,faces,threads,bct_build_ms,near_entries,product1_ms,product2_ms,product3_ms,product4_ms,\
                 vertex_tree_ms,precond_build_ms,precond1_ms,precond2_ms"
            ),
            "--bench" => {}
            "--mesh" => {
                mesh = argv[i + 1].clone();
                i += 1;
            }
            "--threads" => {
                threads = argv[i + 1].split(',').map(|t| t.parse().expect("threads")).collect();
                i += 1;
            }
            "--reps" => {
                reps = argv[i + 1].parse().expect("reps");
                i += 1;
            }
            other => panic!("unknown argument {other}"),
        }
        i += 1;
    }
    let (s, p) = load(&mesh);
    let ex = Exponents::default();
    let nd = 3 * s.n_vertices();
    let xs: Vec<Vec<f64>> = (0..4)
        .map(|q| (0..nd).map(|k| ((k * (q + 3)) as f64 * 0.37).sin()).collect())
        .collect();
    let set = ConstraintSet::new(
        &s,
        &p,
        &[Constraint::Barycenter { component: None, target: None }],
        BarycenterMode::Global,
    )
    .expect("constraints");
    let cq = set.jacobian(&p, &set.positional_rows()).expect("jacobian");
    for &t in &threads {
        let pool = Arc::new(make_pool(Some(t)).expect("pool"));
        let opts = ClusterTreeOptions::default();
        let mut tree = None;
        let build = median_ms(reps.min(3), || {
            tree = Some(FaceMetricTree::build(&pool, &s, &p, &ex, &opts).expect("tree"));
        });
        let tree = tree.expect("built");
        let mut prod = [0.0; 4];
        for (k, pk) in prod.iter_mut().enumerate() {
            let xr: Vec<&[f64]> = xs[..=k].iter().map(Vec::as_slice).collect();
            let mut ys: Vec<Vec<f64>> = vec![vec![0.0; nd]; k + 1];
            *pk = median_ms(reps, || {
                let mut yr: Vec<&mut [f64]> = ys.iter_mut().map(Vec::as_mut_slice).collect();
                tree.apply_batch(&pool, 1.0, 1.0, None, &xr, &mut yr).expect("product");
            });
        }
        let vt = median_ms(reps.min(3), || {
            VertexFractionalTree::for_surface(&pool, &s, &p, 2.0 - ex.s(), &opts).expect("vertex tree");
        });
        let mut saddle = None;
        let pb = median_ms(1, || {
            saddle = Some(
                hierarchical_saddle(
                    pool.clone(),
                    &s,
                    &p,
                    &ex,
                    &SobolevWeights::default(),
                    &MetricExtras::default(),
                    &cq,
                    &HierarchicalOptions::default(),
                    None,
                )
                .expect("saddle"),
            );
        });
        let saddle = saddle.expect("saddle");
        let n = saddle.operator().dim();
        let rs: Vec<Vec<f64>> = (0..2).map(|q| (0..n).map(|k| ((k * (q + 5)) as f64 * 0.21).cos()).collect()).collect();
        let mut pre = [0.0; 2];
        for (k, pk) in pre.iter_mut().enumerate() {
            let rr: Vec<&[f64]> = rs[..=k].iter().map(Vec::as_slice).collect();
            let mut zs: Vec<Vec<f64>> = vec![vec![0.0; n]; k + 1];
            *pk = median_ms(reps, || {
                let mut zr: Vec<&mut [f64]> = zs.iter_mut().map(Vec::as_mut_slice).collect();
                saddle.preconditioner().apply_batch(&rr, &mut zr).expect("preconditioner");
            });
        }
        let _ = saddle.operator().apply(&rs[0], &mut vec![0.0; n]);
        println!(
            "{mesh},{},{t},{build:.1},{},{:.2},{:.2},{:.2},{:.2},{vt:.1},{pb:.1},{:.2},{:.2}",
            s.n_faces(),
            tree.stats().near_entries,
            prod[0],
            prod[1],
            prod[2],
            prod[3],
            pre[0],
            pre[1]
        );
    }
}
