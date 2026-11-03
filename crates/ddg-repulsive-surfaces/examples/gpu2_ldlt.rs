//! GPU-2 measurements of the preconditioner's `L′` factor (`NOTES-GPU2.md`): structure of the
//! simplicial `LDLᵀ` (elimination tree height, chains, long rows), factor and solve times of the
//! simplicial factor against `ddg_core`'s (faer), and their difference.
//!
//! ```text
//! gpu2_ldlt <mesh> [--pins K]       <mesh>: tunnel8:FACES | fig12:I | obj:PATH
//! ```

use std::sync::Arc;
use std::time::Instant;

use ddg_core::SparseMatrix;
use ddg_core::sparse::{EliminationOrder, Method, Parallelism, Symbolic, fill_reducing_ordering};
use ddg_repulsive_surfaces::constraints::{BarycenterMode, Constraint, ConstraintSet};
use ddg_repulsive_surfaces::energy::exact::make_pool;
use ddg_repulsive_surfaces::generate::shapes::{Knot, fig12, knotted_tunnel};
use ddg_repulsive_surfaces::solve::ldlt::SimplicialLdlt;
use ddg_repulsive_surfaces::solve::precond::{PreconditionerOptions, SparsePreconditioner};
use ddg_repulsive_surfaces::{Exponents, Surface, read_obj};

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let spec = args.first().cloned().unwrap_or_else(|| "fig12:3".into());
    let pins: usize = args
        .iter()
        .position(|a| a == "--pins")
        .map(|i| args[i + 1].parse().unwrap())
        .unwrap_or(0);
    let parts: Vec<&str> = spec.split(':').collect();
    let (s, p) = match parts[0] {
        "tunnel8" => {
            let m = knotted_tunnel(Knot::FigureEight, parts[1].parse().unwrap()).unwrap().mesh;
            (Surface::new(m.positions.len(), m.faces).unwrap(), m.positions)
        }
        "fig12" => {
            let m = fig12(parts[1].parse().unwrap()).unwrap().mesh;
            (Surface::new(m.positions.len(), m.faces).unwrap(), m.positions)
        }
        _ => read_obj(parts[1..].join(":")).unwrap(),
    };
    let pool = Arc::new(make_pool(Some(4)).unwrap());
    let mut cons = vec![Constraint::Barycenter { component: None, target: None }];
    if pins > 0 {
        let step = s.n_vertices() / pins;
        cons.push(Constraint::Pins {
            vertices: (0..pins).map(|k| (k * step) as u32).collect(),
            targets: None,
        });
    }
    let set = ConstraintSet::new(&s, &p, &cons, BarycenterMode::PerComponent).unwrap();
    let cq = set.jacobian(&p, &set.positional_rows()).unwrap();
    let prec = SparsePreconditioner::new(
        pool,
        &s,
        &p,
        &Exponents::default(),
        &cq,
        &PreconditionerOptions::default(),
        None,
        None,
    )
    .unwrap();
    let ls = prec.laplacian_saddle();
    let a: &SparseMatrix = ls.matrix();
    let n = a.nrows();
    let nv = s.n_vertices();
    println!("mesh {spec}: {} faces, {nv} vertices, saddle n = {n}, per-axis {}", s.n_faces(), ls.is_per_axis());
    // The order of LaplacianSaddle: AMD of the Laplacian block, then the multipliers.
    let kb: Vec<u32> = (0..nv as u32).collect();
    let k = a.submatrix(&kb, &kb).unwrap();
    let t0 = Instant::now();
    let amd = fill_reducing_ordering(&k).unwrap();
    let t_amd = t0.elapsed().as_secs_f64();
    let order: Vec<u32> = amd.iter().copied().chain(nv as u32..n as u32).collect();
    let t0 = Instant::now();
    let f = SimplicialLdlt::factor(a, &order).unwrap();
    let t_own = t0.elapsed().as_secs_f64();
    let t0 = Instant::now();
    let sym = Symbolic::analyze_with_order(a, Method::Ldlt, EliminationOrder::Custom(&order)).unwrap();
    let fa = sym.factor(a, Parallelism::Sequential).unwrap();
    let t_faer = t0.elapsed().as_secs_f64();
    println!(
        "amd {:.1} ms, own factor {:.1} ms, faer analyze+factor {:.1} ms; faer {:?}",
        1e3 * t_amd,
        1e3 * t_own,
        1e3 * t_faer,
        sym.info()
    );
    println!("{:?}", f.stats());
    // The GPU schedule.
    let sched = f.schedule().unwrap();
    let in_sub: usize = sched.sub_start.iter().zip(&sched.sub_end).map(|(a, b)| (b - a) as usize).sum();
    println!(
        "schedule: {} leaf subtrees ({} nodes), {} chains, {} forward / {} backward levels, slots {:?}",
        sched.sub_start.len(),
        in_sub,
        sched.chain_start.len(),
        sched.fwd.level_ptr.len() - 1,
        sched.bwd.level_ptr.len() - 1,
        sched.n_slots
    );
    let sw = &sched.fwd;
    let mut crit = 0;
    for l in 0..sw.level_ptr.len() - 1 {
        let ch = &sw.chains[sw.level_ptr[l]..sw.level_ptr[l + 1]];
        let longest = ch.iter().map(|&c| (sched.chain_end[c as usize] - sched.chain_start[c as usize]) as usize).max().unwrap_or(0);
        crit += longest;
        println!("  forward level {l}: {} chains, longest {longest}, block {}, slots {}", ch.len(), sw.level_block[l], sw.slot_ptr[l + 1] - sw.slot_ptr[l]);
    }
    println!("  sum of the longest chain per level: {crit}");
    // Solve times, 3 columns (one vector), own vs faer, and the difference.
    let m = 3;
    let b: Vec<f64> = (0..n * m).map(|i| ((i * 7919) % 1000) as f64 / 500.0 - 1.0).collect();
    let mut xo = b.clone();
    let reps = 20;
    let t0 = Instant::now();
    for _ in 0..reps {
        xo.copy_from_slice(&b);
        f.solve_rows(&mut xo, m).unwrap();
    }
    let t_os = t0.elapsed().as_secs_f64() / reps as f64;
    // faer: column-major, original numbering.
    let mut bc = vec![0.0; n * m];
    for (new, &old) in f.perm().iter().enumerate() {
        for c in 0..m {
            bc[c * n + old as usize] = b[new * m + c];
        }
    }
    let mut xf = bc.clone();
    let t0 = Instant::now();
    for _ in 0..reps {
        xf.copy_from_slice(&bc);
        fa.solve_many_in_place(&mut xf, m, Parallelism::Sequential).unwrap();
    }
    let t_fs = t0.elapsed().as_secs_f64() / reps as f64;
    let mut maxd: f64 = 0.0;
    let mut maxx: f64 = 0.0;
    for (new, &old) in f.perm().iter().enumerate() {
        for c in 0..m {
            maxd = maxd.max((xo[new * m + c] - xf[c * n + old as usize]).abs());
            maxx = maxx.max(xf[c * n + old as usize].abs());
        }
    }
    println!(
        "solve (3 columns): own {:.2} ms, faer {:.2} ms; max |own - faer| / max |x| = {:.2e}",
        1e3 * t_os,
        1e3 * t_fs,
        maxd / maxx
    );
}
