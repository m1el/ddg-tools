//! Times the sparse direct solvers on matrices from a real mesh.
//!
//! ```text
//! cargo run --release -p ddg-core --example factor_timings -- [mesh.obj] [repeats]
//! ```
//!
//! With no argument the largest dataset mesh (`brain.obj`) is used. Two families of systems are
//! timed, for several thread counts:
//!
//! * the cotan Laplacian `L`, made positive definite either by a diagonal shift (`L + 1e-8 I`)
//!   or by pinning one vertex (deleting its row and column), factorized with `LLᵀ`;
//! * the `2n × 2n` symmetric indefinite system `[[−M, L], [L, D]]`, with `M` the
//!   area-normalized lumped mass and `D` a diagonal with about 300 positive entries, factorized
//!   with `LDLᵀ`, `LBLᵀ` and `LU`.
//!
//! Each row reports the symbolic analysis time, the best and the first numeric factorization
//! time, the best solve time for one right-hand side, and the relative residual of that solve.

use std::time::Instant;

use ddg_core::geometry::Geometry;
use ddg_core::io::read_obj;
use ddg_core::operators::{cotan_laplacian, lumped_mass_normalized};
use ddg_core::sparse::{Method, Parallelism, SparseMatrix, Symbolic, relative_residual};
use ddg_core::topology::TriMesh;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let path = args.next().unwrap_or_else(|| {
        concat!(env!("CARGO_MANIFEST_DIR"), "/../../data/optimal-cones/optimal-cones-data/data/brain.obj").to_string()
    });
    let repeats: usize = args.next().map(|s| s.parse()).transpose()?.unwrap_or(3);

    let t = Instant::now();
    let obj = read_obj(&path)?;
    let read_time = t.elapsed().as_secs_f64();
    let faces = obj.triangles().ok_or("not a triangle mesh")?;
    let t = Instant::now();
    let mesh = TriMesh::new(obj.positions.len(), faces)?;
    let topology_time = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let geometry = Geometry::from_positions(&mesh, &obj.positions)?;
    let geometry_time = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let laplacian = cotan_laplacian(&mesh, &geometry)?;
    let mass = lumped_mass_normalized(&mesh, &geometry)?;
    let operator_time = t.elapsed().as_secs_f64();
    let n = mesh.n_vertices();

    println!("mesh: {path}");
    println!(
        "V = {n}, F = {}, E = {}, boundary loops = {}, nnz(L) = {}",
        mesh.n_faces(),
        mesh.n_edges(),
        mesh.boundary_loops().len(),
        laplacian.nnz()
    );
    println!(
        "setup: read OBJ {:.1} ms, topology {:.1} ms, geometry {:.1} ms, Laplacian + mass {:.1} ms",
        1e3 * read_time,
        1e3 * topology_time,
        1e3 * geometry_time,
        1e3 * operator_time
    );
    println!(
        "threads available to rayon: {}, numeric repeats per row: {repeats}",
        std::thread::available_parallelism().map_or(0, |n| n.get())
    );
    println!();

    // L + shift * I.
    let mut shifted = laplacian.clone();
    shifted.add_to_diagonal(&vec![1e-8; n])?;
    // L with vertex 0 pinned.
    let keep: Vec<u32> = (1..n as u32).collect();
    let pinned = laplacian.submatrix(&keep, &keep)?;
    // [[-M, L], [L, D]], D has ~300 positive entries; every diagonal entry is stored.
    let stride = (n / 300).max(1);
    let mut triplets: Vec<(u32, u32, f64)> = Vec::with_capacity(2 * laplacian.nnz() + 2 * n);
    let mut n_positive = 0;
    for i in 0..n {
        triplets.push((i as u32, i as u32, -mass[i]));
        let d = if i % stride == 0 { 1.0 } else { 0.0 };
        n_positive += usize::from(d > 0.0);
        triplets.push(((n + i) as u32, (n + i) as u32, d));
    }
    for (i, j, v) in laplacian.triplets() {
        triplets.push((n as u32 + i, j, v));
        triplets.push((i, n as u32 + j, v));
    }
    let saddle = SparseMatrix::from_triplets(2 * n, 2 * n, &triplets)?;
    println!("saddle system: {} x {}, nnz = {}, positive entries in D = {n_positive}", 2 * n, 2 * n, saddle.nnz());
    println!();

    let thread_counts = [
        Parallelism::Sequential,
        Parallelism::Threads(2),
        Parallelism::Threads(4),
        Parallelism::Threads(8),
        Parallelism::Threads(16),
        Parallelism::Threads(0),
    ];

    println!(
        "| {:<22} | {:<6} | {:>7} | {:>11} | {:>12} | {:>12} | {:>10} | {:>11} | {:>9} |",
        "system", "method", "threads", "symbolic ms", "numeric ms", "(first) ms", "solve ms", "factor nnz", "residual"
    );
    println!("|{:-<24}|{:-<8}|{:-<9}|{:-<13}|{:-<14}|{:-<14}|{:-<12}|{:-<13}|{:-<11}|", "", "", "", "", "", "", "", "", "");
    let cases: [(&str, &SparseMatrix, Method); 5] = [
        ("L + 1e-8 I", &shifted, Method::Llt),
        ("L, vertex 0 pinned", &pinned, Method::Llt),
        ("[[-M, L], [L, D]]", &saddle, Method::Ldlt),
        ("[[-M, L], [L, D]]", &saddle, Method::Lblt),
        ("[[-M, L], [L, D]]", &saddle, Method::Lu),
    ];
    for (name, matrix, method) in cases {
        let dim = matrix.nrows();
        let x_true: Vec<f64> = (0..dim).map(|i| (0.001 * i as f64).sin() + 0.5).collect();
        let mut b = vec![0.0; dim];
        matrix.mul_vec(&x_true, &mut b)?;

        let t = Instant::now();
        let symbolic = match Symbolic::analyze(matrix, method) {
            Ok(s) => s,
            Err(e) => {
                println!("| {name:<22} | {method:<6?} | symbolic analysis FAILED: {e}");
                continue;
            }
        };
        let symbolic_time = t.elapsed().as_secs_f64();
        let info = symbolic.info();
        let factor_nnz = info.factor_values.map_or("n/a".to_string(), |v| v.to_string());

        for par in thread_counts {
            let threads = match par {
                Parallelism::Sequential => "1".to_string(),
                Parallelism::Threads(0) => "all".to_string(),
                Parallelism::Threads(k) => k.to_string(),
            };
            let t = Instant::now();
            let mut factorization = match symbolic.factor(matrix, par) {
                Ok(f) => f,
                Err(e) => {
                    println!(
                        "| {name:<22} | {method:<6?} | {threads:>7} | {:>11.1} | numeric factorization FAILED after {:.1} ms: {e}",
                        1e3 * symbolic_time,
                        1e3 * t.elapsed().as_secs_f64()
                    );
                    continue;
                }
            };
            let first = t.elapsed().as_secs_f64();
            let mut best = first;
            let mut failed = None;
            for _ in 1..repeats {
                let t = Instant::now();
                if let Err(e) = factorization.refactor(matrix, par) {
                    failed = Some(e);
                    break;
                }
                best = best.min(t.elapsed().as_secs_f64());
            }
            if let Some(e) = failed {
                println!("| {name:<22} | {method:<6?} | {threads:>7} | refactorization FAILED: {e}");
                continue;
            }
            let mut solve_best = f64::INFINITY;
            let mut x = b.clone();
            let mut solve_error = None;
            for _ in 0..repeats.max(1) {
                x.copy_from_slice(&b);
                let t = Instant::now();
                if let Err(e) = factorization.solve_in_place(&mut x, par) {
                    solve_error = Some(e);
                    break;
                }
                solve_best = solve_best.min(t.elapsed().as_secs_f64());
            }
            let residual = match solve_error {
                Some(e) => format!("FAILED: {e}"),
                None => format!("{:.1e}", relative_residual(matrix, &x, &b)?),
            };
            println!(
                "| {name:<22} | {:<6} | {threads:>7} | {:>11.1} | {:>12.1} | {:>12.1} | {:>10.2} | {:>11} | {:>9} |",
                format!("{method:?}"),
                1e3 * symbolic_time,
                1e3 * best,
                1e3 * first,
                1e3 * solve_best,
                factor_nnz,
                residual
            );
        }
    }
    Ok(())
}
