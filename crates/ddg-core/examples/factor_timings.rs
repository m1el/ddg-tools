//! Times the sparse direct solvers on matrices from a real mesh.
//!
//! ```text
//! cargo run --release -p ddg-core --example factor_timings -- \
//!     [mesh.obj] [--repeats N] [--threads 1,2,4,8] [--only TEXT]
//! ```
//!
//! With no mesh argument the largest dataset mesh (`brain.obj`) is used. `--threads` lists the
//! pool sizes to try (1 means sequential, no pool). `--only` keeps the rows whose
//! "system / method" label contains the text.
//!
//! Systems (with `L` the cotan Laplacian, `M` the area-normalized lumped mass, and `D` a
//! diagonal with about 300 entries equal to 1 and zeros elsewhere):
//!
//! * `L + 1e-8 I` and `L` with one vertex pinned, by `LLᵀ`;
//! * the `2n × 2n` symmetric indefinite `K = [[−M, L], [L, D]]`, by `LDLᵀ`, `LBLᵀ` and `LU` with
//!   the default ordering, and by `LDLᵀ` with a paired elimination order that cannot meet a
//!   zero pivot;
//! * its quasi-definite regularization `[[−M, L], [L, D + 1e-8 I]]`, by `LDLᵀ`;
//! * its `n × n` positive definite Schur complement `L M⁻¹ L + D`, by `LLᵀ`.
//!
//! Each row reports the symbolic analysis time, the best and the first numeric factorization
//! time, the best solve time for one right-hand side, the relative residual of that solve, and
//! the residual after iterative refinement (at most 5 steps, target 1e-13).

use std::time::Instant;

use ddg_core::geometry::Geometry;
use ddg_core::io::read_obj;
use ddg_core::operators::{cotan_laplacian, lumped_mass_normalized};
use ddg_core::sparse::{
    EliminationOrder, Method, Parallelism, RefineOptions, SparseMatrix, Symbolic, ThreadPool,
    fill_reducing_ordering, relative_residual,
};
use ddg_core::topology::TriMesh;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../data/optimal-cones/optimal-cones-data/data/brain.obj"
    )
    .to_string();
    let mut repeats = 3usize;
    let mut threads: Vec<usize> = vec![1, 2, 4, 8, 16, 32];
    let mut only: Option<String> = None;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--repeats" => repeats = args.next().ok_or("--repeats needs a value")?.parse()?,
            "--threads" => {
                threads = args
                    .next()
                    .ok_or("--threads needs a value")?
                    .split(',')
                    .map(str::parse)
                    .collect::<Result<_, _>>()?
            }
            "--only" => only = Some(args.next().ok_or("--only needs a value")?),
            _ => path = arg,
        }
    }

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
        "logical CPUs: {}, numeric repeats per row: {repeats}",
        std::thread::available_parallelism().map_or(0, |n| n.get())
    );

    // L + shift * I.
    let mut shifted = laplacian.clone();
    shifted.add_to_diagonal(&vec![1e-8; n])?;
    // L with vertex 0 pinned.
    let keep: Vec<u32> = (1..n as u32).collect();
    let pinned = laplacian.submatrix(&keep, &keep)?;
    // K = [[-M, L], [L, D]]; every diagonal entry of D is stored.
    let stride = (n / 300).max(1);
    let d: Vec<f64> = (0..n)
        .map(|i| if i % stride == 0 { 1.0 } else { 0.0 })
        .collect();
    let n_positive = d.iter().filter(|&&x| x > 0.0).count();
    let mut triplets: Vec<(u32, u32, f64)> = Vec::with_capacity(2 * laplacian.nnz() + 2 * n);
    for i in 0..n {
        triplets.push((i as u32, i as u32, -mass[i]));
        triplets.push(((n + i) as u32, (n + i) as u32, d[i]));
    }
    for (i, j, v) in laplacian.triplets() {
        triplets.push((n as u32 + i, j, v));
        triplets.push((i, n as u32 + j, v));
    }
    let saddle = SparseMatrix::from_triplets(2 * n, 2 * n, &triplets)?;
    // Quasi-definite regularization: D + delta * I.
    let mut quasi_definite = saddle.clone();
    let mut delta = vec![0.0; 2 * n];
    delta[n..].fill(1e-8);
    quasi_definite.add_to_diagonal(&delta)?;
    // Schur complement S = L M^-1 L + D.
    let t = Instant::now();
    let mut triplets: Vec<(u32, u32, f64)> = Vec::new();
    for (k, &m) in mass.iter().enumerate() {
        let column = laplacian.col_ptr()[k]..laplacian.col_ptr()[k + 1];
        for a in column.clone() {
            for b in column.clone() {
                triplets.push((
                    laplacian.row_idx()[a],
                    laplacian.row_idx()[b],
                    laplacian.values()[a] * laplacian.values()[b] / m,
                ));
            }
        }
    }
    for (i, &x) in d.iter().enumerate() {
        triplets.push((i as u32, i as u32, x));
    }
    let schur = SparseMatrix::from_triplets(n, n, &triplets)?;
    drop(triplets);
    let schur_time = t.elapsed().as_secs_f64();
    println!(
        "K = [[-M, L], [L, D]]: {} x {}, nnz = {}, positive entries in D = {n_positive}",
        2 * n,
        2 * n,
        saddle.nnz()
    );
    println!(
        "Schur complement L M^-1 L + D: {n} x {n}, nnz = {}, assembled in {:.1} ms",
        schur.nnz(),
        1e3 * schur_time
    );
    println!();

    // Elimination order for K that pairs unknown i of the first block with unknown n + i of
    // the second, the pairs following the fill-reducing order of L.
    let t = Instant::now();
    let paired: Vec<u32> = fill_reducing_ordering(&laplacian)?
        .iter()
        .flat_map(|&i| [i, n as u32 + i])
        .collect();
    let paired_ms = 1e3 * t.elapsed().as_secs_f64();
    println!("paired elimination order for K (AMD of L, then interleave): {paired_ms:.1} ms");
    println!();

    let pools: Vec<(usize, Option<ThreadPool>)> = threads
        .iter()
        .map(|&k| {
            Ok((
                k,
                if k <= 1 {
                    None
                } else {
                    Some(ThreadPool::new(k)?)
                },
            ))
        })
        .collect::<Result<_, ddg_core::sparse::SparseError>>()?;

    println!(
        "| {:<26} | {:<6} | {:>7} | {:>11} | {:>10} | {:>10} | {:>8} | {:>10} | {:<9} |",
        "system",
        "method",
        "threads",
        "symbolic ms",
        "numeric ms",
        "(first) ms",
        "solve ms",
        "factor nnz",
        "residual | after refinement"
    );
    println!(
        "|{:-<28}|{:-<8}|{:-<9}|{:-<13}|{:-<12}|{:-<12}|{:-<10}|{:-<12}|{:-<11}|",
        "", "", "", "", "", "", "", "", ""
    );
    let default = EliminationOrder::FillReducing;
    let cases: [(&str, &SparseMatrix, Method, EliminationOrder<'_>); 8] = [
        ("L + 1e-8 I", &shifted, Method::Llt, default),
        ("L, vertex 0 pinned", &pinned, Method::Llt, default),
        ("K = [[-M, L], [L, D]]", &saddle, Method::Ldlt, default),
        (
            "K, paired order",
            &saddle,
            Method::Ldlt,
            EliminationOrder::Custom(&paired),
        ),
        ("K = [[-M, L], [L, D]]", &saddle, Method::Lblt, default),
        ("K = [[-M, L], [L, D]]", &saddle, Method::Lu, default),
        (
            "K + diag(0, 1e-8 I)",
            &quasi_definite,
            Method::Ldlt,
            default,
        ),
        ("Schur L M^-1 L + D", &schur, Method::Llt, default),
    ];
    for (name, matrix, method, ordering) in cases {
        let method_name = format!("{method:?}");
        if let Some(filter) = &only
            && !format!("{name} / {method_name}").contains(filter.as_str())
        {
            continue;
        }
        let dim = matrix.nrows();
        let x_true: Vec<f64> = (0..dim).map(|i| (0.001 * i as f64).sin() + 0.5).collect();
        let mut b = vec![0.0; dim];
        matrix.mul_vec(&x_true, &mut b)?;

        let t = Instant::now();
        let symbolic = match Symbolic::analyze_with_order(matrix, method, ordering) {
            Ok(s) => s,
            Err(e) => {
                println!("| {name:<26} | {method_name:<6} | symbolic analysis FAILED: {e}");
                continue;
            }
        };
        let symbolic_ms = 1e3 * t.elapsed().as_secs_f64();
        let factor_nnz = symbolic
            .info()
            .factor_values
            .map_or("n/a".to_string(), |v| v.to_string());

        for (k, pool) in &pools {
            let par = match pool {
                None => Parallelism::Sequential,
                Some(pool) => Parallelism::Pool(pool),
            };
            let prefix =
                format!("| {name:<26} | {method_name:<6} | {k:>7} | {symbolic_ms:>11.1} |");
            let t = Instant::now();
            let mut factorization = match symbolic.factor(matrix, par) {
                Ok(f) => f,
                Err(e) => {
                    println!(
                        "{prefix} numeric FAILED after {:.1} ms: {e}",
                        1e3 * t.elapsed().as_secs_f64()
                    );
                    continue;
                }
            };
            let first = t.elapsed().as_secs_f64();
            let mut best = first;
            let mut failed = None;
            // Do not repeat runs that take longer than a few seconds.
            let numeric_repeats = if first > 5.0 { 1 } else { repeats };
            for _ in 1..numeric_repeats {
                let t = Instant::now();
                if let Err(e) = factorization.refactor(matrix, par) {
                    failed = Some(e);
                    break;
                }
                best = best.min(t.elapsed().as_secs_f64());
            }
            if let Some(e) = failed {
                println!("{prefix} refactorization FAILED: {e}");
                continue;
            }
            let mut solve_best = f64::INFINITY;
            let mut x = b.clone();
            let mut solve_error = None;
            for _ in 0..repeats.max(1) {
                x.copy_from_slice(&b);
                let t = Instant::now();
                let result = factorization.solve_in_place(&mut x, par);
                solve_best = solve_best.min(t.elapsed().as_secs_f64());
                if let Err(e) = result {
                    solve_error = Some(e);
                    break;
                }
            }
            let residual = match solve_error {
                Some(e) => format!("solve FAILED: {e}"),
                None => {
                    let plain = relative_residual(matrix, &x, &b)?;
                    // Iterative refinement: how many extra solves until 1e-13?
                    let options = RefineOptions {
                        max_steps: 5,
                        tolerance: 1e-13,
                    };
                    let report = factorization.solve_refined(matrix, &b, &mut x, options, par)?;
                    format!(
                        "{plain:.1e} | {:.1e} after {} refinement steps",
                        report.residual, report.steps
                    )
                }
            };
            println!(
                "{prefix} {:>10.1} | {:>10.1} | {:>8.2} | {:>10} | {:<9} |",
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
