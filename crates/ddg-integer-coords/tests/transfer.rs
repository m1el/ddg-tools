//! Attribute transfer (WP7; SPEC §6, §8; paper Sec. 4.3, Eq. 11-14).
//!
//! States: WP1's converted signpost states (`tests/common`) of the standard shapes, spot and
//! camelbad, in the 7 scenarios (seed 11). Run with `--nocapture` to see the tables.

#![allow(clippy::needless_range_loop)]

mod common;

use std::sync::OnceLock;

use common::{SCENARIOS, Scenario, Shape, all_shapes, convert, state};
use ddg_integer_coords::{IntegerTriangulation, Transfer};

fn shapes() -> &'static [Shape] {
    static SHAPES: OnceLock<Vec<Shape>> = OnceLock::new();
    SHAPES.get_or_init(all_shapes)
}

fn par_map<T: Sync, R: Send>(items: &[T], f: impl Fn(&T) -> R + Sync) -> Vec<R> {
    let threads = std::thread::available_parallelism()
        .map_or(4, |n| n.get())
        .min(16);
    let next = std::sync::atomic::AtomicUsize::new(0);
    let out: Vec<std::sync::Mutex<Option<R>>> =
        items.iter().map(|_| std::sync::Mutex::new(None)).collect();
    std::thread::scope(|scope| {
        for _ in 0..threads {
            scope.spawn(|| {
                loop {
                    let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    if i >= items.len() {
                        break;
                    }
                    *out[i].lock().unwrap() = Some(f(&items[i]));
                }
            });
        }
    });
    out.into_iter()
        .map(|m| m.into_inner().unwrap().unwrap())
        .collect()
}

/// The Galerkin mass matrix of the input mesh from its own geometry: `A/12 (1 + δ)` per face,
/// as dense rows (maps), independent of the common subdivision.
fn input_mass(shape: &Shape) -> std::collections::HashMap<(u32, u32), f64> {
    let areas = shape.geometry.face_areas();
    let mut m = std::collections::HashMap::new();
    for (f, face) in shape.mesh.faces().iter().enumerate() {
        for a in 0..3 {
            for b in 0..3 {
                let w = areas[f] / 12.0 * if a == b { 2.0 } else { 1.0 };
                *m.entry((face[a], face[b])).or_insert(0.0) += w;
            }
        }
    }
    m
}

fn smooth(p: [f64; 3], scale: f64) -> f64 {
    let k = 1.0 / scale;
    (k * p[0]).sin() + 0.5 * (1.3 * k * p[1] + 0.4).cos() + 0.25 * (0.7 * k * p[2]).sin()
}

fn affine(p: [f64; 3]) -> f64 {
    0.7 + 1.3 * p[0] - 0.4 * p[1] + 2.1 * p[2]
}

#[derive(Debug, Default)]
struct Row {
    name: String,
    n0: usize,
    n1: usize,
    ns: usize,
    mass_err: f64,
    area_err: f64,
    const_err: f64,
    lin_err: f64,
    residual: f64,
    /// transfer / copy distances (to P1 f), T1 -> T0, and T0 -> T1.
    to_input: (f64, f64),
    to_intrinsic: (f64, f64),
    /// errors against the exact function, transfer / copy (T1 -> T0).
    exact: (f64, f64),
    stale_refused: bool,
    failure: Option<String>,
}

fn run(shape: &Shape, scenario: Scenario, tri: &IntegerTriangulation) -> Result<Row, String> {
    let e = |x: ddg_integer_coords::IntegerError| x.to_string();
    let sub = tri.common_subdivision().map_err(e)?;
    let tr = Transfer::new(tri, &sub).map_err(e)?;
    let (n0, n1, ns) = (
        tr.n_input_vertices(),
        tr.n_intrinsic_vertices(),
        tr.mass_subdivision().nrows(),
    );
    let mut row = Row {
        name: format!("{} / {scenario:?}", shape.name),
        n0,
        n1,
        ns,
        ..Row::default()
    };
    // Mass identity.
    let reference = input_mass(shape);
    let scale = shape.geometry.total_area();
    row.area_err = ((tr.total_area() - scale) / scale).abs();
    let diag = tr.mass_input().diagonal();
    let mut worst = 0.0f64;
    let mut seen = std::collections::HashSet::new();
    for (i, j, v) in tr.mass_input().triplets() {
        let r = reference.get(&(i, j)).copied().unwrap_or(0.0);
        seen.insert((i, j));
        let s = diag[i as usize].max(diag[j as usize]);
        worst = worst.max((v - r).abs() / s);
    }
    for (&(i, j), &r) in &reference {
        if !seen.contains(&(i, j)) {
            worst = worst.max(r.abs() / diag[i as usize].max(diag[j as usize]));
        }
    }
    row.mass_err = worst;
    // Constants.
    let ones1 = vec![1.0; n1];
    let ones0 = vec![1.0; n0];
    let a = tr.to_input(&ones1).map_err(e)?;
    let b = tr.to_intrinsic(&ones0).map_err(e)?;
    let c0 = a.iter().map(|x| (x - 1.0).abs()).fold(0.0, f64::max);
    let c1 = b.iter().map(|x| (x - 1.0).abs()).fold(0.0, f64::max);
    row.const_err = c0.max(c1);
    // Positions on S; V1 positions are the first n1.
    let pos = tr.subdivision_positions(&shape.positions).map_err(e)?;
    // Linear functions (flat shapes).
    if shape.flat {
        let f: Vec<f64> = pos[..n1].iter().map(|&p| affine(p)).collect();
        let g: Vec<f64> = shape.positions.iter().map(|&p| affine(p)).collect();
        let fh = tr.to_input(&f).map_err(e)?;
        let gh = tr.to_intrinsic(&g).map_err(e)?;
        let mut w = 0.0f64;
        for v in 0..n0 {
            w = w.max((fh[v] - g[v]).abs());
        }
        for v in 0..n1 {
            w = w.max((gh[v] - f[v]).abs());
        }
        row.lin_err = w;
    }
    // Smooth function: transfer vs copy.
    let h = shape.mean_edge_length();
    let size = (shape.geometry.total_area()).sqrt().max(h);
    let u = |p| smooth(p, size / 2.0);
    let f: Vec<f64> = pos[..n1].iter().map(|&p| u(p)).collect();
    let g: Vec<f64> = shape.positions.iter().map(|&p| u(p)).collect();
    let dist =
        |xs: &[f64], from: &ddg_core::SparseMatrix, y: &[f64], to: &ddg_core::SparseMatrix| {
            let mut d = vec![0.0; ns];
            let mut t = vec![0.0; ns];
            from.mul_vec(xs, &mut d).unwrap();
            to.mul_vec(y, &mut t).unwrap();
            for (a, b) in d.iter_mut().zip(&t) {
                *a -= b;
            }
            tr.mass_norm(&d).unwrap()
        };
    let (fh, rep) = tr.to_input_report(&f).map_err(e)?;
    let (gh, rep2) = tr.to_intrinsic_report(&g).map_err(e)?;
    row.residual = rep.residual.max(rep2.residual);
    let fc = tr.copy_to_input(&f).map_err(e)?;
    let gc = tr.copy_to_intrinsic(&g).map_err(e)?;
    row.to_input = (
        dist(&f, tr.p1(), &fh, tr.p0()),
        dist(&f, tr.p1(), &fc, tr.p0()),
    );
    row.to_intrinsic = (
        dist(&g, tr.p0(), &gh, tr.p1()),
        dist(&g, tr.p0(), &gc, tr.p1()),
    );
    if (rep.l2_distance - row.to_input.0).abs() > 1e-9 * (1.0 + row.to_input.0) {
        return Err(format!(
            "report l2 {} vs {}",
            rep.l2_distance, row.to_input.0
        ));
    }
    // Against the exact function, on T0.
    let mut on_s = vec![0.0; ns];
    tr.p0().mul_vec(&fh, &mut on_s).unwrap();
    let et = tr.l2_error_to_function(&pos, &on_s, u).map_err(e)?;
    tr.p0().mul_vec(&fc, &mut on_s).unwrap();
    let ec = tr.l2_error_to_function(&pos, &on_s, u).map_err(e)?;
    row.exact = (et, ec);
    // Stale subdivision.
    let mut t2 = tri.clone();
    let mut refused = false;
    for edge in 0..t2.complex().n_edges() {
        if t2.flip_edge(edge).is_ok() {
            refused = Transfer::new(&t2, &sub).is_err() && !tr.is_current_for(&t2);
            break;
        }
    }
    row.stale_refused = refused || t2.fingerprint() == tri.fingerprint();
    Ok(row)
}

#[test]
fn transfer_on_all_states() {
    let shapes = shapes();
    let mut specs = Vec::new();
    for s in 0..shapes.len() {
        for sc in SCENARIOS {
            specs.push((s, sc));
        }
    }
    let rows = par_map(&specs, |&(s, sc)| {
        let t = state(&shapes[s], sc, 11);
        let (tri, _) = convert(&t).ok()?;
        Some(match run(&shapes[s], sc, &tri) {
            Ok(r) => r,
            Err(m) => Row {
                name: format!("{} / {sc:?}", shapes[s].name),
                failure: Some(m),
                ..Row::default()
            },
        })
    });
    let rows: Vec<Row> = rows.into_iter().flatten().collect();
    println!(
        "{:<44} {:>6} {:>6} {:>7} {:>9} {:>9} {:>9} {:>9} {:>9}  {:>10} {:>10} {:>10} {:>10}",
        "state",
        "V0",
        "V1",
        "VS",
        "M0 err",
        "area",
        "const",
        "linear",
        "resid",
        "L2 T1>T0",
        "copy",
        "L2 T0>T1",
        "copy"
    );
    let mut worst_mass = 0.0f64;
    let mut worst_resid = 0.0f64;
    let mut fails = Vec::new();
    for r in &rows {
        if let Some(m) = &r.failure {
            fails.push(format!("{}: {m}", r.name));
            continue;
        }
        println!(
            "{:<44} {:>6} {:>6} {:>7} {:>9.2e} {:>9.2e} {:>9.2e} {:>9.2e} {:>9.2e}  {:>10.3e} {:>10.3e} {:>10.3e} {:>10.3e}",
            r.name,
            r.n0,
            r.n1,
            r.ns,
            r.mass_err,
            r.area_err,
            r.const_err,
            r.lin_err,
            r.residual,
            r.to_input.0,
            r.to_input.1,
            r.to_intrinsic.0,
            r.to_intrinsic.1
        );
        worst_mass = worst_mass.max(r.mass_err);
        worst_resid = worst_resid.max(r.residual);
        // Pieces tile the surface.
        if r.area_err > 1e-9 {
            fails.push(format!("{}: area {:e}", r.name, r.area_err));
        }
        if r.mass_err > 1e-8 {
            fails.push(format!("{}: mass identity {:e}", r.name, r.mass_err));
        }
        if r.const_err > 1e-8 {
            fails.push(format!("{}: constants {:e}", r.name, r.const_err));
        }
        if r.lin_err > 1e-7 {
            fails.push(format!("{}: linear {:e}", r.name, r.lin_err));
        }
        if r.residual > 1e-10 {
            fails.push(format!("{}: residual {:e}", r.name, r.residual));
        }
        if r.to_input.0 > r.to_input.1 * (1.0 + 1e-9) + 1e-12 {
            fails.push(format!("{}: transfer worse than copy T1>T0", r.name));
        }
        if r.to_intrinsic.0 > r.to_intrinsic.1 * (1.0 + 1e-9) + 1e-12 {
            fails.push(format!("{}: transfer worse than copy T0>T1", r.name));
        }
        if !r.stale_refused {
            fails.push(format!("{}: stale subdivision accepted", r.name));
        }
    }
    println!("exact-function errors (transfer / copy), T1 > T0:");
    for r in &rows {
        if r.failure.is_none() {
            println!("  {:<44} {:>10.3e} {:>10.3e}", r.name, r.exact.0, r.exact.1);
        }
    }
    println!(
        "worst relative M0 entry error {worst_mass:.3e}; worst normal-equation residual {worst_resid:.3e}; {} states",
        rows.len()
    );
    assert!(!rows.is_empty());
    assert!(fails.is_empty(), "{}", fails.join("\n"));
}

#[test]
fn bad_input_is_an_error() {
    let shape = &shapes()[0];
    let t = state(shape, Scenario::Delaunay, 11);
    let (tri, _) = convert(&t).unwrap();
    let tr = Transfer::from_state(&tri).unwrap();
    let (n0, n1) = (tr.n_input_vertices(), tr.n_intrinsic_vertices());
    assert!(tr.to_input(&vec![1.0; n1 + 1]).is_err());
    assert!(tr.to_input(&[]).is_err());
    assert!(tr.to_intrinsic(&vec![1.0; n0.saturating_sub(1)]).is_err());
    let mut f = vec![1.0; n1];
    f[0] = f64::NAN;
    assert!(tr.to_input(&f).is_err());
    f[0] = f64::INFINITY;
    assert!(tr.copy_to_input(&f).is_err());
    let mut g = vec![0.0; n0];
    g[n0 - 1] = f64::NAN;
    assert!(tr.to_intrinsic(&g).is_err());
    assert!(tr.copy_to_intrinsic(&g).is_err());
    assert!(tr.mass_norm(&[1.0]).is_err());
    assert!(tr.subdivision_positions(&[[0.0; 3]]).is_err());
    assert!(tr.l2_error_to_function(&[], &[], |_| 0.0).is_err());
    // Zero in, zero out.
    let z = tr.to_input(&vec![0.0; n1]).unwrap();
    assert!(z.iter().all(|&x| x == 0.0));
}
