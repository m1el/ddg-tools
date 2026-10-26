//! Golden fixtures from the authors' program (`SPEC.md` §13.7, §14.7, WP15).
//!
//! The fixtures in `tests/fixtures/reference/*.fx` were written by `rcurves_dump`, a driver that
//! links the reference (Yu, Schumacher, Crane; MIT; commit pinned in `SPEC.md` §17c, batch-mode
//! patch plus `dump_fixtures.patch`, both next to the build in
//! `/var/data/discrete-geometry/RepulsiveCurves-build/`) without changing any of its algorithms.
//! Curves: `trefoil40`, `freedman100` (every 6th vertex of the 592-vertex Freedman-He-Wang curve),
//! `freedman600` (the full curve; Gram matrices only as a digest: diagonal, row sums and three
//! rows), `open30` (an open curve, pins at 5 and 12, tangent pin at 20), `junction` (a tripod, a
//! loop attached at an arm end, a separate polyline; pins 3 and 14, tangent pin 7) and
//! `smallloops` (a 4-loop, a 5-loop and a 12-loop). Exponents (3, 6) and (2, 4.5). Numbers are
//! printed with 17 digits.
//!
//! # What is compared, with the tolerances
//!
//! * **Dense Sobolev Gram matrix** (`SobolevCurves::SobolevGramMatrix`, which assembles
//!   `B + B0` with the midpoint rule): our `A = B + B0` with `GramQuadrature::Midpoint`, and `B`
//!   and `B0` separately (the driver re-runs the reference's own two per-pair routines). Relative
//!   to the largest entry of the matrix: `GRAM_TOL = 1e-12` (measured 4e-16 to 5e-15; the
//!   two programs sum the pairs in different orders). This confirms D2 and the claim of
//!   `SPEC.md` §13.7: the reference uses *edge* tangents in its Gram matrix, so D1 (vertex
//!   tangents in its energy) does not reach the metric. The trapezoid variant is checked to
//!   *differ* from the reference (it is not what the reference computes).
//! * **Constraint Jacobians**, absolute `JAC_TOL = 1e-13` (entries are O(1); measured 5e-16): edge lengths, total
//!   length and point rows. The reference differentiates the lengths `l`, we differentiate the
//!   residuals `l0 - l`, so our rows are exactly the negatives of the reference's. Point rows
//!   (`gamma_v - x_v`) are identical. Column layout is the same (interleaved, `3 v + coordinate`),
//!   the edge orientation is the one of the input.
//! * **Coarsening, one level** (`PolyCurveNetwork::Coarsen` against `multigrid::coarsen` with the
//!   pins as `Point` and `Tangent(Vertex)` constraints): the set of surviving vertices, the two
//!   black neighbours of every white vertex, the prolongation weights (1, and 1/2 + 1/2) and the
//!   set of coarse edges, all as sets of fine vertex indices. Exact (integers and 0.5, 1).
//!
//! # What is not comparable, and why
//!
//! * **Energies and differentials** (D1): the reference sums over vertex pairs with vertex
//!   tangents and masses `sum l / deg`, we sum edge pairs (paper Eq. 17-18). They are reported
//!   side by side in the README, never matched.
//! * **Barycenter Jacobian** (D6): the reference uses constant dual-length weights
//!   `(sum_{I in v} l_I / deg v) / L` (and its `L` is the total length); we differentiate the
//!   length-weighted edge-midpoint barycenter exactly. The test checks the reference's rows
//!   against that formula, and that the weights agree with ours at degree-2 vertices up to the
//!   `(x_I - xbar) T_I^T` term, which only we have; it asserts no equality of the two matrices.
//! * **Tangent, surface constraints** (D7, D8) and penalties (D19, D20): other formulations or
//!   reference bugs; not dumped.
//! * **Coarsening on closed loops with fewer than 6 vertices** (D17): the reference coarsens them
//!   (a 4-loop gets two parallel coarse edges, a 5-loop becomes a triangle), we keep them black.
//!   The test asserts exactly this difference on `smallloops`.
//! * **Coarsening beyond level 1**: the reference numbers coarse vertices in BFS order and the
//!   BFS of the next level starts at coarse vertex 0 of each component, so the second level
//!   starts elsewhere than ours; also the coarse positions differ (`J+ gamma` against kept fine
//!   positions, D17). Only the topology of level 1, which depends on nothing but the fine
//!   network and the pins, is compared. The coarse vertex *numbering* differs as well, hence the
//!   comparison through fine indices.
//! * **Block cluster tree / BVH statistics**: the reference's admissibility is spatial only
//!   with different leaf sizes (D15, D16), so block counts have no counterpart. Our BCT product
//!   is tested against our own dense matrix in `tests/bct.rs`; its agreement with the
//!   *reference's* dense matrix follows from the equality tested here.

#![allow(clippy::needless_range_loop)]

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::PathBuf;

use ddg_repulsive_curves::constraints::{Constraint, ConstraintSet, TangentSite};
use ddg_repulsive_curves::multigrid::{ProlongationRow, coarsen, forced_black};
use ddg_repulsive_curves::{CurveNetwork, Exponents, GramQuadrature, Metric, gram};

const GRAM_TOL: f64 = 1e-12;
const JAC_TOL: f64 = 1e-13;

type P = [f64; 3];

/// rows, columns, entries
type Jac = (usize, usize, Vec<(usize, usize, f64)>);

#[derive(Default)]
struct Fixture {
    pos: Vec<P>,
    edges: Vec<[u32; 2]>,
    pins: Vec<u32>,
    tpins: Vec<u32>,
    digest: bool,
    /// name -> (rows, cols, entries)
    jac: HashMap<String, Jac>,
    /// packed upper triangles (full fixtures)
    sym: HashMap<String, (usize, Vec<f64>)>,
    gdiag: HashMap<String, Vec<f64>>,
    grows: HashMap<String, Vec<f64>>,
    grow: HashMap<String, Vec<(usize, Vec<f64>)>>,
    /// fine row, coarse col, weight
    prolong: Vec<(usize, usize, f64)>,
    n_coarse: usize,
    coarse_edges: Vec<[usize; 2]>,
}

fn load(name: &str) -> Fixture {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/reference")
        .join(format!("{name}.fx"));
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path:?}: {e}"));
    let mut f = Fixture::default();
    let mut cur: Option<String> = None;
    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let t: Vec<&str> = line.split_whitespace().collect();
        let num = |i: usize| -> f64 { t[i].parse().unwrap() };
        let idx = |i: usize| -> usize { t[i].parse().unwrap() };
        match t[0] {
            "n" => {}
            "v" => f.pos.push([num(1), num(2), num(3)]),
            "e" => f.edges.push([idx(1) as u32, idx(2) as u32]),
            "pin" => f.pins.push(idx(1) as u32),
            "tpin" => f.tpins.push(idx(1) as u32),
            "digest" => f.digest = idx(1) == 1,
            "jac" => {
                cur = Some(t[1].to_string());
                f.jac.insert(t[1].to_string(), (idx(2), idx(3), Vec::new()));
            }
            "t" => {
                let k = cur.clone().unwrap();
                f.jac.get_mut(&k).unwrap().2.push((idx(1), idx(2), num(3)));
            }
            "coarse" => f.n_coarse = idx(1),
            "P" => f.prolong.push((idx(1), idx(2), num(3))),
            "ce" => f.coarse_edges.push([idx(1), idx(2)]),
            "cv" => {}
            "sym" => {
                let n = idx(2);
                f.sym.insert(
                    t[1].to_string(),
                    (n, t[3..].iter().map(|s| s.parse().unwrap()).collect()),
                );
            }
            "gdiag" => {
                f.gdiag.insert(
                    t[1].to_string(),
                    t[3..].iter().map(|s| s.parse().unwrap()).collect(),
                );
            }
            "grows" => {
                f.grows.insert(
                    t[1].to_string(),
                    t[3..].iter().map(|s| s.parse().unwrap()).collect(),
                );
            }
            "grow" => {
                f.grow
                    .entry(t[1].to_string())
                    .or_default()
                    .push((idx(2), t[4..].iter().map(|s| s.parse().unwrap()).collect()));
            }
            other => panic!("unknown record {other}"),
        }
    }
    f
}

fn network(f: &Fixture) -> CurveNetwork {
    let net = CurveNetwork::new(f.pos.len(), f.edges.clone()).unwrap();
    assert_eq!(
        net.edges(),
        f.edges.as_slice(),
        "edge order and orientation"
    );
    net
}

fn max_rel(ours: impl Iterator<Item = (f64, f64)>) -> f64 {
    let (mut err, mut scale) = (0.0f64, 0.0f64);
    for (a, b) in ours {
        err = err.max((a - b).abs());
        scale = scale.max(b.abs());
    }
    err / scale
}

/// Unpacks the fixture's upper triangle.
fn full(n: usize, packed: &[f64]) -> Vec<f64> {
    let mut m = vec![0.0; n * n];
    let mut k = 0;
    for i in 0..n {
        for j in i..n {
            m[i * n + j] = packed[k];
            m[j * n + i] = packed[k];
            k += 1;
        }
    }
    assert_eq!(k, packed.len());
    m
}

const AB: [(f64, f64, &str); 2] = [(3.0, 6.0, "3_6"), (2.0, 4.5, "2_4.5")];

fn check_gram(name: &str) {
    let f = load(name);
    let net = network(&f);
    let n = f.pos.len();
    for (alpha, beta, tag) in AB {
        let ex = Exponents::new(alpha, beta).unwrap();
        let g = gram(
            &net,
            &f.pos,
            ex,
            &Metric::Sobolev {
                quadrature: GramQuadrature::Midpoint,
            },
            Some(1),
        )
        .unwrap();
        for (label, ours) in [("A", &g.a), ("B", &g.b), ("B0", &g.b0)] {
            let key = format!("{label}_{tag}");
            let rel = if !f.digest {
                let (m, packed) = &f.sym[&key];
                assert_eq!(*m, n);
                let r = full(n, packed);
                max_rel((0..n * n).map(|k| (ours.get(k / n, k % n), r[k])))
            } else {
                let diag = &f.gdiag[&key];
                let sums = &f.grows[&key];
                let mut pairs: Vec<(f64, f64)> = Vec::new();
                for i in 0..n {
                    pairs.push((ours.get(i, i), diag[i]));
                    pairs.push((ours.row(i).iter().sum::<f64>(), sums[i]));
                }
                for (r, row) in &f.grow[&key] {
                    for j in 0..n {
                        pairs.push((ours.get(*r, j), row[j]));
                    }
                }
                max_rel(pairs.into_iter())
            };
            eprintln!("{name} {key}: relative difference {rel:.2e}");
            assert!(rel < GRAM_TOL, "{name} {key}: {rel:.3e}");
        }
        // The trapezoid rule of Eq. 20-21 is not what the reference computes.
        if !f.digest {
            let t = gram(
                &net,
                &f.pos,
                ex,
                &Metric::Sobolev {
                    quadrature: GramQuadrature::Trapezoid,
                },
                Some(1),
            )
            .unwrap();
            let (_, packed) = &f.sym[&format!("A_{tag}")];
            let r = full(n, packed);
            let rel = max_rel((0..n * n).map(|k| (t.a.get(k / n, k % n), r[k])));
            eprintln!("{name} A_{tag}: trapezoid against the reference {rel:.2e}");
            assert!(
                rel > 1e-5,
                "trapezoid unexpectedly equal to the reference: {rel:e}"
            );
        }
    }
}

#[test]
fn gram_trefoil40() {
    check_gram("trefoil40");
}

#[test]
fn gram_freedman100() {
    check_gram("freedman100");
}

#[test]
fn gram_open30() {
    check_gram("open30");
}

#[test]
fn gram_junction() {
    check_gram("junction");
}

/// The full 592-vertex Freedman-He-Wang curve; the reference's matrix is digested.
#[test]
fn gram_freedman600_digest() {
    check_gram("freedman600");
}

fn dense_jacobian(
    set: &ConstraintSet,
    net: &CurveNetwork,
    pos: &[P],
) -> BTreeMap<(usize, usize), f64> {
    let c = set.jacobian(net, pos).unwrap();
    let mut m = BTreeMap::new();
    for (r, col, v) in c.triplets() {
        *m.entry((r as usize, col as usize)).or_insert(0.0) += v;
    }
    m
}

fn compare_jacobian(what: &str, reference: &Jac, ours: &BTreeMap<(usize, usize), f64>, sign: f64) {
    let mut theirs: BTreeMap<(usize, usize), f64> = BTreeMap::new();
    for &(r, c, v) in &reference.2 {
        *theirs.entry((r, c)).or_insert(0.0) += sign * v;
    }
    let keys: BTreeSet<_> = theirs.keys().chain(ours.keys()).copied().collect();
    let mut worst = 0.0f64;
    for k in keys {
        let a = ours.get(&k).copied().unwrap_or(0.0);
        let b = theirs.get(&k).copied().unwrap_or(0.0);
        worst = worst.max((a - b).abs());
    }
    eprintln!("{what}: max abs difference {worst:.2e}");
    assert!(worst < JAC_TOL, "{what}: {worst:e}");
}

fn check_jacobians(name: &str) {
    let f = load(name);
    let net = network(&f);
    let n = f.pos.len();
    // Edge lengths: ours are the rows of l0 - l, the reference's those of l.
    let set = ConstraintSet::new(
        &net,
        &f.pos,
        &[Constraint::EdgeLengths {
            edges: None,
            targets: None,
        }],
        None,
    )
    .unwrap();
    assert_eq!(set.n_active_rows(), f.edges.len());
    let r = &f.jac["edgelengths"];
    assert_eq!((r.0, r.1), (f.edges.len(), 3 * n));
    compare_jacobian("edge lengths", r, &dense_jacobian(&set, &net, &f.pos), -1.0);
    // Total length.
    let set = ConstraintSet::new(
        &net,
        &f.pos,
        &[Constraint::TotalLength { target: None }],
        None,
    )
    .unwrap();
    compare_jacobian(
        "total length",
        &f.jac["totallength"],
        &dense_jacobian(&set, &net, &f.pos),
        -1.0,
    );
    // Points.
    if !f.pins.is_empty() {
        let cs: Vec<Constraint> = f
            .pins
            .iter()
            .map(|&v| Constraint::Point {
                vertex: v,
                target: None,
            })
            .collect();
        let set = ConstraintSet::new(&net, &f.pos, &cs, None).unwrap();
        compare_jacobian(
            "points",
            &f.jac["points"],
            &dense_jacobian(&set, &net, &f.pos),
            1.0,
        );
    }
}

#[test]
fn jacobians_trefoil40() {
    check_jacobians("trefoil40");
}

#[test]
fn jacobians_freedman100() {
    check_jacobians("freedman100");
}

#[test]
fn jacobians_open30() {
    check_jacobians("open30");
}

#[test]
fn jacobians_junction() {
    check_jacobians("junction");
}

#[test]
fn jacobians_freedman600() {
    check_jacobians("freedman600");
}

/// D6: the reference's barycenter rows are constant dual-length weights. We check the fixture
/// against that formula and against our first term at degree-2 vertices; the matrices are not
/// equal (ours has the extra `(x_I - xbar) T_I^T` term and edge weights `l_I / 2`).
#[test]
fn barycenter_is_dual_length_weights() {
    for name in ["trefoil40", "open30", "junction"] {
        let f = load(name);
        let net = network(&f);
        let n = f.pos.len();
        let geom = net.edge_geometry(&f.pos);
        let total: f64 = geom.length.iter().sum();
        let (_, _, entries) = &f.jac["barycenter"];
        let mut w = vec![[0.0; 3]; n];
        for &(r, c, v) in entries {
            assert_eq!(c % 3, r, "barycenter rows are diagonal blocks");
            w[c / 3][r] += v;
        }
        let mut max_dev_deg2 = 0.0f64;
        let mut max_dev_other = 0.0f64;
        for v in 0..n as u32 {
            let inc = net.incident_edges(v);
            let dual: f64 =
                inc.iter().map(|&e| geom.length[e as usize]).sum::<f64>() / inc.len() as f64;
            let expect = dual / total;
            for r in 0..3 {
                assert!(
                    (w[v as usize][r] - expect).abs() < 1e-14,
                    "{name} vertex {v}"
                );
            }
            let ours: f64 = inc
                .iter()
                .map(|&e| 0.5 * geom.length[e as usize])
                .sum::<f64>()
                / total;
            let dev = (expect - ours).abs();
            if inc.len() == 2 {
                max_dev_deg2 = max_dev_deg2.max(dev);
            } else {
                max_dev_other = max_dev_other.max(dev);
            }
        }
        eprintln!(
            "{name}: reference weight minus our first term: degree 2 {max_dev_deg2:.1e}, other degrees {max_dev_other:.1e}"
        );
        assert!(max_dev_deg2 < 1e-14);
    }
}

/// Coarsening topology of one level, as sets of fine indices.
struct RefLevel {
    /// fine vertex -> Black(fine id of itself) or White(two fine ids)
    kept: Vec<bool>,
    white: Vec<Option<[u32; 2]>>,
    edges: Vec<[u32; 2]>,
}

fn ref_level(f: &Fixture) -> RefLevel {
    let n = f.pos.len();
    let mut cols_of: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    for &(r, c, w) in &f.prolong {
        cols_of[r].push((c, w));
    }
    let mut coarse_to_fine = vec![u32::MAX; f.n_coarse];
    let mut kept = vec![false; n];
    for v in 0..n {
        if cols_of[v].len() == 1 {
            assert_eq!(cols_of[v][0].1, 1.0);
            kept[v] = true;
            coarse_to_fine[cols_of[v][0].0] = v as u32;
        } else {
            assert_eq!(cols_of[v].len(), 2);
            assert!(cols_of[v].iter().all(|&(_, w)| w == 0.5));
        }
    }
    let white = (0..n)
        .map(|v| {
            (!kept[v]).then(|| {
                let mut p = [
                    coarse_to_fine[cols_of[v][0].0],
                    coarse_to_fine[cols_of[v][1].0],
                ];
                p.sort_unstable();
                p
            })
        })
        .collect();
    let mut edges: Vec<[u32; 2]> = f
        .coarse_edges
        .iter()
        .map(|&[a, b]| {
            let mut p = [coarse_to_fine[a], coarse_to_fine[b]];
            p.sort_unstable();
            p
        })
        .collect();
    edges.sort_unstable();
    RefLevel { kept, white, edges }
}

fn pin_constraints(f: &Fixture) -> Vec<Constraint> {
    let mut cs: Vec<Constraint> = f
        .pins
        .iter()
        .map(|&v| Constraint::Point {
            vertex: v,
            target: None,
        })
        .collect();
    cs.extend(f.tpins.iter().map(|&v| Constraint::Tangent {
        site: TangentSite::Vertex(v),
        direction: None,
    }));
    cs
}

/// Vertices that lie in closed components of fewer than 6 vertices (not coarsened by us, D17).
fn small_loop_vertices(net: &CurveNetwork) -> BTreeSet<u32> {
    net.components()
        .iter()
        .filter(|c| c.closed && c.vertices.len() < 6)
        .flat_map(|c| c.vertices.iter().copied())
        .collect()
}

fn check_coarsening(name: &str) {
    let f = load(name);
    let net = network(&f);
    let rl = ref_level(&f);
    let cs = pin_constraints(&f);
    let set = (!cs.is_empty()).then(|| ConstraintSet::new(&net, &f.pos, &cs, None).unwrap());
    let forced = forced_black(f.pos.len(), set.as_ref());
    let ours = coarsen(&net, &forced).unwrap().expect("coarsenable");
    let small = small_loop_vertices(&net);
    let n = f.pos.len();
    let c2f = ours.coarse_to_fine();
    for v in 0..n {
        if small.contains(&(v as u32)) {
            assert!(ours.fine_to_coarse()[v].is_some(), "small loops stay black");
            continue;
        }
        let ours_black = ours.fine_to_coarse()[v].is_some();
        assert_eq!(ours_black, rl.kept[v], "{name}: vertex {v} kept");
        if let ProlongationRow::White(a, b) = ours.rows()[v] {
            let mut p = [c2f[a as usize], c2f[b as usize]];
            p.sort_unstable();
            assert_eq!(
                Some(p),
                rl.white[v],
                "{name}: neighbours of white vertex {v}"
            );
        }
    }
    // Coarse edges, restricted to edges not in small loops.
    let in_small = |e: &[u32; 2]| small.contains(&e[0]);
    let mut ours_edges: Vec<[u32; 2]> = ours
        .coarse()
        .edges()
        .iter()
        .map(|&[a, b]| {
            let mut p = [c2f[a as usize], c2f[b as usize]];
            p.sort_unstable();
            p
        })
        .filter(|e| !in_small(e))
        .collect();
    ours_edges.sort_unstable();
    let ref_edges: Vec<[u32; 2]> = rl.edges.iter().copied().filter(|e| !in_small(e)).collect();
    assert_eq!(ours_edges, ref_edges, "{name}: coarse edges");
    let kept = rl.kept.iter().filter(|&&k| k).count();
    eprintln!(
        "{name}: {} of {} vertices kept, {} coarse edges agree",
        kept,
        n,
        ref_edges.len()
    );
}

#[test]
fn coarsening_trefoil40() {
    check_coarsening("trefoil40");
}

#[test]
fn coarsening_freedman100() {
    check_coarsening("freedman100");
}

#[test]
fn coarsening_freedman600() {
    check_coarsening("freedman600");
}

#[test]
fn coarsening_open30_with_pins() {
    check_coarsening("open30");
}

#[test]
fn coarsening_junction_with_pins() {
    check_coarsening("junction");
}

/// D17: loops of fewer than 6 vertices. The reference coarsens them, we do not.
#[test]
fn coarsening_small_loops_differ_as_documented() {
    let f = load("smallloops");
    let net = network(&f);
    let rl = ref_level(&f);
    // Vertices 0..4 are the 4-loop, 4..9 the 5-loop, 9..21 the 12-loop.
    let ours = coarsen(&net, &forced_black(f.pos.len(), None))
        .unwrap()
        .expect("the 12-loop coarsens");
    for v in 0..9 {
        assert!(ours.fine_to_coarse()[v].is_some(), "ours: vertex {v} stays");
    }
    let ref_kept_small = (0..9).filter(|&v| rl.kept[v]).count();
    // 4-loop: vertices 0 and 2 survive; 5-loop: three of five.
    assert_eq!(ref_kept_small, 2 + 3);
    let small_edges: Vec<_> = rl.edges.iter().filter(|e| e[0] < 4).collect();
    assert_eq!(small_edges.len(), 2);
    assert_eq!(
        small_edges[0], small_edges[1],
        "the reference makes a double edge"
    );
    // The 12-loop agrees.
    for v in 9..21 {
        assert_eq!(ours.fine_to_coarse()[v].is_some(), rl.kept[v], "vertex {v}");
    }
}
