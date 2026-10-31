//! Meshes and states shared by the integration tests.
#![allow(dead_code)]

mod shapes;

pub use shapes::*;

use ddg_core::generators::SplitMix64;
use ddg_integer_coords::{ConversionReport, IntegerError, IntegerTriangulation};
use ddg_intrinsic::{DelaunayOptions, EdgeTraceOptions, IntrinsicTriangulation};

/// The standard shapes of `ddg-intrinsic`'s tests, then spot and camelbad if present.
pub fn all_shapes() -> Vec<Shape> {
    let mut shapes = standard_shapes();
    shapes.extend(load_data("models/spot/spot/spot_triangulated.obj"));
    shapes.extend(load_data(
        "optimal-cones/optimal-cones-data/data/camelbad.obj",
    ));
    shapes
}

/// How a test state is made from the input triangulation with the signpost representation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Scenario {
    /// The input itself.
    Initial,
    /// Many random flips.
    RandomFlips,
    /// Delaunay flipping.
    Delaunay,
    /// Random vertices inserted in faces, then Delaunay flipping.
    Insertions,
    /// Random edges split, then Delaunay flipping.
    EdgeSplits,
    /// Edges along input edges split (inserted vertices on input edges), then random flips:
    /// curves of input edges pass through inserted vertices (SPEC I17).
    SplitInputEdgesThenFlips,
    /// Everything interleaved, then Delaunay flipping.
    Mixed,
}

pub const SCENARIOS: [Scenario; 7] = [
    Scenario::Initial,
    Scenario::RandomFlips,
    Scenario::Delaunay,
    Scenario::Insertions,
    Scenario::EdgeSplits,
    Scenario::SplitInputEdgesThenFlips,
    Scenario::Mixed,
];

/// Random barycentric coordinates, each at least 0.05.
pub fn inner_bary(rng: &mut SplitMix64) -> [f64; 3] {
    let (mut a, mut b) = (rng.next_f64(), rng.next_f64());
    if a + b > 1.0 {
        (a, b) = (1.0 - a, 1.0 - b);
    }
    let shrink = |x: f64| 0.05 + 0.85 * x;
    [shrink(a), shrink(b), shrink(1.0 - a - b)]
}

fn random_flips(t: &mut IntrinsicTriangulation, rng: &mut SplitMix64, n: usize) -> usize {
    let mut done = 0;
    for _ in 0..n {
        let e = (rng.next_u64() % t.n_edges() as u64) as usize;
        if t.flip_edge(e).is_ok() {
            done += 1;
        }
    }
    done
}

fn random_insertions(t: &mut IntrinsicTriangulation, rng: &mut SplitMix64, n: usize) {
    for _ in 0..n {
        let f = (rng.next_u64() % t.n_faces() as u64) as usize;
        let _ = t.insert_vertex(f, inner_bary(rng));
    }
}

fn random_splits(t: &mut IntrinsicTriangulation, rng: &mut SplitMix64, n: usize, on_input: bool) {
    for _ in 0..n {
        let e = (rng.next_u64() % t.n_edges() as u64) as usize;
        if on_input && t.edge_on_input(e).is_none() {
            continue;
        }
        let _ = t.split_edge(e, rng.uniform(0.2, 0.8));
    }
}

/// A signpost triangulation of `shape` in the given scenario. The amount of work grows with the
/// mesh and is capped so that the tests stay fast in debug builds.
pub fn state(shape: &Shape, scenario: Scenario, seed: u64) -> IntrinsicTriangulation {
    let mut t = shape.triangulation();
    let mut rng = SplitMix64::new(seed);
    let ne = t.n_edges();
    let nf = t.n_faces();
    let delaunay = |t: &mut IntrinsicTriangulation| {
        t.flip_to_delaunay(&DelaunayOptions::default())
            .unwrap_or_else(|e| panic!("{}: {e}", shape.name));
    };
    match scenario {
        Scenario::Initial => {}
        Scenario::RandomFlips => {
            random_flips(&mut t, &mut rng, (10 * ne).min(40_000));
        }
        Scenario::Delaunay => delaunay(&mut t),
        Scenario::Insertions => {
            random_insertions(&mut t, &mut rng, (nf / 4).clamp(3, 150));
            delaunay(&mut t);
        }
        Scenario::EdgeSplits => {
            random_splits(&mut t, &mut rng, (ne / 6).clamp(3, 150), false);
            delaunay(&mut t);
        }
        Scenario::SplitInputEdgesThenFlips => {
            random_splits(&mut t, &mut rng, (ne / 6).clamp(3, 150), true);
            let n = (10 * t.n_edges()).min(40_000);
            random_flips(&mut t, &mut rng, n);
        }
        Scenario::Mixed => {
            for _ in 0..3 {
                random_flips(&mut t, &mut rng, (2 * ne).min(5_000));
                random_insertions(&mut t, &mut rng, (nf / 20).clamp(2, 40));
                random_splits(&mut t, &mut rng, (ne / 30).clamp(2, 40), false);
                random_splits(&mut t, &mut rng, (ne / 30).clamp(2, 40), true);
            }
            delaunay(&mut t);
        }
    }
    t
}

/// The integer coordinates of a state, or the error of the conversion.
pub fn convert(
    t: &IntrinsicTriangulation,
) -> Result<(IntegerTriangulation, ConversionReport), IntegerError> {
    IntegerTriangulation::from_intrinsic(t, &EdgeTraceOptions::default())
}

/// The seeds used per scenario.
pub const SEEDS: [u64; 2] = [11, 2021];
