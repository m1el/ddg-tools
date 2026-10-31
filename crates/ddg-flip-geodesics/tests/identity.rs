//! Representation identity (SPEC §1.3, §7 layer 2 and test 8): FlipOut on the signposts of
//! `ddg_intrinsic` and on the integer coordinates of `ddg_integer_coords` makes the same
//! flips in the same order and returns the same edge paths with bitwise-identical lengths;
//! the integer triangulation certifies at the end of every run (F25). Work package F2 adds
//! the paths, loops and networks; F4 extends this file with Bézier curves and the iCDT.

mod common;

use common::{Shape, data_shapes, generated_shapes, random_pairs, vpath, waist_tube};
use ddg_flip_geodesics::{EdgePath, GeodesicNetwork, Segment, ShortenOptions, dijkstra_path};
use ddg_integer_coords::IntegerTriangulation;
use ddg_intrinsic::{DelaunayOptions, IntrinsicTriangulation, Triangulation};

/// Everything a run produces that both representations must agree on, as bits.
#[derive(Debug, PartialEq)]
struct Trace {
    steps: Vec<String>,
    report: String,
    paths: Vec<EdgePath>,
    faces: Vec<[u32; 3]>,
    lengths: Vec<u64>,
    length: u64,
}

fn trace<T: Triangulation + Clone>(tri: T, paths: &[EdgePath]) -> (Trace, T) {
    let mut net = GeodesicNetwork::new(tri, paths, &[], None).unwrap();
    let mut steps = Vec::new();
    let report = net
        .shorten_with(&ShortenOptions::default(), |step, _| {
            // Per iteration: the joint, side, angle bits, flips, outcome and running length.
            steps.push(format!(
                "{:?} {:?} {:x} {} {} {:?} {:x}",
                step.joint,
                step.side,
                step.angle.to_bits(),
                step.flips,
                step.refused_flips,
                step.shortening,
                step.length.to_bits()
            ));
        })
        .unwrap();
    let t = Trace {
        steps,
        report: format!("{report:?}"),
        paths: net.paths(),
        faces: net.triangulation().complex().faces().to_vec(),
        lengths: net
            .triangulation()
            .edge_lengths()
            .iter()
            .map(|x| x.to_bits())
            .collect(),
        length: net.length().to_bits(),
    };
    (t, net.into_triangulation())
}

/// Both representations of `shape`, Delaunay-flipped first if asked (the same flips: the
/// shared `flip_to_delaunay`).
fn both(shape: &Shape, delaunay: bool) -> (IntrinsicTriangulation, IntegerTriangulation) {
    let mut signposts = shape.triangulation();
    let mut integers = IntegerTriangulation::new(&shape.mesh, &shape.geometry)
        .unwrap_or_else(|e| panic!("{}: {e}", shape.name));
    if delaunay {
        signposts
            .flip_to_delaunay(&DelaunayOptions::default())
            .unwrap();
        integers
            .flip_to_delaunay(&DelaunayOptions::default())
            .unwrap();
        assert_eq!(signposts.complex().faces(), integers.complex().faces());
    }
    (signposts, integers)
}

fn compare(shape: &Shape, delaunay: bool, inputs: &[Vec<EdgePath>]) -> usize {
    let (signposts, integers) = both(shape, delaunay);
    let mut iterations = 0;
    for paths in inputs {
        if ddg_flip_geodesics::Network::new(signposts.complex(), paths, &[], None).is_err() {
            continue;
        }
        let (a, _) = trace(signposts.clone(), paths);
        let (b, int_after) = trace(integers.clone(), paths);
        assert!(
            a == b,
            "{}: the representations differ:\n{a:?}\n{b:?}",
            shape.name
        );
        int_after
            .certify()
            .unwrap_or_else(|e| panic!("{}: integer certificate: {e}", shape.name));
        iterations += a.steps.len();
    }
    iterations
}

#[test]
fn signposts_and_integers_shorten_identically() {
    let mut shapes = generated_shapes();
    shapes.extend(data_shapes());
    let mut total = 0;
    for shape in &shapes {
        let n = shape.mesh.n_vertices();
        for delaunay in [false, true] {
            let (tri, _) = both(shape, delaunay);
            let mut inputs: Vec<Vec<EdgePath>> = random_pairs(n, 6, 99)
                .into_iter()
                .map(|(a, b)| vec![dijkstra_path(&tri, a, b).unwrap()])
                .collect();
            // Loops along self-edges (App. B moves, contractions).
            for &s in &tri.irregularities().self_edges {
                inputs.push(vec![EdgePath::closed(vec![Segment::new(s as usize, true)])]);
            }
            let steps = compare(shape, delaunay, &inputs);
            eprintln!(
                "{} (delaunay {delaunay}): {} inputs, {steps} iterations identical",
                shape.name,
                inputs.len()
            );
            total += steps;
        }
    }
    // Loops and two parallel loops on a tube with a waist.
    let (n, w) = (10, 3);
    let tube = waist_tube(n, w);
    let tri = tube.triangulation();
    let ring = |j: usize| vpath(&tri, &(0..n).map(|k| j * n + k).collect::<Vec<_>>(), true);
    total += compare(
        &tube,
        false,
        &[vec![ring(w + 2)], vec![ring(w - 2), ring(w + 2)]],
    );
    assert!(total > 50, "only {total} iterations compared");
}
