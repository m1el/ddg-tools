//! SPEC §5 item 8: against the authors' program, the same intrinsic Delaunay triangulation on
//! the same inputs, compared as sets of edges with their lengths.
//!
//! The reference outputs are in `data/intrinsic-reference/<input>/flip/` (not part of the
//! repository; made by `scripts/intrinsic-reference.sh`, format in
//! `scripts/intrinsic-reference/FORMAT.md`). An input that is absent is skipped with a notice.
//!
//! What is compared, per input, for the triangulation flipped with the reference's test
//! (cotan weight below `−1e-6`, [`DelaunayCriterion::REFERENCE`]) and with this crate's default
//! (angle sum above `π + 1e-10`):
//!
//! * the multiset of halfedges `(min vertex, max vertex, length)`: every interior edge twice,
//!   every boundary edge once; lengths to `1e-8` relative;
//! * the number of self-edges and of vertex pairs joined by several edges;
//! * the cotan Laplacian against the reference's `laplace.spmat`, when the triangulations
//!   are the same;
//! * both sets of lengths against a third measurement: the length in space of the path that
//!   this crate traces for the edge over the input mesh.
//!
//! Measured with the reference's test: the same triangulation on 14 of the 15 inputs, with
//! the same numbers of self-edges and double edges; on `blub` 68 halfedges differ, within its
//! 43 ties. With the default test `bob` differs too, in 24 halfedges at edges whose weight is
//! within `1e-6` of zero, which the reference does not flip. Lengths agree to `4e-12` relative
//! or better, except on `camelshotnoise` (input angles of 0.003°), where 12 edges differ by up
//! to `1.9e-9`; there the lengths of this crate agree with the traced paths to `2e-12` and the
//! reference's do not, so that difference is the reference's rounding.
//!
//! Without cocircular ties the Delaunay triangulation is unique and the two must be equal. A
//! tie is a quadrilateral both of whose diagonals are Delaunay; each program keeps the
//! diagonal it meets, so up to 4 halfedges (two per list) may differ for every tie.
//!
//! Run with `cargo test -p ddg-intrinsic --test reference -- --nocapture` for the table.

mod common;

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use common::{Shape, path_length};
use ddg_core::io::read_obj;
use ddg_intrinsic::{
    DelaunayCriterion, DelaunayOptions, DelaunayStatus, EdgeTraceOptions, IntrinsicTriangulation,
};

const INPUTS: [&str; 15] = [
    "spot",
    "bob",
    "blub",
    "bunny",
    "camel",
    "camelbad",
    "camelbadnonuniform",
    "camel_noise1",
    "camelshotnoise",
    "bump",
    "icosphere",
    "icosphere_jitter",
    "thin_cylinder",
    "long_torus",
    "sheared_grid",
];

fn reference_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/intrinsic-reference")
}

/// A `# dense R C` file: `R` rows of `C` numbers.
fn read_dense(path: &Path) -> Vec<Vec<f64>> {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap().split_whitespace().collect();
    assert_eq!(&header[..2], ["#", "dense"], "{}", path.display());
    let (rows, cols): (usize, usize) = (header[2].parse().unwrap(), header[3].parse().unwrap());
    let data: Vec<Vec<f64>> = lines
        .filter(|l| !l.trim().is_empty())
        .map(|l| l.split_whitespace().map(|x| x.parse().unwrap()).collect())
        .collect();
    assert_eq!(data.len(), rows, "{}", path.display());
    assert!(
        data.iter().all(|r: &Vec<f64>| r.len() == cols),
        "{}",
        path.display()
    );
    data
}

/// A `# sparse R C` file: `row col value`, 1-based.
fn read_sparse(path: &Path) -> (usize, Vec<(usize, usize, f64)>) {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap().split_whitespace().collect();
    assert_eq!(&header[..2], ["#", "sparse"], "{}", path.display());
    let n: usize = header[2].parse().unwrap();
    let entries = lines
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let mut it = l.split_whitespace();
            let i: usize = it.next().unwrap().parse().unwrap();
            let j: usize = it.next().unwrap().parse().unwrap();
            (i - 1, j - 1, it.next().unwrap().parse().unwrap())
        })
        .collect();
    (n, entries)
}

/// The integer after `"key":` in a flat JSON text, or after `key ` in a `key value` text.
fn number(text: &str, key: &str) -> usize {
    let at = text
        .find(&format!("\"{key}\":"))
        .map(|i| i + key.len() + 3)
        .or_else(|| {
            text.lines()
                .position(|l| l.split_whitespace().next() == Some(key))
                .map(|n| text.lines().take(n).map(|l| l.len() + 1).sum::<usize>() + key.len())
        })
        .unwrap_or_else(|| panic!("no {key}"));
    text[at..]
        .trim_start()
        .split(|c: char| !c.is_ascii_digit())
        .next()
        .unwrap()
        .parse()
        .unwrap()
}

type Halfedges = Vec<(usize, usize, f64)>;

/// The halfedges of the reference triangulation as `(min, max, length)`.
fn reference_halfedges(inds: &[Vec<f64>], lengths: &[Vec<f64>]) -> Halfedges {
    let mut out = Vec::with_capacity(3 * inds.len());
    for (face, l) in inds.iter().zip(lengths) {
        for k in 0..3 {
            let (a, b) = (face[k] as usize, face[(k + 1) % 3] as usize);
            out.push((a.min(b), a.max(b), l[k]));
        }
    }
    out
}

fn our_halfedges(tri: &IntrinsicTriangulation) -> Halfedges {
    (0..tri.n_halfedges())
        .map(|h| {
            let (a, b) = (tri.from_vertex(h), tri.to_vertex(h));
            (a.min(b), a.max(b), tri.edge_length(tri.edge(h)))
        })
        .collect()
}

/// The number of halfedges of either list that have no partner in the other (same vertex pair,
/// length within `tolerance` relative), and the largest relative length difference among the
/// matched ones.
fn compare(ours: &Halfedges, theirs: &Halfedges, tolerance: f64) -> (usize, f64) {
    let mut by_pair: HashMap<(usize, usize), Vec<f64>> = HashMap::new();
    for &(a, b, l) in theirs {
        by_pair.entry((a, b)).or_default().push(l);
    }
    let mut unmatched = 0;
    let mut largest = 0.0f64;
    for &(a, b, l) in ours {
        let partner = by_pair.get_mut(&(a, b)).and_then(|lengths| {
            let (k, _) = lengths
                .iter()
                .enumerate()
                .min_by(|x, y| (x.1 - l).abs().total_cmp(&(y.1 - l).abs()))?;
            ((lengths[k] - l).abs() <= tolerance * l).then(|| lengths.swap_remove(k))
        });
        match partner {
            Some(m) => largest = largest.max((m - l).abs() / l),
            None => unmatched += 1,
        }
    }
    (
        unmatched + by_pair.values().map(Vec::len).sum::<usize>(),
        largest,
    )
}

fn irregular_counts(tri: &IntrinsicTriangulation) -> (usize, usize) {
    let irregular = tri.irregularities();
    let mut pairs: Vec<[usize; 2]> = irregular
        .multi_edges
        .iter()
        .map(|&e| {
            let [a, b] = tri.edge_vertices(e as usize);
            [a.min(b), a.max(b)]
        })
        .collect();
    pairs.sort_unstable();
    pairs.dedup();
    (irregular.self_edges.len(), pairs.len())
}

/// The largest difference between our Laplacian and the reference's, over the largest entry.
fn laplacian_difference(tri: &IntrinsicTriangulation, path: &Path) -> f64 {
    let (n, entries) = read_sparse(path);
    assert_eq!(n, tri.n_vertices());
    let ours = tri.cotan_laplacian();
    let mut scale = 0.0f64;
    let mut largest = 0.0f64;
    let mut seen: HashMap<(usize, usize), f64> = HashMap::new();
    for &(i, j, value) in &entries {
        *seen.entry((i, j)).or_default() += value;
    }
    for (&(i, j), &value) in &seen {
        scale = scale.max(value.abs());
        largest = largest.max((ours.get(i, j) - value).abs());
    }
    // Entries we store that the reference does not write are zeros.
    for (i, j, value) in ours.triplets() {
        if !seen.contains_key(&(i as usize, j as usize)) {
            largest = largest.max(value.abs());
        }
    }
    largest / scale
}

#[test]
fn the_delaunay_triangulation_is_the_reference_s() {
    let root = reference_dir();
    if !root.is_dir() {
        eprintln!(
            "SKIPPED: reference directory {} not found; nothing was checked",
            root.display()
        );
        return;
    }
    println!(
        "{:<20} {:>7} | {:>6} {:>6} {:>6} | {:>5} {:>5} | {:>9} {:>9} {:>9} | {:>9} {:>5} | {:>9} {:>9}",
        "input",
        "edges",
        "flips",
        "(ref)",
        "(dflt)",
        "self",
        "pairs",
        "unmatched",
        "max dlen",
        "Laplacian",
        "unm. dflt",
        "ties",
        "ours-path",
        "ref-path"
    );
    let mut compared = 0;
    for name in INPUTS {
        let dir = root.join(name);
        let obj_path = dir.join("input.obj");
        if !obj_path.is_file() {
            eprintln!(
                "SKIPPED: {} not found; nothing was checked for it",
                obj_path.display()
            );
            continue;
        }
        let obj = read_obj(&obj_path).unwrap_or_else(|e| panic!("{e}"));
        let faces = obj
            .triangles()
            .unwrap_or_else(|| panic!("{name}: not triangles"));
        let shape = Shape::new(name, obj.positions, faces, false);
        let flip = dir.join("flip");
        let inds = read_dense(&flip.join("faceInds.dmat"));
        let lengths = read_dense(&flip.join("faceLengths.dmat"));
        let theirs = reference_halfedges(&inds, &lengths);
        let analysis = std::fs::read_to_string(flip.join("analysis.json")).unwrap();
        let stats = std::fs::read_to_string(flip.join("stats.txt")).unwrap();
        let reference_flips = number(&stats, "flips");

        // With the reference's test.
        let mut tri = shape.triangulation();
        let report = tri
            .flip_to_delaunay(&DelaunayOptions {
                criterion: DelaunayCriterion::REFERENCE,
                ..DelaunayOptions::default()
            })
            .unwrap();
        assert_eq!(
            report.status,
            DelaunayStatus::Delaunay,
            "{name}: {report:?}"
        );
        assert_eq!(inds.len(), tri.n_faces(), "{name}");
        let (unmatched, max_difference) = compare(&our_halfedges(&tri), &theirs, 1e-8);
        let (self_edges, pairs) = irregular_counts(&tri);
        // Edges whose weight is within the reference's tolerance of zero: the ties of its test.
        let ties = report.n_ties;

        // With the default test.
        let mut ours = shape.triangulation();
        let default_report = ours.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        assert_eq!(default_report.status, DelaunayStatus::Delaunay, "{name}");
        let (unmatched_default, _) = compare(&our_halfedges(&ours), &theirs, 1e-8);
        let default_ties = (0..ours.n_edges())
            .filter(|&e| !ours.is_boundary_edge(e) && ours.edge_cotan_weight(e).abs() <= 1e-6)
            .count();

        let laplacian = if unmatched == 0 {
            laplacian_difference(&tri, &flip.join("laplace.spmat"))
        } else {
            f64::NAN
        };
        // A third measurement of every edge: the length in space of its traced path.
        let traces = tri
            .trace_edges(&EdgeTraceOptions::default())
            .unwrap_or_else(|e| panic!("{name}: {e}"));
        let mut reference_lengths: HashMap<(usize, usize), Vec<f64>> = HashMap::new();
        for &(a, b, l) in &theirs {
            reference_lengths.entry((a, b)).or_default().push(l);
        }
        let (mut ours_off, mut reference_off) = (0.0f64, 0.0f64);
        for e in 0..tri.n_edges() {
            let [a, b] = tri.edge_vertices(e);
            let length = tri.edge_length(e);
            let path = path_length(
                &shape,
                tri.vertex_location(a),
                traces.crossings(e),
                tri.vertex_location(b),
            );
            ours_off = ours_off.max((length - path).abs() / length);
            if let Some(candidates) = reference_lengths.get(&(a.min(b), a.max(b))) {
                let nearest = candidates
                    .iter()
                    .map(|&r| (r - path).abs() / length)
                    .fold(f64::INFINITY, f64::min);
                reference_off = reference_off.max(nearest);
            }
        }
        println!(
            "{:<20} {:>7} | {:>6} {:>6} {:>6} | {:>5} {:>5} | {:>9} {:>9.1e} {:>9.1e} | {:>9} {:>5} | {:>9.1e} {:>9.1e}",
            name,
            tri.n_edges(),
            reference_flips,
            report.n_flips,
            default_report.n_flips,
            self_edges,
            pairs,
            unmatched,
            max_difference,
            laplacian,
            unmatched_default,
            default_ties,
            ours_off,
            reference_off,
        );
        assert!(
            ours_off < 1e-10,
            "{name}: an edge and its traced path differ by {ours_off}"
        );
        // The same triangulation, up to ties.
        assert!(
            unmatched <= 4 * ties,
            "{name}: {unmatched} halfedges differ from the reference, {ties} ties"
        );
        assert!(
            unmatched_default <= 4 * default_ties.max(ties),
            "{name}: {unmatched_default} halfedges differ with the default test, {default_ties} ties"
        );
        if ties == 0 {
            assert_eq!(unmatched, 0, "{name}");
            assert_eq!(
                (self_edges, pairs),
                (
                    number(&analysis, "self_edges"),
                    number(&analysis, "double_edge_pairs")
                ),
                "{name}"
            );
            // The lengths agree to the rounding of the flips (see the module documentation
            // for camelshotnoise); the Laplacian to the rounding of the cotangents of the
            // thinnest triangles (the reference's own check of its file reached 1.4e-7
            // absolute on entries of 1e5 for camelshotnoise).
            assert!(
                max_difference < 5e-9,
                "{name}: lengths differ by {max_difference}"
            );
            if name != "camelshotnoise" {
                assert!(
                    max_difference < 1e-11,
                    "{name}: lengths differ by {max_difference}"
                );
            }
            assert!(laplacian < 1e-6, "{name}: Laplacians differ by {laplacian}");
        }
        compared += 1;
    }
    println!("{compared} inputs compared");
}
