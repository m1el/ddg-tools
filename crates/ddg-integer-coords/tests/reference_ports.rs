//! The reference's own tests of the integer coordinates (geometry-central `ab116f1`,
//! `test/src/intrinsic_triangulation_test.cpp`, MIT, see `LICENSES/`), ported to this crate on
//! the reference's assets `fox.ply` (closed) and `cat_head.obj` (12 boundary edges), SPEC 10.3.
//!
//! Ported: `IntegerFlip`, `DelaunayTriangulationsAgree`, `IntegerTrace`, `IntegerRefine`,
//! `IntegerCommonSubdivision`, `CommonSubdivisionCompareIntegerSignpost`, `FunctionTransfer`.
//! Not ported: the two trace comparisons (`IntegerEdgeTraceAgreesWithBulk` and the signpost
//! variant: the bulk extraction of this crate is tested against single traces in
//! `tests/subdivision.rs`) and `CommonSubdivisionGeometry` (the demo takes its masses from
//! extrinsic positions; here they come from lengths).
//!
//! As in the reference's integer constructor, the lengths are mollified by `1e-5` (and the
//! signposts, to compare, are made from the same mollified lengths).

use ddg_core::io::read_obj;
use ddg_core::mollify::mollify_edge_lengths;
use ddg_core::{Geometry, TriMesh};
use ddg_integer_coords::{
    InitOptions, IntegerRefineOptions, IntegerTriangulation, PAPER_MOLLIFY_EPSILON, Transfer,
};
use ddg_intrinsic::{
    DelaunayCriterion, DelaunayOptions, EdgeTraceOptions, IntrinsicTriangulation, Triangulation,
};
use std::path::PathBuf;

struct Asset {
    name: &'static str,
    mesh: TriMesh,
    lengths: Vec<f64>,
    boundary: bool,
}

fn asset_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/assets")
        .join(name)
}

/// A binary little-endian PLY with float `x y z` and a `list uchar int` of faces.
fn read_ply(name: &str) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let bytes = std::fs::read(asset_path(name)).expect("asset");
    let end = bytes
        .windows(10)
        .position(|w| w == b"end_header")
        .expect("header");
    let header = String::from_utf8_lossy(&bytes[..end]).to_string();
    let count = |what: &str| -> usize {
        header
            .lines()
            .find_map(|l| l.strip_prefix(&format!("element {what} ")))
            .and_then(|n| n.trim().parse().ok())
            .expect("count")
    };
    let (nv, nf) = (count("vertex"), count("face"));
    let mut at = end + 11;
    let f32_at = |at: &mut usize| {
        let v = f32::from_le_bytes(bytes[*at..*at + 4].try_into().unwrap());
        *at += 4;
        f64::from(v)
    };
    let positions: Vec<[f64; 3]> = (0..nv)
        .map(|_| [f32_at(&mut at), f32_at(&mut at), f32_at(&mut at)])
        .collect();
    let mut faces = Vec::new();
    for _ in 0..nf {
        let k = bytes[at] as usize;
        at += 1;
        let mut c = Vec::new();
        for _ in 0..k {
            c.push(i32::from_le_bytes(bytes[at..at + 4].try_into().unwrap()) as u32);
            at += 4;
        }
        for t in 1..k - 1 {
            faces.push([c[0], c[t], c[t + 1]]);
        }
    }
    (positions, faces)
}

fn assets() -> Vec<Asset> {
    let (p, f) = read_ply("fox.ply");
    let mut out = vec![("fox.ply", p, f)];
    let obj = read_obj(asset_path("cat_head.obj")).expect("cat_head.obj");
    let mut faces = Vec::new();
    for i in 0..obj.n_faces() {
        let c = &obj.corner_verts[obj.face_corners(i)];
        for k in 1..c.len() - 1 {
            faces.push([c[0], c[k], c[k + 1]]);
        }
    }
    out.push(("cat_head.obj", obj.positions, faces));
    out.into_iter()
        .map(|(name, p, f)| {
            let mesh = TriMesh::new(p.len(), f).expect("manifold");
            let geometry = Geometry::from_positions(&mesh, &p).expect("geometry");
            let lengths = geometry.edge_lengths().to_vec();
            Asset {
                name,
                boundary: mesh.n_boundary_edges() > 0,
                mesh,
                lengths,
            }
        })
        .collect()
}

fn integer(a: &Asset) -> IntegerTriangulation {
    IntegerTriangulation::with_options(
        &a.mesh,
        &a.lengths,
        &InitOptions {
            mollify: Some(PAPER_MOLLIFY_EPSILON),
        },
    )
    .expect("init")
    .0
}

fn signpost(a: &Asset) -> IntrinsicTriangulation {
    let m = mollify_edge_lengths(&a.mesh, &a.lengths, PAPER_MOLLIFY_EPSILON).unwrap();
    let g = Geometry::from_edge_lengths(&a.mesh, m.edge_lengths).unwrap();
    IntrinsicTriangulation::new(&a.mesh, &g).unwrap()
}

#[test]
fn the_assets_are_those_of_the_reference() {
    let a = assets();
    assert_eq!(
        (a[0].mesh.n_vertices(), a[0].mesh.n_faces(), a[0].boundary),
        (313, 622, false)
    );
    assert_eq!(
        (
            a[1].mesh.n_vertices(),
            a[1].mesh.n_faces(),
            a[1].mesh.n_boundary_edges()
        ),
        (131, 248, 12)
    );
}

/// `IntegerFlip`.
#[test]
fn integer_flip_makes_it_delaunay() {
    for a in assets() {
        let mut t = integer(&a);
        assert!(!t.is_delaunay(&DelaunayCriterion::default()), "{}", a.name);
        let (r, _) = t
            .flip_to_delaunay_certified(&DelaunayOptions::default())
            .unwrap();
        assert_eq!(r.n_not_delaunay, 0, "{}", a.name);
        assert!(t.is_delaunay(&DelaunayCriterion::default()), "{}", a.name);
    }
}

/// `DelaunayTriangulationsAgree`: the lengths by edge index agree to 1e-5.
#[test]
fn integer_and_signpost_delaunay_triangulations_agree() {
    for a in assets() {
        let (mut t, mut s) = (integer(&a), signpost(&a));
        t.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        s.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        assert_eq!(t.complex().n_edges(), s.n_edges());
        let worst = (0..s.n_edges())
            .map(|e| (t.edge_lengths()[e] - s.edge_lengths()[e]).abs())
            .fold(0.0, f64::max);
        assert!(worst < 1e-5, "{}: {worst}", a.name);
        println!("{}: largest length difference {worst:e}", a.name);
    }
}

/// `IntegerTrace`: the trace of every intrinsic edge has `max(0, n) + 2` points.
#[test]
fn integer_trace_has_n_plus_two_points() {
    for a in assets() {
        let mut t = integer(&a);
        t.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        for e in 0..t.complex().n_edges() {
            let trace = t
                .extract_intrinsic_edge(e, &EdgeTraceOptions::default())
                .unwrap();
            let n = t.normal_coordinate(e).max(0) as usize;
            assert_eq!(trace.crossings.len() + 2, n + 2, "{} edge {e}", a.name);
        }
    }
}

/// `IntegerRefine`: Delaunay, more vertices, minimum angle 25 degrees (the reference asserts it
/// over every corner; with a boundary that holds here where the angle sum allows, see the
/// printed values).
#[test]
fn integer_refine_reaches_25_degrees() {
    for a in assets() {
        let mut t = integer(&a);
        let n0 = t.complex().n_vertices();
        let r = t.refine(&IntegerRefineOptions::default()).unwrap();
        assert!(t.is_delaunay(&DelaunayCriterion::default()), "{}", a.name);
        assert!(t.complex().n_vertices() > n0, "{}", a.name);
        println!(
            "{}: {} -> {} vertices, min angle {:.4}, min over valid faces {:.4}, below bound {}",
            a.name,
            n0,
            t.complex().n_vertices(),
            r.refine.min_angle.to_degrees(),
            r.refine.min_angle_valid.to_degrees(),
            r.refine.below_bound.len()
        );
        assert!(
            r.refine.min_angle_valid.to_degrees() >= 25.0 - 1e-9,
            "{}: {}",
            a.name,
            r.refine.min_angle_valid.to_degrees()
        );
        if !a.boundary {
            assert!(r.refine.min_angle.to_degrees() >= 25.0 - 1e-9, "{}", a.name);
        }
        assert!(r.certificate.is_some());
    }
}

/// `IntegerCommonSubdivision`: `V_S = V1 + sum n+`, and the pieces form a surface of the
/// input's Euler characteristic with `V - E + F` counted from the polygons.
#[test]
fn integer_common_subdivision_counts() {
    for a in assets() {
        let mut t = integer(&a);
        t.refine(&IntegerRefineOptions::default()).unwrap();
        let s = t.common_subdivision().unwrap();
        let (v, f) = (s.overlay.vertices.len(), s.overlay.n_faces());
        assert_eq!(
            v as u64,
            t.complex().n_vertices() as u64 + t.total_crossings() as u64
        );
        let mut edges = std::collections::HashSet::new();
        for i in 0..f {
            let c = s.overlay.face(i);
            for k in 0..c.len() {
                let (p, q) = (c[k], c[(k + 1) % c.len()]);
                edges.insert((p.min(q), p.max(q)));
            }
        }
        assert_eq!(
            v as i64 - edges.len() as i64 + f as i64,
            a.mesh.euler_characteristic(),
            "{}",
            a.name
        );
        let (tris, _) = s.overlay.triangles();
        let m = TriMesh::new(v, tris).expect("a manifold mesh");
        assert_eq!(m.n_boundary_edges() > 0, a.boundary, "{}", a.name);
    }
}

/// `CommonSubdivisionCompareIntegerSignpost`: after flips only, equal V/E/F.
#[test]
fn common_subdivisions_of_both_representations_agree_after_flips() {
    for a in assets() {
        let (mut t, mut s) = (integer(&a), signpost(&a));
        t.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        s.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        let ci = t.common_subdivision().unwrap().overlay;
        let cs = s
            .overlay(&s.trace_edges(&EdgeTraceOptions::default()).unwrap())
            .unwrap();
        let count = |o: &ddg_intrinsic::Overlay| {
            let mut e = std::collections::HashSet::new();
            for i in 0..o.n_faces() {
                let c = o.face(i);
                for k in 0..c.len() {
                    let (p, q) = (c[k], c[(k + 1) % c.len()]);
                    e.insert((p.min(q), p.max(q)));
                }
            }
            (o.vertices.len(), e.len(), o.n_faces())
        };
        assert_eq!(count(&ci), count(&cs), "{}", a.name);
    }
}

/// `FunctionTransfer`: after refinement the `L2` transfer is no worse than the pointwise one in
/// the norm of `M_S`, and satisfies its normal equations `M0 x = P0^T M_S P1 b`.
#[test]
fn function_transfer_l2_beats_pointwise() {
    let a = &assets()[0];
    let mut t = integer(a);
    t.refine(&IntegerRefineOptions::default()).unwrap();
    let tr = Transfer::from_state(&t).unwrap();
    let mut rng = ddg_core::generators::SplitMix64::new(7);
    let b: Vec<f64> = (0..t.complex().n_vertices())
        .map(|_| rng.uniform(-1.0, 1.0))
        .collect();
    let pointwise = tr.copy_to_input(&b).unwrap();
    let l2 = tr.to_input(&b).unwrap();
    let ns = tr.mass_subdivision().nrows();
    let apply = |m: &ddg_core::sparse::SparseMatrix, x: &[f64]| {
        let mut y = vec![0.0; m.nrows()];
        m.mul_vec(x, &mut y).unwrap();
        y
    };
    let truth = apply(tr.p1(), &b);
    let err = |x: &[f64]| {
        let r: Vec<f64> = apply(tr.p0(), x)
            .iter()
            .zip(&truth)
            .map(|(p, q)| p - q)
            .collect();
        let mr = apply(tr.mass_subdivision(), &r);
        assert_eq!(r.len(), ns);
        r.iter().zip(&mr).map(|(p, q)| p * q).sum::<f64>()
    };
    let (e_point, e_l2) = (err(&pointwise), err(&l2));
    println!("pointwise {e_point:e}, L2 {e_l2:e}");
    assert!(e_l2 <= e_point);
    // Normal equations: P0^T M_S (P0 x - P1 b) = 0.
    let r: Vec<f64> = apply(tr.p0(), &l2)
        .iter()
        .zip(&truth)
        .map(|(p, q)| p - q)
        .collect();
    let g = apply(&tr.p0().transpose(), &apply(tr.mass_subdivision(), &r));
    let norm = g.iter().map(|x| x * x).sum::<f64>().sqrt();
    assert!(norm < 1e-6, "{norm}");
}
