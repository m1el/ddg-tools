//! Failures found by the Thingi10k evaluation (WP10). The meshes are `data/thingi10k/converted/`
//! (produced by `examples/thingi_convert.rs`); the tests skip with a notice if they are absent.
//! Run in release mode.
//!
//! **Refinement did not terminate on 176193 and 190253** (IC-SPEC I30; both ~16 000 faces,
//! closed, minimum input angle ~1e-5 degrees). Chew's rules, 25 degrees, mollified 1e-5: the
//! authors' demo (`--backend integer --refineDelaunay --refineMaxInsertions=0`) finishes in
//! 0.6 s / 1.4 s with 17 871 / 37 078 vertices and minimum angle 25.04 / 25.01 degrees. The
//! shared driver `ddg_intrinsic::refine` kept inserting (162 790 insertions at 10 per vertex,
//! minimum angle 6.7e-6 degrees) because Delaunay flips that the shape limit of IT-SPEC N5
//! refused (new triangles with shape numbers ~2e-7, made by the mollified near-duplicate
//! vertices) left edges that are not Delaunay; the circumcentres of their triangles landed next
//! to existing vertices (a quarter of all insertions at less than the circumradius, down to
//! 6e-4 of it) and every such insertion made new short edges. Refinement now makes those flips
//! (`ddg_intrinsic::refine::DELAUNAY_FLIP_SHAPE`, as the reference's `flipEdgeIfNotDelaunay`):
//! 1595 / 21 695 insertions, 25.04 / 25.01 degrees.

use ddg_core::{TriMesh, io::read_obj};
use ddg_integer_coords::{
    InitOptions, IntegerRefineOptions, IntegerTriangulation, PAPER_MOLLIFY_EPSILON,
};
use ddg_intrinsic::{InsertionLimit, RefineOptions, RefineStatus};

fn load(id: &str) -> Option<(TriMesh, Vec<f64>)> {
    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../data/thingi10k/converted")
        .join(format!("{id}.obj"));
    if !path.is_file() {
        eprintln!("SKIPPED: {} not found", path.display());
        return None;
    }
    let obj = read_obj(&path).ok()?;
    let mesh = TriMesh::new(obj.positions.len(), obj.triangles()?).ok()?;
    let lengths = mesh
        .edges()
        .iter()
        .map(|&[a, b]| {
            let (p, q) = (obj.positions[a as usize], obj.positions[b as usize]);
            ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt()
        })
        .collect();
    Some((mesh, lengths))
}

/// `flip_first`: also calls `flip_to_delaunay` before `refine`, as the evaluation does.
fn refine_to_25(id: &str, flip_first: bool) {
    let Some((mesh, lengths)) = load(id) else {
        return;
    };
    let (mut t, _) = IntegerTriangulation::with_options(
        &mesh,
        &lengths,
        &InitOptions {
            mollify: Some(PAPER_MOLLIFY_EPSILON),
        },
    )
    .unwrap();
    if flip_first {
        t.flip_to_delaunay_certified(&ddg_intrinsic::DelaunayOptions::default())
            .unwrap();
    }
    // The demo needs under 25 000 insertions on both meshes; allow 60 000.
    let options = IntegerRefineOptions {
        refine: RefineOptions {
            max_insertions: InsertionLimit::Vertices(60_000),
            ..RefineOptions::default()
        },
        ..IntegerRefineOptions::default()
    };
    let r = t.refine(&options).unwrap();
    eprintln!(
        "mesh {id}: {} inserted, min angle {:.4}, {} Delaunay flips below the shape limit",
        r.refine.n_inserted,
        r.refine.min_angle_valid.to_degrees(),
        r.refine.n_thin_flips
    );
    assert_eq!(
        r.refine.status,
        RefineStatus::Done,
        "mesh {id}: insertion limit reached with {} inserted, minimum angle {} degrees",
        r.refine.n_inserted,
        r.refine.min_angle_valid.to_degrees()
    );
    assert!(r.refine.min_angle_valid.to_degrees() >= 25.0 - 1e-9);
    assert_eq!(r.refine.n_not_delaunay, 0, "mesh {id}");
    assert!(
        r.refine.n_thin_flips > 0,
        "mesh {id}: the case of I30 was not met"
    );
    assert!(r.certificate.is_some(), "mesh {id}: not certified");
}

#[test]
fn refinement_of_thingi_176193_terminates() {
    refine_to_25("176193", true);
}

#[test]
fn refinement_of_thingi_190253_terminates() {
    refine_to_25("190253", false);
}

#[test]
fn refinement_of_thingi_190253_after_idt_terminates() {
    refine_to_25("190253", true);
}
