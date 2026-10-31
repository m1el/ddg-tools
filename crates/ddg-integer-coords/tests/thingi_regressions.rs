//! Failures found by the Thingi10k evaluation (WP10). The meshes are `data/thingi10k/converted/`
//! (produced by `examples/thingi_convert.rs`); the tests skip with a notice if they are absent and
//! are `#[ignore]`d while they fail. Run with `-- --ignored` (release mode).
//!
//! **Refinement does not terminate on 176193 and 190253** (both ~16 000 faces, closed, minimum
//! input angle ~1e-5 degrees). Chew's rules, 25 degrees, mollified 1e-5: the authors' demo
//! (`--backend integer --refineDelaunay --refineMaxInsertions=0`) finishes in 0.6 s / 1.4 s with
//! 17 871 / 37 078 vertices and minimum angle 25.04 / 25.01 degrees. Here the driver keeps
//! inserting (162 790 insertions at 10 per vertex, minimum angle over valid faces still
//! 6.7e-6 degrees, 18 671 corners below the bound, none of them for an angle sum or a needle
//! vertex) and runs past the 10-minute cap without a limit. The signpost representation does
//! exactly the same (equal numbers of insertions, 8 trace failures), so the cause is in the
//! shared driver `ddg_intrinsic::refine`, not in the integers.

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
        "mesh {id}: {} inserted, min angle {:.4}",
        r.refine.n_inserted,
        r.refine.min_angle_valid.to_degrees()
    );
    assert_eq!(
        r.refine.status,
        RefineStatus::Done,
        "mesh {id}: insertion limit reached with {} inserted, minimum angle {} degrees",
        r.refine.n_inserted,
        r.refine.min_angle_valid.to_degrees()
    );
    assert!(r.refine.min_angle_valid.to_degrees() >= 25.0 - 1e-9);
}

#[test]
#[ignore = "WP10 finding: refinement does not terminate (see the module documentation)"]
fn refinement_of_thingi_176193_terminates() {
    refine_to_25("176193", true);
}

#[test]
#[ignore = "WP10 finding: refinement does not terminate (see the module documentation)"]
fn refinement_of_thingi_190253_terminates() {
    refine_to_25("190253", false);
}

#[test]
#[ignore = "WP10 finding: refinement does not terminate (see the module documentation)"]
fn refinement_of_thingi_190253_after_idt_terminates() {
    refine_to_25("190253", true);
}
