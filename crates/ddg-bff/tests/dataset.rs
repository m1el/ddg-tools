//! The optimal-cones dataset: real output of the reference BFF.
//!
//! Every `<name>.obj` stores 3D positions and, per corner, the UVs that the reference
//! produced for the cones in `<name>-cones.txt`. Two comparisons:
//!
//! 1. **Same cut, same map.** With the cut recovered from the stored UVs, our layout must
//!    agree with the stored one up to a similarity per chart.
//! 2. **Own cut, same metric.** With our own cut the layout differs but its edge lengths must
//!    not, up to one global scale.
//!
//! How well they can agree is limited by the files, not by the port: positions, angles and
//! UVs are stored with six digits, and the stored maps are not exactly what the public
//! reference code computes from the published cones (the reference, run today, misses each
//! stored map by as much as this crate does; see the README). The bounds below are a few
//! times the measured values. The port itself agrees with the reference C++ code to
//! rounding error when both get the same input (`examples/reference_diff.rs`,
//! `reference/run.sh`).
//!
//! The dataset is not part of the repository; without it the tests print a notice and pass.
//! A debug build checks a few small meshes, a release build all 23:
//! `cargo test --release -p ddg-bff --test dataset -- --nocapture`.

mod common;

use common::{
    compare_layouts, dataset_names, length_ratio_quantiles, load, stored_cone_angles, total_length,
};
use ddg_bff::{
    Cone, ConeFlattenInput, ConeFlattenOptions, ConeFlattening, CutSource, GaussBonnetPolicy,
    ScaleFactorOptions, component_topology, flatten_with_cones,
};

/// The meshes to check: all of them in a release build, a few small ones otherwise.
fn names() -> Option<Vec<String>> {
    let Some(mut names) = dataset_names() else {
        eprintln!(
            "SKIPPED: dataset directory {} not found; nothing was checked",
            common::dataset_dir().display()
        );
        return None;
    };
    assert_eq!(names.len(), 23, "expected the 23 meshes of the dataset");
    if cfg!(debug_assertions) {
        names.retain(|n| ["airplane", "axes", "bump", "hand"].contains(&n.as_str()));
    }
    Some(names)
}

fn flatten(data: &common::DatasetMesh, cones: &[Cone], cut: CutSource<'_>) -> ConeFlattening {
    flatten_with_cones(
        &ConeFlattenInput {
            mesh: &data.mesh,
            geometry: &data.geometry,
            cones,
            scale_factors: None,
            cut,
            anchors: &[],
        },
        &ConeFlattenOptions {
            // The published angles of the closed meshes miss Gauss-Bonnet by up to 3e-3.
            scale_factors: ScaleFactorOptions {
                gauss_bonnet: GaussBonnetPolicy::SpreadOverVertices,
                ..ScaleFactorOptions::default()
            },
            ..ConeFlattenOptions::default()
        },
    )
    .unwrap_or_else(|e| panic!("{}: {e}", data.name))
}

/// Upper bound on the root mean square distance to the stored map, relative to the diagonal
/// of its bounding box, after the best similarity. Measured values are in the README.
fn same_cut_bound(name: &str) -> f64 {
    match name {
        // Agreement to the precision of the file.
        "axes" | "bump" | "camel_noise2" | "face" => 1e-5,
        // The stored map has a cone angle that is not the published one (see below).
        "airplane" => 1.5e-2,
        "camelshotnoise" => 1.2e-2,
        // The stored map is not reproduced on any cut, by this crate or by the reference; it
        // turns by 0.105 more than the published cones allow (summed over the vertices of its
        // cut that are not cones). The stored cut itself is a valid cut through the published
        // cones: the mirror image of the one the reference chooses today.
        "body" => 2e-2,
        // The stored map does not have a constant scale along the original boundary (the
        // ratio of lengths there varies by 0.37 to 1.11 and 0.83 to 1.18; ours and the
        // reference's by a few percent and a few tenths of a percent): it was computed with
        // other boundary data. Which, is not known.
        "camel_noise1" => 5e-3,
        "rhino" => 2e-3,
        // Closed; its published cone angles miss Gauss-Bonnet by 2.7e-3.
        "brain" => 3e-4,
        // One cut edge leads to a cone whose published angle is 0.0 and cannot be seen in
        // the stored UVs; the planner completes the cut.
        "cow00001" => 5e-4,
        _ => 1e-4,
    }
}

#[test]
fn same_cut_gives_the_stored_map() {
    let Some(names) = names() else { return };
    println!("{:<20} {:>9} {:>9} {:>9}", "mesh", "rms", "max", "bound");
    for name in &names {
        let data = load(name);
        let (cut, added) = data.completed_reference_cut();
        if name != "cow00001" {
            assert_eq!(added, 0, "{name}: the stored cut is complete");
        }
        let flat = flatten(&data, &data.cones, CutSource::Edges(&cut));
        assert_eq!(flat.n_charts(), 1, "{name}");
        assert!(component_topology(&flat.cut.mesh)[0].is_disk());
        let agreement = compare_layouts(
            &flat.face_chart,
            flat.n_charts(),
            &flat.corner_uvs,
            &data.reference_uvs,
        );
        let bound = same_cut_bound(name);
        println!(
            "{:<20} {:9.1e} {:9.1e} {:9.1e}",
            name, agreement.rms, agreement.max, bound
        );
        assert!(
            agreement.rms < bound,
            "{name}: rms {} >= {bound}",
            agreement.rms
        );
        // Seamless, and the scale factors imply the prescribed cone angles (on a closed mesh
        // up to the share of the Gauss-Bonnet discrepancy that every vertex receives).
        assert!(flat.seams.max_length_mismatch < 1e-10, "{name}");
        for cone in &flat.cones {
            assert!(
                (cone.implied_angle - cone.angle).abs() < 1e-6,
                "{name}: {cone:?}"
            );
        }
    }
}

/// In two meshes the stored map realizes, at one or two cones, an angle that is not the one in
/// the cone file. With the angle of the stored map the layouts agree again.
#[test]
fn stored_maps_with_other_cone_angles_than_published() {
    let Some(names) = names() else { return };
    for (name, vertex, difference) in [
        ("airplane", 2145u32, -0.110),
        ("camelshotnoise", 854, -0.259),
    ] {
        if !names.iter().any(|n| n == name) {
            continue;
        }
        let data = load(name);
        let stored = stored_cone_angles(&data);
        let mut cones = data.cones.clone();
        let mut changed = Vec::new();
        for (cone, &angle) in cones.iter_mut().zip(&stored) {
            if (angle - cone.angle).abs() > 1e-4 {
                changed.push((cone.vertex, angle - cone.angle));
                cone.angle = angle;
            } else {
                // Everywhere else the stored map has the published angle, to its six digits.
                assert!(
                    (angle - cone.angle).abs() < 6e-6,
                    "{name}: {cone:?} vs {angle}"
                );
            }
        }
        assert!(
            changed
                .iter()
                .any(|&(v, d)| v == vertex && (d - difference).abs() < 2e-3),
            "{name}: {changed:?}"
        );
        let cut = CutSource::Edges(&data.reference_cut);
        let residual = |cones: &[Cone]| {
            let flat = flatten(&data, cones, cut);
            compare_layouts(&flat.face_chart, 1, &flat.corner_uvs, &data.reference_uvs).rms
        };
        let (published, corrected) = (residual(&data.cones), residual(&cones));
        println!(
            "{name}: rms {published:.1e} with the published angles, {corrected:.1e} with the stored ones"
        );
        assert!(
            published > 3e-3 && corrected < 3e-4,
            "{name}: {published} -> {corrected}"
        );
    }
}

#[test]
fn own_cut_gives_the_same_metric() {
    let Some(names) = names() else { return };
    println!(
        "{:<20} {:>7} {:>7} {:>7} {:>8} {:>8} {:>8} {:>8}",
        "mesh", "refE", "ownE", "len/ref", "min", "p1", "p99", "max"
    );
    for name in &names {
        let data = load(name);
        let flat = flatten(
            &data,
            &data.cones,
            CutSource::Plan {
                may_cut: None,
                must_cut: &[],
            },
        );
        let plan = flat.plan.as_ref().unwrap();
        // One disk, all cones on the cut, seamless.
        assert_eq!(flat.n_charts(), 1, "{name}");
        assert!(component_topology(&flat.cut.mesh)[0].is_disk(), "{name}");
        assert!(flat.cones.iter().all(|c| c.n_wedges >= 1), "{name}");
        assert!(flat.seams.max_length_mismatch < 1e-10, "{name}");
        assert!(flat.max_interior_curvature < 1e-6, "{name}");
        assert_eq!(flat.charts[0].closure.dropped_directions, 0, "{name}");
        // The planner's cut is not longer than the reference's.
        let ratio = total_length(&data.geometry, &plan.cut_edges)
            / total_length(&data.geometry, &data.reference_cut);
        assert!(ratio < 1.02, "{name}: cut length ratio {ratio}");
        // Edge lengths agree with the stored map up to one global factor.
        let q = length_ratio_quantiles(&flat.corner_uvs, &data.reference_uvs);
        println!(
            "{:<20} {:>7} {:>7} {:>7.3} {:>8.4} {:>8.4} {:>8.4} {:>8.4}",
            name,
            data.reference_cut.len(),
            plan.cut_edges.len(),
            ratio,
            q[0],
            q[1],
            q[3],
            q[4]
        );
        // The meshes whose cones are all positive are flattened without any overlap: no
        // flipped triangle, no crossing of the boundary, and no sample point covered twice.
        let positive = data.cones.iter().all(|c| c.angle > 0.0);
        if positive {
            assert!(
                [
                    "axes", "bump", "bunny", "cow1", "cow5e-1", "cow75e-2", "face"
                ]
                .contains(&name.as_str()),
                "{name}"
            );
            assert!(flat.charts[0].is_injective(), "{name}");
            let coverage = flat.charts[0].overlap.coverage.expect("estimated");
            assert_eq!(coverage.overlapped_area, 0.0, "{name}");
        }
        // 98 % of all edges within these factors; wider where the stored map is known to
        // come from other data, or has many flipped triangles.
        let (low, high) = match name.as_str() {
            "body" => (0.9, 1.4),
            "airplane" | "camelshotnoise" | "camel_noise1" => (0.85, 1.2),
            "axes" | "bump" | "bunny" | "chair" | "cow00001" | "cow01" | "cow5e-1" | "face"
            | "mannequin" | "tree" | "brain" => (0.99, 1.01),
            _ => (0.95, 1.06),
        };
        assert!(q[1] > low && q[3] < high, "{name}: {q:?}");
    }
}
