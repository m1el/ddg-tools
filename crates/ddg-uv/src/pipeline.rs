//! The pipeline: the stages in order, the final checks and the report.

use std::time::Instant;

use ddg_core::io::Cone;

use crate::cones::{ConeStage, solve_cones};
use crate::error::{UvError, internal};
use crate::flatten::{FlatStage, flatten};
use crate::input::MeshInput;
use crate::measure::Accumulator;
use crate::options::UnwrapOptions;
use crate::pack::pack_in_place;
use crate::prepare::{Prepared, prepare};
use crate::report::{ChartReport, ConeStats, Report, SeamLengths, SteeringReport, Timings};

/// A UV map of a polygon mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct UvMap {
    /// One UV per polygon corner, in corner order. Packed into the unit square unless
    /// [`UnwrapOptions::packing`] is `None`.
    pub uvs: Vec<[f64; 2]>,
    /// The chart of every polygon. A chart is a connected piece of the surface after the seams
    /// were cut; charts are numbered in order of their lowest polygon.
    pub polygon_charts: Vec<u32>,
    /// The cones that were used, as vertices of the caller's mesh, ascending.
    pub cones: Vec<Cone>,
    /// The edges the surface was cut along, as vertex pairs `[min, max]`, ascending, each
    /// once: the seams (and what was added to one-edge seams), automatic slits, and the cuts
    /// through the cones and around handles. Every one is an edge of a polygon.
    pub cut_edges: Vec<[u32; 2]>,
    /// What happened.
    pub report: Report,
}

/// Unwraps a polygon mesh: see the [crate documentation](crate) for the stages.
///
/// # Errors
/// [`UvError`]; `kind()` says whether the arguments are malformed, the mesh cannot be
/// unwrapped as it is, or the solve failed. Nothing panics and the input is not modified.
pub fn unwrap(input: &MeshInput<'_>, options: &UnwrapOptions) -> Result<UvMap, UvError> {
    let start = Instant::now();
    let prepared = prepare(input, options)?;
    let cones = solve_cones(&prepared, options)?;
    let flat = flatten(&prepared, &cones, options)?;
    let mut map = assemble(input, &prepared, &cones, &flat, options)?;
    map.report.timings.total = start.elapsed().as_secs_f64();
    Ok(map)
}

/// Stage 4 and the final checks: collects the UVs per polygon corner, verifies them, packs
/// the charts and writes the report. `input` must be the input the stages were computed from.
///
/// The report's total time is the sum of the stage times; [`unwrap`] replaces it by the
/// measured total.
///
/// # Errors
/// [`UvError::Internal`] if a final check fails.
pub fn assemble(
    input: &MeshInput<'_>,
    prepared: &Prepared,
    cones: &ConeStage,
    flat: &FlatStage,
    options: &UnwrapOptions,
) -> Result<UvMap, UvError> {
    let start = Instant::now();
    let n_corners = prepared.n_corners;
    if input.n_corners() != n_corners
        || input.n_polygons() != prepared.n_polygons
        || input.positions.len() != prepared.n_input_vertices
        || cones.pieces.len() != prepared.pieces.len()
        || flat.pieces.len() != prepared.pieces.len()
    {
        return Err(internal(
            "the stages passed to assemble were not computed from this input",
        ));
    }

    // One UV per polygon corner; all triangle corners of a polygon corner must agree.
    let mut uvs = vec![[f64::NAN; 2]; n_corners];
    let mut corner_set = vec![false; n_corners];
    let mut polygon_charts = vec![u32::MAX; prepared.n_polygons];
    let mut all = Accumulator::default();
    let mut cone_list: Vec<Cone> = Vec::new();
    let mut cut_edges: Vec<[u32; 2]> = prepared.seam_cut_edges.clone();
    let mut charts: Vec<ChartReport> = Vec::with_capacity(prepared.pieces.len());
    let (mut overlapped, mut covered, mut any_coverage) = (0.0, 0.0, false);
    for (chart, piece) in prepared.pieces.iter().enumerate() {
        let (piece_cones, piece_flat) = (&cones.pieces[chart], &flat.pieces[chart]);
        let on = piece_cones.slit_piece.as_ref().unwrap_or(piece);
        if piece_flat.corner_uvs.len() != 3 * on.mesh.n_faces() {
            return Err(internal("a flattening has the wrong number of corners"));
        }
        for (f, &t) in on.face_triangle.iter().enumerate() {
            let polygon = prepared.triangle_polygon[t as usize];
            polygon_charts[polygon as usize] = chart as u32;
            let q = &piece_flat.corner_uvs[3 * f..3 * f + 3];
            for (&uv, &corner) in q.iter().zip(&prepared.triangle_corners[t as usize]) {
                let c = corner as usize;
                if !(uv[0].is_finite() && uv[1].is_finite()) {
                    return Err(internal(format!(
                        "the UV ({}, {}) of corner {c} (face {polygon}) is not finite",
                        uv[0], uv[1]
                    )));
                }
                if corner_set[c] && uvs[c] != uv {
                    return Err(internal(format!(
                        "corner {c} of face {polygon} received two different UVs, \
                         ({}, {}) and ({}, {}): a cut runs through the face",
                        uvs[c][0], uvs[c][1], uv[0], uv[1]
                    )));
                }
                uvs[c] = uv;
                corner_set[c] = true;
            }
            all.add(
                on.geometry.face_areas()[f],
                0.5 * ddg_bff::quality::signed_double_area(q[0], q[1], q[2]),
            );
        }
        for cone in &piece_cones.cones {
            cone_list.push(Cone {
                vertex: on.vertex_input[cone.vertex as usize],
                angle: cone.angle,
            });
        }
        cut_edges.extend_from_slice(&piece_cones.slit_edges);
        cut_edges.extend_from_slice(&piece_flat.cut_edges);
        if let Some((o, c)) = piece_flat.coverage {
            overlapped += o;
            covered += c;
            any_coverage = true;
        }
        charts.push(ChartReport {
            first_polygon: piece.first_polygon,
            n_polygons: piece.n_polygons,
            n_triangles: piece.mesh.n_faces(),
            genus: piece.genus(),
            n_boundary_loops: piece.mesh.boundary_loops().len(),
            boundary_optimized: piece_cones.boundary_optimized && piece_flat.boundary_values_used,
            closed: piece_cones.closed,
            positive_cones_only: piece_cones.positive_only,
            solve: piece_cones.solve.clone(),
            spikes: piece_cones.spikes.clone(),
            cones: ConeStats::of(piece_cones.cones.iter().map(|c| c.angle)),
            flatten: piece_flat.report,
        });
    }
    if let Some(c) = corner_set.iter().position(|&set| !set) {
        return Err(internal(format!("corner {c} received no UV")));
    }
    if let Some(p) = polygon_charts.iter().position(|&c| c == u32::MAX) {
        return Err(internal(format!("face {p} is in no chart")));
    }
    cone_list.sort_by_key(|c| c.vertex);
    if let Some(w) = cone_list.windows(2).find(|w| w[0].vertex == w[1].vertex) {
        return Err(internal(format!(
            "vertex {} carries two cones",
            w[0].vertex
        )));
    }
    cut_edges.sort_unstable();
    cut_edges.dedup();
    if let Some(&[a, b]) = cut_edges
        .iter()
        .find(|e| prepared.polygon_edges.binary_search(e).is_err())
    {
        return Err(internal(format!(
            "the cut edge ({a}, {b}) is not an edge of a face"
        )));
    }
    if let Some(&[a, b]) = prepared
        .seam_cut_edges
        .iter()
        .find(|e| cut_edges.binary_search(e).is_err())
    {
        return Err(internal(format!("the seam ({a}, {b}) was not cut")));
    }
    let seam_lengths = seam_lengths(input, &uvs, &cut_edges, &prepared.seam_cut_edges);
    let seconds_verify = start.elapsed().as_secs_f64();

    // Packing.
    let pack_start = Instant::now();
    let packing = options.packing.map(|pack| {
        let mut corner_chart = vec![0u32; n_corners];
        for (p, w) in input.polygon_offsets.windows(2).enumerate() {
            corner_chart[w[0] as usize..w[1] as usize].fill(polygon_charts[p]);
        }
        pack_in_place(&mut uvs, &corner_chart, input.polygon_offsets, &pack)
    });
    if let Some(c) = uvs.iter().position(|q| !(q[0].is_finite() && q[1].is_finite())) {
        return Err(internal(format!("packing made the UV of corner {c} non-finite")));
    }
    let seconds_pack = pack_start.elapsed().as_secs_f64();

    // Steering.
    let mut steering = SteeringReport::default();
    for (v, &role) in prepared.steering_role.iter().enumerate() {
        let applied = prepared.vertex_used[v] && cones.steering_applied[v];
        match (role, applied) {
            (1, true) => steering.n_pinned_applied += 1,
            (1, false) => steering.n_pinned_without_effect += 1,
            (2, true) => steering.n_prescribed_applied += 1,
            (2, false) => steering.prescribed_not_applied.push(v as u32),
            (4, true) => steering.n_forbidden_applied += 1,
            (4, false) => steering.n_forbidden_without_effect += 1,
            _ => {}
        }
    }
    let mut warnings = prepared.warnings.clone();
    warnings.extend(cones.warnings.iter().cloned());
    warnings.extend(flat.warnings.iter().cloned());
    if let Some(&first) = steering.prescribed_not_applied.first() {
        warnings.push(format!(
            "{} prescribed cones could not be imposed because their vertex lies on a boundary, \
             a seam or a slit, or is used by no face (the first is vertex {first})",
            steering.prescribed_not_applied.len()
        ));
    }

    // Results that are formally a map but that nobody should use without looking.
    for (chart, piece_cones) in charts.iter().zip(&cones.pieces) {
        let at_bound = piece_cones
            .cones
            .iter()
            .filter(|c| c.angle >= options.max_cone_angle * (1.0 - 1e-9))
            .count();
        if at_bound > 0 {
            warnings.push(format!(
                "{at_bound} cones of the chart that contains face {} sit at the upper bound \
                 {:.4} of the cone angle (max_cone_angle): the surface has spikes there, or a \
                 wide cluster of cones was merged into one; expect strong distortion, mirrored \
                 triangles or overlap around them",
                chart.first_polygon, options.max_cone_angle
            ));
        }
        let flipped = chart.flatten.flipped_triangles;
        if 100 * flipped > chart.n_triangles {
            warnings.push(format!(
                "{flipped} of the {} triangles of the chart that contains face {} are mirrored \
                 in the UV map",
                chart.n_triangles, chart.first_polygon
            ));
        }
        if let Some(fraction) = chart.flatten.overlapped_fraction
            && fraction > 0.01
        {
            warnings.push(format!(
                "the chart that contains face {} overlaps itself: {:.1} % of its footprint is \
                 covered more than once",
                chart.first_polygon,
                100.0 * fraction
            ));
        }
    }

    let timings = Timings {
        validate: prepared.seconds_validate,
        seams: prepared.seconds_seams,
        cones: cones.seconds,
        cone_solver: cones.seconds_solver,
        spike_check: cones.seconds_spikes,
        flatten: flat.seconds,
        pack: seconds_pack,
        verify: seconds_verify,
        total: prepared.seconds_validate
            + prepared.seconds_seams
            + cones.seconds
            + flat.seconds
            + seconds_pack
            + seconds_verify,
    };
    let max_overlapped_fraction = charts
        .iter()
        .filter_map(|c| c.flatten.overlapped_fraction)
        .fold(None, |best: Option<f64>, x| Some(best.map_or(x, |b| b.max(x))));
    let report = Report {
        n_vertices: prepared.n_input_vertices,
        n_polygons: prepared.n_polygons,
        n_corners,
        n_triangles: input.triangles.len(),
        validation: prepared.validation,
        steering,
        cones: ConeStats::of(cone_list.iter().map(|c| c.angle)),
        n_cut_edges: cut_edges.len(),
        area_distortion: all.finish(),
        flipped_triangles: charts.iter().map(|c| c.flatten.flipped_triangles).sum(),
        degenerate_triangles: charts.iter().map(|c| c.flatten.degenerate_triangles).sum(),
        max_overlapped_fraction,
        overlapped_fraction: (any_coverage && covered > 0.0).then(|| overlapped / covered),
        charts,
        seam_lengths,
        packing,
        timings,
        warnings,
    };
    Ok(UvMap {
        uvs,
        polygon_charts,
        cones: cone_list,
        cut_edges,
        report,
    })
}

/// The mismatch of the lengths of the two images of every cut edge that has a polygon on each
/// side; see [`SeamLengths`]. `uvs` are per polygon corner, `cut_edges` and `seams` sorted
/// vertex pairs `[min, max]`.
fn seam_lengths(
    input: &MeshInput<'_>,
    uvs: &[[f64; 2]],
    cut_edges: &[[u32; 2]],
    seams: &[[u32; 2]],
) -> SeamLengths {
    // (edge, length of one image); each edge of a manifold mesh has at most two.
    let mut images: Vec<([u32; 2], f64)> = Vec::new();
    for w in input.polygon_offsets.windows(2) {
        let (first, end) = (w[0] as usize, w[1] as usize);
        for c in first..end {
            let d = if c + 1 == end { first } else { c + 1 };
            let (a, b) = (input.corner_vertices[c], input.corner_vertices[d]);
            let key = [a.min(b), a.max(b)];
            if cut_edges.binary_search(&key).is_ok() {
                let (p, q) = (uvs[c], uvs[d]);
                images.push((key, (p[0] - q[0]).hypot(p[1] - q[1])));
            }
        }
    }
    images.sort_by(|x, y| x.0.cmp(&y.0).then(x.1.total_cmp(&y.1)));
    let mut out = SeamLengths::default();
    let mut seam_sum = 0.0;
    for pair in images.windows(2).filter(|w| w[0].0 == w[1].0) {
        let (edge, (la, lb)) = (pair[0].0, (pair[0].1, pair[1].1));
        let largest = la.max(lb);
        let mismatch = if largest > 0.0 { (la - lb).abs() / largest } else { 0.0 };
        if seams.binary_search(&edge).is_ok() {
            out.n_seam_edges += 1;
            seam_sum += mismatch;
            if out.worst_seam_edge.is_none() || mismatch > out.max_seam_mismatch {
                out.max_seam_mismatch = mismatch;
                out.worst_seam_edge = Some(edge);
            }
        } else {
            out.n_other_cut_edges += 1;
            out.max_other_mismatch = out.max_other_mismatch.max(mismatch);
        }
    }
    if out.n_seam_edges > 0 {
        out.mean_seam_mismatch = seam_sum / out.n_seam_edges as f64;
    }
    out
}
