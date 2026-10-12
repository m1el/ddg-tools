//! Stage 3a: the cones of every piece.

use std::f64::consts::TAU;
use std::time::Instant;

use ddg_core::io::Cone;
use ddg_core::operators::angle_defects;
use ddg_optimal_cones::{ConeError, ConeProblem, SolveOptions, Status, solve};

use crate::error::{UvError, internal};
use crate::options::{ClosedSurfaces, ConeSigns, MergeClusters, SlitPlacement, UnwrapOptions};
use crate::prepare::{Piece, Prepared};
use crate::report::{ClosedHandling, ConeSolveReport, MergeOutcome};
use crate::slit::{choose_slit, distances};

/// The cones of one piece.
#[derive(Debug, Clone)]
pub struct PieceCones {
    /// The piece with its automatic slit, if one was cut; the cones and the scale factors then
    /// refer to this mesh instead of the prepared piece.
    pub slit_piece: Option<Piece>,
    /// The edges of the slit as vertex pairs of the caller, `[min, max]`.
    pub slit_edges: Vec<[u32; 2]>,
    /// The cones, as vertices of the piece's mesh, ascending.
    pub cones: Vec<Cone>,
    /// The log scale factor of these cones per vertex of the piece, or `None` if the cone
    /// solver did not run (the flattener then solves for it).
    pub log_scale: Option<Vec<f64>>,
    /// How the piece was handled if it is closed.
    pub closed: ClosedHandling,
    /// Whether cone angles were restricted to positive values.
    pub positive_only: bool,
    /// The cone solve, if it ran.
    pub solve: Option<ConeSolveReport>,
}

/// The result of stage 3a.
#[derive(Debug, Clone)]
pub struct ConeStage {
    /// One entry per piece of the [`Prepared`] mesh.
    pub pieces: Vec<PieceCones>,
    /// Per vertex of the caller: whether its steering was applied (it is an interior vertex of
    /// the piece the cones were computed on).
    pub(crate) steering_applied: Vec<bool>,
    pub(crate) warnings: Vec<String>,
    pub(crate) seconds: f64,
    pub(crate) seconds_solver: f64,
}

/// What one solve on one mesh gave.
struct Solved {
    cones: Vec<Cone>,
    log_scale: Option<Vec<f64>>,
    report: Option<ConeSolveReport>,
    positive_only: bool,
    /// Number of free sites the pipeline added.
    n_free_sites: usize,
    n_pinned_sites: usize,
}

/// Runs stage 3a on every piece. See [`crate::unwrap`].
///
/// # Errors
/// [`UvError::ClosedSurface`], [`UvError::NegativeConesRequired`],
/// [`UvError::InfeasibleSteering`], [`UvError::ConeSolver`], [`UvError::ConesNotCertified`],
/// [`UvError::UnrealizableCone`].
pub fn solve_cones(prepared: &Prepared, options: &UnwrapOptions) -> Result<ConeStage, UvError> {
    let start = Instant::now();
    options.validate()?;
    let mut stage = ConeStage {
        pieces: Vec::with_capacity(prepared.pieces.len()),
        steering_applied: vec![false; prepared.n_input_vertices],
        warnings: Vec::new(),
        seconds: 0.0,
        seconds_solver: 0.0,
    };
    for piece in &prepared.pieces {
        let cones = piece_cones(prepared, piece, options, &mut stage)?;
        stage.pieces.push(cones);
    }
    stage.seconds_solver = stage
        .pieces
        .iter()
        .filter_map(|p| p.solve.as_ref())
        .map(|s| s.seconds)
        .sum();
    stage.seconds = start.elapsed().as_secs_f64();
    Ok(stage)
}

fn piece_cones(
    prepared: &Prepared,
    piece: &Piece,
    options: &UnwrapOptions,
    stage: &mut ConeStage,
) -> Result<PieceCones, UvError> {
    let closed = piece.is_closed();
    let genus = piece.genus();
    if closed && options.closed_surfaces == ClosedSurfaces::Refuse {
        return Err(UvError::ClosedSurface {
            polygon: piece.first_polygon,
            n_polygons: piece.n_polygons,
            genus,
        });
    }
    let sphere = closed && genus == 0;
    let slit_of = |placement: SlitPlacement, relative_length: f64| {
        // Keep the slit away from prescribed cones: a vertex on the slit is on the boundary.
        let avoid: Vec<bool> = piece
            .vertex_input
            .iter()
            .map(|&v| prepared.steering_role[v as usize] & 2 != 0)
            .collect();
        let slit = choose_slit(piece, &avoid, placement, relative_length)?;
        let cut = piece.cut(&slit.edges)?;
        let mut edges: Vec<[u32; 2]> = slit
            .edges
            .iter()
            .map(|&e| {
                let [a, b] = piece.mesh.edge_vertices(e as usize);
                let (a, b) = (piece.vertex_input[a], piece.vertex_input[b]);
                [a.min(b), a.max(b)]
            })
            .collect();
        edges.sort_unstable();
        Ok::<_, UvError>((cut, edges, slit.length))
    };

    let mut warnings = Vec::new();
    let mut result = match (sphere, options.closed_surfaces) {
        (true, ClosedSurfaces::Slit {
            placement,
            relative_length,
        }) => {
            let (cut, edges, length) = slit_of(placement, relative_length)?;
            let solved = solve_on(prepared, &cut, options, None, &mut warnings)?;
            finish(
                solved,
                Some(cut),
                edges.clone(),
                ClosedHandling::Slit {
                    n_edges: edges.len(),
                    length,
                },
            )
        }
        (true, ClosedSurfaces::Auto {
            max_cones,
            max_cone_fraction,
            fallback_sites,
        }) => {
            let mut raw_warnings = Vec::new();
            let raw = solve_on(prepared, piece, options, None, &mut raw_warnings)?;
            let sparse = raw.cones.len() as f64
                <= max_cone_fraction * piece.mesh.n_vertices() as f64;
            if raw.cones.len() <= max_cones || sparse {
                warnings.append(&mut raw_warnings);
                finish(raw, None, Vec::new(), ClosedHandling::RawOptimum)
            } else {
                let n_raw = raw.cones.len();
                let raw_seconds = raw.report.as_ref().map_or(0.0, |r| r.seconds);
                let sites = Some((fallback_sites, options.max_cone_angle, true));
                let mut solved = solve_on(prepared, piece, options, sites, &mut warnings)?;
                if let Some(report) = &mut solved.report {
                    report.seconds += raw_seconds;
                }
                let handling = ClosedHandling::FreeSites {
                    n_automatic: solved.n_free_sites,
                    n_pinned: solved.n_pinned_sites,
                    fallback_from_cones: Some(n_raw),
                };
                finish(solved, None, Vec::new(), handling)
            }
        }
        (true, ClosedSurfaces::FreeSites {
            count,
            max_angle,
            forbid_elsewhere,
        }) => {
            let sites = Some((count, max_angle, forbid_elsewhere));
            let solved = solve_on(prepared, piece, options, sites, &mut warnings)?;
            let handling = ClosedHandling::FreeSites {
                n_automatic: solved.n_free_sites,
                n_pinned: solved.n_pinned_sites,
                fallback_from_cones: None,
            };
            finish(solved, None, Vec::new(), handling)
        }
        (true, _) | (false, _) => {
            let solved = solve_on(prepared, piece, options, None, &mut warnings)?;
            let handling = if closed {
                ClosedHandling::RawOptimum
            } else {
                ClosedHandling::NotClosed
            };
            finish(solved, None, Vec::new(), handling)
        }
    };

    // Which steering took effect, and whether every cone can exist.
    let on = result.slit_piece.as_ref().unwrap_or(piece);
    for v in 0..on.mesh.n_vertices() {
        if !on.mesh.is_boundary_vertex(v) {
            stage.steering_applied[on.vertex_input[v] as usize] = true;
        }
    }
    if let Some(cone) = result.cones.iter().find(|c| c.angle >= TAU) {
        return Err(UvError::UnrealizableCone {
            vertex: on.vertex_input[cone.vertex as usize],
            angle: cone.angle,
        });
    }
    result.cones.sort_by_key(|c| c.vertex);
    stage.warnings.append(&mut warnings);
    Ok(result)
}

fn finish(
    solved: Solved,
    slit_piece: Option<Piece>,
    slit_edges: Vec<[u32; 2]>,
    closed: ClosedHandling,
) -> PieceCones {
    PieceCones {
        slit_piece,
        slit_edges,
        cones: solved.cones,
        log_scale: solved.log_scale,
        closed,
        positive_only: solved.positive_only,
        solve: solved.report,
    }
}

/// Solves the cone problem of one mesh (a piece, or a piece with its slit). `free_sites` is
/// `(count, angle bound, forbid cones elsewhere)` if the policy asks for free cone sites.
fn solve_on(
    prepared: &Prepared,
    piece: &Piece,
    options: &UnwrapOptions,
    free_sites: Option<(usize, f64, bool)>,
    warnings: &mut Vec<String>,
) -> Result<Solved, UvError> {
    let mesh = &piece.mesh;
    let n = mesh.n_vertices();
    let closed = piece.is_closed();
    let genus = piece.genus();
    let positive_only = match options.cone_signs {
        ConeSigns::Any => false,
        ConeSigns::PositiveOnly => {
            if closed && genus >= 2 {
                return Err(UvError::NegativeConesRequired {
                    polygon: piece.first_polygon,
                    genus,
                });
            }
            true
        }
        ConeSigns::PositiveWhereFeasible => {
            if closed && genus >= 2 {
                warnings.push(format!(
                    "the piece of the surface that contains face {} is closed with genus \
                     {genus}: its cone angles must sum to a negative number, so negative cones \
                     are allowed on it; expect overlap, or mark seams",
                    piece.first_polygon
                ));
                false
            } else {
                true
            }
        }
    };
    let mut solved = Solved {
        cones: Vec::new(),
        log_scale: None,
        report: None,
        positive_only,
        n_free_sites: 0,
        n_pinned_sites: 0,
    };
    let interior: Vec<bool> = (0..n).map(|v| !mesh.is_boundary_vertex(v)).collect();
    if !interior.iter().any(|&i| i) {
        return Ok(solved);
    }

    // Steering, at interior vertices only.
    let role = |v: usize| prepared.steering_role[piece.vertex_input[v] as usize];
    let mut costs = vec![1.0; n];
    let mut any_cost = false;
    let mut prescribed = Vec::new();
    for v in (0..n).filter(|&v| interior[v]) {
        match role(v) {
            1 => {
                costs[v] = 0.0;
                any_cost = true;
                solved.n_pinned_sites += 1;
            }
            2 => prescribed.push(Cone {
                vertex: v as u32,
                angle: prepared.prescribed_angle[piece.vertex_input[v] as usize],
            }),
            4 => {
                costs[v] = f64::INFINITY;
                any_cost = true;
            }
            _ => {}
        }
    }
    // A closed piece of genus 1 with positive cones: the angles sum to zero, so there is no
    // cone, and nothing to solve.
    if closed && genus == 1 && positive_only && prescribed.is_empty() {
        return Ok(solved);
    }
    let mut upper = vec![options.max_cone_angle; n];
    if let Some((count, max_angle, forbid_elsewhere)) = free_sites {
        // The pinned vertices are sites already; more are added by farthest-point sampling,
        // starting at the vertex of largest absolute angle defect.
        let mut is_site: Vec<bool> = (0..n).map(|v| interior[v] && role(v) == 1).collect();
        let mut nearest = vec![f64::INFINITY; n];
        let nothing_blocked = vec![false; n];
        for v in (0..n).filter(|&v| is_site[v]) {
            for (d, x) in nearest.iter_mut().zip(distances(piece, &nothing_blocked, v)) {
                *d = d.min(x);
            }
        }
        let defects = angle_defects(mesh, &piece.geometry)
            .map_err(|e| internal(format!("angle defects of a piece: {e}")))?;
        while solved.n_pinned_sites + solved.n_free_sites < count {
            let candidates = (0..n).filter(|&v| interior[v] && role(v) == 0 && !is_site[v]);
            let mut site: Option<usize> = None;
            for v in candidates {
                let better = match site {
                    None => true,
                    // The first site: largest curvature. Later ones: farthest from all sites.
                    Some(s) if nearest[s].is_infinite() && nearest[v].is_infinite() => {
                        defects[v].abs() > defects[s].abs()
                    }
                    Some(s) => nearest[v] > nearest[s],
                };
                if better {
                    site = Some(v);
                }
            }
            let Some(site) = site else { break };
            is_site[site] = true;
            solved.n_free_sites += 1;
            for (d, x) in nearest.iter_mut().zip(distances(piece, &nothing_blocked, site)) {
                *d = d.min(x);
            }
        }
        for v in (0..n).filter(|&v| interior[v]) {
            if is_site[v] {
                costs[v] = 0.0;
                upper[v] = max_angle;
                any_cost = true;
            } else if forbid_elsewhere && role(v) != 2 {
                costs[v] = f64::INFINITY;
                any_cost = true;
            }
        }
    }

    let (merge, max_cluster_size) = match options.merge {
        MergeClusters::Never => (false, 0),
        MergeClusters::Always { max_cluster_size }
        | MergeClusters::Guarded {
            max_cluster_size, ..
        } => (true, max_cluster_size),
    };
    let problem = ConeProblem {
        mesh,
        geometry: &piece.geometry,
        lambda: options.lambda,
        distortion_weights: None,
        cone_costs: any_cost.then(|| costs.clone()),
        lower_bounds: positive_only.then(|| vec![0.0; n]),
        upper_bounds: Some(upper),
        prescribed_cones: prescribed,
    };
    let defaults = SolveOptions::default();
    let cone_threshold = defaults.cone_threshold;
    let solve_options = SolveOptions {
        merge_cone_clusters: merge,
        max_merge_cluster_size: if merge {
            max_cluster_size
        } else {
            defaults.max_merge_cluster_size
        },
        ..defaults
    };
    let solution = solve(&problem, &solve_options).map_err(|e| match e {
        ConeError::Infeasible { .. } => UvError::InfeasibleSteering {
            polygon: piece.first_polygon,
            message: e.to_string(),
        },
        other => UvError::ConeSolver {
            polygon: piece.first_polygon,
            message: other.to_string(),
        },
    })?;

    let certificate = &solution.report.certificate;
    let status = solution.report.status;
    if status != Status::Converged {
        let message = format!(
            "status {status:?}{}",
            solution
                .report
                .message
                .as_ref()
                .map_or(String::new(), |m| format!(" ({m})"))
        );
        if options.require_certified_cones {
            return Err(UvError::ConesNotCertified {
                polygon: piece.first_polygon,
                message,
            });
        }
        warnings.push(format!(
            "the cones of the piece of the surface that contains face {} are not a certified \
             optimum: {message}",
            piece.first_polygon
        ));
    }

    // The merge guard.
    let before = solution.distortion_optimal;
    let after = solution.distortion;
    let anything_merged = !solution.extraction.merged.is_empty();
    let outcome = if !merge {
        MergeOutcome::NotRequested
    } else if !anything_merged {
        MergeOutcome::NothingToMerge
    } else {
        match options.merge {
            MergeClusters::Guarded {
                max_factor,
                max_increase,
                ..
            } if !(after <= (max_factor * before).max(before + max_increase)) => {
                MergeOutcome::Rejected
            }
            _ => MergeOutcome::Merged,
        }
    };
    let n_cones_optimal = solution
        .measure
        .iter()
        .filter(|m| m.abs() > cone_threshold)
        .count();
    if outcome == MergeOutcome::Rejected {
        warnings.push(format!(
            "merging the cone clusters of the piece of the surface that contains face {} would \
             raise the distortion from {before:.4e} to {after:.4e}; the {n_cones_optimal} \
             unmerged cones are used instead of the {} merged ones",
            piece.first_polygon,
            solution.cones.len()
        ));
        solved.cones = solution
            .measure
            .iter()
            .enumerate()
            .filter(|(_, m)| m.abs() > cone_threshold)
            .map(|(v, &angle)| Cone {
                vertex: v as u32,
                angle,
            })
            .collect();
        solved.log_scale = Some(solution.log_scale_optimal.clone());
    } else {
        solved.cones = solution.cones.clone();
        solved.log_scale = Some(solution.log_scale.clone());
    }
    let linear = &solution.report.linear_solver;
    solved.report = Some(ConeSolveReport {
        status,
        max_residual: certificate
            .state_residual
            .max(certificate.adjoint_residual)
            .max(certificate.cone_residual)
            .max(certificate.bound_violation)
            .max(certificate.dual_infeasibility)
            .max(certificate.gauss_bonnet_residual),
        duality_gap: certificate.duality_gap,
        complementarity_gap: certificate.complementarity_gap,
        gap_uncertainty: certificate.gap_uncertainty,
        objective: solution.objective,
        newton_iterations: [
            solution.report.regularized_iterations,
            solution.report.exact_iterations,
            solution.report.safeguard_rounds,
        ],
        factorizations: linear.ldlt_factorizations + linear.lu_factorizations,
        n_cones_optimal,
        n_clusters: solution.extraction.n_clusters,
        merge: outcome,
        distortion_optimal: before,
        distortion_merged: after,
        distortion: if outcome == MergeOutcome::Rejected {
            before
        } else {
            after
        },
        before_reoptimization: None,
        seconds: solution.report.timings.total,
        message: solution.report.message.clone(),
    });

    // Optionally: the best angles for these cone positions, without the penalty.
    if options.reoptimize_angles && !solved.cones.is_empty() {
        let mut free = vec![f64::INFINITY; n];
        for cone in &solved.cones {
            free[cone.vertex as usize] = 0.0;
        }
        for v in (0..n).filter(|&v| interior[v] && costs[v] == 0.0) {
            free[v] = 0.0;
        }
        let again = ConeProblem {
            cone_costs: Some(free),
            ..problem.clone()
        };
        match solve(&again, &SolveOptions::default()) {
            Ok(second) if second.report.status == Status::Converged => {
                let report = solved.report.as_mut().expect("set above");
                report.before_reoptimization = Some((report.distortion, solved.cones.len()));
                report.distortion = second.distortion;
                report.seconds += second.report.timings.total;
                report.factorizations += second.report.linear_solver.ldlt_factorizations
                    + second.report.linear_solver.lu_factorizations;
                solved.cones = second.cones;
                solved.log_scale = Some(second.log_scale);
            }
            Ok(second) => warnings.push(format!(
                "re-optimizing the cone angles of the piece of the surface that contains face \
                 {} did not converge (status {:?}); the angles of the penalized solve are used",
                piece.first_polygon, second.report.status
            )),
            Err(e) => warnings.push(format!(
                "re-optimizing the cone angles of the piece of the surface that contains face \
                 {} failed ({e}); the angles of the penalized solve are used",
                piece.first_polygon
            )),
        }
    }
    Ok(solved)
}
