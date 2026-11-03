//! The objective `F = w_TP ℰ + Σ w_i P_i` of one step with frozen Barnes–Hut lists
//! (`SPEC.md` §3.4, §5.3, §7.2): the energy's tree and every tangent-point obstacle's
//! [`CrossEvaluator`] are built at the step's base positions and only refit for the line-search
//! trials, so every trial evaluates the same function whose exact gradient drove the direction.

use std::sync::Arc;

use crate::energy::barnes_hut::BhTree;
use crate::energy::{BhOptions, Energy, EnergyEval, exact};
use crate::error::RsError;
use crate::geometry::v3::axpy;
use crate::geometry::{FaceGeometry, vertex_areas};
use crate::kernel::Exponents;
use crate::metric::FaceDiagonal;
use crate::obstacle::{CrossEvaluator, ObstacleKind};
use crate::penalties::{Penalty, coulomb_energy, evaluate_one, willmore_and_gradient};
use crate::surface::Surface;

/// The tangent-point part.
enum TpEval {
    /// Not the tangent-point energy (Coulomb, Willmore).
    None,
    /// All pairs.
    Exact,
    /// Barnes–Hut with frozen lists.
    Bh(Box<BhTree>),
}

/// One evaluation of the objective.
#[derive(Debug, Clone)]
pub(crate) struct Evaluation {
    /// `ℰ` (unweighted: tangent-point, Coulomb or Willmore).
    pub energy: f64,
    /// `w_i P_i` per penalty.
    pub penalties: Vec<f64>,
    /// `F = w_TP ℰ + Σ w_i P_i`.
    pub total: f64,
    /// `dF` (when asked for).
    pub differential: Option<Vec<[f64; 3]>>,
}

/// The objective of one step (module docs).
pub(crate) struct Objective<'a> {
    pool: &'a rayon::ThreadPool,
    surface: &'a Surface,
    exponents: Exponents,
    energy: Energy,
    eval: EnergyEval,
    tp_weight: f64,
    tp: TpEval,
    penalties: &'a [Penalty],
    /// Per penalty: the frozen cross evaluator of a mesh or point-cloud obstacle.
    cross: Vec<Option<CrossEvaluator>>,
    /// Seconds spent building trees and lists.
    pub build_seconds: f64,
}

/// Errors that mean "this configuration has no finite objective" (a rejected trial, §2.3).
pub(crate) fn is_geometric(e: &RsError) -> bool {
    matches!(
        e,
        RsError::DegenerateFace { .. }
            | RsError::NonFiniteEnergy { .. }
            | RsError::NonFinitePosition { .. }
    )
}

impl<'a> Objective<'a> {
    /// Builds the trees and lists at `positions`.
    #[allow(clippy::too_many_arguments)]
    pub fn build(
        pool: &'a rayon::ThreadPool,
        surface: &'a Surface,
        positions: &[[f64; 3]],
        exponents: Exponents,
        energy: Energy,
        eval: EnergyEval,
        tp_weight: f64,
        penalties: &'a [Penalty],
        gpu: Option<&Arc<crate::gpu::GpuContext>>,
    ) -> Result<Objective<'a>, RsError> {
        let t0 = std::time::Instant::now();
        let tp = match (energy, eval) {
            (Energy::TangentPoint, EnergyEval::Exact) => TpEval::Exact,
            (Energy::TangentPoint, EnergyEval::BarnesHut(o)) => {
                let mut t = BhTree::build(pool, surface, positions, o)?;
                if let Some(g) = gpu {
                    t.attach_gpu(g)?;
                }
                TpEval::Bh(Box::new(t))
            }
            _ => TpEval::None,
        };
        let bh = match eval {
            EnergyEval::Exact => BhOptions {
                theta: 0.0,
                ..BhOptions::default()
            },
            EnergyEval::BarnesHut(o) => o,
        };
        let mut cross = Vec::with_capacity(penalties.len());
        for p in penalties {
            cross.push(match p {
                Penalty::MeshObstacle { obstacle, kind, .. } => Some(CrossEvaluator::build(
                    pool,
                    surface,
                    positions,
                    Arc::clone(obstacle),
                    *kind,
                    &exponents,
                    bh,
                )?),
                Penalty::PointCloudObstacle { points, .. } => Some(CrossEvaluator::build(
                    pool,
                    surface,
                    positions,
                    Arc::clone(points),
                    ObstacleKind::TangentPoint,
                    &exponents,
                    bh,
                )?),
                _ => None,
            });
        }
        Ok(Objective {
            pool,
            surface,
            exponents,
            energy,
            eval,
            tp_weight,
            tp,
            penalties,
            cross,
            build_seconds: t0.elapsed().as_secs_f64(),
        })
    }

    /// The Barnes–Hut tree of the energy (for the hierarchical back end to reuse its face tree).
    pub fn bh_tree(&self) -> Option<&BhTree> {
        match &self.tp {
            TpEval::Bh(t) => Some(t),
            _ => None,
        }
    }

    /// `F` (and `dF` if `gradient`) at `positions`. The trees are refit to `positions` (they
    /// were built at the base positions; call with the base positions first). Geometric
    /// failures (degenerate faces, coincident barycenters) are errors here; the line search
    /// turns them into rejected trials.
    pub fn evaluate(
        &mut self,
        positions: &[[f64; 3]],
        gradient: bool,
    ) -> Result<Evaluation, RsError> {
        let pool = self.pool;
        let surface = self.surface;
        let nv = surface.n_vertices();
        let (energy, mut diff): (f64, Option<Vec<[f64; 3]>>) = match (&mut self.tp, self.energy) {
            (TpEval::Exact, _) => {
                let geom = FaceGeometry::new(surface, positions)?;
                if gradient {
                    let (e, g) =
                        exact::energy_and_differential(pool, surface, positions, &geom, &self.exponents)?;
                    (e, Some(g))
                } else {
                    (exact::energy_value(pool, &geom, &self.exponents), None)
                }
            }
            (TpEval::Bh(tree), _) => {
                tree.refit(positions)?;
                if gradient {
                    let ev = tree.energy_and_differential(pool, &self.exponents)?;
                    (ev.energy, Some(ev.differential))
                } else {
                    (tree.try_energy(pool, &self.exponents)?, None)
                }
            }
            (TpEval::None, Energy::Coulomb { alpha }) => {
                let (e, g) = coulomb_energy(pool, surface, positions, alpha)?;
                (e, gradient.then_some(g))
            }
            (TpEval::None, Energy::Willmore) => {
                let (e, g) = willmore_and_gradient(surface, positions)?;
                (e, gradient.then_some(g))
            }
            (TpEval::None, Energy::TangentPoint) => unreachable!("built with a tangent-point tree"),
        };
        if let Some(d) = diff.as_mut() {
            if self.tp_weight != 1.0 {
                for g in d.iter_mut() {
                    *g = crate::geometry::v3::scale(*g, self.tp_weight);
                }
            }
            if d.len() != nv {
                return Err(RsError::Internal {
                    what: "objective",
                    detail: format!("differential of length {} for {nv} vertices", d.len()),
                });
            }
        }
        let mut total = self.tp_weight * energy;
        let mut values = Vec::with_capacity(self.penalties.len());
        for (i, p) in self.penalties.iter().enumerate() {
            let (v, g) = match self.cross[i].as_mut() {
                Some(ev) => {
                    ev.refit(positions)?;
                    let w = p.weight();
                    if gradient {
                        let (e, mut g) = ev.energy_and_differential(pool);
                        let wv = w * e;
                        if wv.is_finite() && g.iter().flatten().all(|x| x.is_finite()) {
                            for x in g.iter_mut() {
                                *x = crate::geometry::v3::scale(*x, w);
                            }
                            (wv, Some(g))
                        } else {
                            (wv, None)
                        }
                    } else {
                        (w * ev.energy(pool), None)
                    }
                }
                None => {
                    let (v, g) =
                        evaluate_one(pool, surface, positions, p, &self.exponents, &self.eval)?;
                    (v, gradient.then_some(g))
                }
            };
            if let (Some(d), Some(g)) = (diff.as_mut(), g.as_ref()) {
                for (a, b) in d.iter_mut().zip(g) {
                    *a = axpy(*a, 1.0, *b);
                }
            }
            total += v;
            values.push(v);
        }
        Ok(Evaluation {
            energy,
            penalties: values,
            total,
            differential: diff,
        })
    }

    /// `O_diag` of the tangent-point mesh obstacles (`SPEC.md` §4.4, D23), summed, at the base
    /// positions (`None` without such obstacles). `low_order: false` drops its `B⁰` part.
    /// `point_clouds` includes the point-cloud obstacles (zero obstacle normals), which the
    /// reference does (`reference_compat`); off by default (SPEC §5.3: the preconditioner does
    /// not contain the term and GMRES needs 1600 iterations per step on the Reconstruction scene).
    pub fn obstacle_diagonal(
        &self,
        positions: &[[f64; 3]],
        low_order: bool,
        point_clouds: bool,
    ) -> Result<Option<FaceDiagonal>, RsError> {
        let mut out: Option<FaceDiagonal> = None;
        let mut areas: Option<Vec<f64>> = None;
        for (i, p) in self.penalties.iter().enumerate() {
            let (weight, Some(ev)) = (match p {
                Penalty::MeshObstacle {
                    kind: ObstacleKind::TangentPoint,
                    weight,
                    ..
                }
                | Penalty::PointCloudObstacle { weight, .. }
                    if point_clouds || matches!(p, Penalty::MeshObstacle { .. }) =>
                {
                    weight
                }
                _ => continue,
            }, self.cross[i].as_ref())
            else {
                continue;
            };
            let a = match &areas {
                Some(a) => a,
                None => {
                    areas = Some(FaceGeometry::new(self.surface, positions)?.areas);
                    areas.as_ref().expect("set")
                }
            };
            let d = ev.metric_diagonal(self.pool, self.exponents.s(), *weight)?;
            let o = out.get_or_insert_with(|| FaceDiagonal {
                hi: vec![0.0; a.len()],
                lo: if low_order { vec![0.0; a.len()] } else { Vec::new() },
            });
            // `ObstacleDiagonal` stores `c_S` with `Vᵀ diag(c) V`; `FaceDiagonal` stores `d_S`
            // with `2 Vᵀ diag(d/a) V`: `d_S = c_S a_S / 2`.
            for s in 0..a.len() {
                o.hi[s] += 0.5 * d.high[s] * a[s];
                if low_order {
                    o.lo[s] += 0.5 * d.low[s] * a[s];
                }
            }
        }
        Ok(out)
    }

    /// The total Willmore penalty weight (the weight of the `W M⁻¹ W` term of the metric, D7).
    pub fn willmore_weight(&self) -> f64 {
        self.penalties
            .iter()
            .map(|p| match p {
                Penalty::Willmore { weight } => *weight,
                _ => 0.0,
            })
            .sum()
    }
}

/// The area-weighted mean of `g` over each component (`SPEC.md` §5.1 `BarycenterShift`).
pub(crate) fn component_means(
    surface: &Surface,
    positions: &[[f64; 3]],
    g: &[[f64; 3]],
) -> Result<Vec<[f64; 3]>, RsError> {
    let geom = FaceGeometry::new(surface, positions)?;
    let va = vertex_areas(surface, &geom.areas);
    let nc = surface.n_components();
    let comp = surface.vertex_components();
    let mut sum = vec![[0.0; 3]; nc];
    let mut w = vec![0.0; nc];
    for v in 0..surface.n_vertices() {
        let c = comp[v] as usize;
        sum[c] = axpy(sum[c], va[v], g[v]);
        w[c] += va[v];
    }
    Ok(sum
        .iter()
        .zip(&w)
        .map(|(s, &w)| s.map(|x| x / w))
        .collect())
}
