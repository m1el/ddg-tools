//! `optimize` on a mesh: the generic optimizer applied to the developability energy (SPEC §4, §7.1).

use ddg_core::topology::TriMesh;

use super::core::{EvalError, Objective, OptimizeError, Problem, minimize};
use super::options::{Control, FlatProgress, OptimizeOptions, OptimizeReport};
use crate::energy::{EnergyOptions, evaluate};
use crate::error::DevError;

/// Progress of [`optimize`], passed to the callback after every accepted iteration.
#[derive(Debug, Clone, Copy)]
pub struct Progress<'a> {
    /// Round (set by the flow driver; 0 for a plain `optimize`).
    pub round: usize,
    /// Accepted iterations so far.
    pub iteration: usize,
    /// Current energy.
    pub energy: f64,
    /// `ℓ̄ ‖g‖_∞`.
    pub gradient_norm: f64,
    /// Current positions.
    pub positions: &'a [[f64; 3]],
}

/// Result of [`optimize`].
#[derive(Debug, Clone, PartialEq)]
pub struct Optimized {
    /// The last accepted positions.
    pub positions: Vec<[f64; 3]>,
    /// Report (SPEC §4.7).
    pub report: OptimizeReport,
}

struct MeshObjective<'a> {
    mesh: &'a TriMesh,
    energy: &'a EnergyOptions,
    near_ties: usize,
    min_face_angle: f64,
    evaluated: bool,
}

/// Relative slack of the admissibility check at the starting point (rounding in a subdivision).
pub(crate) const START_SLACK: f64 = 1e-6;

/// The admissibility rule of D31: every face has a smallest corner angle of at least `floor`
/// (`floor = 0`: no check). `slack` (relative, `0` for a trial point) lets a starting point that
/// was accepted at the floor survive rounding, e.g. in a subdivision. Lowest face id first.
pub(crate) fn check_admissible(
    mesh: &TriMesh,
    pos: &[[f64; 3]],
    floor: f64,
    slack: f64,
) -> Result<(), DevError> {
    use crate::geometry::v3::{cross, norm, sub};
    if floor <= 0.0 {
        return Ok(());
    }
    let sin_floor = floor.sin() * (1.0 - slack);
    for (f, t) in mesh.faces().iter().enumerate() {
        let p = t.map(|v| pos[v as usize]);
        let nl = norm(cross(sub(p[1], p[0]), sub(p[2], p[0])));
        // the smallest angle is opposite the shortest edge: sin = |n| / (product of the two longer edges)
        let mut l = [norm(sub(p[1], p[0])), norm(sub(p[2], p[1])), norm(sub(p[0], p[2]))];
        l.sort_by(f64::total_cmp);
        let sin = nl / (l[1] * l[2]);
        if sin.is_nan() || sin < sin_floor {
            let a = crate::geometry::corner_angles(p).into_iter().fold(f64::INFINITY, f64::min);
            return Err(DevError::NearDegenerateFace { face: f as u32, min_angle: a, limit: floor });
        }
    }
    Ok(())
}

/// Faces whose smallest angle is below this multiple of the floor block their vertices (D31).
const BLOCK_FACTOR: f64 = 2.0;

impl Objective for MeshObjective<'_> {
    type Error = DevError;
    fn blocked(&mut self, x: &[f64], blocked: &mut [bool]) {
        blocked.fill(false);
        if self.min_face_angle <= 0.0 {
            return;
        }
        use crate::geometry::v3::{cross, norm, sub};
        let (pos, _) = x.as_chunks::<3>();
        let sin_block = (BLOCK_FACTOR * self.min_face_angle).sin();
        for t in self.mesh.faces() {
            let p = t.map(|v| pos[v as usize]);
            let nl = norm(cross(sub(p[1], p[0]), sub(p[2], p[0])));
            let mut l = [norm(sub(p[1], p[0])), norm(sub(p[2], p[1])), norm(sub(p[0], p[2]))];
            l.sort_by(f64::total_cmp);
            let sin = nl / (l[1] * l[2]);
            if sin.is_nan() || sin < sin_block {
                for &v in t {
                    blocked[v as usize] = true;
                }
            }
        }
    }
    fn evaluate(&mut self, x: &[f64], g: &mut [f64]) -> Result<f64, EvalError<DevError>> {
        let (pos, _) = x.as_chunks::<3>();
        let slack = if self.evaluated { 0.0 } else { START_SLACK };
        self.evaluated = true;
        check_admissible(self.mesh, pos, self.min_face_angle, slack).map_err(EvalError::Inadmissible)?;
        match evaluate(self.mesh, pos, self.energy) {
            Ok(ev) => {
                for (o, v) in g.chunks_exact_mut(3).zip(&ev.gradient) {
                    o.copy_from_slice(v);
                }
                self.near_ties = ev.near_ties;
                Ok(ev.energy)
            }
            // an undefined energy (D12) or a collapsed face at a trial point is +infinity
            Err(
                e @ (DevError::UndefinedEnergy { .. }
                | DevError::DegenerateFace { .. }
                | DevError::NonFinitePosition { .. }),
            ) => Err(EvalError::Undefined(e)),
            Err(e) => Err(EvalError::Fatal(e)),
        }
    }
    fn near_ties(&self) -> usize {
        self.near_ties
    }
    /// Every per-vertex energy is a function of unit normals and angles whose terms are at most
    /// about `2π` (the trace of the angle-weighted covariance), so its absolute rounding error is
    /// a few `2π ε`; summed over the vertices (D34).
    fn energy_noise(&self) -> f64 {
        super::core::NOISE_ULPS * f64::EPSILON * std::f64::consts::TAU * self.mesh.n_vertices() as f64
    }
}

fn mean_edge_length(mesh: &TriMesh, pos: &[[f64; 3]]) -> f64 {
    let sum: f64 = mesh
        .edges()
        .iter()
        .map(|&[a, b]| {
            let (p, q) = (pos[a as usize], pos[b as usize]);
            ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt()
        })
        .sum();
    sum / mesh.n_edges() as f64
}

/// Minimizes the developability energy over the vertex positions (SPEC §4).
///
/// Lengths are measured in units of `ℓ̄`, the mean edge length of `positions` (D17); positions are
/// never rescaled. Pinned vertices (`opts.pinned`) do not move (D18). The callback runs after each
/// accepted iteration; returning [`Control::Break`] ends the run with the last accepted state
/// (`Cancelled`). An energy that is undefined at a trial point of a line search counts as `+∞`
/// (D12); at the starting point it is an error.
///
/// Nothing is modified: the result is returned.
///
/// # Errors
/// Invalid options, wrong lengths, non-finite positions, pinned ids out of range, a degenerate
/// face or undefined energy at the starting point ([`evaluate`]'s errors).
pub fn optimize(
    mesh: &TriMesh,
    positions: &[[f64; 3]],
    energy: &EnergyOptions,
    opts: &OptimizeOptions,
    progress: Option<&mut dyn FnMut(&Progress) -> Control>,
) -> Result<Optimized, DevError> {
    opts.validate().map_err(|e| DevError::InvalidOption {
        name: e.name,
        value: e.value,
        requirement: e.requirement,
    })?;
    let nv = mesh.n_vertices();
    if positions.len() != nv {
        return Err(DevError::LengthMismatch { what: "positions", expected: nv, got: positions.len() });
    }
    if let Some(v) = positions.iter().position(|p| p.iter().any(|c| !c.is_finite())) {
        return Err(DevError::NonFinitePosition { vertex: v as u32 });
    }
    if let Some(&p) = opts.pinned.iter().find(|&&p| p as usize >= nv) {
        return Err(DevError::VertexOutOfRange { what: "pinned", vertex: p, n_vertices: nv });
    }
    let ell = mean_edge_length(mesh, positions);
    if !(ell.is_finite() && ell > 0.0) {
        return Err(DevError::DegenerateFace { face: 0 });
    }
    let x0: Vec<f64> = positions.iter().flatten().copied().collect();
    let mut obj = MeshObjective {
        mesh,
        energy,
        near_ties: 0,
        min_face_angle: opts.min_face_angle,
        evaluated: false,
    };
    let mut adapter = progress.map(|cb| {
        move |p: &FlatProgress| {
            let (pos, _) = p.x.as_chunks::<3>();
            cb(&Progress {
                round: 0,
                iteration: p.iteration,
                energy: p.energy,
                gradient_norm: p.gradient_norm,
                positions: pos,
            })
        }
    });
    let cb: Option<&mut dyn FnMut(&FlatProgress) -> Control> = adapter.as_mut().map(|f| f as _);
    let out = minimize(&mut obj, &x0, &Problem { length_scale: ell, group: 3 }, opts, cb).map_err(|e| match e {
        OptimizeError::Options(e) => DevError::InvalidOption { name: e.name, value: e.value, requirement: e.requirement },
        OptimizeError::PinnedOutOfRange { id, n_groups } => {
            DevError::VertexOutOfRange { what: "pinned", vertex: id, n_vertices: n_groups }
        }
        OptimizeError::UndefinedAtStart(e) | OptimizeError::Objective(e) => e,
        OptimizeError::NonFiniteStart => DevError::UndefinedEnergy {
            vertex: 0,
            reason: "the energy or its gradient is not finite at the starting point",
        },
        e @ (OptimizeError::BadLayout { .. } | OptimizeError::BadLengthScale(_)) => DevError::Internal(e.to_string()),
    })?;
    let (pos, _) = out.x.as_chunks::<3>();
    Ok(Optimized { positions: pos.to_vec(), report: out.report })
}
