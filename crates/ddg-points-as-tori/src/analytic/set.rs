//! The 70 analytic test instances (SPEC.md §11): the 14 shapes of App. E × 5 instances, seeds
//! 0..69, recorded in `data/points-as-tori/analytic.json`.
//!
//! Instance `seed` has shape `PAPER_KINDS[seed / 5]`; its own generator `Rng::new(seed)` draws, in
//! this order, the profile parameters (see [`Shape2d::random`]), a coin (revolution if `< 0.5`,
//! else extrusion) and the offset `o` or half height `h` uniform in `[0.1, 0.3]`.

use super::sdf2d::{Kind, PAPER_KINDS, Shape2d};
use super::solid::{Solid, SolidKind};
use crate::PatError;
use crate::generate::rng::Rng;

/// Number of instances per shape.
pub const PER_SHAPE: usize = 5;
/// Number of instances in the set.
pub const SET_SIZE: usize = PAPER_KINDS.len() * PER_SHAPE;

/// One instance of the set.
#[derive(Debug, Clone)]
pub struct Instance {
    /// Seed (and index) of the instance.
    pub seed: u64,
    /// The 3D shape.
    pub solid: Solid,
}

impl Instance {
    /// The instance with this seed (any `u64`; the set uses 0..70; the shape class is
    /// `PAPER_KINDS[(seed / 5) % 14]`).
    pub fn from_seed(seed: u64) -> Result<Instance, PatError> {
        let kind = PAPER_KINDS[(seed as usize / PER_SHAPE) % PAPER_KINDS.len()];
        Instance::with_kind(kind, seed)
    }

    /// An instance of the given class drawn from `Rng::new(seed)`.
    pub fn with_kind(kind: Kind, seed: u64) -> Result<Instance, PatError> {
        let mut rng = Rng::new(seed);
        let profile = Shape2d::random(kind, &mut rng);
        let revolve = rng.f64() < 0.5;
        let x = rng.uniform(0.1, 0.3);
        let sk = if revolve { SolidKind::Revolution { o: x } } else { SolidKind::Extrusion { h: x } };
        Ok(Instance { seed, solid: Solid::new(profile, sk)? })
    }

    /// `"revolution"` or `"extrusion"`.
    pub fn solid_name(&self) -> &'static str {
        match self.solid.kind() {
            SolidKind::Revolution { .. } => "revolution",
            SolidKind::Extrusion { .. } => "extrusion",
        }
    }

    /// A human-readable name such as `pie-revolution-17`.
    pub fn name(&self) -> String {
        format!("{}-{}-{}", self.solid.profile().kind().name(), self.solid_name(), self.seed)
    }

    fn json(&self) -> String {
        let (key, val) = match self.solid.kind() {
            SolidKind::Revolution { o } => ("o", o),
            SolidKind::Extrusion { h } => ("h", h),
        };
        let params: Vec<String> = self.solid.profile().params().iter().map(|p| format!("{p}")).collect();
        format!(
            "    {{\"seed\": {}, \"shape\": \"{}\", \"solid\": \"{}\", \"{}\": {}, \"params\": [{}]}}",
            self.seed,
            self.solid.profile().kind().name(),
            self.solid_name(),
            key,
            val,
            params.join(", ")
        )
    }
}

/// The 70 instances, seeds 0..69.
pub fn analytic_set() -> Result<Vec<Instance>, PatError> {
    (0..SET_SIZE as u64).map(Instance::from_seed).collect()
}

/// The JSON record of the set (parameters printed with shortest round-trip precision).
pub fn analytic_json(set: &[Instance]) -> String {
    let rows: Vec<String> = set.iter().map(Instance::json).collect();
    format!(
        "{{\n  \"description\": \"Points as Tori analytic test set (SPEC 11): 14 App. E shapes x 5 instances, seed = index; shape = PAPER_KINDS[seed/5]; profile params drawn by Shape2d::random from Rng::new(seed), then a coin (revolution if < 0.5) and o or h in [0.1, 0.3]. Revolution about y: phi2(|(x,z)| - o, y); extrusion to z in [-h, h]. Params are listed in the order of Shape2d::new.\",\n  \"count\": {},\n  \"instances\": [\n{}\n  ]\n}}\n",
        set.len(),
        rows.join(",\n")
    )
}

/// Write `analytic.json` (creating the directory).
pub fn write_analytic_json(path: &std::path::Path) -> Result<(), PatError> {
    let set = analytic_set()?;
    if let Some(dir) = path.parent() {
        std::fs::create_dir_all(dir).map_err(|e| PatError::Format { reason: format!("{}: {e}", dir.display()) })?;
    }
    std::fs::write(path, analytic_json(&set)).map_err(|e| PatError::Format { reason: format!("{}: {e}", path.display()) })
}
