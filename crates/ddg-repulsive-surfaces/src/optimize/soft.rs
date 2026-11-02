//! The soft problem of the comparison optimizers (`SPEC.md` §14.2).

use crate::constraints::Constraint;
use crate::error::RsError;
use crate::geometry::FaceGeometry;
use crate::penalties::Penalty;
use crate::penalties::simple::volume_and_gradient;
use crate::step::FlowOptions;
use crate::surface::Surface;

/// The weight of the area penalty of the reference's `scene_soft.txt` and the paper's
/// comparisons.
pub const SOFT_AREA_WEIGHT: f64 = 100.0;

/// Turns a hard-constraint problem into the soft one the comparison optimizers need: every
/// `Area` constraint becomes `AreaDeviation { weight: 100 }` and every `Volume` constraint a
/// `VolumeDeviation { weight: 100 }` (targets: the constraint's, or the initial value);
/// positional constraints are kept. Pin schedules are re-indexed.
///
/// # Errors
/// [`RsError::Unsupported`] for an area/volume constraint of a single component (the penalties
/// are global) or for a growth schedule on a removed constraint.
pub fn soft_problem(
    surface: &Surface,
    positions: &[[f64; 3]],
    options: &FlowOptions,
) -> Result<FlowOptions, RsError> {
    let mut o = options.clone();
    o.constraints.clear();
    let mut map: Vec<Option<usize>> = Vec::new();
    for c in &options.constraints {
        match c {
            Constraint::Area { component, target } => {
                if component.is_some() {
                    return Err(RsError::Unsupported {
                        what: "a per-component area constraint has no global area penalty",
                    });
                }
                let t = match target {
                    Some(t) => *t,
                    None => FaceGeometry::new(surface, positions)?.areas.iter().sum(),
                };
                o.penalties.push(Penalty::AreaDeviation { weight: SOFT_AREA_WEIGHT, target: t });
                map.push(None);
            }
            Constraint::Volume { component, target } => {
                if component.is_some() {
                    return Err(RsError::Unsupported {
                        what: "a per-component volume constraint has no global volume penalty",
                    });
                }
                let t = target.unwrap_or_else(|| volume_and_gradient(surface, positions).0);
                o.penalties.push(Penalty::VolumeDeviation { weight: SOFT_AREA_WEIGHT, target: t });
                map.push(None);
            }
            c => {
                map.push(Some(o.constraints.len()));
                o.constraints.push(c.clone());
            }
        }
    }
    if !o.schedules.is_empty() {
        return Err(RsError::Unsupported {
            what: "growth schedules on a soft problem (the area/volume constraints became penalties)",
        });
    }
    for s in &mut o.pin_schedules {
        match map.get(s.constraint) {
            Some(Some(j)) => s.constraint = *j,
            _ => {
                return Err(RsError::ConstraintOutOfRange {
                    constraint: s.constraint,
                    index: s.constraint as u32,
                });
            }
        }
    }
    Ok(o)
}
