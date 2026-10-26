//! Field alignment penalty (`SPEC.md` §5.3, D20).

use std::fmt;
use std::sync::Arc;

use super::{check_args, check_weight};
use crate::error::RcError;
use crate::kernel::v3::{cross, dot, scale, sub};
use crate::network::CurveNetwork;
use crate::step::PenaltyTerm;
use crate::subdivide::SubdivisionMap;
use crate::surfaces::{VectorField, jacobian_transpose_apply};

/// `w Σ_I ℓ_I |T_I × X(x_I)|²` with `x_I` the edge midpoint: it is zero where every edge is
/// tangent to the field. `X` is not normalized (the paper assumes a unit field).
///
/// With `c = T_I × X_I` and `J_X` the Jacobian of `X`, the derivative with respect to the
/// endpoint `i_a` (sign `s_a = +1` for `i₂`, `−1` for `i₁`) is
/// `w (s_a |c|² T_I + 2 s_a P_{T_I}(X_I × c) + ℓ_I J_X(x_I)ᵀ(c × T_I))`,
/// `P_T = I − T Tᵀ`. The reference's value omits the factor `ℓ_I` that its gradient contains
/// (D20); here both are the same function.
///
/// [formula: rc-penalty-field]
#[derive(Clone)]
pub struct FieldAlignment {
    field: Arc<dyn VectorField>,
    weight: f64,
}

impl fmt::Debug for FieldAlignment {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("FieldAlignment")
            .field("weight", &self.weight)
            .finish_non_exhaustive()
    }
}

impl FieldAlignment {
    /// The penalty for `field` with weight `w >= 0`.
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for a negative or non-finite weight.
    pub fn new(field: Arc<dyn VectorField>, weight: f64) -> Result<Self, RcError> {
        check_weight("FieldAlignment.weight", weight, false)?;
        Ok(FieldAlignment { field, weight })
    }

    /// The weight.
    pub fn weight(&self) -> f64 {
        self.weight
    }
}

impl PenaltyTerm for FieldAlignment {
    fn name(&self) -> &str {
        "FieldAlignment"
    }

    /// [formula: rc-penalty-field]
    fn evaluate(
        &self,
        network: &CurveNetwork,
        positions: &[[f64; 3]],
        mut gradient: Option<&mut [[f64; 3]]>,
    ) -> Result<f64, RcError> {
        check_args(network, positions, &gradient)?;
        let geo = network.edge_geometry(positions);
        let mut total = 0.0;
        for (e, &[a, b]) in network.edges().iter().enumerate() {
            let (l, t, x) = (geo.length[e], geo.tangent[e], geo.midpoint[e]);
            let xv = self.field.value(x);
            let c = cross(t, xv);
            let c2 = dot(c, c);
            total += l * c2;
            if let Some(g) = gradient.as_deref_mut() {
                let u = cross(xv, c);
                let pu = sub(u, scale(t, dot(t, u)));
                let jt = jacobian_transpose_apply(&self.field.jacobian(x), cross(c, t));
                for k in 0..3 {
                    let common = l * jt[k];
                    let s = c2 * t[k] + 2.0 * pu[k];
                    g[b as usize][k] += self.weight * (s + common);
                    g[a as usize][k] += self.weight * (-s + common);
                }
            }
        }
        Ok(self.weight * total)
    }

    fn subdivided(&self, _map: &SubdivisionMap) -> Result<Arc<dyn PenaltyTerm>, RcError> {
        Ok(Arc::new(self.clone()))
    }
}
