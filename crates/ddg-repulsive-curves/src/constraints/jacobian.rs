//! Constraint values and Jacobians (`SPEC.md` §5.1).
//!
//! With `s_{I,v} = +1` if `v = i₂`, `−1` if `v = i₁`, `L = Σ ℓ_I`, `x̄` the current barycenter:
//!
//! | constraint | `Φ` | `∂Φ/∂γ_v` |
//! |---|---|---|
//! | `Barycenter` | `(Σ ℓ_I x_I)/L − x₀` | `(1/L) Σ_{I ∋ v} [(ℓ_I/2) I₃ + s_{I,v} (x_I − x̄) T_Iᵀ]` |
//! | `TotalLength` | `L⁰ − Σ ℓ_I` | `−Σ_{I ∋ v} s_{I,v} T_Iᵀ` |
//! | `EdgeLengths` | `ℓ⁰_I − ℓ_I` | `i₁: +T_Iᵀ`, `i₂: −T_Iᵀ` |
//! | `Point` | `γ_v − x_v` | `I₃` |
//! | `Surface` | `f(γ_v)` | `∇f(γ_v)ᵀ` |
//! | `Tangent` | `(e₁·s, e₂·s)` | `e_kᵀ ∂s/∂γ_w`, with `∂T/∂γ_{i₂} = P_T/ℓ = −∂T/∂γ_{i₁}` |
//!
//! The barycenter Jacobian is the exact derivative of the length-weighted mean (D6; the
//! reference uses dual-length weights `Σℓ/deg` held constant). The tangent rows are the two
//! components orthogonal to `X` (D7; the paper's three rows `T_I − X` have rank 2, and the
//! reference's vertex block has the wrong sign).

use super::{ConstraintSet, ResolvedConstraint, TangentSite};
use crate::kernel::v3::{add, dot, norm, scale, sub};
use crate::network::{CurveNetwork, EdgeGeometry};

/// Derivatives `∂s/∂γ_w` of a tangent site: `(w, 3×3 block [r][c])`.
pub type TangentBlocks = Vec<(u32, [[f64; 3]; 3])>;

/// Jacobian entries `(row, vertex, coordinate, value)` in 3-D coordinates.
pub type JacobianEntries = Vec<(usize, u32, usize, f64)>;

/// The barycenter `x̄ = (Σ_I ℓ_I x_I)/L` and `L = Σ_I ℓ_I` (sums in edge order).
///
/// [formula: rc-constraint-barycenter]
pub fn barycenter(network: &CurveNetwork, geom: &EdgeGeometry) -> ([f64; 3], f64) {
    let mut s = [0.0; 3];
    let mut l = 0.0;
    for i in 0..network.n_edges() {
        s = add(s, scale(geom.midpoint[i], geom.length[i]));
        l += geom.length[i];
    }
    (scale(s, 1.0 / l), l)
}

/// `∂x̄/∂γ_v` as a 3×3 block `[r][c] = ∂x̄_r/∂γ_{v,c}`:
/// `(1/L) Σ_{I ∋ v} [(ℓ_I/2) δ_rc + s_{I,v} (x_I − x̄)_r T_{I,c}]`.
///
/// [formula: rc-constraint-barycenter]
pub fn barycenter_block(
    network: &CurveNetwork,
    geom: &EdgeGeometry,
    xbar: [f64; 3],
    total_length: f64,
    v: u32,
) -> [[f64; 3]; 3] {
    let mut j = [[0.0; 3]; 3];
    for &e in network.incident_edges(v) {
        let i = e as usize;
        let s = if network.edges()[i][1] == v {
            1.0
        } else {
            -1.0
        };
        let dx = sub(geom.midpoint[i], xbar);
        let t = geom.tangent[i];
        for r in 0..3 {
            j[r][r] += 0.5 * geom.length[i];
            for c in 0..3 {
                j[r][c] += s * dx[r] * t[c];
            }
        }
    }
    for row in &mut j {
        for x in row.iter_mut() {
            *x /= total_length;
        }
    }
    j
}

/// `∂Φ/∂γ_v = −Σ_{I ∋ v} s_{I,v} T_I` for `Φ = L⁰ − Σ ℓ_I`.
///
/// [formula: rc-constraint-total-length]
pub fn total_length_gradient(network: &CurveNetwork, geom: &EdgeGeometry, v: u32) -> [f64; 3] {
    let mut g = [0.0; 3];
    for &e in network.incident_edges(v) {
        let i = e as usize;
        let s = if network.edges()[i][1] == v {
            1.0
        } else {
            -1.0
        };
        g = sub(g, scale(geom.tangent[i], s));
    }
    g
}

/// `(∂Φ/∂γ_{i₁}, ∂Φ/∂γ_{i₂}) = (T_I, −T_I)` for `Φ = ℓ⁰_I − ℓ_I`.
///
/// [formula: rc-constraint-edge-length]
pub fn edge_length_gradient(geom: &EdgeGeometry, edge: usize) -> ([f64; 3], [f64; 3]) {
    let t = geom.tangent[edge];
    (t, scale(t, -1.0))
}

/// `P_u / ℓ` for a unit vector `u`: the derivative of `(q − p)/|q − p|` with respect to `q`.
fn unit_derivative(u: [f64; 3], l: f64) -> [[f64; 3]; 3] {
    let mut p = [[0.0; 3]; 3];
    for r in 0..3 {
        for c in 0..3 {
            p[r][c] = ((r == c) as u8 as f64 - u[r] * u[c]) / l;
        }
    }
    p
}

fn neg(m: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    m.map(|r| r.map(|x| -x))
}

fn unit(p: [f64; 3], q: [f64; 3]) -> ([f64; 3], f64) {
    let d = sub(q, p);
    let l = norm(d);
    (scale(d, 1.0 / l), l)
}

/// The tangent `s` of a site (D7) and its derivatives `∂s/∂γ_w` (3×3 blocks `[r][c]`) for
/// every vertex `w` it depends on:
///
/// * `Edge(I)`: `s = T_I`, `∂s/∂γ_{i₂} = P_T/ℓ_I = −∂s/∂γ_{i₁}`;
/// * `Vertex(v)`, degree 1, edge to `w`: `s = (γ_w − γ_v)/|γ_w − γ_v|` (oriented away from `v`);
/// * `Vertex(v)`, degree 2, `a` the other end of the first incident edge (edge order) and `b`
///   of the second: `s = (γ_v − γ_a)/|γ_v − γ_a| + (γ_b − γ_v)/|γ_b − γ_v|`.
///
/// The site must be valid (checked by [`ConstraintSet::new`]).
///
/// [formula: rc-constraint-tangent]
pub fn tangent_vector(
    network: &CurveNetwork,
    positions: &[[f64; 3]],
    site: TangentSite,
) -> ([f64; 3], TangentBlocks) {
    let edges = network.edges();
    let other = |e: u32, v: u32| {
        let [a, b] = edges[e as usize];
        if a == v { b } else { a }
    };
    match site {
        TangentSite::Edge(e) => {
            let [a, b] = edges[e as usize];
            let (t, l) = unit(positions[a as usize], positions[b as usize]);
            let p = unit_derivative(t, l);
            (t, vec![(a, neg(p)), (b, p)])
        }
        TangentSite::Vertex(v) => {
            let inc = network.incident_edges(v);
            if inc.len() == 1 {
                let w = other(inc[0], v);
                let (t, l) = unit(positions[v as usize], positions[w as usize]);
                let p = unit_derivative(t, l);
                (t, vec![(v, neg(p)), (w, p)])
            } else {
                let a = other(inc[0], v);
                let b = other(inc[1], v);
                let (u1, l1) = unit(positions[a as usize], positions[v as usize]);
                let (u2, l2) = unit(positions[v as usize], positions[b as usize]);
                let p1 = unit_derivative(u1, l1);
                let p2 = unit_derivative(u2, l2);
                let mut pv = [[0.0; 3]; 3];
                for r in 0..3 {
                    for c in 0..3 {
                        pv[r][c] = p1[r][c] - p2[r][c];
                    }
                }
                (add(u1, u2), vec![(a, neg(p1)), (v, pv), (b, p2)])
            }
        }
    }
}

/// A deterministic orthonormal pair `(e₁, e₂)` completing `X/|X|` to an orthonormal frame
/// (Duff, Burgess, Christensen, Hery, Kensler, Liani, Villemin, *Building an Orthonormal
/// Basis, Revisited*, JCGT 6(1), 2017). Then `(e₁·s, e₂·s) = 0 ⇔ s ∥ X`.
///
/// [formula: rc-tangent-frame]
pub fn tangent_frame(x: [f64; 3]) -> [[f64; 3]; 2] {
    let n = scale(x, 1.0 / norm(x));
    let sign = 1.0f64.copysign(n[2]);
    let a = -1.0 / (sign + n[2]);
    let b = n[0] * n[1] * a;
    [
        [1.0 + sign * n[0] * n[0] * a, sign * b, -sign * n[0]],
        [b, sign + n[1] * n[1] * a, -n[1]],
    ]
}

/// `Φ` on all rows of `set` (dropped ones included) and, if asked, the Jacobian entries
/// `(row, vertex, coordinate, value)` in 3-D coordinates (the caller drops the held coordinate
/// in planar mode). Rows in the order of [`ConstraintSet::rows`].
///
/// [formula: rc-constraint-barycenter]
/// [formula: rc-constraint-total-length]
/// [formula: rc-constraint-edge-length]
/// [formula: rc-constraint-point]
/// [formula: rc-constraint-surface]
/// [formula: rc-constraint-tangent]
pub(crate) fn evaluate(
    set: &ConstraintSet,
    network: &CurveNetwork,
    positions: &[[f64; 3]],
    with_jacobian: bool,
) -> (Vec<f64>, JacobianEntries) {
    let geom = network.edge_geometry(positions);
    let coords = set.layout.coords();
    let n = network.n_vertices();
    let edges = network.edges();
    let mut values = Vec::with_capacity(set.rows.len());
    let mut jac = Vec::new();
    let mut need_bary = None;
    for c in &set.resolved {
        let row0 = values.len();
        match c {
            ResolvedConstraint::Barycenter { target } => {
                let (xbar, l) = *need_bary.get_or_insert_with(|| barycenter(network, &geom));
                for &r in coords {
                    values.push(xbar[r] - target[r]);
                }
                if with_jacobian {
                    for v in 0..n as u32 {
                        let blk = barycenter_block(network, &geom, xbar, l, v);
                        for (k, &r) in coords.iter().enumerate() {
                            for c in 0..3 {
                                jac.push((row0 + k, v, c, blk[r][c]));
                            }
                        }
                    }
                }
            }
            ResolvedConstraint::TotalLength { target } => {
                let l: f64 = geom.length.iter().sum();
                values.push(target - l);
                if with_jacobian {
                    for v in 0..n as u32 {
                        let g = total_length_gradient(network, &geom, v);
                        for c in 0..3 {
                            jac.push((row0, v, c, g[c]));
                        }
                    }
                }
            }
            ResolvedConstraint::EdgeLengths {
                edges: list,
                targets,
            } => {
                for (k, (&e, &t)) in list.iter().zip(targets).enumerate() {
                    values.push(t - geom.length[e as usize]);
                    if with_jacobian {
                        let (g1, g2) = edge_length_gradient(&geom, e as usize);
                        let [a, b] = edges[e as usize];
                        for c in 0..3 {
                            jac.push((row0 + k, a, c, g1[c]));
                            jac.push((row0 + k, b, c, g2[c]));
                        }
                    }
                }
            }
            ResolvedConstraint::Point { vertex, target } => {
                let p = positions[*vertex as usize];
                for (k, &r) in coords.iter().enumerate() {
                    values.push(p[r] - target[r]);
                    if with_jacobian {
                        jac.push((row0 + k, *vertex, r, 1.0));
                    }
                }
            }
            ResolvedConstraint::Surface { vertices, surface } => {
                for (k, &v) in vertices.iter().enumerate() {
                    let p = positions[v as usize];
                    values.push(surface.value(p));
                    if with_jacobian {
                        let g = surface.gradient(p);
                        for c in 0..3 {
                            jac.push((row0 + k, v, c, g[c]));
                        }
                    }
                }
            }
            ResolvedConstraint::Tangent {
                site,
                frame,
                n_rows,
                ..
            } => {
                let (s, blocks) = tangent_vector(network, positions, *site);
                for (k, e) in frame.iter().enumerate().take(*n_rows) {
                    values.push(dot(*e, s));
                    if with_jacobian {
                        for (w, blk) in &blocks {
                            for c in 0..3 {
                                let val = e[0] * blk[0][c] + e[1] * blk[1][c] + e[2] * blk[2][c];
                                jac.push((row0 + k, *w, c, val));
                            }
                        }
                    }
                }
            }
        }
    }
    (values, jac)
}
