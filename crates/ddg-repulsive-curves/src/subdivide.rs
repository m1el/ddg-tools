//! Edge subdivision with the constraint mapping of `SPEC.md` §5.4 (D13).
//!
//! Every edge is split at its midpoint. Numbering: old vertices keep their indices, the
//! midpoint of edge `I` is vertex `|V| + I`; child `(i₁, m)` keeps index `I`, child `(m, i₂)` is
//! `|E| + I` (orientation kept). Constraints map as follows: `EdgeLengths` targets are halved on
//! both children; `TotalLength`, `Barycenter`, `Point` are unchanged (splitting at midpoints keeps
//! the total length and the length-weighted barycenter); `Tangent(Edge(I))` applies to both
//! children; `Tangent(Vertex(v))` is unchanged (the unit edge vectors at `v` keep their
//! directions); `Surface` rows are added for the midpoints of edges whose two endpoints are
//! constrained to the same surface. All constraints still hold after the split except `Surface`
//! on the new vertices, which the next re-target projection of the session fixes.
//!
//! [formula: rc-subdivision-constraints]

use crate::constraints::{Constraint, ResolvedConstraint, TangentSite};
use crate::error::RcError;
use crate::network::CurveNetwork;

/// The index map of one subdivision.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SubdivisionMap {
    /// `|V|` before the split.
    pub old_vertices: usize,
    /// `|E|` before the split.
    pub old_edges: usize,
}

impl SubdivisionMap {
    /// The midpoint vertex of old edge `edge`: `|V| + I`.
    pub fn midpoint(&self, edge: u32) -> u32 {
        (self.old_vertices + edge as usize) as u32
    }

    /// The two children of old edge `edge`: `[I, |E| + I]` (`(i₁, m)` and `(m, i₂)`).
    pub fn children(&self, edge: u32) -> [u32; 2] {
        [edge, (self.old_edges + edge as usize) as u32]
    }

    /// Maps per-edge data (e.g. rest lengths) to the children, each child getting `f(value)`.
    pub fn map_edge_values(&self, values: &[f64], f: impl Fn(f64) -> f64) -> Vec<f64> {
        let mut out: Vec<f64> = values.iter().map(|&x| f(x)).collect();
        out.extend(values.iter().map(|&x| f(x)));
        out
    }
}

/// Splits every edge of `network` at its midpoint.
///
/// # Errors
/// [`RcError::PositionCountMismatch`]; [`RcError::InvalidOption`] if `|V| + |E|` does not fit
/// the `u32` indices.
pub fn subdivide(
    network: &CurveNetwork,
    positions: &[[f64; 3]],
) -> Result<(CurveNetwork, Vec<[f64; 3]>, SubdivisionMap), RcError> {
    let n = network.n_vertices();
    let m = network.n_edges();
    if positions.len() != n {
        return Err(RcError::PositionCountMismatch {
            expected: n,
            found: positions.len(),
        });
    }
    if n + m > u32::MAX as usize || 2 * m > u32::MAX as usize {
        return Err(RcError::InvalidOption {
            name: "subdivision",
            value: format!("{n} vertices, {m} edges"),
            requirement: "|V| + |E| and 2|E| must fit 32-bit indices",
        });
    }
    let map = SubdivisionMap {
        old_vertices: n,
        old_edges: m,
    };
    let edges = network.edges();
    let mut pos = Vec::with_capacity(n + m);
    pos.extend_from_slice(positions);
    for &[a, b] in edges {
        let (p, q) = (positions[a as usize], positions[b as usize]);
        pos.push([
            0.5 * (p[0] + q[0]),
            0.5 * (p[1] + q[1]),
            0.5 * (p[2] + q[2]),
        ]);
    }
    let mut new_edges = Vec::with_capacity(2 * m);
    for (i, &[a, _]) in edges.iter().enumerate() {
        new_edges.push([a, map.midpoint(i as u32)]);
    }
    for (i, &[_, b]) in edges.iter().enumerate() {
        new_edges.push([map.midpoint(i as u32), b]);
    }
    let net = CurveNetwork::new(n + m, new_edges)?;
    Ok((net, pos, map))
}

/// The constraints of the subdivided network (explicit targets), and for every old constraint
/// the indices of the new constraints it became (`Tangent(Edge(I))` becomes two).
///
/// [formula: rc-subdivision-constraints]
pub fn subdivide_constraints(
    constraints: &[ResolvedConstraint],
    old_network: &CurveNetwork,
    map: &SubdivisionMap,
) -> (Vec<Constraint>, Vec<Vec<usize>>) {
    let mut out = Vec::with_capacity(constraints.len() + 1);
    let mut index = Vec::with_capacity(constraints.len());
    for c in constraints {
        let first = out.len();
        match c {
            ResolvedConstraint::Barycenter { target } => out.push(Constraint::Barycenter {
                target: Some(*target),
            }),
            ResolvedConstraint::TotalLength { target } => out.push(Constraint::TotalLength {
                target: Some(*target),
            }),
            ResolvedConstraint::EdgeLengths { edges, targets } => {
                let mut e: Vec<u32> = edges.iter().map(|&i| map.children(i)[0]).collect();
                e.extend(edges.iter().map(|&i| map.children(i)[1]));
                let t = map.map_edge_values(targets, |x| 0.5 * x);
                out.push(Constraint::EdgeLengths {
                    edges: Some(e),
                    targets: Some(t),
                });
            }
            ResolvedConstraint::Point { vertex, target } => out.push(Constraint::Point {
                vertex: *vertex,
                target: Some(*target),
            }),
            ResolvedConstraint::Surface { vertices, surface } => {
                let mut on = vec![false; map.old_vertices];
                for &v in vertices {
                    if let Some(x) = on.get_mut(v as usize) {
                        *x = true;
                    }
                }
                let mut vs = vertices.clone();
                for (i, &[a, b]) in old_network.edges().iter().enumerate() {
                    if on[a as usize] && on[b as usize] {
                        vs.push(map.midpoint(i as u32));
                    }
                }
                out.push(Constraint::Surface {
                    vertices: Some(vs),
                    surface: surface.clone(),
                });
            }
            ResolvedConstraint::Tangent {
                site, direction, ..
            } => match site {
                TangentSite::Edge(i) => {
                    for child in map.children(*i) {
                        out.push(Constraint::Tangent {
                            site: TangentSite::Edge(child),
                            direction: Some(*direction),
                        });
                    }
                }
                TangentSite::Vertex(v) => out.push(Constraint::Tangent {
                    site: TangentSite::Vertex(*v),
                    direction: Some(*direction),
                }),
            },
        }
        index.push((first..out.len()).collect());
    }
    (out, index)
}
