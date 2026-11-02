//! Mesh and point-cloud obstacles (`SPEC.md` §5.3, D23).
//!
//! An obstacle is a static set of weighted points `X_T` with weights `a_T`, optionally with unit
//! normals `N_T` ([`StaticSet`]): the centroids, areas and normals of a triangle soup (the
//! curves crate's [`Obstacle`]), or the points of a cloud with weight 1.
//!
//! # Energies
//!
//! With the mesh faces `S` (`a_S`, `X_S`, `N_S`) and the obstacle elements `T`:
//!
//! ```text
//! TangentPoint (triangles):  Σ_{S,T} [ K(X_S, N_S; X_T) + K(X_T, N_T; X_S) ] a_S a_T      (two-sided)
//! TangentPoint (point cloud): Σ_{S,p} K(X_S, N_S; y_p) a_S                                (one-sided)
//! Coulomb { power }:          Σ_{S,T} |X_S − X_T|^{−power} a_S a_T                         (as printed)
//! ```
//!
//! evaluated by Barnes–Hut with the admissibility of `SPEC.md` §3.4 (`max(r_ℓ, r_C) < θ d_box`):
//! the first sum over mesh leaves `ℓ` against nodes `C` of the static tree (a lumped entry
//! `K(X_S, N_S; X_C) a_S A_C`, an inadmissible leaf pair exactly), the second over static leaves
//! `λ` against nodes `D` of the mesh face tree (`K(X_T, N_T; X_D) a_T A_D`). `θ = 0` is the exact
//! double sum. The gradient is the exact gradient of the approximation for frozen lists
//! (including the variation of the mesh node aggregates `A_D`, `X_D`, as `rs-bh-aggregate-adjoint`).
//!
//! # The metric term
//!
//! [`CrossEvaluator::metric_diagonal`] returns the data of `O_diag` (`SPEC.md` §4.4).
//!
//! [formula: rs-mesh-obstacle]

use std::sync::Arc;

use ddg_repulsive_curves::obstacle::Obstacle;
use ddg_repulsive_curves::spatial::{Aabb, Bvh, Element};
use rayon::prelude::*;

use crate::energy::BhOptions;
use crate::error::RsError;
use crate::geometry::v3::{axpy, dot, scale, sub};
use crate::geometry::{FaceGeometry, corners, face_pull_back, face_quantities};
use crate::kernel::{Exponents, KernelEval, kernel_of_v, kernel_with_gradient_of_v};
use crate::spatial::{Cluster, FaceTree, bh_admissible, farthest_corner_distance};
use crate::surface::Surface;

/// Which interaction an obstacle exerts (`SPEC.md` §5.3).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ObstacleKind {
    /// The two-sided tangent-point cross energy (the reference; keeps the barrier property),
    /// with the exponents of the energy.
    TangentPoint,
    /// The printed Coulomb-like potential `Σ |X_S − X_T|^{−power} a_S a_T`; no metric term.
    Coulomb {
        /// The exponent `p > 0`.
        power: f64,
    },
}

#[derive(Debug, Clone, Copy)]
struct StaticElem {
    point: [f64; 3],
    bounds: Aabb<3>,
    mass: f64,
}

impl Element<3> for StaticElem {
    fn point(&self) -> [f64; 3] {
        self.point
    }
    fn bounds(&self) -> Aabb<3> {
        self.bounds
    }
    fn mass(&self) -> f64 {
        self.mass
    }
}

/// A static weighted point set with a BVH: the obstacle side of the cross energies.
#[derive(Debug, Clone)]
pub struct StaticSet {
    points: Vec<[f64; 3]>,
    normals: Option<Vec<[f64; 3]>>,
    weights: Vec<f64>,
    bvh: Bvh<3>,
    clusters: Vec<Cluster>,
}

impl StaticSet {
    /// The centroids, normals and areas of the triangles of `obstacle`.
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] for `leaf_size = 0`; [`RsError::Curves`] if the BVH fails.
    pub fn from_obstacle(obstacle: &Obstacle, leaf_size: usize) -> Result<StaticSet, RsError> {
        let v = obstacle.vertices();
        let mut points = Vec::new();
        let mut normals = Vec::new();
        let mut elems = Vec::new();
        for (t, tri) in obstacle.triangle_indices().iter().enumerate() {
            let p = tri.map(|i| v[i as usize]);
            let (_, n, x) = face_quantities(p);
            let mut b = Aabb::<3>::empty();
            for c in p {
                b.grow(c);
            }
            points.push(x);
            normals.push(n);
            elems.push(StaticElem {
                point: x,
                bounds: b,
                mass: obstacle.areas()[t],
            });
        }
        Self::finish(points, Some(normals), elems, leaf_size)
    }

    /// A point cloud (weight 1 each, no normals).
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] for an empty cloud, a non-finite point or `leaf_size = 0`.
    pub fn from_points(points: &[[f64; 3]], leaf_size: usize) -> Result<StaticSet, RsError> {
        if points.is_empty() {
            return Err(RsError::InvalidOption {
                name: "points",
                value: "[]".into(),
                requirement: "at least one point",
            });
        }
        if let Some(p) = points.iter().find(|p| !p.iter().all(|x| x.is_finite())) {
            return Err(RsError::InvalidOption {
                name: "points",
                value: format!("{p:?}"),
                requirement: "finite coordinates",
            });
        }
        let elems = points
            .iter()
            .map(|&p| StaticElem {
                point: p,
                bounds: Aabb::from_point(p),
                mass: 1.0,
            })
            .collect();
        Self::finish(points.to_vec(), None, elems, leaf_size)
    }

    fn finish(
        points: Vec<[f64; 3]>,
        normals: Option<Vec<[f64; 3]>>,
        elems: Vec<StaticElem>,
        leaf_size: usize,
    ) -> Result<StaticSet, RsError> {
        if leaf_size == 0 {
            return Err(RsError::InvalidOption {
                name: "leaf_size",
                value: "0".into(),
                requirement: "at least 1",
            });
        }
        let weights: Vec<f64> = elems.iter().map(|e| e.mass).collect();
        let bvh = Bvh::build(&elems, leaf_size)?;
        let sums = bvh.fold_up(
            |idx| {
                let mut s = [0.0; 7];
                if let Some(nr) = &normals {
                    for &i in idx {
                        let (a, n) = (weights[i as usize], nr[i as usize]);
                        s[0] += a;
                        s[1] += a * n[0] * n[0];
                        s[2] += a * n[0] * n[1];
                        s[3] += a * n[0] * n[2];
                        s[4] += a * n[1] * n[1];
                        s[5] += a * n[1] * n[2];
                        s[6] += a * n[2] * n[2];
                    }
                }
                s
            },
            |l, r| std::array::from_fn(|k| l[k] + r[k]),
        );
        let clusters = bvh
            .nodes()
            .iter()
            .zip(sums)
            .map(|(nd, s)| {
                let centroid = nd.centroid();
                Cluster {
                    area: nd.mass,
                    centroid,
                    projector: std::array::from_fn(|k| s[k + 1] / s[0]),
                    bbox: nd.bbox,
                    radius: farthest_corner_distance(&nd.bbox, centroid),
                }
            })
            .collect();
        Ok(StaticSet {
            points,
            normals,
            weights,
            bvh,
            clusters,
        })
    }

    /// Number of elements (triangles or points).
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// Always `false` (construction rejects empty sets).
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    /// `true` for a triangle obstacle (normals available).
    pub fn has_normals(&self) -> bool {
        self.normals.is_some()
    }

    /// The element points `X_T`.
    pub fn points(&self) -> &[[f64; 3]] {
        &self.points
    }

    /// The weights `a_T` (areas, or 1).
    pub fn weights(&self) -> &[f64] {
        &self.weights
    }
}

#[derive(Debug, Clone, Copy)]
enum Pair {
    Tp(Exponents),
    Coulomb(f64),
}

impl Pair {
    #[inline]
    fn value(&self, v: [f64; 3], n: [f64; 3]) -> f64 {
        match self {
            Pair::Tp(e) => kernel_of_v(e, v, n),
            Pair::Coulomb(p) => dot(v, v).powf(-0.5 * p),
        }
    }
    #[inline]
    fn eval(&self, v: [f64; 3], n: [f64; 3]) -> KernelEval {
        match self {
            Pair::Tp(e) => kernel_with_gradient_of_v(e, v, n),
            Pair::Coulomb(p) => {
                let r2 = dot(v, v);
                let value = r2.powf(-0.5 * p);
                KernelEval {
                    value,
                    d_v: scale(v, -p * value / r2),
                    d_n: [0.0; 3],
                }
            }
        }
    }
}

/// The data of `O_diag` for one obstacle (`SPEC.md` §4.4): for each mesh face `S`,
/// `c_S = (2 w / a_S) Σ_T a_T h(S, T)`, with `h_hi = ρ^{−2s}` and
/// `h_lo = ½ (K₂(S,T) + K₂(T,S)) ρ^{−2s}` (`ρ = |X_S − X_T|`, `K₂(S,T) = (N_S·v)² ρ^{−4}`;
/// lumped nodes use `vᵀP̄_C v` for `K₂(T,S)`). The metric adds
/// `Σ_S a_S² high_S ⟨D_f u(S), D_f v(S)⟩` (that is `Vᵀ diag(high) V`) and
/// `Σ_S a_S² low_S ū_S v̄_S` (`Uᵀ diag(low) U`) to `uᵀ A v`.
#[derive(Debug, Clone, PartialEq)]
pub struct ObstacleDiagonal {
    /// `c_S` for the `B` part (`V` channels).
    pub high: Vec<f64>,
    /// `c_S` for the `B⁰` part (`U` channel).
    pub low: Vec<f64>,
}

struct LeafOut {
    energy: f64,
    own: Vec<(f64, [f64; 3], [f64; 3])>,
}

struct StaticLeafOut {
    energy: f64,
    /// Per lumped entry: `(∂/∂A_D, ∂/∂X_D)`.
    node: Vec<(f64, [f64; 3])>,
    /// Per exact block, per mesh face of the block: `(∂/∂a_S, ∂/∂X_S)`.
    partner: Vec<(f64, [f64; 3])>,
}

/// The evaluator of one obstacle: mesh face tree, static tree and the frozen interaction
/// lists. Build once per step ([`CrossEvaluator::build`]), [`refit`](CrossEvaluator::refit) for
/// line-search trials.
#[derive(Debug, Clone)]
pub struct CrossEvaluator {
    pair: Pair,
    set: Arc<StaticSet>,
    opts: BhOptions,
    faces: Vec<[u32; 3]>,
    n_vertices: usize,
    positions: Vec<[f64; 3]>,
    geom: FaceGeometry,
    tree: FaceTree,
    /// Static leaves (node indices) in increasing order.
    static_leaves: Vec<u32>,
    /// Per mesh leaf: lumped static nodes / exact static leaves.
    l1_lumped: Vec<Vec<u32>>,
    l1_exact: Vec<Vec<u32>>,
    /// Per static leaf (two-sided only): lumped mesh nodes / exact mesh leaves.
    l2_lumped: Vec<Vec<u32>>,
    l2_exact: Vec<Vec<u32>>,
}

impl CrossEvaluator {
    /// Builds the face tree and the interaction lists at `positions`.
    ///
    /// `exponents` are used by [`ObstacleKind::TangentPoint`] only. A Coulomb obstacle needs
    /// a finite positive power; a tangent-point obstacle on a point cloud is one-sided.
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] for bad options; the errors of
    /// [`Surface::validate_positions`].
    pub fn build(
        pool: &rayon::ThreadPool,
        surface: &Surface,
        positions: &[[f64; 3]],
        set: Arc<StaticSet>,
        kind: ObstacleKind,
        exponents: &Exponents,
        opts: BhOptions,
    ) -> Result<CrossEvaluator, RsError> {
        opts.validate()?;
        let pair = match kind {
            ObstacleKind::TangentPoint => Pair::Tp(*exponents),
            ObstacleKind::Coulomb { power } => {
                if !(power.is_finite() && power > 0.0) {
                    return Err(RsError::InvalidOption {
                        name: "power",
                        value: power.to_string(),
                        requirement: "a finite positive number",
                    });
                }
                Pair::Coulomb(power)
            }
        };
        let geom = FaceGeometry::new(surface, positions)?;
        let tree = FaceTree::build(surface, positions, &geom, opts.leaf_size)?;
        let two_sided = matches!(pair, Pair::Tp(_)) && set.normals.is_some();
        let theta = opts.theta;
        let snodes = set.bvh.nodes();
        let static_leaves: Vec<u32> = (0..snodes.len())
            .filter(|&i| snodes[i].is_leaf())
            .map(|i| i as u32)
            .collect();
        let mesh_leaves = tree.leaves();
        let l1: Vec<(Vec<u32>, Vec<u32>)> = pool.install(|| {
            mesh_leaves
                .par_iter()
                .map(|&l| {
                    let cl = &tree.clusters()[l as usize];
                    let (mut lump, mut exact) = (Vec::new(), Vec::new());
                    let mut stack = vec![0usize];
                    while let Some(c) = stack.pop() {
                        if bh_admissible(cl, &set.clusters[c], theta) {
                            lump.push(c as u32);
                        } else if snodes[c].is_leaf() {
                            exact.push(c as u32);
                        } else {
                            stack.push(snodes[c].right as usize);
                            stack.push(c + 1);
                        }
                    }
                    (lump, exact)
                })
                .collect()
        });
        let (l1_lumped, l1_exact) = l1.into_iter().unzip();
        let (l2_lumped, l2_exact) = if two_sided {
            let mnodes = tree.bvh().nodes();
            let l2: Vec<(Vec<u32>, Vec<u32>)> = pool.install(|| {
                static_leaves
                    .par_iter()
                    .map(|&l| {
                        let cl = &set.clusters[l as usize];
                        let (mut lump, mut exact) = (Vec::new(), Vec::new());
                        let mut stack = vec![0usize];
                        while let Some(c) = stack.pop() {
                            if bh_admissible(cl, &tree.clusters()[c], theta) {
                                lump.push(c as u32);
                            } else if mnodes[c].is_leaf() {
                                exact.push(c as u32);
                            } else {
                                stack.push(mnodes[c].right as usize);
                                stack.push(c + 1);
                            }
                        }
                        (lump, exact)
                    })
                    .collect()
            });
            l2.into_iter().unzip()
        } else {
            (Vec::new(), Vec::new())
        };
        Ok(CrossEvaluator {
            pair,
            set,
            opts,
            faces: surface.faces().to_vec(),
            n_vertices: surface.n_vertices(),
            positions: positions.to_vec(),
            geom,
            tree,
            static_leaves,
            l1_lumped,
            l1_exact,
            l2_lumped,
            l2_exact,
        })
    }

    /// Moves the evaluator to `positions` keeping tree topology, boxes, radii and the lists.
    ///
    /// # Errors
    /// [`RsError::PositionCountMismatch`], [`RsError::NonFinitePosition`].
    pub fn refit(&mut self, positions: &[[f64; 3]]) -> Result<(), RsError> {
        if positions.len() != self.n_vertices {
            return Err(RsError::PositionCountMismatch {
                expected: self.n_vertices,
                found: positions.len(),
            });
        }
        if let Some(v) = positions
            .iter()
            .position(|p| !p.iter().all(|x| x.is_finite()))
        {
            return Err(RsError::NonFinitePosition { vertex: v as u32 });
        }
        self.positions.clear();
        self.positions.extend_from_slice(positions);
        for (f, t) in self.faces.iter().enumerate() {
            let (a, n, x) = face_quantities(corners(positions, t));
            self.geom.areas[f] = a;
            self.geom.normals[f] = n;
            self.geom.barycenters[f] = x;
        }
        self.tree.refit_faces(&self.faces, positions, &self.geom)
    }

    /// The options of the build.
    pub fn options(&self) -> &BhOptions {
        &self.opts
    }

    /// Number of (leaf, node) interactions lumped, both sums.
    pub fn n_lumped(&self) -> u64 {
        let c = |v: &[Vec<u32>]| v.iter().map(|x| x.len() as u64).sum::<u64>();
        c(&self.l1_lumped) + c(&self.l2_lumped)
    }

    fn two_sided(&self) -> bool {
        !self.l2_lumped.is_empty()
    }

    fn mesh_leaf_faces(&self, k: usize) -> &[u32] {
        self.tree.faces_of(self.tree.leaves()[k] as usize)
    }

    fn static_leaf_elems(&self, k: usize) -> &[u32] {
        self.set.bvh.elements_of(self.static_leaves[k] as usize)
    }

    fn leaf1_energy(&self, k: usize) -> f64 {
        let g = &self.geom;
        let mut energy = 0.0;
        for &s in self.mesh_leaf_faces(k) {
            let s = s as usize;
            let (xs, ns) = (g.barycenters[s], g.normals[s]);
            let mut row = 0.0;
            for &c in &self.l1_lumped[k] {
                let cl = &self.set.clusters[c as usize];
                row += cl.area * self.pair.value(sub(xs, cl.centroid), ns);
            }
            for &c in &self.l1_exact[k] {
                for &t in self.set.bvh.elements_of(c as usize) {
                    let t = t as usize;
                    row += self.set.weights[t]
                        * self.pair.value(sub(xs, self.set.points[t]), ns);
                }
            }
            energy += g.areas[s] * row;
        }
        energy
    }

    fn leaf2_energy(&self, k: usize) -> f64 {
        let g = &self.geom;
        let normals = self.set.normals.as_ref().expect("two-sided");
        let mut energy = 0.0;
        for &t in self.static_leaf_elems(k) {
            let t = t as usize;
            let (xt, nt) = (self.set.points[t], normals[t]);
            let mut row = 0.0;
            for &d in &self.l2_lumped[k] {
                let cl = &self.tree.clusters()[d as usize];
                row += cl.area * self.pair.value(sub(xt, cl.centroid), nt);
            }
            for &d in &self.l2_exact[k] {
                for &s in self.tree.faces_of(d as usize) {
                    let s = s as usize;
                    row += g.areas[s] * self.pair.value(sub(xt, g.barycenters[s]), nt);
                }
            }
            energy += self.set.weights[t] * row;
        }
        energy
    }

    /// The approximate cross energy at the current positions on the frozen lists. May be
    /// non-finite (coincident points): the caller decides.
    ///
    /// [formula: rs-mesh-obstacle]
    pub fn energy(&self, pool: &rayon::ThreadPool) -> f64 {
        let n1 = self.l1_lumped.len();
        let n2 = self.l2_lumped.len();
        let (e1, e2): (Vec<f64>, Vec<f64>) = pool.install(|| {
            (
                (0..n1).into_par_iter().map(|k| self.leaf1_energy(k)).collect(),
                (0..n2).into_par_iter().map(|k| self.leaf2_energy(k)).collect(),
            )
        });
        e1.iter().sum::<f64>() + e2.iter().sum::<f64>()
    }

    fn leaf1_gradient(&self, k: usize) -> LeafOut {
        let g = &self.geom;
        let faces = self.mesh_leaf_faces(k);
        let mut out = LeafOut {
            energy: 0.0,
            own: vec![(0.0, [0.0; 3], [0.0; 3]); faces.len()],
        };
        for (i, &s) in faces.iter().enumerate() {
            let s = s as usize;
            let (xs, ns, a_s) = (g.barycenters[s], g.normals[s], g.areas[s]);
            let mut row = 0.0;
            let mut dx = [0.0; 3];
            let mut dn = [0.0; 3];
            for &c in &self.l1_lumped[k] {
                let cl = &self.set.clusters[c as usize];
                let kv = self.pair.eval(sub(xs, cl.centroid), ns);
                row += cl.area * kv.value;
                dx = axpy(dx, cl.area, kv.d_v);
                dn = axpy(dn, cl.area, kv.d_n);
            }
            for &c in &self.l1_exact[k] {
                for &t in self.set.bvh.elements_of(c as usize) {
                    let t = t as usize;
                    let w = self.set.weights[t];
                    let kv = self.pair.eval(sub(xs, self.set.points[t]), ns);
                    row += w * kv.value;
                    dx = axpy(dx, w, kv.d_v);
                    dn = axpy(dn, w, kv.d_n);
                }
            }
            out.energy += a_s * row;
            let o = &mut out.own[i];
            o.0 = row;
            o.1 = scale(dx, a_s);
            o.2 = scale(dn, a_s);
        }
        out
    }

    fn leaf2_gradient(&self, k: usize) -> StaticLeafOut {
        let g = &self.geom;
        let normals = self.set.normals.as_ref().expect("two-sided");
        let mut out = StaticLeafOut {
            energy: 0.0,
            node: vec![(0.0, [0.0; 3]); self.l2_lumped[k].len()],
            partner: Vec::new(),
        };
        let n_partner: usize = self.l2_exact[k]
            .iter()
            .map(|&d| self.tree.faces_of(d as usize).len())
            .sum();
        out.partner = vec![(0.0, [0.0; 3]); n_partner];
        for &t in self.static_leaf_elems(k) {
            let t = t as usize;
            let (xt, nt, a_t) = (self.set.points[t], normals[t], self.set.weights[t]);
            for (j, &d) in self.l2_lumped[k].iter().enumerate() {
                let cl = &self.tree.clusters()[d as usize];
                let kv = self.pair.eval(sub(xt, cl.centroid), nt);
                out.energy += a_t * cl.area * kv.value;
                let o = &mut out.node[j];
                o.0 += a_t * kv.value;
                o.1 = axpy(o.1, -a_t * cl.area, kv.d_v);
            }
            let mut p = 0;
            for &d in &self.l2_exact[k] {
                for &s in self.tree.faces_of(d as usize) {
                    let s = s as usize;
                    let kv = self.pair.eval(sub(xt, g.barycenters[s]), nt);
                    out.energy += a_t * g.areas[s] * kv.value;
                    let o = &mut out.partner[p];
                    o.0 += a_t * kv.value;
                    o.1 = axpy(o.1, -a_t * g.areas[s], kv.d_v);
                    p += 1;
                }
            }
        }
        out
    }

    /// The energy and its exact differential `∂Ẽ/∂f_v` on the frozen lists.
    ///
    /// [formula: rs-mesh-obstacle]
    /// [formula: rs-bh-aggregate-adjoint]
    pub fn energy_and_differential(&self, pool: &rayon::ThreadPool) -> (f64, Vec<[f64; 3]>) {
        let nf = self.faces.len();
        let n1 = self.l1_lumped.len();
        let n2 = self.l2_lumped.len();
        let (o1, o2): (Vec<LeafOut>, Vec<StaticLeafOut>) = pool.install(|| {
            (
                (0..n1).into_par_iter().map(|k| self.leaf1_gradient(k)).collect(),
                (0..n2).into_par_iter().map(|k| self.leaf2_gradient(k)).collect(),
            )
        });
        let mut energy = 0.0;
        let mut d_area = vec![0.0; nf];
        let mut d_bary = vec![[0.0; 3]; nf];
        let mut d_normal = vec![[0.0; 3]; nf];
        for (k, o) in o1.iter().enumerate() {
            energy += o.energy;
            for (i, &s) in self.mesh_leaf_faces(k).iter().enumerate() {
                let s = s as usize;
                d_area[s] += o.own[i].0;
                d_bary[s] = axpy(d_bary[s], 1.0, o.own[i].1);
                d_normal[s] = axpy(d_normal[s], 1.0, o.own[i].2);
            }
        }
        let nn = self.tree.clusters().len();
        let mut node_c = vec![0.0; nn];
        let mut node_x = vec![[0.0; 3]; nn];
        for (k, o) in o2.iter().enumerate() {
            energy += o.energy;
            for (j, &d) in self.l2_lumped[k].iter().enumerate() {
                let cl = &self.tree.clusters()[d as usize];
                let (ea, ex) = o.node[j];
                let inv = 1.0 / cl.area;
                // ∂/∂a_S = ea + ex·(X_S − X_D)/A_D ; ∂/∂X_S = a_S ex/A_D.
                node_c[d as usize] += ea - dot(ex, cl.centroid) * inv;
                node_x[d as usize] = axpy(node_x[d as usize], inv, ex);
            }
            let mut p = 0;
            for &d in &self.l2_exact[k] {
                for &s in self.tree.faces_of(d as usize) {
                    let s = s as usize;
                    d_area[s] += o.partner[p].0;
                    d_bary[s] = axpy(d_bary[s], 1.0, o.partner[p].1);
                    p += 1;
                }
            }
        }
        if self.two_sided() {
            let nodes = self.tree.bvh().nodes();
            for n in 0..nn {
                if !nodes[n].is_leaf() {
                    for c in [n + 1, nodes[n].right as usize] {
                        node_c[c] += node_c[n];
                        node_x[c] = axpy(node_x[c], 1.0, node_x[n]);
                    }
                }
            }
            for s in 0..nf {
                let l = self.tree.leaf_of_face(s);
                d_area[s] += node_c[l] + dot(node_x[l], self.geom.barycenters[s]);
                d_bary[s] = axpy(d_bary[s], self.geom.areas[s], node_x[l]);
            }
        }
        let mut out = vec![[0.0; 3]; self.n_vertices];
        for (f, t) in self.faces.iter().enumerate() {
            let gr = face_pull_back(
                corners(&self.positions, t),
                d_area[f],
                d_bary[f],
                d_normal[f],
            );
            for c in 0..3 {
                let v = t[c] as usize;
                out[v] = axpy(out[v], 1.0, gr[c]);
            }
        }
        (energy, out)
    }

    /// The data of `O_diag` (`SPEC.md` §4.4) for the exponent `s` of the energy, with the
    /// obstacle's weight `w` included, from the first interaction lists (same `θ`).
    ///
    /// # Errors
    /// [`RsError::ObstacleError`] for a point cloud (no normals) or a Coulomb obstacle.
    ///
    /// [formula: rs-obstacle-diagonal]
    pub fn metric_diagonal(
        &self,
        pool: &rayon::ThreadPool,
        s: f64,
        weight: f64,
    ) -> Result<ObstacleDiagonal, RsError> {
        let normals = match (&self.pair, &self.set.normals) {
            (Pair::Tp(_), Some(n)) => n,
            _ => {
                return Err(RsError::ObstacleError {
                    obstacle: 0,
                    message: "the metric term needs a tangent-point obstacle with triangles".into(),
                });
            }
        };
        let g = &self.geom;
        let rows: Vec<Vec<(f64, f64)>> = pool.install(|| {
            (0..self.l1_lumped.len())
                .into_par_iter()
                .map(|k| {
                    self.mesh_leaf_faces(k)
                        .iter()
                        .map(|&f| {
                            let f = f as usize;
                            let (xs, ns) = (g.barycenters[f], g.normals[f]);
                            let (mut hi, mut lo) = (0.0, 0.0);
                            for &c in &self.l1_lumped[k] {
                                let cl = &self.set.clusters[c as usize];
                                let v = sub(xs, cl.centroid);
                                let r2 = dot(v, v);
                                let rs = r2.powf(-s);
                                let k2 = 0.5
                                    * (dot(ns, v).powi(2) + cl.projector_quadratic_form(v))
                                    / (r2 * r2);
                                hi += cl.area * rs;
                                lo += cl.area * rs * k2;
                            }
                            for &c in &self.l1_exact[k] {
                                for &t in self.set.bvh.elements_of(c as usize) {
                                    let t = t as usize;
                                    let v = sub(xs, self.set.points[t]);
                                    let r2 = dot(v, v);
                                    let rs = r2.powf(-s);
                                    let k2 = 0.5
                                        * (dot(ns, v).powi(2) + dot(normals[t], v).powi(2))
                                        / (r2 * r2);
                                    hi += self.set.weights[t] * rs;
                                    lo += self.set.weights[t] * rs * k2;
                                }
                            }
                            let c = 2.0 * weight / g.areas[f];
                            (c * hi, c * lo)
                        })
                        .collect()
                })
                .collect()
        });
        let nf = self.faces.len();
        let mut high = vec![0.0; nf];
        let mut low = vec![0.0; nf];
        for (k, row) in rows.iter().enumerate() {
            for (i, &f) in self.mesh_leaf_faces(k).iter().enumerate() {
                high[f as usize] = row[i].0;
                low[f as usize] = row[i].1;
            }
        }
        Ok(ObstacleDiagonal { high, low })
    }
}

/// One-shot `O_diag` data: builds the evaluator and calls
/// [`CrossEvaluator::metric_diagonal`].
///
/// # Errors
/// As [`CrossEvaluator::build`] and [`CrossEvaluator::metric_diagonal`].
pub fn obstacle_diagonal(
    pool: &rayon::ThreadPool,
    surface: &Surface,
    positions: &[[f64; 3]],
    set: Arc<StaticSet>,
    exponents: &Exponents,
    opts: BhOptions,
    weight: f64,
) -> Result<ObstacleDiagonal, RsError> {
    let ev = CrossEvaluator::build(
        pool,
        surface,
        positions,
        set,
        ObstacleKind::TangentPoint,
        exponents,
        opts,
    )?;
    ev.metric_diagonal(pool, exponents.s(), weight)
}
