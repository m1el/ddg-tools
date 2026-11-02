//! The face BVH and its cluster data for Barnes–Hut and the block cluster tree (`SPEC.md` §3.4,
//! §8.1).
//!
//! The hierarchy is the curves crate's generic [`Bvh<3>`] over the faces: element point the
//! barycenter `X_S`, element bounds the box of the triangle's three vertices, mass the area
//! `a_S` (leaf size 8 by default; its split rule is that crate's, D5). On top of it every node
//! `I` carries a [`Cluster`]:
//!
//! ```text
//! A_I = Σ_{S∈I} a_S,   X_I = Σ a_S X_S / A_I,   P̄_I = Σ a_S N_S N_Sᵀ / A_I,
//! box_I = AABB of the vertices of I's triangles,   r_I = |X_I − farthest corner of box_I|
//! ```
//!
//! `A_I` and `X_I` are the BVH's own node mass and moment (`Σ a`, `Σ a X`, summed bottom-up
//! by [`Bvh::refit`]); `P̄_I` comes from [`Bvh::fold_up`]. The box and the radius are those of
//! the **build** positions and are not changed by [`FaceTree::refit`] (the admissibility of a
//! step is frozen during its line search, paper Sec. 4.1, `ref:
//! optimized_cluster_tree.cpp:1313–1369`); a refit refreshes `A_I`, `X_I` and `P̄_I`.
//!
//! The Barnes–Hut admissibility (`SPEC.md` §3.4, D5; `ref: tpe_barnes_hut_0.cpp:81–100`) of
//! node `C` for leaf `ℓ` is [`bh_admissible`]: `max(r_ℓ, r_C) < θ · d_box(ℓ, C)` with
//! [`box_distance_sq`]. The block cluster tree (WP6) uses `≤` on the same quantities.

pub mod block_tree;

use ddg_repulsive_curves::spatial::{Aabb, Bvh, Element};

use crate::error::RsError;
use crate::geometry::FaceGeometry;
use crate::surface::Surface;

/// A face as an element of the [`Bvh<3>`]: barycenter, vertex box, area.
#[derive(Debug, Clone, Copy)]
pub struct FaceElement {
    /// `X_S`.
    pub barycenter: [f64; 3],
    /// The box of the triangle's three vertices.
    pub bounds: Aabb<3>,
    /// `a_S`.
    pub area: f64,
}

impl Element<3> for FaceElement {
    fn point(&self) -> [f64; 3] {
        self.barycenter
    }
    fn bounds(&self) -> Aabb<3> {
        self.bounds
    }
    fn mass(&self) -> f64 {
        self.area
    }
}

/// The elements of all faces, in face order.
pub fn face_elements(
    faces: &[[u32; 3]],
    positions: &[[f64; 3]],
    geom: &FaceGeometry,
) -> Vec<FaceElement> {
    faces
        .iter()
        .enumerate()
        .map(|(f, t)| {
            let mut b = Aabb::<3>::empty();
            for &v in t {
                b.grow(positions[v as usize]);
            }
            FaceElement {
                barycenter: geom.barycenters[f],
                bounds: b,
                area: geom.areas[f],
            }
        })
        .collect()
}

/// The data of one node `I` of a [`FaceTree`] (`SPEC.md` §3.4).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cluster {
    /// `A_I = Σ a_S` (current positions).
    pub area: f64,
    /// `X_I = Σ a_S X_S / A_I` (current positions).
    pub centroid: [f64; 3],
    /// `P̄_I = Σ a_S N_S N_Sᵀ / A_I` (current positions), the upper triangle
    /// `[xx, xy, xz, yy, yz, zz]` of the symmetric matrix (`rs-projector-average`: `vᵀP̄v` is the
    /// area average of `(N_S·v)²`).
    pub projector: [f64; 6],
    /// The box of the vertices of the node's triangles at the **build** positions.
    pub bbox: Aabb<3>,
    /// `r_I`: distance from the build-time `X_I` to the farthest corner of [`Cluster::bbox`].
    pub radius: f64,
}

impl Cluster {
    /// `vᵀ P̄_I v`.
    #[inline]
    pub fn projector_quadratic_form(&self, v: [f64; 3]) -> f64 {
        let p = &self.projector;
        p[0] * v[0] * v[0]
            + p[3] * v[1] * v[1]
            + p[5] * v[2] * v[2]
            + 2.0 * (p[1] * v[0] * v[1] + p[2] * v[0] * v[2] + p[4] * v[1] * v[2])
    }
}

/// Squared distance between two boxes: `d² = Σ_k max(0, max(lo_a, lo_b) − min(hi_a, hi_b))²`
/// (`ref: tpe_barnes_hut_0.cpp:87–96`). Zero when the boxes overlap or touch; a lower bound of
/// the squared distance between any two points of the boxes.
///
/// [formula: rs-box-distance]
#[inline]
pub fn box_distance_sq(a: &Aabb<3>, b: &Aabb<3>) -> f64 {
    let mut d2 = 0.0;
    for k in 0..3 {
        let d = (a.lo[k].max(b.lo[k]) - a.hi[k].min(b.hi[k])).max(0.0);
        d2 += d * d;
    }
    d2
}

/// Distance from `x` to the farthest corner of `b`: `sqrt(Σ_k max(hi_k − x_k, x_k − lo_k)²)`
/// (`ref: optimized_cluster_tree.cpp:532–540`).
#[inline]
pub fn farthest_corner_distance(b: &Aabb<3>, x: [f64; 3]) -> f64 {
    let mut r2 = 0.0;
    for k in 0..3 {
        let d = (b.hi[k] - x[k]).max(x[k] - b.lo[k]);
        r2 += d * d;
    }
    r2.sqrt()
}

/// The Barnes–Hut admissibility of `SPEC.md` §3.4 (strict): `max(r_ℓ, r_C) < θ · d_box(ℓ, C)`,
/// evaluated as `max(r_ℓ, r_C)² < θ² d²` (`ref: tpe_barnes_hut_0.cpp:100`). Boxes that overlap
/// or touch (`d = 0`) are never admissible; `θ = 0` admits nothing.
///
/// [formula: rs-bh-energy]
/// [formula: rs-box-distance]
#[inline]
pub fn bh_admissible(leaf: &Cluster, c: &Cluster, theta: f64) -> bool {
    let h = leaf.radius.max(c.radius);
    h * h < theta * theta * box_distance_sq(&leaf.bbox, &c.bbox)
}

/// The face BVH with its [`Cluster`] data (`SPEC.md` §3.4, §8.1). Built once per step,
/// [refitted](FaceTree::refit) for line-search trials.
#[derive(Debug, Clone)]
pub struct FaceTree {
    bvh: Bvh<3>,
    clusters: Vec<Cluster>,
    /// Node indices of the leaves, increasing.
    leaves: Vec<u32>,
    /// For every node, its position in `leaves` (`u32::MAX` for internal nodes).
    leaf_ordinal: Vec<u32>,
    /// For every face, the node index of its leaf.
    leaf_of_face: Vec<u32>,
}

impl FaceTree {
    /// Builds the hierarchy over the faces of `surface` at `positions` (`geom` must be the face
    /// geometry of `positions`, e.g. [`FaceGeometry::new`], which validates them).
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] for `leaf_size = 0`; [`RsError::PositionCountMismatch`];
    /// [`RsError::Curves`] if the BVH cannot be built (an empty surface).
    pub fn build(
        surface: &Surface,
        positions: &[[f64; 3]],
        geom: &FaceGeometry,
        leaf_size: usize,
    ) -> Result<FaceTree, RsError> {
        if leaf_size == 0 {
            return Err(RsError::InvalidOption {
                name: "leaf_size",
                value: "0".into(),
                requirement: "at least 1",
            });
        }
        check_counts(surface.faces(), surface.n_vertices(), positions, geom)?;
        let elems = face_elements(surface.faces(), positions, geom);
        let bvh = Bvh::build(&elems, leaf_size)?;
        let nn = bvh.nodes().len();
        let mut leaves = Vec::new();
        let mut leaf_ordinal = vec![u32::MAX; nn];
        let mut leaf_of_face = vec![0u32; surface.n_faces()];
        for (n, node) in bvh.nodes().iter().enumerate() {
            if node.is_leaf() {
                leaf_ordinal[n] = leaves.len() as u32;
                leaves.push(n as u32);
                for &f in bvh.elements_of(n) {
                    leaf_of_face[f as usize] = n as u32;
                }
            }
        }
        let mut tree = FaceTree {
            clusters: Vec::with_capacity(nn),
            bvh,
            leaves,
            leaf_ordinal,
            leaf_of_face,
        };
        let projectors = tree.projectors(geom);
        tree.clusters = tree
            .bvh
            .nodes()
            .iter()
            .zip(projectors)
            .map(|(nd, projector)| {
                let centroid = nd.centroid();
                Cluster {
                    area: nd.mass,
                    centroid,
                    projector,
                    bbox: nd.bbox,
                    radius: farthest_corner_distance(&nd.bbox, centroid),
                }
            })
            .collect();
        Ok(tree)
    }

    /// Moves the tree to `positions` (with their face geometry `geom`), keeping the topology,
    /// the boxes and the radii of the build: refreshes `A_I`, `X_I`, `P̄_I` (`SPEC.md` §3.4,
    /// §7.2). Degenerate faces are not an error here; they give non-finite node data.
    ///
    /// # Errors
    /// [`RsError::PositionCountMismatch`] (also for a `geom` of another face count).
    pub fn refit(
        &mut self,
        surface: &Surface,
        positions: &[[f64; 3]],
        geom: &FaceGeometry,
    ) -> Result<(), RsError> {
        check_counts(surface.faces(), surface.n_vertices(), positions, geom)?;
        self.refit_faces(surface.faces(), positions, geom)
    }

    /// As [`FaceTree::refit`] with the face list of the surface (`positions` and `geom` must
    /// match it; vertex indices are not checked again).
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] for a face list or `geom` of another face count than the
    /// build.
    pub fn refit_faces(
        &mut self,
        faces: &[[u32; 3]],
        positions: &[[f64; 3]],
        geom: &FaceGeometry,
    ) -> Result<(), RsError> {
        if faces.len() != self.bvh.n_elements() || geom.areas.len() != faces.len() {
            return Err(RsError::InvalidOption {
                name: "faces",
                value: format!("{} faces", faces.len()),
                requirement: "the face count of the surface the tree was built for",
            });
        }
        let elems = face_elements(faces, positions, geom);
        self.bvh.refit(&elems)?;
        let projectors = self.projectors(geom);
        for ((c, nd), p) in self
            .clusters
            .iter_mut()
            .zip(self.bvh.nodes())
            .zip(projectors)
        {
            c.area = nd.mass;
            c.centroid = nd.centroid();
            c.projector = p;
        }
        Ok(())
    }

    /// `P̄_I` of every node by [`Bvh::fold_up`] (sums of `a N Nᵀ`, divided by `A_I`).
    fn projectors(&self, geom: &FaceGeometry) -> Vec<[f64; 6]> {
        let sums = self.bvh.fold_up(
            |idx| {
                let mut s = [0.0; 7];
                for &f in idx {
                    let a = geom.areas[f as usize];
                    let n = geom.normals[f as usize];
                    s[0] += a;
                    s[1] += a * n[0] * n[0];
                    s[2] += a * n[0] * n[1];
                    s[3] += a * n[0] * n[2];
                    s[4] += a * n[1] * n[1];
                    s[5] += a * n[1] * n[2];
                    s[6] += a * n[2] * n[2];
                }
                s
            },
            |l, r| std::array::from_fn(|k| l[k] + r[k]),
        );
        sums.iter()
            .map(|s| std::array::from_fn(|k| s[k + 1] / s[0]))
            .collect()
    }

    /// The underlying hierarchy (nodes in pre-order, node 0 the root; its node boxes, masses and
    /// moments are those of the last build or refit).
    pub fn bvh(&self) -> &Bvh<3> {
        &self.bvh
    }

    /// The cluster data, indexed like [`Bvh::nodes`].
    pub fn clusters(&self) -> &[Cluster] {
        &self.clusters
    }

    /// The node indices of the leaves, in increasing order.
    pub fn leaves(&self) -> &[u32] {
        &self.leaves
    }

    /// The position of node `n` in [`FaceTree::leaves`], `None` for an internal node.
    pub fn leaf_ordinal(&self, n: usize) -> Option<usize> {
        let k = self.leaf_ordinal[n];
        (k != u32::MAX).then_some(k as usize)
    }

    /// The node index of the leaf that holds face `f`.
    pub fn leaf_of_face(&self, f: usize) -> usize {
        self.leaf_of_face[f] as usize
    }

    /// The faces of node `n` (a contiguous range of [`Bvh::order`]).
    pub fn faces_of(&self, n: usize) -> &[u32] {
        self.bvh.elements_of(n)
    }

    /// Number of faces.
    pub fn n_faces(&self) -> usize {
        self.bvh.n_elements()
    }
}

fn check_counts(
    faces: &[[u32; 3]],
    n_vertices: usize,
    positions: &[[f64; 3]],
    geom: &FaceGeometry,
) -> Result<(), RsError> {
    if positions.len() != n_vertices {
        return Err(RsError::PositionCountMismatch {
            expected: n_vertices,
            found: positions.len(),
        });
    }
    let nf = faces.len();
    if geom.areas.len() != nf || geom.normals.len() != nf || geom.barycenters.len() != nf {
        return Err(RsError::InvalidOption {
            name: "geometry",
            value: format!("face data of {} faces", geom.areas.len()),
            requirement: "the face geometry of the surface (one entry per face)",
        });
    }
    Ok(())
}
