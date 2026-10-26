//! A bounding volume hierarchy over weighted points in `R^D` (`SPEC.md` §8.1, §16, D15).
//!
//! The tree is generic in the dimension `D` and in the element type (trait [`Element`]); there is
//! no curve type in any signature, so the later surface papers can reuse it. For curves the
//! elements are the edges' tangent-points `(x_I, T_I) ∈ R⁶` with mass `ℓ_I`.
//!
//! * **Build** (top down, App. B.1.1): the split axis cycles through the `D` coordinates with the
//!   depth, starting at coordinate 0; along the elements sorted by that coordinate, the split
//!   position minimizes the sum of the squared diameters of the two children's boxes (prefix and
//!   suffix boxes, `O(n)` per node after the sort), restricted to positions that leave at least a
//!   quarter of the elements on each side ([`DEFAULT_MIN_SPLIT_FRACTION`]; the unrestricted
//!   rule peels single outliers and gives a degenerate tree on curves, SPEC §17d). The reference
//!   (`src/spatial/tpe_bvh.cpp`, MIT, Copyright (c) 2019 Christopher Yu) splits at the position
//!   minimizing squared 1-D widths along the axis; this one follows the paper (D15). Leaves hold
//!   at most `leaf_size` elements. The build is sequential and deterministic (ties broken by
//!   element index), and iterative (no recursion, so degenerate unbalanced trees cannot overflow
//!   the stack).
//! * **Layout**: nodes in pre-order, the left child of node `i` is `i + 1`; the elements of a node
//!   are the contiguous range `start..end` of [`Bvh::order`].
//! * **Node data** [formula: rc-bvh-lumped]: the box of the element bounds, the mass `Σ m` and the
//!   first moment `Σ m p` (so that the mass-weighted centroid `Σ m p / Σ m` of a parent is the
//!   combination of its children's, with the moments summed and not the means).
//! * **Refit** keeps the topology and recomputes the node data from the new elements (D26).
//!   [`Bvh::fold_up`] runs a user-defined bottom-up aggregation over the same topology for extra
//!   node payloads.

use crate::error::RcError;
use crate::spatial::aabb::Aabb;

/// The default of [`Bvh::build`]: each child of a split gets at least a quarter of the elements.
pub const DEFAULT_MIN_SPLIT_FRACTION: f64 = 0.25;

/// An element of a [`Bvh`]: a weighted point with a bounding box.
pub trait Element<const D: usize> {
    /// The point in `R^D` used to split (for curves: `(x_I, T_I)`).
    fn point(&self) -> [f64; D];
    /// A box that contains the element (default: the degenerate box at [`Element::point`]).
    fn bounds(&self) -> Aabb<D> {
        Aabb::from_point(self.point())
    }
    /// The mass of the element (for curves: the edge length `ℓ_I`).
    fn mass(&self) -> f64;
}

/// A node of the hierarchy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Node<const D: usize> {
    /// Union of the element bounds.
    pub bbox: Aabb<D>,
    /// `Σ m`.
    pub mass: f64,
    /// `Σ m p` (first moment of the element points).
    pub moment: [f64; D],
    /// First position of the node's elements in [`Bvh::order`].
    pub start: u32,
    /// One past the last position.
    pub end: u32,
    /// Index of the right child; `0` for a leaf (the root is never a child). The left child is
    /// the node after this one.
    pub right: u32,
}

impl<const D: usize> Node<D> {
    /// `true` for a leaf.
    #[inline]
    pub fn is_leaf(&self) -> bool {
        self.right == 0
    }

    /// Number of elements below the node.
    #[inline]
    pub fn len(&self) -> usize {
        (self.end - self.start) as usize
    }

    /// Always `false`: nodes are never empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.end == self.start
    }

    /// The mass-weighted centroid `Σ m p / Σ m` (NaN if the mass is 0).
    #[inline]
    pub fn centroid(&self) -> [f64; D] {
        let mut c = self.moment;
        for x in &mut c {
            *x /= self.mass;
        }
        c
    }
}

/// A bounding volume hierarchy over `n` elements of `R^D`.
#[derive(Debug, Clone, PartialEq)]
pub struct Bvh<const D: usize> {
    nodes: Vec<Node<D>>,
    order: Vec<u32>,
    leaf_size: usize,
}

fn opt_err(name: &'static str, value: String, requirement: &'static str) -> RcError {
    RcError::InvalidOption {
        name,
        value,
        requirement,
    }
}

impl<const D: usize> Bvh<D> {
    /// Builds the hierarchy with the default balance bound [`DEFAULT_MIN_SPLIT_FRACTION`].
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for `leaf_size = 0` or an empty element list.
    pub fn build<E: Element<D>>(elements: &[E], leaf_size: usize) -> Result<Self, RcError> {
        Self::build_with(elements, leaf_size, DEFAULT_MIN_SPLIT_FRACTION)
    }

    /// Builds the hierarchy; every split leaves at least `⌈min_fraction · n⌉` elements (and at
    /// least one) on each side of a node with `n` elements. `min_fraction = 0` is the paper's
    /// rule as printed (the sum of squared diameters over all split positions), which on curves
    /// peels single outliers off the node and builds a tree of depth `O(n)` (D15, SPEC §17d).
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for `leaf_size = 0`, an empty element list, or `min_fraction`
    /// outside `[0, ½]`.
    pub fn build_with<E: Element<D>>(
        elements: &[E],
        leaf_size: usize,
        min_fraction: f64,
    ) -> Result<Self, RcError> {
        if !(0.0..=0.5).contains(&min_fraction) {
            return Err(opt_err(
                "min_split_fraction",
                min_fraction.to_string(),
                "a number in [0, 0.5]",
            ));
        }
        if leaf_size == 0 {
            return Err(opt_err(
                "leaf_size",
                "0".into(),
                "a leaf holds at least 1 element",
            ));
        }
        if elements.is_empty() {
            return Err(opt_err("elements", "[]".into(), "at least 1 element"));
        }
        let n = elements.len();
        let pts: Vec<[f64; D]> = elements.iter().map(|e| e.point()).collect();
        let bnds: Vec<Aabb<D>> = elements.iter().map(|e| e.bounds()).collect();
        let mut order: Vec<u32> = (0..n as u32).collect();
        let mut nodes: Vec<Node<D>> = Vec::with_capacity(2 * n / leaf_size.max(1) + 1);
        let mut suffix: Vec<Aabb<D>> = Vec::new();
        // (start, end, depth, parent whose `right` to patch (u32::MAX: none))
        let mut stack = vec![(0usize, n, 0usize, u32::MAX)];
        while let Some((start, end, depth, patch)) = stack.pop() {
            let id = nodes.len() as u32;
            if patch != u32::MAX {
                nodes[patch as usize].right = id;
            }
            let slice = &mut order[start..end];
            nodes.push(Node {
                bbox: Aabb::empty(),
                mass: 0.0,
                moment: [0.0; D],
                start: start as u32,
                end: end as u32,
                right: 0,
            });
            let cnt = end - start;
            if cnt <= leaf_size {
                continue;
            }
            let axis = depth % D;
            slice.sort_unstable_by(|&a, &b| {
                pts[a as usize][axis]
                    .total_cmp(&pts[b as usize][axis])
                    .then(a.cmp(&b))
            });
            // suffix[k] = box of slice[k..]
            suffix.clear();
            suffix.resize(cnt, Aabb::empty());
            let mut acc = Aabb::empty();
            for k in (0..cnt).rev() {
                acc = acc.union(&bnds[slice[k] as usize]);
                suffix[k] = acc;
            }
            let mut prefix = Aabb::empty();
            let mut best = (f64::INFINITY, cnt / 2);
            let lo = ((min_fraction * cnt as f64).ceil() as usize).clamp(1, cnt / 2);
            for k in 1..=cnt - lo {
                prefix = prefix.union(&bnds[slice[k - 1] as usize]);
                if k < lo {
                    continue;
                }
                let cost = prefix.diameter_sq() + suffix[k].diameter_sq();
                // A NaN or infinite cost never wins; the median split is the fallback.
                if cost < best.0 {
                    best = (cost, k);
                }
            }
            let split = start + best.1;
            stack.push((split, end, depth + 1, id));
            stack.push((start, split, depth + 1, u32::MAX));
        }
        let mut bvh = Bvh {
            nodes,
            order,
            leaf_size,
        };
        // Node data comes from the same bottom-up pass as a refit, so that a refit on unchanged
        // elements reproduces the build bit for bit.
        bvh.refit(elements)?;
        Ok(bvh)
    }

    /// Recomputes the node data from the (moved) elements, keeping the topology (D26).
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] if the number of elements differs from the build.
    pub fn refit<E: Element<D>>(&mut self, elements: &[E]) -> Result<(), RcError> {
        if elements.len() != self.order.len() {
            return Err(opt_err(
                "elements",
                elements.len().to_string(),
                "the number of elements must equal the number at build time",
            ));
        }
        let pts: Vec<[f64; D]> = elements.iter().map(|e| e.point()).collect();
        let bnds: Vec<Aabb<D>> = elements.iter().map(|e| e.bounds()).collect();
        for i in (0..self.nodes.len()).rev() {
            let nd = self.nodes[i];
            let (bbox, mass, moment) = if nd.is_leaf() {
                summarize(
                    &self.order[nd.start as usize..nd.end as usize],
                    &pts,
                    &bnds,
                    elements,
                )
            } else {
                let (l, r) = (&self.nodes[i + 1], &self.nodes[nd.right as usize]);
                let mut mom = l.moment;
                for k in 0..D {
                    mom[k] += r.moment[k];
                }
                (l.bbox.union(&r.bbox), l.mass + r.mass, mom)
            };
            let nd = &mut self.nodes[i];
            nd.bbox = bbox;
            nd.mass = mass;
            nd.moment = moment;
        }
        Ok(())
    }

    /// Bottom-up aggregation of a node payload `S`: `leaf` maps the element indices of a leaf,
    /// `merge` combines the payloads of the left and right child. The result is indexed like
    /// [`Bvh::nodes`]. Sequential, so deterministic.
    pub fn fold_up<S>(&self, leaf: impl Fn(&[u32]) -> S, merge: impl Fn(&S, &S) -> S) -> Vec<S> {
        let mut out: Vec<Option<S>> = (0..self.nodes.len()).map(|_| None).collect();
        for i in (0..self.nodes.len()).rev() {
            let nd = &self.nodes[i];
            let s = if nd.is_leaf() {
                leaf(&self.order[nd.start as usize..nd.end as usize])
            } else {
                match (&out[i + 1], &out[nd.right as usize]) {
                    (Some(a), Some(b)) => merge(a, b),
                    _ => unreachable!("children have larger indices than their parent"),
                }
            };
            out[i] = Some(s);
        }
        out.into_iter()
            .map(|s| s.expect("every node visited"))
            .collect()
    }

    /// The nodes in pre-order; node 0 is the root.
    pub fn nodes(&self) -> &[Node<D>] {
        &self.nodes
    }

    /// The permutation of the element indices: node `N` owns `order[N.start..N.end]`.
    pub fn order(&self) -> &[u32] {
        &self.order
    }

    /// The element indices of node `i`.
    pub fn elements_of(&self, i: usize) -> &[u32] {
        let nd = &self.nodes[i];
        &self.order[nd.start as usize..nd.end as usize]
    }

    /// The leaf size the tree was built with.
    pub fn leaf_size(&self) -> usize {
        self.leaf_size
    }

    /// Number of elements.
    pub fn n_elements(&self) -> usize {
        self.order.len()
    }

    /// Index of the left child of an internal node.
    #[inline]
    pub fn left(&self, i: usize) -> usize {
        i + 1
    }
}

/// Box, mass and moment of the elements `idx`, accumulated in the order of `idx`.
fn summarize<const D: usize, E: Element<D>>(
    idx: &[u32],
    pts: &[[f64; D]],
    bnds: &[Aabb<D>],
    elements: &[E],
) -> (Aabb<D>, f64, [f64; D]) {
    let mut bbox = Aabb::empty();
    let mut mass = 0.0;
    let mut mom = [0.0; D];
    for &e in idx {
        let m = elements[e as usize].mass();
        bbox = bbox.union(&bnds[e as usize]);
        mass += m;
        for k in 0..D {
            mom[k] += m * pts[e as usize][k];
        }
    }
    (bbox, mass, mom)
}
