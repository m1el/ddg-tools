//! Test-mesh generators for the experiments and applications of the paper (`SPEC.md` §14.4, WP12).
//!
//! Every generator is deterministic, takes sizes as parameters and returns a validated [`Mesh`]
//! (positions and triangles; [`Mesh::surface`] gives the [`Surface`]). Faces are oriented
//! outward (counter-clockwise seen from outside) on the closed meshes.
//!
//! | module | content |
//! |---|---|
//! | [`primitives`] | box, capsule, sphere, torus, punctured torus, slab with spheres (Fig. 20, 25), soup with holes (Fig. 23), coarse "H" (Fig. 22) |
//! | [`field`] | distance fields: tubes around polylines and networks, smooth union / difference, boxes, spheres |
//! | [`tetra`] | marching tetrahedra on a grid (manifold by construction), clamping and projected smoothing; [`tetra::mesh_with_faces`] hits a face count |
//! | [`shapes`] | genus-g surfaces (Fig. 15), handcuffs, tunnels (Fig. 12), tube around a trivalent graph (for Fig. 16) |
//! | [`perturbed`] | the smooth perturbed torus of Sec. 8.1 (Fig. 10) |
//! | [`check`] | self-intersection test, linking number of closed polylines |
//!
//! # Not done
//! * The knotted genus-2 surfaces of Ishii et al. 2012 (Fig. 16): the table's diagrams are not
//!   available offline and were not transcribed (never invented); [`shapes::graph_tube`] meshes
//!   any trivalent spatial graph given as polylines, so entries can be added from the paper.
//! * No remeshing: the meshes come out of marching tetrahedra with clamped vertices and
//!   projected smoothing ([`tetra`]); callers who need better triangles use the WP4 remesher.

pub mod check;
pub mod field;
pub mod perturbed;
pub mod primitives;
pub mod shapes;
pub mod tetra;

use crate::error::RsError;
use crate::surface::Surface;

/// A triangle mesh: positions and faces (validated by the generators).
#[derive(Debug, Clone, PartialEq)]
pub struct Mesh {
    /// Vertex positions.
    pub positions: Vec<[f64; 3]>,
    /// Triangles.
    pub faces: Vec<[u32; 3]>,
}

impl Mesh {
    /// Validates topology and positions.
    ///
    /// # Errors
    /// The errors of [`Surface::new`] and [`Surface::validate_positions`].
    pub fn validated(self) -> Result<Mesh, RsError> {
        let s = Surface::new(self.positions.len(), self.faces.clone())?;
        s.validate_positions(&self.positions)?;
        Ok(self)
    }

    /// The validated topology.
    ///
    /// # Errors
    /// As [`Surface::new`].
    pub fn surface(&self) -> Result<Surface, RsError> {
        Surface::new(self.positions.len(), self.faces.clone())
    }

    /// Disjoint union with `other` (components are numbered by lowest face index: `self` first).
    pub fn merged(mut self, other: &Mesh) -> Mesh {
        let off = self.positions.len() as u32;
        self.positions.extend_from_slice(&other.positions);
        self.faces
            .extend(other.faces.iter().map(|f| f.map(|i| i + off)));
        self
    }

    /// Translates and scales (`p ↦ s p + t`).
    pub fn transformed(mut self, scale: f64, t: [f64; 3]) -> Mesh {
        for p in &mut self.positions {
            for k in 0..3 {
                p[k] = scale * p[k] + t[k];
            }
        }
        self
    }

    /// Total area.
    pub fn area(&self) -> f64 {
        use crate::geometry::v3::{cross, norm, sub};
        self.faces
            .iter()
            .map(|f| {
                let [a, b, c] = f.map(|i| self.positions[i as usize]);
                0.5 * norm(cross(sub(b, a), sub(c, a)))
            })
            .sum()
    }

    /// Signed volume (closed, consistently oriented meshes).
    pub fn volume(&self) -> f64 {
        use crate::geometry::v3::{cross, dot};
        self.faces
            .iter()
            .map(|f| {
                let [a, b, c] = f.map(|i| self.positions[i as usize]);
                dot(a, cross(b, c)) / 6.0
            })
            .sum()
    }
}

pub(crate) fn bad(name: &'static str, value: impl ToString, requirement: &'static str) -> RsError {
    RsError::InvalidOption {
        name,
        value: value.to_string(),
        requirement,
    }
}

/// The SplitMix64 generator (reproducible across platforms).
#[derive(Debug, Clone)]
pub struct Rng(u64);

impl Rng {
    /// Seeded generator.
    pub fn new(seed: u64) -> Rng {
        Rng(seed)
    }
    /// 64 random bits.
    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    /// Uniform in `[0, 1)`.
    pub fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    /// Uniform in `[lo, hi)`.
    pub fn uniform(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }
}
