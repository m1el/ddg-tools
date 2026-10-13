//! Handing an intrinsic triangulation to code that takes a `TriMesh` with a `Geometry`
//! (SPEC §6), and carrying a cut found there back to the input mesh.
//!
//! * **Export.** A triangulation that is not a simplicial complex (a self-edge, several edges
//!   between two vertices, a face that uses a vertex twice) is made one by splitting the
//!   offending edges at their midpoints, which is an intrinsic operation: the surface does not
//!   change. The result is the triangulation itself, handed out as a `TriMesh` with the same
//!   vertex, face and corner ids, a `Geometry` from the intrinsic edge lengths, and the place
//!   of every vertex on the input mesh. The number of splits is reported.
//! * **The way back.** A cut on the exported mesh consists of points on its edges, points
//!   inside its faces where three regions meet, and straight segments across its faces between
//!   them. [`SimplicialExport::edge_point`] turns a point of an exported edge into a point of
//!   the triangulation; [`IntrinsicTriangulation::segments_on_input`] then gives every node
//!   its place on the input mesh and every segment its pieces, each inside one input face.
//! * **State.** Values on the vertices of an exported mesh (the level-set functions of the cut
//!   flow) belong to that export. [`SimplicialExport::fingerprint`] identifies it: a hash of
//!   its connectivity, its edge lengths and the places of its vertices. Refinement and export
//!   are deterministic, so the same input and options give the same fingerprint again;
//!   [`SimplicialExport::check_state`] refuses a state made for another one.

use ddg_core::{Geometry, TriMesh};

use crate::delaunay::{DelaunayOptions, DelaunayStatus};
use crate::error::IntrinsicError;
use crate::point::SurfacePoint;
use crate::triangulation::IntrinsicTriangulation;

/// Options of [`IntrinsicTriangulation::export_simplicial`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ExportOptions {
    /// Whether, and with which test, the triangulation is flipped back to Delaunay after
    /// edges were split (a split can leave neighbouring edges that are not Delaunay). The
    /// splits and the flips alternate until the triangulation is a simplicial complex that
    /// needs no flip. `None`: no flips, the result may fail the Delaunay test at edges next to
    /// the splits.
    pub delaunay: Option<DelaunayOptions>,
    /// The export fails when more than this many edges would have to be split.
    pub max_splits: usize,
}

impl Default for ExportOptions {
    /// Flips with the default Delaunay options; at most 100 000 splits.
    fn default() -> Self {
        ExportOptions {
            delaunay: Some(DelaunayOptions::default()),
            max_splits: 100_000,
        }
    }
}

/// An intrinsic triangulation as a simplicial mesh. Returned by
/// [`IntrinsicTriangulation::export_simplicial`].
#[derive(Debug, Clone, PartialEq)]
pub struct SimplicialExport {
    /// The connectivity. Vertices, faces, corners and halfedges have the ids of the
    /// triangulation; edges are numbered by the `TriMesh` (sorted vertex pairs).
    pub mesh: TriMesh,
    /// The geometry, from the intrinsic edge lengths.
    pub geometry: Geometry,
    /// The place of every vertex on the input mesh. The first vertices are the input
    /// vertices; those inserted by refinement and by the export follow.
    pub locations: Vec<SurfacePoint>,
    /// Number of edges that were split to make the triangulation a simplicial complex.
    pub n_splits: usize,
    /// Number of flips made after splits.
    pub n_flips: usize,
    /// Whether every interior edge passes the Delaunay test of the options (`true` without
    /// test if none was asked for and nothing was split).
    pub is_delaunay: bool,
    /// What identifies this export (see the [module documentation](self)).
    pub fingerprint: u64,
    /// Per edge of `mesh`: the edge of the triangulation.
    edge: Vec<u32>,
    /// Per edge of `mesh`: whether its lower-numbered vertex is the tail of that edge.
    forward: Vec<bool>,
}

impl SimplicialExport {
    /// The edge of the triangulation that is edge `mesh_edge` of [`mesh`](Self::mesh).
    pub fn triangulation_edge(&self, mesh_edge: usize) -> usize {
        self.edge[mesh_edge] as usize
    }

    /// The point at parameter `t` of edge `mesh_edge` of [`mesh`](Self::mesh), measured from
    /// its lower-numbered vertex (as the cut flow gives its cut points), as a point of the
    /// triangulation.
    pub fn edge_point(&self, mesh_edge: usize, t: f64) -> SurfacePoint {
        SurfacePoint::Edge {
            edge: self.edge[mesh_edge],
            t: if self.forward[mesh_edge] { t } else { 1.0 - t },
        }
    }

    /// Checks that a state with `n_values` values per vertex function, made for the export
    /// with the given fingerprint, belongs to this export.
    ///
    /// # Errors
    /// [`IntrinsicError::StaleState`] for another fingerprint and
    /// [`IntrinsicError::IndexOutOfRange`] for a wrong number of values.
    pub fn check_state(&self, fingerprint: u64, n_values: usize) -> Result<(), IntrinsicError> {
        if fingerprint != self.fingerprint {
            return Err(IntrinsicError::StaleState {
                expected: self.fingerprint,
                got: fingerprint,
            });
        }
        if n_values != self.mesh.n_vertices() {
            return Err(IntrinsicError::IndexOutOfRange {
                what: "number of state values",
                index: n_values,
                bound: self.mesh.n_vertices() + 1,
            });
        }
        Ok(())
    }
}

/// FNV-1a over 64-bit words.
struct Hasher(u64);

impl Hasher {
    fn new() -> Hasher {
        Hasher(0xcbf2_9ce4_8422_2325)
    }
    fn word(&mut self, w: u64) {
        for byte in w.to_le_bytes() {
            self.0 ^= u64::from(byte);
            self.0 = self.0.wrapping_mul(0x0000_0100_0000_01b3);
        }
    }
}

impl IntrinsicTriangulation {
    /// A hash of the triangulation: its faces, its edge lengths (per halfedge) and the places
    /// of its vertices on the input mesh. Two triangulations made by the same operations from
    /// the same input have the same fingerprint; any flip, insertion or split changes it
    /// (up to hash collisions, which are not excluded).
    pub fn fingerprint(&self) -> u64 {
        let mut hasher = Hasher::new();
        hasher.word(self.n_vertices() as u64);
        hasher.word(self.n_faces() as u64);
        for (f, face) in self.faces.iter().enumerate() {
            for (k, &vertex) in face.iter().enumerate() {
                hasher.word(u64::from(vertex));
                hasher.word(self.edge_lengths[self.halfedge_edge[3 * f + k] as usize].to_bits());
            }
        }
        for location in &self.locations {
            match *location {
                SurfacePoint::Vertex(v) => {
                    hasher.word(0);
                    hasher.word(u64::from(v));
                }
                SurfacePoint::Edge { edge, t } => {
                    hasher.word(1);
                    hasher.word(u64::from(edge));
                    hasher.word(t.to_bits());
                }
                SurfacePoint::Face { face, bary } => {
                    hasher.word(2);
                    hasher.word(u64::from(face));
                    for b in bary {
                        hasher.word(b.to_bits());
                    }
                }
            }
        }
        hasher.0
    }

    /// Makes the triangulation a simplicial complex by splitting edges, and hands it out as a
    /// `TriMesh` with a `Geometry` (SPEC §6).
    ///
    /// Every edge from a vertex to itself is split at its midpoint, and of every group of
    /// edges between the same two vertices all but one. This is repeated (a split can make a
    /// new group) until nothing is left, then, if `options.delaunay` is set, the triangulation
    /// is flipped to Delaunay, which can undo it; splits and flips alternate until both hold.
    /// The triangulation is changed by this; the export is its state afterwards, and the way
    /// back ([`segments_on_input`](Self::segments_on_input)) must be asked of the
    /// triangulation in that state.
    ///
    /// # Errors
    /// [`IntrinsicError::ExportFailed`] if more than `options.max_splits` splits would be
    /// needed, if an edge that must be split cannot be (its end point has degree one and no
    /// other edge can be split first), or if two faces have the same three vertices (a closed
    /// component of two faces).
    pub fn export_simplicial(
        &mut self,
        options: &ExportOptions,
    ) -> Result<SimplicialExport, IntrinsicError> {
        let mut n_splits = 0;
        let mut n_flips = 0;
        let mut is_delaunay = options.delaunay.is_none();
        loop {
            // Split until simplicial.
            loop {
                let irregular = self.irregularities();
                let mut to_split: Vec<u32> = irregular.self_edges.clone();
                // Of each group of parallel edges, all but the first.
                let mut pairs: Vec<([usize; 2], u32)> = irregular
                    .multi_edges
                    .iter()
                    .map(|&e| {
                        let [a, b] = self.edge_vertices(e as usize);
                        ([a.min(b), a.max(b)], e)
                    })
                    .collect();
                pairs.sort_unstable();
                for w in pairs.windows(2) {
                    if w[0].0 == w[1].0 {
                        to_split.push(w[1].1);
                    }
                }
                if to_split.is_empty() {
                    if let Some(&f) = irregular.faces_with_repeated_vertex.first() {
                        return Err(IntrinsicError::ExportFailed {
                            n_splits,
                            reason: format!(
                                "face {f} uses a vertex more than once although no edge is a \
                                 self-edge or one of several between two vertices"
                            ),
                        });
                    }
                    break;
                }
                let mut progress = false;
                let mut last_refusal = String::new();
                for e in to_split {
                    if n_splits >= options.max_splits {
                        return Err(IntrinsicError::ExportFailed {
                            n_splits,
                            reason: format!(
                                "more than {} splits would be needed",
                                options.max_splits
                            ),
                        });
                    }
                    match self.split_edge(e as usize, 0.5) {
                        Ok(_) => {
                            n_splits += 1;
                            progress = true;
                        }
                        Err(IntrinsicError::InsertionRefused { reason, .. }) => {
                            last_refusal = format!("edge {e} cannot be split: {reason}");
                        }
                        Err(other) => return Err(other),
                    }
                }
                if !progress {
                    return Err(IntrinsicError::ExportFailed {
                        n_splits,
                        reason: last_refusal,
                    });
                }
            }
            let Some(delaunay) = &options.delaunay else {
                break;
            };
            let report = self.flip_to_delaunay(delaunay)?;
            n_flips += report.n_flips;
            is_delaunay = report.status == DelaunayStatus::Delaunay;
            if report.n_flips == 0 || self.irregularities().is_empty() {
                break;
            }
        }
        let (mesh, geometry) = self
            .to_tri_mesh()
            .map_err(|e| IntrinsicError::ExportFailed {
                n_splits,
                reason: e.to_string(),
            })?;
        let mut edge = vec![0u32; mesh.n_edges()];
        let mut forward = vec![true; mesh.n_edges()];
        for e in 0..self.n_edges() {
            let [tail, head] = self.edge_vertices(e);
            // The faces are the same, so the halfedge ids are.
            let mesh_edge = mesh.edge(self.edge_halfedge[e] as usize);
            edge[mesh_edge] = e as u32;
            forward[mesh_edge] = tail < head;
        }
        Ok(SimplicialExport {
            mesh,
            geometry,
            locations: self.locations.clone(),
            n_splits,
            n_flips,
            is_delaunay,
            fingerprint: self.fingerprint(),
            edge,
            forward,
        })
    }
}
