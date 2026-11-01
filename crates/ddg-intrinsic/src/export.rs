//! Handing an intrinsic triangulation to code that takes a `TriMesh` with a `Geometry`
//! (SPEC §6), and carrying a cut found there back to the input mesh.
//!
//! * **Export.** A triangulation that is not a simplicial complex (a self-edge, several edges
//!   between two vertices, a face that uses a vertex twice) is made one by splitting the
//!   offending edges at their midpoints, which is an intrinsic operation: the surface does not
//!   change. The result is the triangulation itself, handed out as a `TriMesh` with the same
//!   vertex, face and corner ids, a `Geometry` from the intrinsic edge lengths, and the place
//!   of every vertex on the input mesh. The number of splits is reported, and the smallest
//!   and largest angle of the result: a triangulation that needed no split keeps its angles,
//!   one that needed splits keeps the Delaunay property but not an angle bound (SPEC N15).
//! * **The way back.** A cut on the exported mesh consists of points on its edges, points
//!   inside its faces where three regions meet, and straight segments across its faces between
//!   them. [`SimplicialExport::edge_point`] turns a point of an exported edge into a point of
//!   the triangulation; [`IntrinsicTriangulation::segments_on_input`] then gives every node
//!   its place on the input mesh and every segment its pieces, each inside one input face.
//! * **State.** Values on the vertices of an exported mesh (the level-set functions of the cut
//!   flow) belong to that export. [`SimplicialExport::fingerprint`] identifies it: a hash of
//!   its connectivity, its edge lengths and the places of its vertices (for integer
//!   coordinates also the normal coordinates and roundabouts, see
//!   [`Exportable::state_fingerprint`]). Refinement and export are deterministic, so the same
//!   input and options give the same fingerprint again; [`SimplicialExport::check_state`]
//!   refuses a state made for another one.
//! * **Representations.** The export is written once over the shared interface
//!   ([`export_simplicial`], generic over [`Exportable`]: [`Triangulation`] plus the
//!   refinement, the fingerprint and the way back of the representation), so the signposts and
//!   integer coordinates (Gillespie, Sharp, Crane 2021, crate `ddg-integer-coords`) split and
//!   flip the same edges in the same order wherever their lengths agree.

use ddg_core::{Geometry, TriMesh};

use crate::delaunay::DelaunayOptions;
use crate::edge_trace::EdgeTraceOptions;
use crate::error::IntrinsicError;
use crate::interface::Triangulation;
use crate::point::SurfacePoint;
use crate::query::{PathsOnInput, Segment};
use crate::refine::{RefineOptions, RefineReport};
use crate::triangulation::{IntrinsicTriangulation, irregularities, simplicial_mesh};

/// Options of [`IntrinsicTriangulation::export_simplicial`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ExportOptions {
    /// Whether, and with which options, the triangulation is flipped to Delaunay: before
    /// anything else if it is not Delaunay yet, and again after edges were split (a split can
    /// leave neighbouring edges that are not Delaunay, and a flip can make a new self-edge or
    /// double edge). The splits and the flips alternate until the triangulation is a
    /// simplicial complex that needs no flip. `None`: no flip at all, only the splits; the
    /// triangulation keeps the edges it has, Delaunay or not.
    pub delaunay: Option<DelaunayOptions>,
    /// If given, the triangulation is refined with these options after every round of
    /// splits, for at most [`max_refine_rounds`](Self::max_refine_rounds) rounds: splitting
    /// does not keep an angle bound, refinement restores it, and can in turn make new edges
    /// to split. Where the rounds end with a refinement that leaves a simplicial complex, the
    /// export has the angle bound and carries that refinement's report
    /// ([`SimplicialExport::refinement`]). The refinements stop, and the export ends with
    /// splits and flips as without this option, when two in a row leave no fewer edges to
    /// split than there were before them, or when one leaves more than twice as many (that
    /// one is undone). The rounds end the first way on thin tubes and handles, and the second
    /// way at a vertex whose angle sum is below three times the bound, where no simplicial
    /// complex has the bound (SPEC N15).
    pub refine: Option<RefineOptions>,
    /// The largest number of refinements between rounds of splits.
    pub max_refine_rounds: usize,
    /// The export fails when more than this many edges would have to be split.
    pub max_splits: usize,
}

impl Default for ExportOptions {
    /// Flips with the default Delaunay options; no refinement between the rounds of splits
    /// (at most 6 rounds if one is asked for); at most 100 000 splits.
    fn default() -> Self {
        ExportOptions {
            delaunay: Some(DelaunayOptions::default()),
            refine: None,
            max_refine_rounds: 6,
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
    /// Number of flips made by the export, those of its refinements included.
    pub n_flips: usize,
    /// Number of refinements made between rounds of splits ([`ExportOptions::refine`]); one
    /// that was undone is not counted.
    pub n_refine_rounds: usize,
    /// Number of vertices those refinements inserted.
    pub n_refined: usize,
    /// The report of the last refinement, if the export ended with it: the triangulation was
    /// a simplicial complex after it, so what it says about the angles holds for the export.
    /// `None` if no refinement was asked for, or if the export ended with splits after the
    /// last one (then [`min_angle`](Self::min_angle) says what the splits left).
    pub refinement: Option<RefineReport>,
    /// Whether every interior edge of the result passes the Delaunay test (the test of
    /// [`ExportOptions::delaunay`], or the default test if no flips were asked for). Tested on
    /// the result, not assumed.
    pub is_delaunay: bool,
    /// The smallest corner angle of the exported mesh. **Splitting does not keep an angle
    /// bound:** where edges were split, this can be far below the bound a refinement reached
    /// before (SPEC N15); where nothing was split, the angles are those of the triangulation
    /// as it was.
    pub min_angle: f64,
    /// The largest corner angle of the exported mesh.
    pub max_angle: f64,
    /// What identifies this export (see the [module documentation](self)).
    pub fingerprint: u64,
    /// Splits made off the midpoint of their edge because the midpoint was refused
    /// ([`Exportable::split_for_export`]): integer coordinates refuse a point within rounding
    /// of a crossing with an input edge (Integer Coordinates SPEC I28) and split near it
    /// instead. Always `0` on signposts.
    pub n_off_midpoint_splits: usize,
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
        for (f, face) in self.complex.faces.iter().enumerate() {
            for (k, &vertex) in face.iter().enumerate() {
                hasher.word(u64::from(vertex));
                hasher.word(
                    self.edge_lengths[self.complex.halfedge_edge[3 * f + k] as usize].to_bits(),
                );
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
    /// `TriMesh` with a `Geometry` (SPEC §6): [`export_simplicial`] on this representation.
    ///
    /// Every edge from a vertex to itself is split at its midpoint, and of every group of
    /// edges between the same two vertices all but one. This is repeated (a split can make a
    /// new group) until nothing is left. If `options.delaunay` is set, the triangulation is
    /// then flipped to Delaunay (also when nothing had to be split: with these options the
    /// export of a triangulation that is not Delaunay is its Delaunay triangulation), which can
    /// make new edges to split; splits and flips alternate until both hold.
    ///
    /// The triangulation is changed by this; the export is its state afterwards, and the way
    /// back ([`segments_on_input`](Self::segments_on_input)) must be asked of the
    /// triangulation in that state. If the export fails, the splits and flips made up to then
    /// stay (the error gives the number of splits); the triangulation is still a valid one of
    /// the same surface.
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
        export_simplicial(self, options)
    }
}

/// A representation of an intrinsic triangulation that can be exported as a simplicial mesh
/// ([`export_simplicial`]) and that can carry straight segments of the export back to the
/// input mesh (SPEC §6, N15; Integer Coordinates SPEC I29). Implemented by the signposts
/// ([`IntrinsicTriangulation`]) and by integer coordinates (`ddg-integer-coords`).
///
/// The export itself needs only the shared interface ([`Triangulation`]: splits, Delaunay
/// flips, the connectivity, the lengths, [`Triangulation::vertex_location`]); what differs
/// between representations is the refinement between rounds of splits, what identifies a
/// state, and how a segment is found on the input mesh.
pub trait Exportable: Triangulation + Clone {
    /// Splits edge `e` to make the triangulation simplicial; returns whether the split is off
    /// the midpoint. Provided: [`Triangulation::split_edge`] at `t = ½` (the signposts). A
    /// representation that refuses the midpoint for a reason of its own (integer coordinates:
    /// within rounding of a crossing, SPEC I28) may split at another parameter instead and say
    /// so.
    ///
    /// # Errors
    /// As [`Triangulation::split_edge`]; a refusal ([`IntrinsicError::InsertionRefused`]) is
    /// not fatal for the export (another edge may be split first).
    fn split_for_export(&mut self, e: usize) -> Result<bool, IntrinsicError> {
        self.split_edge(e, 0.5).map(|_| false)
    }

    /// The refinement [`ExportOptions::refine`] asks for. The signposts:
    /// [`IntrinsicTriangulation::refine`]. Every representation of this workspace refines by
    /// the same rules ([`RefineRules::Ruppert`](crate::RefineRules::Ruppert), IT-SPEC N1, N4),
    /// so that the export is the same function of the triangulation in each.
    ///
    /// # Errors
    /// The refinement's.
    fn refine_for_export(
        &mut self,
        options: &RefineOptions,
    ) -> Result<RefineReport, IntrinsicError>;

    /// What identifies the state of the triangulation ([`SimplicialExport::fingerprint`]): a
    /// hash of everything that determines the export and the way back. The signposts:
    /// [`IntrinsicTriangulation::fingerprint`] (faces, lengths per halfedge, places of the
    /// vertices).
    fn state_fingerprint(&self) -> u64;

    /// Straight segments of the triangulation as paths over the input mesh: the contract of
    /// [`IntrinsicTriangulation::segments_on_input`] (the same nodes, segments, output and
    /// errors).
    ///
    /// # Errors
    /// As [`IntrinsicTriangulation::segments_on_input`].
    fn segments_on_input(
        &self,
        nodes: &[SurfacePoint],
        segments: &[Segment],
        options: &EdgeTraceOptions,
    ) -> Result<PathsOnInput, IntrinsicError>;
}

impl Exportable for IntrinsicTriangulation {
    fn refine_for_export(
        &mut self,
        options: &RefineOptions,
    ) -> Result<RefineReport, IntrinsicError> {
        self.refine(options)
    }

    fn state_fingerprint(&self) -> u64 {
        self.fingerprint()
    }

    fn segments_on_input(
        &self,
        nodes: &[SurfacePoint],
        segments: &[Segment],
        options: &EdgeTraceOptions,
    ) -> Result<PathsOnInput, IntrinsicError> {
        IntrinsicTriangulation::segments_on_input(self, nodes, segments, options)
    }
}

/// Makes any representation of an intrinsic triangulation a simplicial complex by splitting
/// edges, and hands it out as a `TriMesh` with a `Geometry` (SPEC §6, N15). What it does is
/// described at [`IntrinsicTriangulation::export_simplicial`], which is this function; it uses
/// only the shared interface ([`Triangulation`]) and the three operations of [`Exportable`],
/// in the same order on every representation. The places of the vertices are
/// [`Triangulation::vertex_location`] (for integer coordinates the stored `q⁰`), the
/// fingerprint is [`Exportable::state_fingerprint`].
///
/// # Errors
/// As [`IntrinsicTriangulation::export_simplicial`]; the representation's errors of its
/// splits, flips and refinement.
pub fn export_simplicial<T: Exportable>(
    tri: &mut T,
    options: &ExportOptions,
) -> Result<SimplicialExport, IntrinsicError> {
    let mut n_splits = 0;
    let mut n_off_midpoint_splits = 0;
    let mut n_flips = 0;
    let mut n_refine_rounds = 0;
    let mut n_refined = 0;
    let mut refinement = None;
    let mut refining = options.refine.is_some();
    let mut stalled = 0;
    let offending = |tri: &T| {
        let irregular = irregularities(tri.complex());
        irregular.self_edges.len() + irregular.multi_edges.len()
    };
    loop {
        let offending_before = offending(tri);
        // Split until simplicial.
        loop {
            let irregular = irregularities(tri.complex());
            let mut to_split: Vec<u32> = irregular.self_edges.clone();
            // Of each group of parallel edges, all but the first.
            let mut pairs: Vec<([usize; 2], u32)> = irregular
                .multi_edges
                .iter()
                .map(|&e| {
                    let [a, b] = tri.complex().edge_vertices(e as usize);
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
                        reason: format!("more than {} splits would be needed", options.max_splits),
                    });
                }
                match tri.split_for_export(e as usize) {
                    Ok(off_midpoint) => {
                        n_splits += 1;
                        n_off_midpoint_splits += usize::from(off_midpoint);
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
        if let Some(refine) = &options.refine
            && refining
            && n_refine_rounds < options.max_refine_rounds
        {
            let before = tri.clone();
            let report = tri.refine_for_export(refine)?;
            let offending_after = offending(tri);
            if offending_after > 2 * offending_before + 16 {
                // Around a spike a refinement brings back the edges that were split and
                // adds vertices that need more: it is undone, and the export ends with
                // splits and flips.
                *tri = before;
                refining = false;
            } else {
                n_refine_rounds += 1;
                n_flips += report.n_flips;
                n_refined += report.n_inserted;
                if offending_after == 0 {
                    refinement = Some(report);
                    break;
                }
                // Two rounds in a row that leave no fewer edges to split: no more.
                stalled = if offending_after >= offending_before {
                    stalled + 1
                } else {
                    0
                };
                refining = stalled < 2;
                continue;
            }
        }
        let Some(delaunay) = &options.delaunay else {
            break;
        };
        let report = tri.flip_to_delaunay(delaunay)?;
        n_flips += report.n_flips;
        if report.n_flips == 0 || irregularities(tri.complex()).is_empty() {
            break;
        }
    }
    let criterion = match (&refinement, &options.refine, &options.delaunay) {
        (Some(_), Some(refine), _) => refine.criterion,
        (_, _, Some(delaunay)) => delaunay.criterion,
        _ => Default::default(),
    };
    let is_delaunay = tri.is_delaunay(&criterion);
    let complex = tri.complex();
    let (mesh, geometry) =
        simplicial_mesh(complex, tri.edge_lengths()).map_err(|e| IntrinsicError::ExportFailed {
            n_splits,
            reason: e.to_string(),
        })?;
    let mut edge = vec![0u32; mesh.n_edges()];
    let mut forward = vec![true; mesh.n_edges()];
    for e in 0..complex.n_edges() {
        let [tail, head] = complex.edge_vertices(e);
        // The faces are the same, so the halfedge ids are.
        let mesh_edge = mesh.edge(complex.edge_halfedge(e));
        edge[mesh_edge] = e as u32;
        forward[mesh_edge] = tail < head;
    }
    let (mut min_angle, mut max_angle) = (f64::INFINITY, 0.0f64);
    for angles in geometry.corner_angles() {
        for &angle in angles {
            min_angle = min_angle.min(angle);
            max_angle = max_angle.max(angle);
        }
    }
    Ok(SimplicialExport {
        mesh,
        geometry,
        locations: (0..complex.n_vertices())
            .map(|v| tri.vertex_location(v))
            .collect(),
        n_splits,
        n_flips,
        n_refine_rounds,
        n_refined,
        refinement,
        is_delaunay,
        min_angle,
        max_angle,
        fingerprint: tri.state_fingerprint(),
        n_off_midpoint_splits,
        edge,
        forward,
    })
}
