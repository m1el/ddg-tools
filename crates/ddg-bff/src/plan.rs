//! Choosing where to cut.
//!
//! [`plan_cut`] selects interior edges such that cutting along them (see [`crate::cut`]) turns
//! every connected piece of the surface into a topological disk with every cone on its
//! boundary. It handles closed surfaces, any number of boundary loops, any genus, edges that
//! must not be cut and edges that must be cut.
//!
//! # Algorithm
//!
//! 1. **Required cuts first.** The mesh is cut along the required edges. They may split it into
//!    several pieces; each piece is treated on its own and becomes one chart. The following
//!    steps never split a piece further.
//! 2. **Terminals.** In each piece, every boundary loop (original boundary or required cut) and
//!    every cone that is still an interior vertex is a terminal. A closed piece of positive
//!    genus without cones gets its lowest-numbered vertex with a cuttable edge as the only
//!    source of step 3.
//! 3. **Shortest-path forest.** One Dijkstra run over the cuttable edges, started from all
//!    terminals at once (a boundary loop counts as a single terminal), assigns every vertex to
//!    its nearest terminal. Neighbouring regions are then joined by their cheapest connections
//!    in order of increasing length `d(a) + ℓ_ab + d(b)` (Kruskal over the terminals). This is
//!    Mehlhorn's approximation of the minimum Steiner tree; it gives a forest `T` in the graph
//!    in which every boundary loop is contracted to a point.
//! 4. **Cotree.** The faces, grouped across the edges that must not be cut, are joined along
//!    the cuttable edges outside `T` by a maximum spanning tree `T*` with the same edge weights.
//!    Cutting every cuttable edge that is not crossed by `T*` leaves the faces glued in the
//!    pattern of a tree, which is a disk (tree–cotree decomposition, Eppstein 2003). The
//!    cuttable edges that are in neither tree — `2g` of them per piece of genus `g`, when `T`
//!    spans the piece — close the handle loops; taking the cotree of maximum weight selects
//!    short loops (Erickson & Whittlesey 2005).
//! 5. **Pruning.** A cut edge that ends at an interior vertex with no other cut edge is glued
//!    back unless that vertex is a cone (zipping a slit keeps a disk a disk). What remains is
//!    the forest joining the terminals plus the handle loops.
//!
//!    A cut that is a single edge between two interior vertices (two adjacent cones on a
//!    sphere, or a required seam of one edge) cannot be represented by an indexed mesh: its
//!    two sides would have the same end points. Such a slit is extended by a cuttable
//!    neighbouring edge, the shortest one that ends at an untouched interior vertex if there
//!    is one ([`CutPlan::n_extension_edges`]).
//! 6. **Verification.** The mesh is actually cut and every chart is checked to be a disk with
//!    all cones on its boundary. Failures are reported with their cause.
//!
//! # Guarantees
//!
//! If [`plan_cut`] returns `Ok`, then cutting along [`CutPlan::cut_edges`] gives a mesh in
//! which every connected component has Euler characteristic 1 and one boundary loop, every cone
//! vertex lies on the boundary, no forbidden edge is cut and every required interior edge is
//! cut. This is verified on the result, not assumed. The result is deterministic: ties are
//! broken by vertex and edge ids.
//!
//! The cut is short but not shortest: finding the shortest cut is NP-hard (Erickson &
//! Har-Peled 2004). The forest of step 3 is within a factor 2 of the shortest forest through
//! the terminals.
//!
//! # References
//!
//! * K. Mehlhorn. *A faster approximation algorithm for the Steiner problem in graphs.*
//!   Information Processing Letters 27 (1988).
//! * D. Eppstein. *Dynamic generators of topologically embedded graphs.* SODA 2003.
//! * J. Erickson, K. Whittlesey. *Greedy optimal homotopy and homology generators.* SODA 2005.
//! * J. Erickson, S. Har-Peled. *Optimally cutting a surface into a disk.* Discrete &
//!   Computational Geometry 31 (2004).
//!
//! The reference implementation of Boundary First Flattening cuts with a similar shortest-path
//! forest (`Cutter.cpp`), computes handle loops from breadth-first trees (`Generators.cpp`) and
//! fills all boundary loops but one instead of cutting between them (`HoleFiller.cpp`).

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use ddg_core::topology::{NONE, TriMesh};

use crate::cut::{
    ComponentTopology, CutMesh, component_topology, cut_mesh_flags, cut_mesh_lenient,
};
use crate::error::{BffError, check_values};

/// What to cut and under which constraints.
#[derive(Debug, Clone, Copy)]
pub struct CutPlanInput<'a> {
    /// The surface.
    pub mesh: &'a TriMesh,
    /// One non-negative length per edge of `mesh`: the cost of cutting that edge.
    pub edge_lengths: &'a [f64],
    /// The vertices that must end up on the boundary of their chart. Repetitions are allowed.
    /// Vertices that already lie on the boundary need no cut.
    pub cone_vertices: &'a [u32],
    /// Per edge of `mesh`: whether the edge may be cut. `None` allows every edge.
    pub may_cut: Option<&'a [bool]>,
    /// Edges that must be cut (edge ids of `mesh`). They can split the surface into several
    /// charts. Boundary edges in this list are ignored.
    pub must_cut: &'a [u32],
}

/// A planned cut.
#[derive(Debug, Clone, PartialEq)]
pub struct CutPlan {
    /// All interior edges to cut, ascending: the required ones and the chosen ones.
    pub cut_edges: Vec<u32>,
    /// Number of required interior edges among them.
    pub n_required: usize,
    /// Total length of the required interior edges.
    pub required_length: f64,
    /// Number of edges the planner added.
    pub n_chosen: usize,
    /// Total length of the edges the planner added.
    pub chosen_length: f64,
    /// Number of chosen edges that close handle loops (edges in neither tree, step 4 of the
    /// algorithm): `2g` summed over the pieces.
    pub n_loop_edges: usize,
    /// Number of chosen edges that were added only to make a one-edge slit representable
    /// (step 5 of the algorithm).
    pub n_extension_edges: usize,
    /// One entry per chart, in chart order (charts are numbered by their lowest face).
    pub charts: Vec<PlannedChart>,
}

/// One chart of a planned cut: a connected piece of the surface after the required cuts.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlannedChart {
    /// Topology of the piece before the planner added cuts.
    pub before: ComponentTopology,
    /// Topology of the chart after all cuts (a disk).
    pub after: ComponentTopology,
    /// Number of cones that were interior vertices of the piece and had to be reached by a cut.
    pub n_interior_cones: usize,
}

/// Plans a cut. See the [module documentation](self) for the algorithm and its guarantees.
///
/// # Errors
/// * [`BffError::LengthMismatch`], [`BffError::NonFiniteInput`], [`BffError::EdgeOutOfRange`],
///   [`BffError::ConeVertexOutOfRange`] for malformed input;
/// * [`BffError::RequiredEdgeForbidden`] if an edge is both required and not cuttable;
/// * [`BffError::SphereNeedsTwoCones`] for a closed piece of genus 0 with fewer than two cones;
/// * [`BffError::ConeEnclosedByUncuttableEdges`] if a cone cannot be reached;
/// * [`BffError::IsolatedCutEdge`] if a one-edge slit cannot be extended;
/// * [`BffError::UncuttableRegionNotDisk`] if the edges that must not be cut make a disk
///   impossible.
pub fn plan_cut(input: &CutPlanInput<'_>) -> Result<CutPlan, BffError> {
    plan_and_cut(input).map(|(plan, _)| plan)
}

/// [`plan_cut`], also returning the mesh cut along the plan (which verification builds anyway).
pub(crate) fn plan_and_cut(input: &CutPlanInput<'_>) -> Result<(CutPlan, CutMesh), BffError> {
    let mesh = input.mesh;
    let n_edges = mesh.n_edges();
    check_values("edge lengths", input.edge_lengths, n_edges)?;
    if let Some(index) = input.edge_lengths.iter().position(|&l| l < 0.0) {
        return Err(BffError::NonFiniteInput {
            what: "edge lengths (must be non-negative)",
            index,
            value: input.edge_lengths[index],
        });
    }
    if let Some(mask) = input.may_cut
        && mask.len() != n_edges
    {
        return Err(BffError::LengthMismatch {
            what: "may_cut mask",
            expected: n_edges,
            got: mask.len(),
        });
    }
    let is_cone = cone_flags(mesh, input.cone_vertices)?;

    // Step 1: the required cuts.
    let mut required = vec![false; n_edges];
    for (position, &e) in input.must_cut.iter().enumerate() {
        if e as usize >= n_edges {
            return Err(BffError::EdgeOutOfRange {
                what: "must_cut",
                position,
                edge: e,
                n_edges,
            });
        }
        if let Some(mask) = input.may_cut
            && !mask[e as usize]
        {
            return Err(BffError::RequiredEdgeForbidden {
                edge: e,
                vertices: mesh.edges()[e as usize],
            });
        }
        if !mesh.is_boundary_edge(e as usize) {
            required[e as usize] = true;
        }
    }
    let n_required = required.iter().filter(|&&r| r).count();
    let required_length: f64 = (0..n_edges)
        .filter(|&e| required[e])
        .map(|e| input.edge_lengths[e])
        .sum();
    // A required seam of a single edge has to be extended before it can be cut.
    let cuttable: Vec<bool> = (0..n_edges)
        .map(|e| !mesh.is_boundary_edge(e) && input.may_cut.is_none_or(|mask| mask[e]))
        .collect();
    let mut extension = vec![false; n_edges];
    let mut n_extension_edges = 0usize;
    if n_required > 0 {
        let before = required.clone();
        n_extension_edges =
            extend_isolated_slits(mesh, &mut required, &cuttable, input.edge_lengths).map_err(
                |e| BffError::IsolatedCutEdge {
                    edge: e,
                    vertices: mesh.edges()[e as usize],
                },
            )?;
        for e in 0..n_edges {
            extension[e] = required[e] && !before[e];
        }
    }
    let has_required = n_required > 0;

    let pieces_owned;
    let (pieces, piece_edge_to_original): (&TriMesh, Option<&[u32]>) = if !has_required {
        (mesh, None)
    } else {
        pieces_owned = cut_mesh_flags(mesh, &required, 0)?;
        (&pieces_owned.mesh, Some(&pieces_owned.edge_to_original))
    };
    let original_edge = |e: usize| piece_edge_to_original.map_or(e, |map| map[e] as usize);

    // Edge data of the mesh of pieces. An edge is usable for cutting if it is interior to its
    // piece and not forbidden.
    let lengths: Vec<f64> = (0..pieces.n_edges())
        .map(|e| input.edge_lengths[original_edge(e)])
        .collect();
    let allowed: Vec<bool> = (0..pieces.n_edges())
        .map(|e| !pieces.is_boundary_edge(e) && cuttable[original_edge(e)])
        .collect();
    // A cone needs a cut only if it is still an interior vertex. An interior vertex of the
    // mesh of pieces has a single copy, which keeps the original id.
    let terminal: Vec<bool> = (0..pieces.n_vertices())
        .map(|v| v < mesh.n_vertices() && is_cone[v] && !pieces.is_boundary_vertex(v))
        .collect();

    let before = component_topology(pieces);
    let mut chosen = plan_pieces(pieces, &before, &lengths, &allowed, &terminal)?;
    n_extension_edges += extend_isolated_slits(pieces, &mut chosen.is_cut, &allowed, &lengths)
        .map_err(|e| BffError::IsolatedCutEdge {
            edge: original_edge(e as usize) as u32,
            vertices: mesh.edges()[original_edge(e as usize)],
        })?;

    let mut is_cut = required;
    let mut n_chosen = 0usize;
    let mut chosen_length = 0.0;
    for e in (0..n_edges).filter(|&e| extension[e]) {
        n_chosen += 1;
        chosen_length += input.edge_lengths[e];
    }
    for (e, _) in chosen.is_cut.iter().enumerate().filter(|&(_, &c)| c) {
        is_cut[original_edge(e)] = true;
        n_chosen += 1;
        chosen_length += lengths[e];
    }

    // Step 6: cut and verify.
    let cut = cut_mesh_flags(mesh, &is_cut, 0)?;
    let after = match check_disk_charts(&cut, &is_cone) {
        Ok(after) => after,
        Err(error) => return Err(diagnose(mesh, &is_cut, input.may_cut, error)),
    };
    debug_assert_eq!(after.len(), before.len());
    let piece_of_vertex = pieces.vertex_components();
    let mut n_interior_cones = vec![0usize; before.len()];
    for (v, _) in terminal.iter().enumerate().filter(|&(_, &t)| t) {
        n_interior_cones[piece_of_vertex[v] as usize] += 1;
    }
    let charts = before
        .iter()
        .zip(&after)
        .zip(&n_interior_cones)
        .map(|((&before, &after), &n_interior_cones)| PlannedChart {
            before,
            after,
            n_interior_cones,
        })
        .collect();
    let plan = CutPlan {
        cut_edges: cut.cut_edges.clone(),
        n_required,
        required_length,
        n_chosen,
        chosen_length,
        n_loop_edges: chosen.n_loop_edges,
        n_extension_edges,
        charts,
    };
    Ok((plan, cut))
}

/// Per vertex: whether it is listed as a cone.
pub(crate) fn cone_flags(mesh: &TriMesh, cone_vertices: &[u32]) -> Result<Vec<bool>, BffError> {
    let mut is_cone = vec![false; mesh.n_vertices()];
    for (position, &v) in cone_vertices.iter().enumerate() {
        if v as usize >= mesh.n_vertices() {
            return Err(BffError::ConeVertexOutOfRange {
                position,
                vertex: v,
                n_vertices: mesh.n_vertices(),
            });
        }
        is_cone[v as usize] = true;
    }
    Ok(is_cone)
}

/// Checks that every component of a cut mesh is a disk and that every flagged original vertex
/// lies on the boundary of the cut mesh. Returns the topology of the charts.
///
/// [formula: genus-from-euler]
pub(crate) fn check_disk_charts(
    cut: &CutMesh,
    is_cone: &[bool],
) -> Result<Vec<ComponentTopology>, BffError> {
    let topology = component_topology(&cut.mesh);
    for (chart, t) in topology.iter().enumerate() {
        if !t.is_disk() {
            return Err(BffError::ChartNotDisk {
                chart: chart as u32,
                first_face: t.first_face,
                n_faces: t.n_faces,
                euler_characteristic: t.euler_characteristic,
                boundary_loops: t.boundary_loops,
                genus: t.genus,
            });
        }
    }
    for (v, &original) in cut.vertex_to_original.iter().enumerate() {
        if is_cone[original as usize] && !cut.mesh.is_boundary_vertex(v) {
            let chart = cut.mesh.face_components()[cut.mesh.outgoing_halfedge(v) / 3];
            return Err(BffError::ConeNotOnCut {
                vertex: original,
                chart,
            });
        }
    }
    Ok(topology)
}

/// Turns a failed verification into the error that names its cause, if the cause is the
/// `may_cut` mask: with every cuttable edge cut, some region is still not a disk.
fn diagnose(
    mesh: &TriMesh,
    is_cut: &[bool],
    may_cut: Option<&[bool]>,
    error: BffError,
) -> BffError {
    let Some(mask) = may_cut else {
        return error;
    };
    let everything: Vec<bool> = (0..mesh.n_edges())
        .map(|e| !mesh.is_boundary_edge(e) && (mask[e] || is_cut[e]))
        .collect();
    let mut slits = Vec::new();
    let Ok(cut) = cut_mesh_lenient(mesh, &everything, &mut slits) else {
        return error;
    };
    let mut topology = component_topology(&cut.mesh);
    // A slit that could not be represented is one more edge and one more boundary loop.
    for &e in &slits {
        let face = mesh.edge_halfedge(e as usize) / 3;
        let t = &mut topology[cut.mesh.face_components()[face] as usize];
        t.n_edges += 1;
        t.euler_characteristic -= 1;
        t.boundary_loops += 1;
    }
    match topology.iter().find(|t| !t.is_disk()) {
        Some(t) => BffError::UncuttableRegionNotDisk {
            first_face: t.first_face,
            n_faces: t.n_faces,
            euler_characteristic: t.euler_characteristic,
            boundary_loops: t.boundary_loops,
        },
        None => error,
    }
}

/// Extends every isolated slit of a cut by one edge: a cut edge whose two end points are
/// interior vertices without another cut edge gets a usable neighbouring edge, preferably the
/// shortest one that ends at an interior vertex without cut edges. Adding such an edge keeps
/// every disk a disk (it is the reverse of pruning). If there is none, the shortest usable
/// neighbouring edge of any kind is taken; it ends on the boundary or on another part of the
/// cut, which joins two boundary loops and cannot disconnect anything. Returns the number of
/// added edges, or the slit that has no usable neighbour at all.
///
/// [formula: cut-leaf-pruning]
fn extend_isolated_slits(
    mesh: &TriMesh,
    is_cut: &mut [bool],
    usable: &[bool],
    lengths: &[f64],
) -> Result<usize, u32> {
    let mut degree = vec![0u32; mesh.n_vertices()];
    for e in (0..mesh.n_edges()).filter(|&e| is_cut[e]) {
        let [a, b] = mesh.edge_vertices(e);
        degree[a] += 1;
        degree[b] += 1;
    }
    let mut added = 0usize;
    for e in 0..mesh.n_edges() {
        if !is_cut[e] {
            continue;
        }
        let [a, b] = mesh.edge_vertices(e);
        let isolated = [a, b]
            .iter()
            .all(|&v| degree[v] == 1 && !mesh.is_boundary_vertex(v));
        if !isolated {
            continue;
        }
        let mut best: Option<(bool, Dist, usize, usize, usize)> = None;
        for x in [a, b] {
            for h in mesh.outgoing_halfedges(x) {
                let (candidate, c) = (mesh.edge(h), mesh.to_vertex(h));
                if usable[candidate] && !is_cut[candidate] {
                    let touched = degree[c] != 0 || mesh.is_boundary_vertex(c);
                    let key = (touched, Dist(lengths[candidate]), candidate, x, c);
                    if best.is_none_or(|b| (key.0, key.1, key.2) < (b.0, b.1, b.2)) {
                        best = Some(key);
                    }
                }
            }
        }
        let Some((_, _, candidate, x, c)) = best else {
            return Err(e as u32);
        };
        is_cut[candidate] = true;
        degree[x] += 1;
        degree[c] += 1;
        added += 1;
    }
    Ok(added)
}

/// A distance with a total order, for the priority queue.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Dist(f64);

impl Eq for Dist {}

impl PartialOrd for Dist {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Dist {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.total_cmp(&other.0)
    }
}

/// Union-find with path halving and union by size.
struct DisjointSets {
    parent: Vec<u32>,
    size: Vec<u32>,
}

impl DisjointSets {
    fn new(n: usize) -> Self {
        DisjointSets {
            parent: (0..n as u32).collect(),
            size: vec![1; n],
        }
    }

    fn find(&mut self, mut x: u32) -> u32 {
        while self.parent[x as usize] != x {
            let grandparent = self.parent[self.parent[x as usize] as usize];
            self.parent[x as usize] = grandparent;
            x = grandparent;
        }
        x
    }

    /// Joins the sets of `a` and `b`; returns whether they were different.
    fn union(&mut self, a: u32, b: u32) -> bool {
        let (mut a, mut b) = (self.find(a), self.find(b));
        if a == b {
            return false;
        }
        if self.size[a as usize] < self.size[b as usize] {
            std::mem::swap(&mut a, &mut b);
        }
        self.parent[b as usize] = a;
        self.size[a as usize] += self.size[b as usize];
        true
    }
}

struct Chosen {
    /// Per edge of the mesh of pieces: whether the planner cuts it.
    is_cut: Vec<bool>,
    n_loop_edges: usize,
}

/// Steps 2 to 5 of the algorithm, for all pieces at once.
///
/// The disk property rests on two facts. (i) If `T` is a set of interior edges without cycles
/// in the graph obtained by contracting every boundary loop to a vertex, then the dual graph
/// without the edges of `T` is connected, so a spanning tree `T*` of the faces avoiding `T`
/// exists; cutting all interior edges not crossed by `T*` glues the faces along a tree, which
/// gives a disk. (ii) Un-cutting an edge at an interior vertex that has no other cut edge
/// glues two consecutive boundary edges of a disk and gives a disk again, unless these two
/// edges are its whole boundary.
///
/// [formula: tree-cotree-disk]
/// [formula: cut-leaf-pruning]
fn plan_pieces(
    mesh: &TriMesh,
    topology: &[ComponentTopology],
    lengths: &[f64],
    allowed: &[bool],
    terminal: &[bool],
) -> Result<Chosen, BffError> {
    let n = mesh.n_vertices();
    let n_edges = mesh.n_edges();
    let piece_of_vertex = mesh.vertex_components();

    // Every terminal cone needs a cuttable edge.
    for v in (0..n).filter(|&v| terminal[v]) {
        if !mesh.outgoing_halfedges(v).any(|h| allowed[mesh.edge(h)]) {
            return Err(BffError::ConeEnclosedByUncuttableEdges { vertex: v as u32 });
        }
    }
    let mut cones_in_piece = vec![0usize; topology.len()];
    for v in (0..n).filter(|&v| terminal[v]) {
        cones_in_piece[piece_of_vertex[v] as usize] += 1;
    }
    for (piece, t) in topology.iter().enumerate() {
        if t.is_closed() && t.genus == 0 && cones_in_piece[piece] < 2 {
            return Err(BffError::SphereNeedsTwoCones {
                component: piece as u32,
                first_face: t.first_face,
                n_cones: cones_in_piece[piece],
            });
        }
    }

    // Step 2: sources. Set ids: boundary loops, then cones, then roots of closed pieces
    // without cones.
    let mut region = vec![NONE; n];
    let mut n_sets = 0u32;
    for boundary in mesh.boundary_loops() {
        for &v in boundary {
            region[v as usize] = n_sets;
        }
        n_sets += 1;
    }
    for v in (0..n).filter(|&v| terminal[v]) {
        region[v] = n_sets;
        n_sets += 1;
    }
    let mut needs_root: Vec<bool> = topology
        .iter()
        .zip(&cones_in_piece)
        .map(|(t, &cones)| t.is_closed() && cones == 0)
        .collect();
    for v in 0..n {
        let piece = piece_of_vertex[v] as usize;
        if needs_root[piece] && mesh.outgoing_halfedges(v).any(|h| allowed[mesh.edge(h)]) {
            needs_root[piece] = false;
            region[v] = n_sets;
            n_sets += 1;
        }
    }

    // Step 3a: multi-source Dijkstra over the cuttable edges.
    let mut dist = vec![f64::INFINITY; n];
    let mut parent_edge = vec![NONE; n];
    let mut heap: BinaryHeap<std::cmp::Reverse<(Dist, u32)>> = BinaryHeap::new();
    for v in (0..n).filter(|&v| region[v] != NONE) {
        dist[v] = 0.0;
        heap.push(std::cmp::Reverse((Dist(0.0), v as u32)));
    }
    while let Some(std::cmp::Reverse((Dist(d), v))) = heap.pop() {
        let v = v as usize;
        if d > dist[v] {
            continue;
        }
        for h in mesh.outgoing_halfedges(v) {
            let e = mesh.edge(h);
            if !allowed[e] {
                continue;
            }
            let w = mesh.to_vertex(h);
            let candidate = d + lengths[e];
            if candidate < dist[w] {
                dist[w] = candidate;
                region[w] = region[v];
                parent_edge[w] = e as u32;
                heap.push(std::cmp::Reverse((Dist(candidate), w as u32)));
            }
        }
    }
    let mut in_tree = vec![false; n_edges];
    for &e in parent_edge.iter().filter(|&&e| e != NONE) {
        in_tree[e as usize] = true;
    }

    // The weight of a cuttable edge: the length of the connection it closes.
    let weight = |e: usize| {
        let [a, b] = mesh.edge_vertices(e);
        dist[a] + lengths[e] + dist[b]
    };

    // Step 3b: join the regions (Kruskal over the terminals).
    let mut bridges: Vec<(Dist, u32)> = (0..n_edges)
        .filter(|&e| allowed[e])
        .filter(|&e| {
            let [a, b] = mesh.edge_vertices(e);
            region[a] != NONE && region[b] != NONE && region[a] != region[b]
        })
        .map(|e| (Dist(weight(e)), e as u32))
        .collect();
    bridges.sort_unstable();
    let mut sets = DisjointSets::new(n_sets as usize);
    for &(_, e) in &bridges {
        let [a, b] = mesh.edge_vertices(e as usize);
        if sets.union(region[a], region[b]) {
            in_tree[e as usize] = true;
        }
    }

    // Step 4: maximum spanning cotree over the groups of faces.
    let mut groups = DisjointSets::new(mesh.n_faces());
    let faces_of = |e: usize| {
        let h = mesh.edge_halfedge(e);
        let t = mesh.twin(h).expect("interior edge");
        ((h / 3) as u32, (t / 3) as u32)
    };
    for e in (0..n_edges).filter(|&e| !allowed[e] && !mesh.is_boundary_edge(e)) {
        let (f, g) = faces_of(e);
        groups.union(f, g);
    }
    let mut dual: Vec<(Dist, u32)> = (0..n_edges)
        .filter(|&e| allowed[e] && !in_tree[e])
        .map(|e| (Dist(weight(e)), e as u32))
        .collect();
    dual.sort_unstable_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
    let mut is_cut = in_tree;
    let mut loop_edge = vec![false; n_edges];
    for &(_, e) in &dual {
        let (f, g) = faces_of(e as usize);
        if !groups.union(f, g) {
            is_cut[e as usize] = true;
            loop_edge[e as usize] = true;
        }
    }

    // Step 5: prune.
    let mut degree = vec![0u32; n];
    for e in (0..n_edges).filter(|&e| is_cut[e]) {
        let [a, b] = mesh.edge_vertices(e);
        degree[a] += 1;
        degree[b] += 1;
    }
    let prunable = |v: usize| !terminal[v] && !mesh.is_boundary_vertex(v);
    let mut queue: std::collections::VecDeque<u32> = (0..n as u32)
        .filter(|&v| degree[v as usize] == 1 && prunable(v as usize))
        .collect();
    while let Some(v) = queue.pop_front() {
        let v = v as usize;
        if degree[v] != 1 {
            continue;
        }
        let Some(h) = mesh.outgoing_halfedges(v).find(|&h| is_cut[mesh.edge(h)]) else {
            continue;
        };
        let w = mesh.to_vertex(h);
        // Keep the last cut edge of a cone: the cone must stay on the boundary.
        if terminal[w] && degree[w] == 1 {
            continue;
        }
        is_cut[mesh.edge(h)] = false;
        degree[v] = 0;
        degree[w] -= 1;
        if degree[w] == 1 && prunable(w) {
            queue.push_back(w as u32);
        }
    }
    let n_loop_edges = (0..n_edges).filter(|&e| is_cut[e] && loop_edge[e]).count();
    Ok(Chosen {
        is_cut,
        n_loop_edges,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use ddg_core::generators;
    use ddg_core::geometry::Geometry;

    fn lengths(g: &generators::GeneratedMesh, mesh: &TriMesh) -> Vec<f64> {
        Geometry::from_positions(mesh, &g.positions)
            .unwrap()
            .edge_lengths()
            .to_vec()
    }

    #[test]
    fn disk_without_cones_needs_no_cut() {
        let g = generators::grid_disk(4, 3).unwrap();
        let mesh = g.topology().unwrap();
        let plan = plan_cut(&CutPlanInput {
            mesh: &mesh,
            edge_lengths: &lengths(&g, &mesh),
            cone_vertices: &[],
            may_cut: None,
            must_cut: &[],
        })
        .unwrap();
        assert!(plan.cut_edges.is_empty());
        assert_eq!(plan.charts.len(), 1);
    }

    #[test]
    fn sphere_with_two_cones_is_cut_along_a_path() {
        let g = generators::icosphere(2).unwrap();
        let mesh = g.topology().unwrap();
        let plan = plan_cut(&CutPlanInput {
            mesh: &mesh,
            edge_lengths: &lengths(&g, &mesh),
            cone_vertices: &[0, 3],
            may_cut: None,
            must_cut: &[],
        })
        .unwrap();
        assert!(!plan.cut_edges.is_empty());
        assert_eq!(plan.n_loop_edges, 0);
        // A path: both cones have one cut edge, every other cut vertex has two.
        let mut degree = vec![0; mesh.n_vertices()];
        for &e in &plan.cut_edges {
            let [a, b] = mesh.edge_vertices(e as usize);
            degree[a] += 1;
            degree[b] += 1;
        }
        assert_eq!(degree[0], 1);
        assert_eq!(degree[3], 1);
        assert!(degree.iter().all(|&d| d <= 2));
        assert!(plan.charts[0].after.is_disk());
    }

    #[test]
    fn sphere_needs_two_cones() {
        let g = generators::icosphere(1).unwrap();
        let mesh = g.topology().unwrap();
        for cones in [&[][..], &[5][..]] {
            let err = plan_cut(&CutPlanInput {
                mesh: &mesh,
                edge_lengths: &lengths(&g, &mesh),
                cone_vertices: cones,
                may_cut: None,
                must_cut: &[],
            })
            .unwrap_err();
            assert!(matches!(err, BffError::SphereNeedsTwoCones { .. }), "{err}");
        }
    }

    #[test]
    fn torus_without_cones_gets_two_loops() {
        let g = generators::torus(12, 8, 2.0, 0.7).unwrap();
        let mesh = g.topology().unwrap();
        let plan = plan_cut(&CutPlanInput {
            mesh: &mesh,
            edge_lengths: &lengths(&g, &mesh),
            cone_vertices: &[],
            may_cut: None,
            must_cut: &[],
        })
        .unwrap();
        assert_eq!(plan.n_loop_edges, 2);
        assert_eq!(plan.charts.len(), 1);
        assert_eq!(plan.charts[0].before.genus, 1);
        assert!(plan.charts[0].after.is_disk());
        // One loop around the tube (8 edges) and one around the axis (12 edges), which meet
        // in a vertex.
        assert_eq!(plan.cut_edges.len(), 20);
    }

    #[test]
    fn one_edge_slits_are_extended() {
        let g = generators::icosphere(1).unwrap();
        let mesh = g.topology().unwrap();
        let l = lengths(&g, &mesh);
        // Two adjacent cones on a sphere: the path between them is one edge.
        let [a, b] = mesh.edge_vertices(7);
        let plan = plan_cut(&CutPlanInput {
            mesh: &mesh,
            edge_lengths: &l,
            cone_vertices: &[a as u32, b as u32],
            may_cut: None,
            must_cut: &[],
        })
        .unwrap();
        assert_eq!(plan.cut_edges.len(), 2);
        assert!(plan.cut_edges.contains(&7));
        assert_eq!(plan.n_extension_edges, 1);
        assert!(plan.charts[0].after.is_disk());
        // A required seam of one edge in a disk.
        let g = generators::grid_disk(6, 6).unwrap();
        let mesh = g.topology().unwrap();
        let l = lengths(&g, &mesh);
        let seam = mesh.find_edge(3 * 7 + 3, 3 * 7 + 4).unwrap() as u32;
        let plan = plan_cut(&CutPlanInput {
            mesh: &mesh,
            edge_lengths: &l,
            cone_vertices: &[],
            may_cut: None,
            must_cut: &[seam],
        })
        .unwrap();
        assert!(plan.cut_edges.contains(&seam));
        assert_eq!(plan.n_required, 1);
        assert_eq!(plan.n_extension_edges, 1);
        // The slit is a second boundary loop, so the planner joins it to the outer boundary.
        assert!(plan.charts[0].after.is_disk());
        assert_eq!(plan.charts[0].before.boundary_loops, 2);
    }

    #[test]
    fn annulus_is_cut_between_its_boundary_loops() {
        let g = generators::annulus(16, 3, 1.0, 2.0).unwrap();
        let mesh = g.topology().unwrap();
        let plan = plan_cut(&CutPlanInput {
            mesh: &mesh,
            edge_lengths: &lengths(&g, &mesh),
            cone_vertices: &[],
            may_cut: None,
            must_cut: &[],
        })
        .unwrap();
        assert_eq!(plan.cut_edges.len(), 3);
        assert!(plan.charts[0].after.is_disk());
    }

    #[test]
    fn forbidden_and_required_edges_are_honoured() {
        let g = generators::icosphere(2).unwrap();
        let mesh = g.topology().unwrap();
        let l = lengths(&g, &mesh);
        let free = plan_cut(&CutPlanInput {
            mesh: &mesh,
            edge_lengths: &l,
            cone_vertices: &[0, 3, 7],
            may_cut: None,
            must_cut: &[],
        })
        .unwrap();
        // Forbid the edges of the unconstrained plan: the new plan must avoid them.
        let mut may = vec![true; mesh.n_edges()];
        for &e in &free.cut_edges {
            may[e as usize] = false;
        }
        // Require an edge far from everything.
        let required = (0..mesh.n_edges() as u32)
            .find(|&e| may[e as usize] && !free.cut_edges.contains(&e))
            .unwrap();
        let plan = plan_cut(&CutPlanInput {
            mesh: &mesh,
            edge_lengths: &l,
            cone_vertices: &[0, 3, 7],
            may_cut: Some(&may),
            must_cut: &[required],
        })
        .unwrap();
        assert!(plan.cut_edges.iter().all(|&e| may[e as usize]));
        assert!(plan.cut_edges.contains(&required));
        assert!(plan.charts.iter().all(|c| c.after.is_disk()));
        // Requiring a forbidden edge is a contradiction.
        let err = plan_cut(&CutPlanInput {
            mesh: &mesh,
            edge_lengths: &l,
            cone_vertices: &[0, 3, 7],
            may_cut: Some(&may),
            must_cut: &[free.cut_edges[0]],
        })
        .unwrap_err();
        assert!(
            matches!(err, BffError::RequiredEdgeForbidden { .. }),
            "{err}"
        );
    }

    #[test]
    fn enclosed_cone_and_non_disk_regions_are_reported() {
        let g = generators::icosphere(1).unwrap();
        let mesh = g.topology().unwrap();
        let l = lengths(&g, &mesh);
        // No edge at vertex 0 may be cut.
        let mut may = vec![true; mesh.n_edges()];
        for h in mesh.outgoing_halfedges(0) {
            may[mesh.edge(h)] = false;
        }
        let err = plan_cut(&CutPlanInput {
            mesh: &mesh,
            edge_lengths: &l,
            cone_vertices: &[0, 5],
            may_cut: Some(&may),
            must_cut: &[],
        })
        .unwrap_err();
        assert_eq!(err, BffError::ConeEnclosedByUncuttableEdges { vertex: 0 });
        // Nothing may be cut: the whole sphere is one uncuttable region.
        let may = vec![false; mesh.n_edges()];
        let err = plan_cut(&CutPlanInput {
            mesh: &mesh,
            edge_lengths: &l,
            cone_vertices: &[],
            may_cut: Some(&may),
            must_cut: &[],
        })
        .unwrap_err();
        assert!(matches!(err, BffError::SphereNeedsTwoCones { .. }), "{err}");
        // A torus on which only one edge may be cut cannot become a disk.
        let g = generators::torus(6, 5, 2.0, 0.7).unwrap();
        let mesh = g.topology().unwrap();
        let mut may = vec![false; mesh.n_edges()];
        may[0] = true;
        let err = plan_cut(&CutPlanInput {
            mesh: &mesh,
            edge_lengths: &lengths(&g, &mesh),
            cone_vertices: &[],
            may_cut: Some(&may),
            must_cut: &[],
        })
        .unwrap_err();
        assert!(
            matches!(err, BffError::UncuttableRegionNotDisk { .. }),
            "{err}"
        );
    }
}
