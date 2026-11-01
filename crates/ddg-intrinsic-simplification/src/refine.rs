//! Intrinsic retriangulation around the coarsening: pre-refinement and post-refinement (SPEC §6;
//! paper Sec. 6.3, Fig. 12 and 13).
//!
//! Both use the workspace's own Delaunay refinement (`ddg-integer-coords`, Gillespie, Sharp,
//! Crane 2021, by default; the signposts of `ddg-intrinsic`, Sharp, Soliman, Crane 2019, as an
//! option; [`Representation::DEFAULT`]), not the reference's, whose refinement is slow and has
//! bugs (reference map R35, R36, R44).
//!
//! # Pre-refinement (SPEC §6.1)
//!
//! [`CoarseMesh::from_intrinsic`] makes a coarse mesh of an intrinsic triangulation `T¹` of the
//! input (generic over [`LocatedTriangulation`], implemented for both representations): its
//! [`Complex`] and edge lengths, the input vertices as `Vertex` points (points `0..|V|`), and
//! the extra [`InputPoint`]s of the input surface located in `T¹` by the representation's
//! `T⁰ → T¹` point query. [`CoarseMesh::from_refined`] is the whole recipe: build `T⁰`, refine
//! it, call `from_intrinsic`.
//!
//! * The vertices that refinement inserted have fine ids `|V|, |V|+1, …` (their ids in `T¹`).
//!   They are ordinary removable vertices (they are not fixed), they are **not** tracked points
//!   (the rows of `P` are the input vertices and the extra points, as before:
//!   [`CoarseMesh::n_fine_vertices`] stays `|V|`), and they are columns of `P` if they stay.
//! * `CoarseMesh::vertices()` lists every current vertex; loop over it (not over
//!   `0..n_fine_vertices`) when removing vertices in a fixed order.
//! * The reference recipe pairs pre-refinement with an area weight of `0.1` (otherwise all the
//!   zero-curvature inserted vertices go first); the default weights are kept here and the
//!   recipe is not imposed.
//! * The one-call convenience "refine, then coarsen" is `from_refined` followed by `simplify`
//!   (WP3); it is wired once `simplify` is merged.
//!
//! # Post-refinement (SPEC §6.2)
//!
//! [`CoarseMesh::post_refine`] treats the coarse mesh (its [`Complex`] and lengths) as the input
//! `T⁰` of a **new** intrinsic triangulation, Delaunay-refines that to a minimum angle
//! (default 25°), and composes the maps: each tracked point `Face(f, b)` of the coarse mesh is a
//! point of `T⁰` and is located in the refined `T¹` by the `T⁰ → T¹` query; vertices of the
//! coarse mesh stay vertices (the input vertices of the refinement are never removed). The
//! result is again a [`CoarseMesh`]:
//!
//! * vertices of the old mesh keep their fine ids; inserted vertices get fresh ids above every
//!   id used so far (`complex_of.len()` of the old mesh, in the order of `T¹`'s vertices);
//! * the fixed mask is carried over, the new vertices are free;
//! * faces, edges and halfedges are the **new** complex's (the old face ids mean nothing any
//!   more; tracked points are located by their new faces);
//! * lengths are those of `T¹`; the angular coordinates are **rebuilt from scratch** (S4:
//!   reference halfedge = the outgoing halfedge of the new complex, `φ = 0` there, the rest
//!   from the corner angles). The new frames are therefore not the old ones: a vertex's
//!   reference direction in the old coarse mesh and in the new are different directions in
//!   general, and a coarse vector field expressed in the old frames has to be re-expressed
//!   (prolongation to the fine mesh uses the frames of the mesh it is computed from and the
//!   fine frames, never both coarse meshes, so `P` and `P^vec` of the new mesh are consistent
//!   by themselves);
//! * the tracked-point order is the old one: the same point ids, and `P` has the same rows.
//!
//! **What is refused.** The refinement crates take a *simplicial* manifold mesh
//! ([`TriMesh`]: distinct vertices in every face, no two faces on the same three vertices, at
//! most one edge between two vertices), while a coarse mesh is a Δ-complex. After heavy
//! coarsening it can have a **self-edge or a face that uses a vertex twice**, **two edges between
//! the same two vertices** (a bigon), or **two faces on the same three vertices**.
//! [`CoarseMesh::post_refine`] checks this first and returns [`RefineError::NotSimplicial`]
//! naming the faces, edges and fine ids; nothing is repaired, subdivided or changed. (The
//! measured frequency is in the README.) Other refusals: a pending removal
//! ([`IemError::RemovalPending`]), and the errors of the refinement and of the point queries,
//! which name the point.
//!
//! **The guarantee.** The refinement reaches the angle bound on every corner except those at
//! vertices whose angle sum is below it, and by Gillespie et al. 2021a (Thm. C.1 of their
//! supplement) on a closed surface every angle sum of the input `≥ 60°` makes `30°` attainable
//! (so the default 25° as well). [`RefinementReport::guarantee_applies`] says whether the input
//! of the refinement was closed with every angle sum `≥ 60°`;
//! [`RefinementReport::min_angle`] is the smallest corner angle of the result, whether or not
//! the guarantee applied.

use std::fmt;
use std::time::Instant;

use ddg_core::topology::NONE;
use ddg_core::{Geometry, TriMesh};
use ddg_integer_coords::{IntegerError, IntegerRefineOptions, IntegerRefineReport, IntegerTriangulation};
use ddg_intrinsic::{
    Complex, IntrinsicError, IntrinsicTriangulation, RefineOptions, RefineReport, RefineRules,
    SurfacePoint, Triangulation,
};
pub use ddg_intrinsic_apps::representation::Representation;
use ddg_intrinsic_apps::representation::InputPoints;

use crate::error::IemError;
use crate::mesh::{Certificate, CoarseMesh, CoarseMeshOptions, InputPoint, Location};
use crate::undo::UndoLog;

/// The tolerance with which a located barycentric coordinate may be negative before it is an
/// error rather than clamped to zero (S16).
pub const LOCATE_CLAMP: f64 = 1e-9;

/// What can go wrong in pre- or post-refinement.
#[derive(Debug)]
pub enum RefineError {
    /// An error of the coarse mesh (invalid extra point, pending removal, certificate).
    Coarse(IemError),
    /// The coarse mesh is a Δ-complex that the refinement crates cannot take as input (see the
    /// [module documentation](self)).
    NotSimplicial {
        /// What is wrong, naming the faces, edges and fine vertex ids.
        what: String,
    },
    /// An error of the integer-coordinate triangulation.
    Integer(IntegerError),
    /// An error of the signpost triangulation or of a point query.
    Intrinsic(IntrinsicError),
    /// A tracked point could not be located in the refined triangulation.
    Locate {
        /// The point id.
        point: usize,
        /// Why.
        what: String,
    },
}

impl fmt::Display for RefineError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            RefineError::Coarse(e) => write!(f, "{e}"),
            RefineError::NotSimplicial { what } => write!(
                f,
                "the coarse mesh is not a simplicial mesh, which the refinement needs as its \
                 input: {what}"
            ),
            RefineError::Integer(e) => write!(f, "integer-coordinate triangulation: {e}"),
            RefineError::Intrinsic(e) => write!(f, "intrinsic triangulation: {e}"),
            RefineError::Locate { point, what } => {
                write!(f, "tracked point {point} could not be located: {what}")
            }
        }
    }
}

impl std::error::Error for RefineError {}

impl From<IemError> for RefineError {
    fn from(e: IemError) -> Self {
        RefineError::Coarse(e)
    }
}

impl From<IntegerError> for RefineError {
    fn from(e: IntegerError) -> Self {
        RefineError::Integer(e)
    }
}

impl From<IntrinsicError> for RefineError {
    fn from(e: IntrinsicError) -> Self {
        RefineError::Intrinsic(e)
    }
}

/// An intrinsic triangulation that can say where points of its input mesh are (`T⁰ → T¹`), in
/// batches. Implemented for the signposts ([`IntrinsicTriangulation`]) and for integer
/// coordinates ([`IntegerTriangulation`], which caches the common subdivision per input face).
pub trait LocatedTriangulation: InputPoints {
    /// [`InputPoints::locate_on_triangulation`] for many points.
    fn locate_points(&self, points: &[SurfacePoint]) -> Vec<Result<SurfacePoint, IntrinsicError>> {
        points
            .iter()
            .map(|&p| self.locate_on_triangulation(p))
            .collect()
    }
}

impl LocatedTriangulation for IntrinsicTriangulation {}

impl LocatedTriangulation for IntegerTriangulation {
    fn locate_points(&self, points: &[SurfacePoint]) -> Vec<Result<SurfacePoint, IntrinsicError>> {
        let mut locator = self.point_locator();
        points
            .iter()
            .map(|&p| locator.point_on_triangulation(p))
            .collect()
    }
}

/// Options of [`CoarseMesh::from_refined`] and [`CoarseMesh::post_refine`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RefinementOptions {
    /// The representation of the intrinsic triangulation (default: integer coordinates).
    pub representation: Representation,
    /// The refinement's options: angle bound (default 25°), circumradius bound, insertion
    /// limit, Delaunay test.
    pub refine: RefineOptions,
    /// The rules; `None`: the representation's default ([`RefineRules::Chew`] for integers,
    /// [`RefineRules::Ruppert`] for signposts).
    pub rules: Option<RefineRules>,
    /// Whether the integer triangulation certifies itself at the end of refinement (its own
    /// default, yes). Ignored for signposts.
    pub certify: bool,
}

impl Default for RefinementOptions {
    fn default() -> Self {
        RefinementOptions {
            representation: Representation::DEFAULT,
            refine: RefineOptions::default(),
            rules: None,
            certify: true,
        }
    }
}

/// What a refinement did.
#[derive(Debug, Clone)]
pub struct RefinementReport {
    /// The representation used.
    pub representation: Representation,
    /// The rules followed.
    pub rules: RefineRules,
    /// Vertices of the input of the refinement (`T⁰`).
    pub n_input_vertices: usize,
    /// Vertices of the refined `T¹`.
    pub n_vertices: usize,
    /// `n_vertices − n_input_vertices`: vertices added (net of the removals of Chew's rules).
    pub n_added: usize,
    /// The shared driver's report (insertions, flips, removals, the angles, what is left below
    /// the bound).
    pub refine: RefineReport,
    /// The integer representation's extras (face splits, snaps, certificate), if used.
    pub integer: Option<IntegerRefineReport>,
    /// The smallest corner angle of the resulting coarse mesh (all corners).
    pub min_angle: f64,
    /// The largest corner angle of the resulting coarse mesh.
    pub max_angle: f64,
    /// The smallest angle sum `Θ` over the vertices of the input of the refinement.
    pub min_input_angle_sum: f64,
    /// Whether the input of the refinement is closed.
    pub closed: bool,
    /// Whether the guarantee of Gillespie et al. 2021a applies: closed and every input angle
    /// sum `≥ 60°`.
    pub guarantee_applies: bool,
    /// Tracked points located by the `T⁰ → T¹` query (not counting vertices that stay vertices).
    pub n_located: usize,
    /// Located barycentric coordinates in `[−1e-9, 0)` that were clamped to zero.
    pub n_clamped: usize,
    /// Located points on a vertex or edge of `T¹`, put in a face (S16).
    pub n_on_boundary_of_face: usize,
    /// The certificate of the resulting coarse mesh (`certify(false)`).
    pub certificate: Certificate,
    /// Seconds spent in the refinement itself (building `T⁰` included).
    pub seconds_refine: f64,
    /// Seconds spent locating the tracked points.
    pub seconds_locate: f64,
}

/// The raw pieces of a coarse mesh (see [`CoarseMesh::assemble`]).
struct Assembly {
    complex: Complex,
    length: Vec<f64>,
    fine_of: Vec<u32>,
    complex_of: Vec<u32>,
    fixed: Vec<bool>,
    n_input: u32,
    points: Vec<Location>,
    n_clamped: u64,
    options: CoarseMeshOptions,
}

impl CoarseMesh {
    /// A coarse mesh from raw pieces: the angular coordinates and angle sums are rebuilt from
    /// the lengths with the outgoing halfedges as reference (S4), the per-face point index from
    /// the points.
    fn assemble(a: Assembly) -> CoarseMesh {
        let Assembly {
            complex,
            length,
            fine_of,
            complex_of,
            fixed,
            n_input,
            points,
            n_clamped,
            options,
        } = a;
        let n = complex.n_vertices();
        let mut face_points: Vec<Vec<u32>> = vec![Vec::new(); complex.n_faces()];
        let mut point_slot = vec![NONE; points.len()];
        for (p, loc) in points.iter().enumerate() {
            if let Location::Face { face, .. } = *loc {
                point_slot[p] = face_points[face as usize].len() as u32;
                face_points[face as usize].push(p as u32);
            }
        }
        let mut mesh = CoarseMesh {
            reference: (0..n)
                .map(|v| complex.outgoing_halfedge(v) as u32)
                .collect(),
            phi: vec![0.0; complex.n_halfedges()],
            angle_sum: vec![0.0; n],
            fine_of,
            complex_of,
            fixed,
            n_input,
            points,
            point_slot,
            face_points,
            options,
            n_clamped,
            undo: UndoLog::default(),
            pending: None,
            vertex_mark: vec![0; n],
            edge_mark: vec![0; complex.n_edges()],
            epoch: 0,
            scratch_fan: Vec::new(),
            complex,
            length,
        };
        for v in 0..n {
            mesh.refresh_vertex(v);
        }
        mesh
    }

    /// The coarse mesh of an intrinsic triangulation `T¹` of the input (pre-refinement, SPEC
    /// §3.1, §6.1; see the [module documentation](self)): the complex and lengths of `tri`; the
    /// vertices `0..tri.n_input_vertices()` are the fine vertices, tracked as `Vertex` points;
    /// the other vertices of `tri` (inserted ones) have fine ids from `|V|` on, are not tracked
    /// and not fixed; each of `points` (a point of the input surface `T⁰`: input face and
    /// barycentric coordinates) is located in `tri` by the `T⁰ → T¹` query and is tracked as
    /// point `|V| + index`.
    ///
    /// # Errors
    /// [`RefineError::Coarse`] with [`IemError::InvalidPoint`] for a point that is not one of the
    /// input, and with [`IemError::Certify`] if the result does not certify (a triangulation
    /// with a face that violates the triangle inequality); [`RefineError::Locate`] naming the
    /// point when the query fails or places it outside its face by more than rounding.
    pub fn from_intrinsic<T: LocatedTriangulation>(
        tri: &T,
        points: &[InputPoint],
        options: CoarseMeshOptions,
    ) -> Result<CoarseMesh, RefineError> {
        let mut stats = LocateStats::default();
        from_intrinsic_with(tri, points, options, &mut stats)
    }

    /// Pre-refinement end to end: builds the intrinsic triangulation of the input
    /// (`options.representation`), Delaunay-refines it, and returns the coarse mesh of the
    /// refined triangulation with the extra `points` located
    /// ([`from_intrinsic`](Self::from_intrinsic)) and the report. The coarse mesh is not
    /// simplified; with a removal order of one's own (or `simplify`) it is coarsened.
    ///
    /// # Errors
    /// [`RefineError::Coarse`] for an invalid mesh, geometry or point, [`RefineError::Integer`]
    /// and [`RefineError::Intrinsic`] for the refinement's, and those of `from_intrinsic`.
    pub fn from_refined(
        mesh: &TriMesh,
        geometry: &Geometry,
        points: &[InputPoint],
        mesh_options: CoarseMeshOptions,
        options: &RefinementOptions,
    ) -> Result<(CoarseMesh, RefinementReport), RefineError> {
        geometry.check_mesh(mesh).map_err(IemError::from)?;
        for (index, p) in points.iter().enumerate() {
            check_input_point(index, p, mesh.n_faces())?;
        }
        let t0 = Instant::now();
        match options.representation {
            Representation::Integer => {
                let mut tri = IntegerTriangulation::new(mesh, geometry)?;
                let report = refine_integer(&mut tri, options)?;
                finish_pre(&tri, points, mesh_options, options, report, t0)
            }
            Representation::Signposts => {
                let mut tri = IntrinsicTriangulation::new(mesh, geometry)?;
                let report = refine_signposts(&mut tri, options)?;
                finish_pre(&tri, points, mesh_options, options, report, t0)
            }
        }
    }

    /// The coarse mesh as the input of a simplicial mesh for the refinement crates: the
    /// [`TriMesh`] of its faces (complex vertex ids) and the lengths of its edges in that
    /// mesh's edge order, if the coarse mesh is simplicial, else the reason it is not.
    ///
    /// # Errors
    /// [`RefineError::NotSimplicial`] (see the [module documentation](self)),
    /// [`RefineError::Coarse`] with [`IemError::RemovalPending`] while a removal is pending.
    pub fn as_simplicial_input(&self) -> Result<(TriMesh, Geometry), RefineError> {
        if self.pending.is_some() {
            return Err(IemError::RemovalPending.into());
        }
        let c = &self.complex;
        let fine = |v: u32| self.fine_of[v as usize];
        let not = |what: String| Err(RefineError::NotSimplicial { what });
        for (f, face) in c.faces().iter().enumerate() {
            let [a, b, d] = *face;
            if a == b || b == d || a == d {
                return not(format!(
                    "face {f} (fine vertices {}, {}, {}) uses a vertex twice (a self-edge or a \
                     loop of the Δ-complex)",
                    fine(a),
                    fine(b),
                    fine(d)
                ));
            }
        }
        let tri_mesh = match TriMesh::new(c.n_vertices(), c.faces().to_vec()) {
            Ok(m) => m,
            Err(e) => {
                return not(format!(
                    "{e} (vertex ids in this message are those of the compacted coarse mesh; \
                     fine id = CoarseMesh::fine_vertex)"
                ));
            }
        };
        for h in 0..c.n_halfedges() {
            if tri_mesh.twin(h) != c.twin(h) {
                let [a, b] = [c.from_vertex(h), c.to_vertex(h)];
                return not(format!(
                    "halfedge {h} (face {}, fine vertices {} to {}) is glued differently in the \
                     simplicial mesh of the same faces: the coarse mesh has two different edges \
                     between the same two vertices (or two faces on the same three vertices)",
                    h / 3,
                    fine(a as u32),
                    fine(b as u32)
                ));
            }
        }
        let mut lengths = vec![0.0; tri_mesh.n_edges()];
        for h in 0..c.n_halfedges() {
            lengths[tri_mesh.edge(h)] = self.length[c.edge(h)];
        }
        let geometry = Geometry::from_edge_lengths(&tri_mesh, lengths).map_err(IemError::from)?;
        Ok((tri_mesh, geometry))
    }

    /// Post-refinement (SPEC §6.2; see the [module documentation](self)): Delaunay refinement of
    /// this coarse mesh as the input of a new intrinsic triangulation, with the tracked points
    /// composed through it. `self` is not changed.
    ///
    /// # Errors
    /// [`RefineError::NotSimplicial`] if the coarse mesh is a Δ-complex that cannot be the input
    /// of the refinement; [`IemError::RemovalPending`]; the refinement's errors; and
    /// [`RefineError::Locate`] naming a tracked point that could not be located.
    pub fn post_refine(
        &self,
        options: &RefinementOptions,
    ) -> Result<(CoarseMesh, RefinementReport), RefineError> {
        let (mesh, geometry) = self.as_simplicial_input()?;
        let t0 = Instant::now();
        match options.representation {
            Representation::Integer => {
                let mut tri = IntegerTriangulation::new(&mesh, &geometry)?;
                let report = refine_integer(&mut tri, options)?;
                self.finish_post(&tri, options, report, t0)
            }
            Representation::Signposts => {
                let mut tri = IntrinsicTriangulation::new(&mesh, &geometry)?;
                let report = refine_signposts(&mut tri, options)?;
                self.finish_post(&tri, options, report, t0)
            }
        }
    }

    fn finish_post<T: LocatedTriangulation>(
        &self,
        tri: &T,
        options: &RefinementOptions,
        refined: Refined,
        t0: Instant,
    ) -> Result<(CoarseMesh, RefinementReport), RefineError> {
        let seconds_refine = t0.elapsed().as_secs_f64();
        let complex = tri.complex().clone();
        let m = self.complex.n_vertices();
        let n_new = complex.n_vertices();
        let n_old_ids = self.complex_of.len();
        // Ids: the coarse vertices keep theirs (the refinement's input vertices come first),
        // inserted vertices get fresh ones.
        let mut fine_of = Vec::with_capacity(n_new);
        fine_of.extend_from_slice(&self.fine_of[..m]);
        fine_of.extend((m..n_new).map(|v| (n_old_ids + v - m) as u32));
        let mut complex_of = vec![NONE; n_old_ids + n_new.saturating_sub(m)];
        for (v, &fine) in fine_of.iter().enumerate() {
            complex_of[fine as usize] = v as u32;
        }
        let mut fixed = self.fixed.clone();
        fixed.resize(complex_of.len(), false);
        // The points: vertices stay, faces are located.
        let mut locations: Vec<Location> = self.points.clone();
        let mut ids = Vec::new();
        let mut queries = Vec::new();
        for (p, loc) in self.points.iter().enumerate() {
            if let Location::Face { face, bary } = *loc {
                ids.push(p);
                queries.push(SurfacePoint::Face { face, bary });
            }
        }
        let t1 = Instant::now();
        let located = tri.locate_points(&queries);
        let seconds_locate = t1.elapsed().as_secs_f64();
        let mut stats = LocateStats::default();
        for ((&p, answer), query) in ids.iter().zip(located).zip(&queries) {
            let sp = answer.map_err(|e| RefineError::Locate {
                point: p,
                what: format!("{e} (query {query})"),
            })?;
            let (face, bary) = place(&complex, sp, &mut stats)
                .map_err(|what| RefineError::Locate { point: p, what })?;
            locations[p] = Location::Face { face, bary };
        }
        let n_clamped = self.n_clamped + stats.clamped as u64;
        let result = CoarseMesh::assemble(Assembly {
            complex,
            length: tri.edge_lengths().to_vec(),
            fine_of,
            complex_of,
            fixed,
            n_input: self.n_input,
            points: locations,
            n_clamped,
            options: self.options,
        });
        let certificate = result.certify(false)?;
        let report = make_report(
            options,
            refined,
            tri.input().vertex_angle_sums(),
            tri.input().mesh().n_boundary_edges() == 0,
            m,
            &result,
            &stats,
            ids.len(),
            certificate,
            seconds_refine,
            seconds_locate,
        );
        Ok((result, report))
    }
}

/// Counters of the point placement.
#[derive(Debug, Default, Clone, Copy)]
struct LocateStats {
    clamped: usize,
    on_boundary: usize,
}

/// What a refinement returned, for the report.
struct Refined {
    rules: RefineRules,
    refine: RefineReport,
    integer: Option<IntegerRefineReport>,
}

fn refine_integer(
    tri: &mut IntegerTriangulation,
    options: &RefinementOptions,
) -> Result<Refined, RefineError> {
    let rules = options.rules.unwrap_or(RefineRules::Chew);
    let report = tri.refine(&IntegerRefineOptions {
        refine: options.refine,
        rules,
        certify: options.certify,
        ..IntegerRefineOptions::default()
    })?;
    Ok(Refined {
        rules,
        refine: report.refine.clone(),
        integer: Some(report),
    })
}

fn refine_signposts(
    tri: &mut IntrinsicTriangulation,
    options: &RefinementOptions,
) -> Result<Refined, RefineError> {
    let rules = options.rules.unwrap_or(RefineRules::Ruppert);
    let report = tri.refine_with_rules(&options.refine, rules)?;
    Ok(Refined {
        rules,
        refine: report,
        integer: None,
    })
}

#[allow(clippy::too_many_arguments)]
fn make_report(
    options: &RefinementOptions,
    refined: Refined,
    input_angle_sums: &[f64],
    closed: bool,
    n_input_vertices: usize,
    result: &CoarseMesh,
    stats: &LocateStats,
    n_located: usize,
    certificate: Certificate,
    seconds_refine: f64,
    seconds_locate: f64,
) -> RefinementReport {
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for h in 0..result.complex.n_halfedges() {
        let a = result.corner_angle(h);
        lo = lo.min(a);
        hi = hi.max(a);
    }
    let min_sum = input_angle_sums.iter().copied().fold(f64::INFINITY, f64::min);
    RefinementReport {
        representation: options.representation,
        rules: refined.rules,
        n_input_vertices,
        n_vertices: result.n_vertices(),
        n_added: result.n_vertices().saturating_sub(n_input_vertices),
        refine: refined.refine,
        integer: refined.integer,
        min_angle: lo,
        max_angle: hi,
        min_input_angle_sum: min_sum,
        closed,
        guarantee_applies: closed && min_sum >= std::f64::consts::FRAC_PI_3,
        n_located,
        n_clamped: stats.clamped,
        n_on_boundary_of_face: stats.on_boundary,
        certificate,
        seconds_refine,
        seconds_locate,
    }
}

fn finish_pre<T: LocatedTriangulation>(
    tri: &T,
    points: &[InputPoint],
    mesh_options: CoarseMeshOptions,
    options: &RefinementOptions,
    refined: Refined,
    t0: Instant,
) -> Result<(CoarseMesh, RefinementReport), RefineError> {
    let seconds_refine = t0.elapsed().as_secs_f64();
    let mut stats = LocateStats::default();
    let t1 = Instant::now();
    let result = from_intrinsic_with(tri, points, mesh_options, &mut stats)?;
    let seconds_locate = t1.elapsed().as_secs_f64();
    let certificate = result.certify(false)?;
    let report = make_report(
        options,
        refined,
        tri.input().vertex_angle_sums(),
        tri.input().mesh().n_boundary_edges() == 0,
        tri.n_input_vertices(),
        &result,
        &stats,
        points.len(),
        certificate,
        seconds_refine,
        seconds_locate,
    );
    Ok((result, report))
}

fn from_intrinsic_with<T: LocatedTriangulation>(
    tri: &T,
    points: &[InputPoint],
    options: CoarseMeshOptions,
    stats: &mut LocateStats,
) -> Result<CoarseMesh, RefineError> {
    let complex = tri.complex().clone();
    let n_total = complex.n_vertices();
    let n_input = tri.n_input_vertices();
    let n_faces_input = tri.input().mesh().n_faces();
    for (index, p) in points.iter().enumerate() {
        check_input_point(index, p, n_faces_input)?;
    }
    let queries: Vec<SurfacePoint> = points
        .iter()
        .map(|p| SurfacePoint::Face {
            face: p.face,
            bary: p.bary,
        })
        .collect();
    let located = tri.locate_points(&queries);
    let mut locations: Vec<Location> = (0..n_input as u32).map(Location::Vertex).collect();
    for (index, answer) in located.into_iter().enumerate() {
        let point = n_input + index;
        let sp = answer.map_err(|e| RefineError::Locate {
            point,
            what: format!("{e} (query {})", queries[index]),
        })?;
        let (face, bary) =
            place(&complex, sp, stats).map_err(|what| RefineError::Locate { point, what })?;
        locations.push(Location::Face { face, bary });
    }
    let mesh = CoarseMesh::assemble(Assembly {
        length: tri.edge_lengths().to_vec(),
        fine_of: (0..n_total as u32).collect(),
        complex_of: (0..n_total as u32).collect(),
        fixed: vec![false; n_total],
        n_input: n_input as u32,
        points: locations,
        n_clamped: stats.clamped as u64,
        options,
        complex,
    });
    mesh.certify(false)?;
    Ok(mesh)
}

/// The validation of an extra input point, as [`CoarseMesh::from_mesh`] does it.
fn check_input_point(index: usize, p: &InputPoint, n_faces: usize) -> Result<(), IemError> {
    let bad = |what: String| Err(IemError::InvalidPoint { index, what });
    if p.face as usize >= n_faces {
        return bad(format!(
            "face {} is out of range (the mesh has {n_faces} faces)",
            p.face
        ));
    }
    if !p.bary.iter().all(|b| b.is_finite() && *b >= 0.0) {
        return bad(format!(
            "barycentric coordinates {:?} are not finite and non-negative",
            p.bary
        ));
    }
    let sum = p.bary[0] + p.bary[1] + p.bary[2];
    if (sum - 1.0).abs() > 1e-12 {
        return bad(format!(
            "barycentric coordinates {:?} sum to {sum}, not one",
            p.bary
        ));
    }
    Ok(())
}

/// A located point as a face and barycentric coordinates of `complex` (S16): a vertex is put
/// in its outgoing halfedge's face, an edge point in the face of the edge's halfedge (the
/// coordinates of its end points), tiny negatives clamped, the coordinates renormalized to
/// sum to one.
fn place(
    complex: &Complex,
    sp: SurfacePoint,
    stats: &mut LocateStats,
) -> Result<(u32, [f64; 3]), String> {
    let (face, mut bary) = match sp {
        SurfacePoint::Vertex(v) => {
            if v as usize >= complex.n_vertices() {
                return Err(format!("located at vertex {v}, which does not exist"));
            }
            let h = complex.outgoing_halfedge(v as usize);
            let mut b = [0.0; 3];
            b[h % 3] = 1.0;
            stats.on_boundary += 1;
            ((h / 3) as u32, b)
        }
        SurfacePoint::Edge { edge, t } => {
            if edge as usize >= complex.n_edges() || !t.is_finite() {
                return Err(format!("located at the invalid point {sp}"));
            }
            let h = complex.edge_halfedge(edge as usize);
            let mut b = [0.0; 3];
            b[h % 3] = 1.0 - t;
            b[((h % 3) + 1) % 3] = t;
            stats.on_boundary += 1;
            ((h / 3) as u32, b)
        }
        SurfacePoint::Face { face, bary } => {
            if face as usize >= complex.n_faces() {
                return Err(format!("located in face {face}, which does not exist"));
            }
            (face, bary)
        }
    };
    for b in bary.iter_mut() {
        if !b.is_finite() || *b < -LOCATE_CLAMP {
            return Err(format!(
                "located at {sp}, with coordinates {bary:?} outside the face by more than \
                 {LOCATE_CLAMP:e}"
            ));
        }
        if *b < 0.0 {
            *b = 0.0;
            stats.clamped += 1;
        }
    }
    let sum = bary[0] + bary[1] + bary[2];
    if !(sum > 0.5 && sum < 1.5) {
        return Err(format!("located at {sp}, whose coordinates sum to {sum}"));
    }
    for b in bary.iter_mut() {
        *b /= sum;
    }
    Ok((face, bary))
}
