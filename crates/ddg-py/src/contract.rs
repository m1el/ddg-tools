//! The data that crosses the Python boundary, as Rust types.
//!
//! [`UnwrapInput`] is the validated form of the arrays Blender hands over; [`UnwrapOptions`] the
//! keyword options; [`UnwrapOutput`] what a pipeline returns. A pipeline is any function
//! `fn(&UnwrapInput, &UnwrapOptions) -> Result<UnwrapOutput, UnwrapError>`
//! (see [`crate::pipeline::run_pipeline`]).
//!
//! # What `UnwrapInput` guarantees
//!
//! After [`UnwrapInput::from_buffers`] succeeds:
//!
//! * every position is finite; there are fewer than `u32::MAX` vertices and corners;
//! * `corner_vertices[c] < n_vertices` for every polygon corner `c` ("loop" in Blender);
//! * `polygon_offsets` has `n_polygons + 1` entries, starts at 0, ends at the number of corners,
//!   and every polygon has at least 3 corners; polygon `p` owns the corners
//!   `polygon_offsets[p]..polygon_offsets[p + 1]`, in winding order;
//! * `triangles` holds corner indices, and the triangles of each polygon are a triangulation of
//!   it: `corners − 2` triangles, wound like the polygon, every side of the polygon in exactly
//!   one triangle and every diagonal in exactly two;
//! * `seam_edges` are vertex pairs `[min, max]` with `min < max`, sorted, without duplicates;
//! * `pinned_vertices` and `forbidden_vertices` are sorted and free of duplicates,
//!   `prescribed_cones` is sorted by vertex with every vertex at most once and every angle finite
//!   and below 2π, and no vertex is in more than one of the three sets.
//!
//! What is **not** checked here, because it is the pipeline's business: whether the mesh is
//! manifold and consistently oriented, whether a seam edge is an edge of the mesh, whether a
//! polygon repeats a vertex, and whether triangles are degenerate.

use crate::error::{InputError, TriangulationProblem, UnwrapError};

/// The decoded but not yet validated arrays of one `unwrap` call, flat as Python passes them.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct MeshBuffers {
    /// Vertex positions, 3 per vertex.
    pub positions: Vec<f64>,
    /// Vertex index of every polygon corner.
    pub corner_vertices: Vec<u32>,
    /// First corner of every polygon (Blender's `MeshPolygon.loop_start`), one per polygon.
    pub polygon_starts: Vec<u32>,
    /// Corner indices of Blender's loop triangles, 3 per triangle.
    pub triangles: Vec<u32>,
    /// Seam edges as vertex pairs, 2 per edge.
    pub seam_edges: Vec<u32>,
    /// Vertices where a cone is free (no cost).
    pub pinned_vertices: Vec<u32>,
    /// Vertices with a prescribed cone angle.
    pub prescribed_vertices: Vec<u32>,
    /// The prescribed cone angles, parallel to `prescribed_vertices`.
    pub prescribed_angles: Vec<f64>,
    /// Vertices where a cone is forbidden.
    pub forbidden_vertices: Vec<u32>,
}

/// A cone singularity: a vertex and its cone angle `µ = 2π − (angle sum in the flat metric)`.
/// Positive is a cone tip (see the math conventions in the repository's `CLAUDE.md`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cone {
    /// The vertex, an index into the caller's vertex array.
    pub vertex: u32,
    /// The cone angle in radians.
    pub angle: f64,
}

/// A polygon mesh with Blender's triangulation and the user's steering, validated as described
/// in the [module documentation](self).
#[derive(Debug, Clone, PartialEq)]
pub struct UnwrapInput {
    positions: Vec<[f64; 3]>,
    corner_vertices: Vec<u32>,
    polygon_offsets: Vec<u32>,
    triangles: Vec<[u32; 3]>,
    seam_edges: Vec<[u32; 2]>,
    pinned_vertices: Vec<u32>,
    prescribed_cones: Vec<Cone>,
    forbidden_vertices: Vec<u32>,
}

fn check_range(
    name: &'static str,
    values: &[u32],
    count: usize,
    of: &'static str,
) -> Result<(), InputError> {
    match values.iter().position(|&v| v as usize >= count) {
        Some(position) => Err(InputError::IndexOutOfRange {
            name,
            position,
            value: values[position],
            count,
            of,
        }),
        None => Ok(()),
    }
}

fn check_finite(name: &'static str, values: &[f64]) -> Result<(), InputError> {
    match values.iter().position(|v| !v.is_finite()) {
        Some(position) => Err(InputError::NonFinite {
            name,
            position,
            value: values[position],
        }),
        None => Ok(()),
    }
}

fn sorted_unique(mut values: Vec<u32>) -> Vec<u32> {
    values.sort_unstable();
    values.dedup();
    values
}

/// First element common to two sorted slices.
fn first_common(a: &[u32], b: &[u32]) -> Option<u32> {
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => return Some(a[i]),
        }
    }
    None
}

/// Turns per-polygon start offsets (Blender's `MeshPolygon.loop_start`, one per polygon) into
/// `n_polygons + 1` offsets by appending the number of corners, after checking that the first
/// polygon starts at 0, that the offsets increase and that every polygon has at least 3 corners.
///
/// # Errors
/// [`InputError::NoPolygons`], [`InputError::FirstPolygonStart`],
/// [`InputError::PolygonTooSmall`], [`InputError::TooLarge`].
pub fn polygon_offsets_from_starts(
    polygon_starts: Vec<u32>,
    n_corners: usize,
) -> Result<Vec<u32>, InputError> {
    if n_corners >= u32::MAX as usize {
        return Err(InputError::TooLarge {
            name: "corner_vertices",
            count: n_corners,
        });
    }
    if polygon_starts.is_empty() {
        return Err(InputError::NoPolygons);
    }
    if polygon_starts[0] != 0 {
        return Err(InputError::FirstPolygonStart {
            start: polygon_starts[0],
        });
    }
    let n_polygons = polygon_starts.len();
    let mut polygon_offsets = polygon_starts;
    polygon_offsets.push(n_corners as u32);
    for p in 0..n_polygons {
        let (start, end) = (polygon_offsets[p], polygon_offsets[p + 1]);
        if end < start || end - start < 3 {
            return Err(InputError::PolygonTooSmall {
                polygon: p,
                start,
                end,
            });
        }
    }
    Ok(polygon_offsets)
}

impl UnwrapInput {
    /// Validates decoded arrays and builds the input. See the [module documentation](self) for
    /// what is checked.
    ///
    /// Set-valued arguments are put in canonical form: seam edges become sorted `[min, max]`
    /// pairs without duplicates, and the pinned and forbidden vertex lists are sorted with
    /// duplicates removed. Nothing else is changed.
    ///
    /// # Errors
    /// An [`InputError`] naming the Python parameter and the offending entry.
    pub fn from_buffers(buffers: MeshBuffers) -> Result<UnwrapInput, InputError> {
        let MeshBuffers {
            positions,
            corner_vertices,
            polygon_starts,
            triangles,
            seam_edges,
            pinned_vertices,
            prescribed_vertices,
            prescribed_angles,
            forbidden_vertices,
        } = buffers;
        const LIMIT: usize = u32::MAX as usize;

        // Positions.
        if !positions.len().is_multiple_of(3) {
            return Err(InputError::LengthNotMultiple {
                name: "positions",
                len: positions.len(),
                width: 3,
                row: "vertex",
            });
        }
        let n_vertices = positions.len() / 3;
        if n_vertices >= LIMIT {
            return Err(InputError::TooLarge {
                name: "positions",
                count: n_vertices,
            });
        }
        check_finite("positions", &positions)?;
        let positions: Vec<[f64; 3]> = positions
            .chunks_exact(3)
            .map(|c| [c[0], c[1], c[2]])
            .collect();

        // Polygons.
        let n_corners = corner_vertices.len();
        if n_corners >= LIMIT {
            return Err(InputError::TooLarge {
                name: "corner_vertices",
                count: n_corners,
            });
        }
        check_range("corner_vertices", &corner_vertices, n_vertices, "vertices")?;
        let polygon_offsets = polygon_offsets_from_starts(polygon_starts, n_corners)?;
        let n_polygons = polygon_offsets.len() - 1;

        // Triangles.
        if !triangles.len().is_multiple_of(3) {
            return Err(InputError::LengthNotMultiple {
                name: "triangles",
                len: triangles.len(),
                width: 3,
                row: "triangle",
            });
        }
        check_range("triangles", &triangles, n_corners, "polygon corners")?;
        let expected = n_corners - 2 * n_polygons;
        if triangles.len() / 3 != expected {
            return Err(InputError::TriangleCount {
                got: triangles.len() / 3,
                expected,
                n_corners,
                n_polygons,
            });
        }
        let triangles: Vec<[u32; 3]> = triangles
            .chunks_exact(3)
            .map(|c| [c[0], c[1], c[2]])
            .collect();
        check_triangulation(&polygon_offsets, &triangles)?;

        // Seams.
        if !seam_edges.len().is_multiple_of(2) {
            return Err(InputError::LengthNotMultiple {
                name: "seam_edges",
                len: seam_edges.len(),
                width: 2,
                row: "edge",
            });
        }
        check_range("seam_edges", &seam_edges, n_vertices, "vertices")?;
        let mut seams: Vec<[u32; 2]> = Vec::with_capacity(seam_edges.len() / 2);
        for (row, pair) in seam_edges.chunks_exact(2).enumerate() {
            let (a, b) = (pair[0], pair[1]);
            if a == b {
                return Err(InputError::DegenerateSeamEdge { row, vertex: a });
            }
            seams.push([a.min(b), a.max(b)]);
        }
        seams.sort_unstable();
        seams.dedup();

        // Steering.
        check_range("pinned_vertices", &pinned_vertices, n_vertices, "vertices")?;
        check_range(
            "forbidden_vertices",
            &forbidden_vertices,
            n_vertices,
            "vertices",
        )?;
        check_range(
            "prescribed_vertices",
            &prescribed_vertices,
            n_vertices,
            "vertices",
        )?;
        if prescribed_vertices.len() != prescribed_angles.len() {
            return Err(InputError::PrescribedLengthMismatch {
                vertices: prescribed_vertices.len(),
                angles: prescribed_angles.len(),
            });
        }
        check_finite("prescribed_angles", &prescribed_angles)?;
        let mut prescribed_cones: Vec<Cone> = prescribed_vertices
            .iter()
            .zip(&prescribed_angles)
            .map(|(&vertex, &angle)| Cone { vertex, angle })
            .collect();
        if let Some(cone) = prescribed_cones
            .iter()
            .find(|c| c.angle >= std::f64::consts::TAU)
        {
            return Err(InputError::PrescribedAngleTooLarge {
                vertex: cone.vertex,
                angle: cone.angle,
            });
        }
        prescribed_cones.sort_by_key(|c| c.vertex);
        if let Some(w) = prescribed_cones
            .windows(2)
            .find(|w| w[0].vertex == w[1].vertex)
        {
            return Err(InputError::DuplicatePrescribed {
                vertex: w[0].vertex,
            });
        }
        let pinned_vertices = sorted_unique(pinned_vertices);
        let forbidden_vertices = sorted_unique(forbidden_vertices);
        let prescribed_sorted: Vec<u32> = prescribed_cones.iter().map(|c| c.vertex).collect();
        for (a, first, b, second) in [
            (&pinned_vertices, "pinned", &forbidden_vertices, "forbidden"),
            (&pinned_vertices, "pinned", &prescribed_sorted, "prescribed"),
            (
                &prescribed_sorted,
                "prescribed",
                &forbidden_vertices,
                "forbidden",
            ),
        ] {
            if let Some(vertex) = first_common(a, b) {
                return Err(InputError::ConflictingSteering {
                    vertex,
                    first,
                    second,
                });
            }
        }

        Ok(UnwrapInput {
            positions,
            corner_vertices,
            polygon_offsets,
            triangles,
            seam_edges: seams,
            pinned_vertices,
            prescribed_cones,
            forbidden_vertices,
        })
    }

    /// Number of vertices.
    #[inline]
    pub fn n_vertices(&self) -> usize {
        self.positions.len()
    }

    /// Number of polygon corners ("loops").
    #[inline]
    pub fn n_corners(&self) -> usize {
        self.corner_vertices.len()
    }

    /// Number of polygons.
    #[inline]
    pub fn n_polygons(&self) -> usize {
        self.polygon_offsets.len() - 1
    }

    /// Number of triangles of the triangulation.
    #[inline]
    pub fn n_triangles(&self) -> usize {
        self.triangles.len()
    }

    /// Vertex positions.
    #[inline]
    pub fn positions(&self) -> &[[f64; 3]] {
        &self.positions
    }

    /// Vertex index of every polygon corner.
    #[inline]
    pub fn corner_vertices(&self) -> &[u32] {
        &self.corner_vertices
    }

    /// Corner offsets of the polygons, `n_polygons + 1` entries: polygon `p` owns the corners
    /// `polygon_offsets[p]..polygon_offsets[p + 1]`.
    #[inline]
    pub fn polygon_offsets(&self) -> &[u32] {
        &self.polygon_offsets
    }

    /// The corners of polygon `p` as an index range into [`Self::corner_vertices`].
    #[inline]
    pub fn polygon_corners(&self, p: usize) -> std::ops::Range<usize> {
        self.polygon_offsets[p] as usize..self.polygon_offsets[p + 1] as usize
    }

    /// The triangulation: three **corner** indices per triangle, wound like their polygon.
    /// `corner_vertices[triangles[t][k]]` is the vertex.
    #[inline]
    pub fn triangles(&self) -> &[[u32; 3]] {
        &self.triangles
    }

    /// For every corner, the polygon it belongs to.
    pub fn corner_polygons(&self) -> Vec<u32> {
        let mut out = vec![0u32; self.n_corners()];
        for p in 0..self.n_polygons() {
            out[self.polygon_corners(p)].fill(p as u32);
        }
        out
    }

    /// Seam edges as sorted, unique vertex pairs `[min, max]`.
    #[inline]
    pub fn seam_edges(&self) -> &[[u32; 2]] {
        &self.seam_edges
    }

    /// Vertices where a cone costs nothing (sorted, unique).
    #[inline]
    pub fn pinned_vertices(&self) -> &[u32] {
        &self.pinned_vertices
    }

    /// Cones with a prescribed angle (sorted by vertex, unique).
    #[inline]
    pub fn prescribed_cones(&self) -> &[Cone] {
        &self.prescribed_cones
    }

    /// Vertices where a cone is forbidden (sorted, unique).
    #[inline]
    pub fn forbidden_vertices(&self) -> &[u32] {
        &self.forbidden_vertices
    }
}

/// Checks that the triangles of every polygon triangulate it. `polygon_offsets` and the corner
/// indices in `triangles` are already range-checked.
fn check_triangulation(polygon_offsets: &[u32], triangles: &[[u32; 3]]) -> Result<(), InputError> {
    let n_polygons = polygon_offsets.len() - 1;
    let n_corners = *polygon_offsets.last().expect("offsets are never empty") as usize;
    let mut corner_polygon = vec![0u32; n_corners];
    for p in 0..n_polygons {
        corner_polygon[polygon_offsets[p] as usize..polygon_offsets[p + 1] as usize].fill(p as u32);
    }
    let size = |p: u32| (polygon_offsets[p as usize + 1] - polygon_offsets[p as usize]) as usize;
    let bad = |polygon: u32, problem: TriangulationProblem| InputError::BadTriangulation {
        polygon,
        n_corners: size(polygon),
        problem,
    };
    // Next corner around the polygon.
    let next = |c: u32, p: u32| {
        if c + 1 == polygon_offsets[p as usize + 1] {
            polygon_offsets[p as usize]
        } else {
            c + 1
        }
    };

    let mut triangle_count = vec![0u32; n_polygons];
    // Per corner: how many triangles contain the polygon side that starts there.
    let mut side_count = vec![0u8; n_corners];
    // Diagonals as (lower corner, higher corner, runs from lower to higher).
    let mut diagonals: Vec<(u32, u32, bool)> = Vec::new();
    for (t, &[a, b, c]) in triangles.iter().enumerate() {
        let polygons = [
            corner_polygon[a as usize],
            corner_polygon[b as usize],
            corner_polygon[c as usize],
        ];
        let p = polygons[0];
        if polygons[1] != p || polygons[2] != p {
            return Err(InputError::TriangleSpansPolygons {
                triangle: t,
                corners: [a, b, c],
                polygons,
            });
        }
        if a == b || b == c || a == c {
            return Err(bad(p, TriangulationProblem::RepeatedCorner { triangle: t }));
        }
        triangle_count[p as usize] += 1;
        for (from, to) in [(a, b), (b, c), (c, a)] {
            if next(from, p) == to {
                side_count[from as usize] = side_count[from as usize].saturating_add(1);
            } else if next(to, p) == from {
                return Err(bad(
                    p,
                    TriangulationProblem::ReversedTriangle { triangle: t },
                ));
            } else {
                diagonals.push((from.min(to), from.max(to), from < to));
            }
        }
    }
    for (p, &count) in triangle_count.iter().enumerate() {
        if count as usize != size(p as u32) - 2 {
            return Err(bad(
                p as u32,
                TriangulationProblem::TriangleCount(count as usize),
            ));
        }
    }
    if let Some(corner) = side_count.iter().position(|&n| n != 1) {
        return Err(bad(
            corner_polygon[corner],
            TriangulationProblem::SideCoverage {
                corner: corner as u32,
                count: side_count[corner] as usize,
            },
        ));
    }
    // Every diagonal must appear exactly twice, once in each direction. After sorting, the two
    // halves of a diagonal are adjacent, `false` (high to low) first.
    diagonals.sort_unstable();
    let mut i = 0;
    while i < diagonals.len() {
        let (lo, hi, forward) = diagonals[i];
        let paired = !forward
            && i + 1 < diagonals.len()
            && diagonals[i + 1] == (lo, hi, true)
            && (i + 2 == diagonals.len() || (diagonals[i + 2].0, diagonals[i + 2].1) != (lo, hi));
        if !paired {
            return Err(bad(
                corner_polygon[lo as usize],
                TriangulationProblem::UnpairedDiagonal { corners: [lo, hi] },
            ));
        }
        i += 2;
    }
    Ok(())
}

/// The keyword options of `ddg_tools.unwrap`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UnwrapOptions {
    /// Trade-off `λ` between area distortion and total cone angle (paper Sec. 6.2). Dimensionless
    /// under the area normalization of `CLAUDE.md`; default 0.5, useful range 0.1 to 1. Larger
    /// means fewer cones.
    pub lambda: f64,
    /// Allow only cones with a positive angle (cone tips, no saddles).
    pub positive_cones_only: bool,
    /// Optimize the scale factor along the boundary of every piece (the mesh's own boundary
    /// and the seams) instead of holding it at zero: lower distortion, but lengths along the
    /// two sides of a seam may differ (`max_seam_length_mismatch` in the report). Default
    /// `false`.
    pub optimize_boundary: bool,
    /// Split charts that overlap themselves automatically (`ddg_uv::SplitOptions`): a cut that
    /// separates the two layers is added to the seams and the charts are flattened again.
    /// Default `false`.
    pub split_overlaps: bool,
    /// With `split_overlaps`: at most this many cuts are kept. Default 32.
    pub split_max: usize,
    /// With `split_overlaps`: a chart is split while more than this fraction of its footprint
    /// is covered twice. Default 0.001.
    pub split_threshold: f64,
    /// The feedback loop of the cone solver (`ddg_uv::UnwrapOptions::cone_feedback`): lower
    /// the distortion weight where the solver's linearized scale factor is far from the one
    /// the flattening realizes (thin spikes) and solve again. Default `true`.
    pub cone_feedback: bool,
    /// With `cone_feedback`: the threshold of the mismatch, in log units
    /// (`ddg_optimal_cones::FeedbackOptions::threshold`). Default 20.
    pub feedback_threshold: f64,
    /// With `cone_feedback`: the most solves, the first included. Default 4.
    pub feedback_max_rounds: usize,
}

impl Default for UnwrapOptions {
    fn default() -> Self {
        let split = ddg_uv::SplitOptions::default();
        let feedback = ddg_uv::FeedbackOptions::default();
        UnwrapOptions {
            lambda: 0.1,
            positive_cones_only: true,
            optimize_boundary: false,
            split_overlaps: false,
            split_max: split.max_splits,
            split_threshold: split.min_overlap_fraction,
            cone_feedback: true,
            feedback_threshold: feedback.threshold,
            feedback_max_rounds: feedback.max_rounds,
        }
    }
}

impl UnwrapOptions {
    /// Checks the options: `lambda` must be finite and not negative.
    ///
    /// # Errors
    /// [`InputError::BadOption`].
    pub fn validate(&self) -> Result<(), InputError> {
        if !(self.lambda.is_finite() && self.lambda >= 0.0) {
            return Err(InputError::BadOption {
                name: "lambda_",
                value: format!("{}", self.lambda),
                requirement: "lambda must be a finite number that is not negative",
            });
        }
        if self.split_overlaps && !(self.split_threshold.is_finite() && (0.0..1.0).contains(&self.split_threshold)) {
            return Err(InputError::BadOption {
                name: "split_threshold",
                value: format!("{}", self.split_threshold),
                requirement: "must be at least 0 and below 1",
            });
        }
        if self.cone_feedback
            && !(self.feedback_threshold.is_finite() && self.feedback_threshold > 0.0)
        {
            return Err(InputError::BadOption {
                name: "feedback_threshold",
                value: format!("{}", self.feedback_threshold),
                requirement: "must be positive and finite",
            });
        }
        if self.cone_feedback && self.feedback_max_rounds == 0 {
            return Err(InputError::BadOption {
                name: "feedback_max_rounds",
                value: format!("{}", self.feedback_max_rounds),
                requirement: "must be at least 1",
            });
        }
        Ok(())
    }
}

/// A plain value in a [`Report`]. Becomes `bool`, `int`, `float`, `str` or `list[str]` in Python.
#[derive(Debug, Clone, PartialEq)]
pub enum ReportValue {
    /// A flag.
    Bool(bool),
    /// A count.
    Int(i64),
    /// A measurement.
    Float(f64),
    /// A name or a status.
    Text(String),
    /// A list of messages.
    TextList(Vec<String>),
}

impl From<bool> for ReportValue {
    fn from(v: bool) -> Self {
        ReportValue::Bool(v)
    }
}
impl From<usize> for ReportValue {
    fn from(v: usize) -> Self {
        ReportValue::Int(v as i64)
    }
}
impl From<i64> for ReportValue {
    fn from(v: i64) -> Self {
        ReportValue::Int(v)
    }
}
impl From<f64> for ReportValue {
    fn from(v: f64) -> Self {
        ReportValue::Float(v)
    }
}
impl From<&str> for ReportValue {
    fn from(v: &str) -> Self {
        ReportValue::Text(v.to_owned())
    }
}
impl From<String> for ReportValue {
    fn from(v: String) -> Self {
        ReportValue::Text(v)
    }
}

/// What happened during an unwrap: an ordered list of named plain values plus warnings.
///
/// Becomes the `report` dict in Python, in insertion order, with the warnings under the key
/// `"warnings"` (always present, possibly empty). Keys starting with `seconds_` are wall-clock
/// times in seconds.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Report {
    entries: Vec<(String, ReportValue)>,
    warnings: Vec<String>,
}

impl Report {
    /// An empty report.
    pub fn new() -> Report {
        Report::default()
    }

    /// Sets `key` to `value`. An existing key keeps its position and gets the new value.
    pub fn set(&mut self, key: &str, value: impl Into<ReportValue>) {
        let value = value.into();
        match self.entries.iter_mut().find(|(k, _)| k == key) {
            Some(entry) => entry.1 = value,
            None => self.entries.push((key.to_owned(), value)),
        }
    }

    /// The value stored under `key`.
    pub fn get(&self, key: &str) -> Option<&ReportValue> {
        self.entries.iter().find(|(k, _)| k == key).map(|(_, v)| v)
    }

    /// Adds a warning: something the caller should know that did not stop the unwrap.
    pub fn warn(&mut self, message: impl Into<String>) {
        self.warnings.push(message.into());
    }

    /// The entries in insertion order.
    pub fn entries(&self) -> &[(String, ReportValue)] {
        &self.entries
    }

    /// The warnings in the order they were raised.
    pub fn warnings(&self) -> &[String] {
        &self.warnings
    }

    /// Appends the entries and warnings of `other`; keys that exist are overwritten in place.
    pub fn extend(&mut self, other: Report) {
        for (key, value) in other.entries {
            self.set(&key, value);
        }
        self.warnings.extend(other.warnings);
    }
}

/// What a pipeline returns.
#[derive(Debug, Clone, PartialEq)]
pub struct UnwrapOutput {
    /// One UV per polygon corner, in corner ("loop") order. Not packed: every chart is in its
    /// own coordinate system, at the scale the flattening produced.
    pub uvs: Vec<[f64; 2]>,
    /// For every polygon, the chart it belongs to, `0..n_charts`. A chart is a piece of the
    /// surface that was flattened as one; different charts overlap in `uvs` until they are
    /// packed ([`crate::pack::pack_charts`]).
    pub polygon_charts: Vec<u32>,
    /// The cones that were used, sorted by vertex.
    pub cones: Vec<Cone>,
    /// The edges the surface was cut along, as vertex pairs `[min, max]`, sorted, each once.
    /// Includes the seam edges of the input that were honoured. Cuts follow polygon edges
    /// only (never a diagonal of the triangulation): with one UV per polygon corner a polygon
    /// cannot be split.
    pub cut_edges: Vec<[u32; 2]>,
    /// Status, counts, timings and warnings of the pipeline.
    pub report: Report,
}

/// Checks that an output fits its input: array lengths, index ranges and finite numbers.
///
/// This guards the binding against a faulty pipeline: Blender must never receive a UV array of
/// the wrong length or with NaNs in it.
///
/// # Errors
/// [`UnwrapError::Solve`] with a message that starts with "internal error".
pub fn check_output(input: &UnwrapInput, output: &UnwrapOutput) -> Result<(), UnwrapError> {
    let fail = |message: String| Err(UnwrapError::Solve(format!("internal error: {message}")));
    if output.uvs.len() != input.n_corners() {
        return fail(format!(
            "the pipeline returned {} UVs for {} polygon corners",
            output.uvs.len(),
            input.n_corners()
        ));
    }
    if let Some(c) = output
        .uvs
        .iter()
        .position(|uv| !(uv[0].is_finite() && uv[1].is_finite()))
    {
        return fail(format!(
            "the pipeline returned a non-finite UV ({}, {}) at polygon corner {c}",
            output.uvs[c][0], output.uvs[c][1]
        ));
    }
    if output.polygon_charts.len() != input.n_polygons() {
        return fail(format!(
            "the pipeline returned {} chart ids for {} polygons",
            output.polygon_charts.len(),
            input.n_polygons()
        ));
    }
    if let Some(cone) = output
        .cones
        .iter()
        .find(|c| c.vertex as usize >= input.n_vertices() || !c.angle.is_finite())
    {
        return fail(format!(
            "the pipeline returned an invalid cone (vertex {}, angle {}) for a mesh with {} vertices",
            cone.vertex,
            cone.angle,
            input.n_vertices()
        ));
    }
    if let Some(&[a, b]) = output
        .cut_edges
        .iter()
        .find(|e| e[0] >= e[1] || e[1] as usize >= input.n_vertices())
    {
        return fail(format!(
            "the pipeline returned an invalid cut edge ({a}, {b}) for a mesh with {} vertices",
            input.n_vertices()
        ));
    }
    Ok(())
}

#[cfg(test)]
pub(crate) mod test_meshes {
    //! Small meshes for the unit tests of this crate.
    use super::*;

    /// Fan triangulation (corner indices) of polygons given by their corner offsets.
    pub fn fan_triangles(polygon_starts: &[u32], n_corners: u32) -> Vec<u32> {
        let mut out = Vec::new();
        for (p, &start) in polygon_starts.iter().enumerate() {
            let end = polygon_starts.get(p + 1).copied().unwrap_or(n_corners);
            for i in start + 1..end - 1 {
                out.extend([start, i, i + 1]);
            }
        }
        out
    }

    /// Buffers of a polygon mesh given as positions and polygons (vertex lists), fan-triangulated.
    pub fn buffers(positions: &[[f64; 3]], polygons: &[&[u32]]) -> MeshBuffers {
        let mut corner_vertices = Vec::new();
        let mut polygon_starts = Vec::new();
        for polygon in polygons {
            polygon_starts.push(corner_vertices.len() as u32);
            corner_vertices.extend_from_slice(polygon);
        }
        let triangles = fan_triangles(&polygon_starts, corner_vertices.len() as u32);
        MeshBuffers {
            positions: positions.iter().flatten().copied().collect(),
            corner_vertices,
            polygon_starts,
            triangles,
            ..MeshBuffers::default()
        }
    }

    /// The unit cube `[0, 1]³` with 6 outward-facing quads. Vertex `i` is at
    /// `(i & 1, (i >> 1) & 1, (i >> 2) & 1)`.
    pub fn cube() -> MeshBuffers {
        let positions: Vec<[f64; 3]> = (0..8u32)
            .map(|i| {
                [
                    f64::from(i & 1),
                    f64::from((i >> 1) & 1),
                    f64::from((i >> 2) & 1),
                ]
            })
            .collect();
        let polygons: [&[u32]; 6] = [
            &[0, 2, 3, 1], // z = 0, facing -z
            &[4, 5, 7, 6], // z = 1, facing +z
            &[0, 1, 5, 4], // y = 0
            &[2, 6, 7, 3], // y = 1
            &[0, 4, 6, 2], // x = 0
            &[1, 3, 7, 5], // x = 1
        ];
        buffers(&positions, &polygons)
    }

    /// A flat `nx × ny` grid of unit quads in the plane `z = 0`, facing `+z`.
    pub fn grid(nx: u32, ny: u32) -> MeshBuffers {
        let mut positions = Vec::new();
        for j in 0..=ny {
            for i in 0..=nx {
                positions.push([f64::from(i), f64::from(j), 0.0]);
            }
        }
        let v = |i: u32, j: u32| j * (nx + 1) + i;
        let quads: Vec<[u32; 4]> = (0..ny)
            .flat_map(|j| (0..nx).map(move |i| (i, j)))
            .map(|(i, j)| [v(i, j), v(i + 1, j), v(i + 1, j + 1), v(i, j + 1)])
            .collect();
        let polygons: Vec<&[u32]> = quads.iter().map(|q| q.as_slice()).collect();
        buffers(&positions, &polygons)
    }

    /// Three triangles sharing the edge (0, 1): non-manifold.
    pub fn three_fins() -> MeshBuffers {
        let positions = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
            [0.5, -1.0, 0.0],
            [0.5, 0.0, 1.0],
        ];
        buffers(&positions, &[&[0, 1, 2], &[1, 0, 3], &[0, 1, 4]])
    }
}

#[cfg(test)]
mod tests {
    use super::test_meshes::*;
    use super::*;

    fn err(buffers: MeshBuffers) -> String {
        UnwrapInput::from_buffers(buffers).unwrap_err().to_string()
    }

    #[test]
    fn cube_is_accepted_and_described() {
        let input = UnwrapInput::from_buffers(cube()).unwrap();
        assert_eq!(input.n_vertices(), 8);
        assert_eq!(input.n_corners(), 24);
        assert_eq!(input.n_polygons(), 6);
        assert_eq!(input.n_triangles(), 12);
        assert_eq!(input.polygon_offsets(), &[0, 4, 8, 12, 16, 20, 24]);
        assert_eq!(input.polygon_corners(2), 8..12);
        assert_eq!(input.corner_polygons()[9], 2);
        assert!(input.seam_edges().is_empty());
    }

    #[test]
    fn positions_are_checked() {
        let mut b = cube();
        b.positions.pop();
        assert_eq!(
            err(b),
            "positions: 23 items is not a multiple of 3 (3 per vertex)"
        );
        let mut b = cube();
        b.positions[7] = f64::NAN;
        assert_eq!(err(b), "positions[7] = NaN is not a finite number");
    }

    #[test]
    fn corners_and_polygons_are_checked() {
        let mut b = cube();
        b.corner_vertices[5] = 8;
        assert_eq!(
            err(b),
            "corner_vertices[5] = 8 is out of range: there are 8 vertices (valid indices 0..=7)"
        );
        let mut b = cube();
        b.polygon_starts.clear();
        assert_eq!(err(b), "polygon_starts is empty: the mesh has no faces");
        let mut b = cube();
        b.polygon_starts[0] = 1;
        assert_eq!(
            err(b),
            "polygon_starts[0] = 1, but the first polygon must start at corner 0"
        );
        let mut b = cube();
        b.polygon_starts[2] = 6;
        assert_eq!(
            err(b),
            "polygon 1 has 2 corners (corners 4..6); a polygon needs at least 3"
        );
        let mut b = cube();
        b.polygon_starts[2] = 3;
        assert_eq!(
            err(b),
            "polygon_starts[1] = 4 is followed by 3: polygon start offsets must increase"
        );
        // The last polygon is measured against the number of corners.
        let mut b = cube();
        b.polygon_starts[5] = 23;
        assert_eq!(
            err(b),
            "polygon 5 has 1 corners (corners 23..24); a polygon needs at least 3"
        );
    }

    #[test]
    fn triangles_are_checked() {
        let mut b = cube();
        b.triangles.truncate(33);
        assert_eq!(
            err(b),
            "triangles: got 11 triangles, expected 12 (24 corners - 2 x 6 polygons)"
        );
        let mut b = cube();
        b.triangles[4] = 24;
        assert_eq!(
            err(b),
            "triangles[4] = 24 is out of range: there are 24 polygon corners (valid indices 0..=23)"
        );
        // Triangle 1 = (0, 2, 3); move one corner into polygon 1.
        let mut b = cube();
        b.triangles[5] = 4;
        assert_eq!(
            err(b),
            "triangles: triangle 1 = corners (0, 2, 4) spans polygons 0, 0, 1; the corners of a \
             triangle must belong to one polygon"
        );
        // Reverse a triangle.
        let mut b = cube();
        b.triangles.swap(1, 2);
        assert_eq!(
            err(b),
            "triangles: polygon 0 (4 corners) is not triangulated correctly: triangle 0 runs \
             against the winding of the polygon"
        );
        // Both triangles of polygon 0 identical: side (2, 3) uncovered, side (0, 1) doubled.
        let mut b = cube();
        b.triangles[3..6].copy_from_slice(&[0, 1, 2]);
        assert_eq!(
            err(b),
            "triangles: polygon 0 (4 corners) is not triangulated correctly: the side starting \
             at corner 0 is covered by 2 triangles instead of 1"
        );
        // Give polygon 0 a triangle of polygon 1: counts per polygon are off.
        let mut b = cube();
        b.triangles[3..6].copy_from_slice(&[4, 5, 6]);
        assert_eq!(
            err(b),
            "triangles: polygon 0 (4 corners) is not triangulated correctly: it has 1 triangles \
             instead of 2"
        );
        let mut b = cube();
        b.triangles[3..6].copy_from_slice(&[0, 0, 3]);
        assert_eq!(
            err(b),
            "triangles: polygon 0 (4 corners) is not triangulated correctly: triangle 1 uses a \
             corner twice"
        );
    }

    #[test]
    fn diagonals_must_pair_up() {
        // A hexagon: fan triangulation uses diagonals (0,2), (0,3), (0,4).
        let positions: Vec<[f64; 3]> = (0..6)
            .map(|i| {
                let a = f64::from(i) * std::f64::consts::FRAC_PI_3;
                [a.cos(), a.sin(), 0.0]
            })
            .collect();
        let hexagon = buffers(&positions, &[&[0, 1, 2, 3, 4, 5]]);
        assert!(UnwrapInput::from_buffers(hexagon.clone()).is_ok());
        // Another valid triangulation: (0,1,2), (2,3,4), (4,5,0), (0,2,4).
        let mut b = hexagon.clone();
        b.triangles = vec![0, 1, 2, 2, 3, 4, 4, 5, 0, 0, 2, 4];
        assert!(UnwrapInput::from_buffers(b).is_ok());
        // Every side is covered once, but the fourth triangle (1, 3, 5) does not use the
        // diagonals the other three leave open.
        let mut b = hexagon.clone();
        b.triangles = vec![0, 1, 2, 2, 3, 4, 4, 5, 0, 1, 3, 5];
        assert_eq!(
            err(b),
            "triangles: polygon 0 (6 corners) is not triangulated correctly: the diagonal \
             between corners 0 and 2 is not shared by exactly two triangles in opposite directions"
        );
    }

    #[test]
    fn seams_are_canonical() {
        let mut b = cube();
        b.seam_edges = vec![3, 1, 0, 1, 1, 3, 1, 0];
        let input = UnwrapInput::from_buffers(b).unwrap();
        assert_eq!(input.seam_edges(), &[[0, 1], [1, 3]]);
        let mut b = cube();
        b.seam_edges = vec![0, 1, 2];
        assert_eq!(
            err(b),
            "seam_edges: 3 items is not a multiple of 2 (2 per edge)"
        );
        let mut b = cube();
        b.seam_edges = vec![0, 1, 5, 5];
        assert_eq!(err(b), "seam_edges: edge 1 joins vertex 5 to itself");
        let mut b = cube();
        b.seam_edges = vec![0, 9];
        assert_eq!(
            err(b),
            "seam_edges[1] = 9 is out of range: there are 8 vertices (valid indices 0..=7)"
        );
    }

    #[test]
    fn steering_is_checked() {
        let mut b = cube();
        b.pinned_vertices = vec![5, 1, 5];
        b.forbidden_vertices = vec![7, 2];
        b.prescribed_vertices = vec![6, 0];
        b.prescribed_angles = vec![1.0, -0.5];
        let input = UnwrapInput::from_buffers(b).unwrap();
        assert_eq!(input.pinned_vertices(), &[1, 5]);
        assert_eq!(input.forbidden_vertices(), &[2, 7]);
        assert_eq!(
            input.prescribed_cones(),
            &[
                Cone {
                    vertex: 0,
                    angle: -0.5
                },
                Cone {
                    vertex: 6,
                    angle: 1.0
                }
            ]
        );

        let mut b = cube();
        b.pinned_vertices = vec![1, 5];
        b.forbidden_vertices = vec![5, 7];
        assert_eq!(
            err(b),
            "vertex 5 is both pinned and forbidden; a vertex can be in only one of the pinned, \
             prescribed and forbidden sets"
        );
        let mut b = cube();
        b.pinned_vertices = vec![8];
        assert_eq!(
            err(b),
            "pinned_vertices[0] = 8 is out of range: there are 8 vertices (valid indices 0..=7)"
        );
        let mut b = cube();
        b.prescribed_vertices = vec![1];
        assert_eq!(
            err(b),
            "prescribed_vertices has 1 entries but prescribed_angles has 0; they must be given \
             together and have equal length"
        );
        let mut b = cube();
        b.prescribed_vertices = vec![1, 1];
        b.prescribed_angles = vec![0.5, 0.5];
        assert_eq!(
            err(b),
            "prescribed_vertices: vertex 1 is listed more than once"
        );
        let mut b = cube();
        b.prescribed_vertices = vec![1];
        b.prescribed_angles = vec![7.0];
        assert_eq!(
            err(b),
            "prescribed_angles: cone angle 7 at vertex 1 is not below 2\u{3c0}"
        );
        let mut b = cube();
        b.prescribed_vertices = vec![1];
        b.prescribed_angles = vec![f64::INFINITY];
        assert_eq!(err(b), "prescribed_angles[0] = inf is not a finite number");
        let mut b = cube();
        b.prescribed_vertices = vec![3];
        b.prescribed_angles = vec![0.5];
        b.forbidden_vertices = vec![3];
        assert!(err(b).starts_with("vertex 3 is both prescribed and forbidden"));
    }

    #[test]
    fn options_are_checked() {
        assert!(UnwrapOptions::default().validate().is_ok());
        assert_eq!(UnwrapOptions::default().lambda, 0.1);
        assert!(UnwrapOptions::default().positive_cones_only);
        assert!(!UnwrapOptions::default().optimize_boundary);
        let bad = UnwrapOptions {
            lambda: -1.0,
            ..UnwrapOptions::default()
        };
        assert_eq!(
            bad.validate().unwrap_err().to_string(),
            "lambda_ = -1: lambda must be a finite number that is not negative"
        );
        let nan = UnwrapOptions {
            lambda: f64::NAN,
            ..UnwrapOptions::default()
        };
        assert!(nan.validate().is_err());
        assert!(UnwrapOptions::default().cone_feedback);
        for (threshold, rounds) in [(0.0, 4), (f64::NAN, 4), (f64::INFINITY, 4), (4.0, 0)] {
            let bad = UnwrapOptions {
                cone_feedback: true,
                feedback_threshold: threshold,
                feedback_max_rounds: rounds,
                ..UnwrapOptions::default()
            };
            let message = bad.validate().unwrap_err().to_string();
            assert!(message.starts_with("feedback_"), "{message}");
            // Not checked while the loop is off.
            let off = UnwrapOptions {
                cone_feedback: false,
                ..bad
            };
            assert!(off.validate().is_ok());
        }
    }

    #[test]
    fn report_keeps_order_and_overwrites_in_place() {
        let mut r = Report::new();
        r.set("a", 1usize);
        r.set("b", "x");
        r.set("a", 2.5);
        r.warn("careful");
        assert_eq!(r.entries().len(), 2);
        assert_eq!(r.entries()[0], ("a".to_owned(), ReportValue::Float(2.5)));
        assert_eq!(r.get("b"), Some(&ReportValue::Text("x".to_owned())));
        assert_eq!(r.warnings(), &["careful".to_owned()]);
        let mut other = Report::new();
        other.set("b", true);
        other.set("c", 3i64);
        other.warn("again");
        r.extend(other);
        assert_eq!(r.get("b"), Some(&ReportValue::Bool(true)));
        assert_eq!(r.entries().len(), 3);
        assert_eq!(r.warnings().len(), 2);
    }
}
