//! The cut expressed on the caller's polygons: what an editor needs to insert it.
//!
//! The flow runs on a triangle mesh; the caller owns polygons. A cut found by the flow is a set
//! of points on the edges of that triangle mesh and straight segments across its triangles.
//! [`CutOnPolygons`] says the same thing in the caller's terms:
//!
//! * **edge points** ([`EdgePoint`]): where the cut crosses an edge of a polygon, as the vertex
//!   pair of the edge, the parameter from the first vertex and the position. On an edge that is
//!   a seam the two sides were cut apart before the flow ran, so they carry separate points
//!   ([`EdgePoint::side`]);
//! * **paths** ([`PolygonPath`]): per polygon, the pieces of the cut inside it. A path goes from
//!   one edge point to another and passes through the points where it crosses a diagonal of the
//!   polygon's triangulation ([`InteriorPoint`]). Two regions never branch. The ends are
//!   [`PathEnd`]s, which can also say "an interior junction point" for the day the flow
//!   handles three or more regions.
//!
//! [`crate::cutting::insert`] turns a polygon mesh and a [`CutOnPolygons`] into a new polygon
//! mesh. Everything here is plain data with public fields, so that a binding can copy it to
//! arrays.

/// Where the cut crosses an edge of a polygon.
#[derive(Debug, Clone, PartialEq)]
pub struct EdgePoint {
    /// The polygon edge as the pair of its vertices, `edge[0] < edge[1]`.
    pub edge: [u32; 2],
    /// The parameter from `edge[0]`: the point is `(1 − t)·x[edge[0]] + t·x[edge[1]]`, with
    /// `0 < t < 1`.
    pub t: f64,
    /// The position in space.
    pub position: [f64; 3],
    /// For an edge that is a seam (the two polygons were cut apart before the flow ran), the
    /// polygon whose side of the seam this point is on. `None` for every other edge: the point
    /// is on the edge itself and is seen from both of its polygons.
    pub side: Option<u32>,
}

/// A point where a path crosses a diagonal of the triangulation of its polygon.
#[derive(Debug, Clone, PartialEq)]
pub struct InteriorPoint {
    /// The diagonal as the pair of **polygon corners** (indices into `corner_vertices`) it
    /// joins, `corners[0] < corners[1]`.
    pub corners: [u32; 2],
    /// The parameter from `corners[0]`, `0 < t < 1`.
    pub t: f64,
    /// The position in space.
    pub position: [f64; 3],
}

/// An end of a path.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PathEnd {
    /// An edge point (an index into [`CutOnPolygons::edge_points`]).
    EdgePoint(u32),
    /// A point inside the polygon where three or more regions meet (an index into
    /// [`CutOnPolygons::junctions`]). Never produced for two regions;
    /// [`crate::cutting::insert`] refuses it for now.
    Junction(u32),
}

/// A junction of three or more regions inside a polygon. Reserved for the flow with `N`
/// regions; [`CutOnPolygons::junctions`] is empty for two.
#[derive(Debug, Clone, PartialEq)]
pub struct JunctionPoint {
    /// The position in space.
    pub position: [f64; 3],
}

/// One piece of the cut inside one polygon.
///
/// The path runs from [`start`](PolygonPath::start) through [`interior`](PolygonPath::interior)
/// to [`end`](PolygonPath::end). It crosses one triangle of the polygon's triangulation between
/// two consecutive points. Its direction is that of the flow's cut: the region with the lower id
/// ([`PolygonPath::regions`]`[0]`) is on the left, seen from the side the polygon's normal
/// points to.
#[derive(Debug, Clone, PartialEq)]
pub struct PolygonPath {
    /// The polygon.
    pub polygon: u32,
    /// Where the path starts.
    pub start: PathEnd,
    /// The diagonals it crosses, in order.
    pub interior: Vec<InteriorPoint>,
    /// Where the path ends.
    pub end: PathEnd,
    /// The regions on its left and on its right.
    pub regions: [u32; 2],
}

/// The cut expressed on the caller's polygons. See the module documentation of `describe` (above) for the meaning.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct CutOnPolygons {
    /// The points where the cut crosses polygon edges, each once.
    pub edge_points: Vec<EdgePoint>,
    /// The paths, ordered by polygon (ascending) and, inside a polygon, by their first edge
    /// point.
    pub paths: Vec<PolygonPath>,
    /// Junction points (none for two regions).
    pub junctions: Vec<JunctionPoint>,
    /// Number of crossings the flow moved away from a vertex into `[δ, 1 − δ]`.
    pub n_clamped_crossings: usize,
}

impl CutOnPolygons {
    /// Whether the cut is empty.
    pub fn is_empty(&self) -> bool {
        self.paths.is_empty()
    }

    /// The paths of polygon `p`, as a slice of [`CutOnPolygons::paths`].
    pub fn polygon_paths(&self, p: u32) -> &[PolygonPath] {
        let lo = self.paths.partition_point(|x| x.polygon < p);
        let hi = self.paths.partition_point(|x| x.polygon <= p);
        &self.paths[lo..hi]
    }
}
