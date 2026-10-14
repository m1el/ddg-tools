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
    /// A point inside the polygon where three regions meet (an index into
    /// [`CutOnPolygons::junctions`]). Never produced for two regions.
    Junction(u32),
}

/// A junction of three regions inside a triangle of a polygon's triangulation (a *triple
/// point* of the cut network). Three paths end here, one through each side of the triangle.
/// [`CutOnPolygons::junctions`] is empty for two regions.
#[derive(Debug, Clone, PartialEq)]
pub struct JunctionPoint {
    /// The position in space: `Σ barycentric[k] · x[corners[k]]`.
    pub position: [f64; 3],
    /// The polygon it lies in.
    pub polygon: u32,
    /// The triangle of the supplied triangulation it lies in (an index into `triangles`).
    pub triangle: u32,
    /// The corners of that triangle (as in `triangles[triangle]`).
    pub corners: [u32; 3],
    /// Its barycentric coordinates with respect to `corners`: positive, summing to 1. They are
    /// the weights an editor gives the three corner positions to place the vertex.
    pub barycentric: [f64; 3],
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

    /// The cut as lines in space, for display: the paths joined at the points where exactly
    /// two of them meet (an edge point shared by two polygons, a bend), so a line runs from a
    /// point where the cut ends or branches to the next one, or is closed. Each line carries
    /// the two regions it separates. A closed line does not repeat its first position.
    pub fn lines(&self) -> Vec<CutLine> {
        let position = |end: PathEnd| match end {
            PathEnd::EdgePoint(i) => self.edge_points[i as usize].position,
            PathEnd::Junction(i) => self.junctions[i as usize].position,
        };
        let key = |end: PathEnd| match end {
            PathEnd::EdgePoint(i) => (0u8, i),
            PathEnd::Junction(i) => (1u8, i),
        };
        let chain = |p: &PolygonPath| -> Vec<[f64; 3]> {
            let mut v = vec![position(p.start)];
            v.extend(p.interior.iter().map(|q| q.position));
            v.push(position(p.end));
            v
        };
        // For every node, the path ends that meet there: (path, at the end).
        let mut at: std::collections::BTreeMap<(u8, u32), Vec<(usize, bool)>> = Default::default();
        for (i, p) in self.paths.iter().enumerate() {
            at.entry(key(p.start)).or_default().push((i, false));
            at.entry(key(p.end)).or_default().push((i, true));
        }
        let mut used = vec![false; self.paths.len()];
        let mut out = Vec::new();
        // Walk from path `first`, entering it at its start (`forward`) or its end.
        let walk = |first: usize, forward: bool, used: &mut Vec<bool>| -> (Vec<[f64; 3]>, bool) {
            let mut line: Vec<[f64; 3]> = Vec::new();
            let (mut current, mut forward) = (first, forward);
            let first_node = key(if forward {
                self.paths[first].start
            } else {
                self.paths[first].end
            });
            let mut closed = false;
            loop {
                used[current] = true;
                let mut c = chain(&self.paths[current]);
                if !forward {
                    c.reverse();
                }
                if line.is_empty() {
                    line.extend(c);
                } else {
                    line.extend(c.into_iter().skip(1));
                }
                let exit = key(if forward {
                    self.paths[current].end
                } else {
                    self.paths[current].start
                });
                let list = &at[&exit];
                if list.len() != 2 {
                    break;
                }
                let Some(&(next, at_end)) = list.iter().find(|&&(p, _)| p != current) else {
                    break;
                };
                if next == first && exit == first_node {
                    closed = true;
                    break;
                }
                if used[next] {
                    break;
                }
                current = next;
                forward = !at_end;
            }
            if closed {
                line.pop();
            }
            (line, closed)
        };
        for i in 0..self.paths.len() {
            if used[i] {
                continue;
            }
            // Start at a path end that is not a pass-through point.
            let start_forward = if at[&key(self.paths[i].start)].len() != 2 {
                Some(true)
            } else if at[&key(self.paths[i].end)].len() != 2 {
                Some(false)
            } else {
                None
            };
            if let Some(f) = start_forward {
                let (positions, closed) = walk(i, f, &mut used);
                out.push(CutLine {
                    positions,
                    closed,
                    regions: self.paths[i].regions,
                });
            }
        }
        for i in 0..self.paths.len() {
            if !used[i] {
                let (positions, closed) = walk(i, true, &mut used);
                out.push(CutLine {
                    positions,
                    closed,
                    regions: self.paths[i].regions,
                });
            }
        }
        out
    }
}

/// A line of the cut for display, from [`CutOnPolygons::lines`].
#[derive(Debug, Clone, PartialEq)]
pub struct CutLine {
    /// The positions, in order.
    pub positions: Vec<[f64; 3]>,
    /// Whether the last position is joined to the first.
    pub closed: bool,
    /// The two regions the line separates, `[lower id, higher id]`.
    pub regions: [u32; 2],
}
