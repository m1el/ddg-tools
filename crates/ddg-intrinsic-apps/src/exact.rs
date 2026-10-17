//! The exact polyhedral geodesic distance from a set of vertices to every vertex of a triangle
//! mesh: the ground truth for the heat method of paper Sec. 5.3 (its Fig. 18 compares against
//! "the exact polyhedral distance").
//!
//! # The algorithm
//!
//! A window propagation in the manner of Mitchell, Mount and Papadimitriou (1987), with the
//! window parametrization of Surazhsky, Surazhsky, Kirsanov, Gortler and Hoppe (2005) and the
//! vertex filter of Xin and Wang (2009) instead of the intersection of overlapping windows.
//! Written for this crate from those papers; no code was taken from an implementation.
//!
//! A *window* is an interval `[b0, b1]` of a halfedge `h`, lit by straight lines from a
//! *pseudo-source* `s` that lies on the other side of the edge in the plane of `face(h)`
//! unfolded: in the frame of `h` (its tail at the origin, its head at `(ℓ, 0)`, `face(h)` above
//! the x-axis) `s = (s_x, s_y)` with `s_y < 0`, and the distance at the point `x` of the
//! interval is `σ + |(x, 0) − s|`, `σ` the distance of `s` from the source. Popping a window
//! (smallest distance first) propagates it across `face(h)`: the rays from `s` through the
//! interval reach the two other sides, which gives at most two child windows on the halfedges
//! across those sides, and, if the opposite vertex `c` is lit, the candidate `σ + |c − s|` for
//! its distance. A geodesic bends only at a vertex whose angle sum is at least `2π` (a saddle,
//! here also a flat vertex) or at a boundary vertex; such a vertex becomes a pseudo-source with
//! its distance once that distance is not improved any more, and lights every face around it.
//!
//! *Filter* (Xin and Wang's theorem, in the form used here). Let a window lie on the edge from
//! `a` to `b` (length `ℓ`) with current vertex distances `D(a)`, `D(b)`, which are lengths of
//! paths that exist. On its interval, `σ + |(x, 0) − s| − x` does not increase with `x`; so if
//! `D(a) + b1 < σ + |(b1, 0) − s|`, every point of the interval is reached by a shorter path
//! through `a` (along the edge to it), every path through the window is longer than some other
//! path, and the window carries no shortest path. The same with `b` and `b0`. Such windows are
//! dropped when popped. So is a window whose smallest distance is not below the largest
//! distance of all vertices once every vertex has been reached: it cannot improve any vertex.
//! A shortest path to a vertex is never dropped by either rule, so the result is exact up to
//! rounding and the thin windows below.
//!
//! *Thin windows.* A child narrower than `min_window_width` times its edge is not created
//! (counted in the report): it can only carry paths that pass within that distance of a vertex,
//! where the neighbouring windows meet. Endpoints of a child that rounding puts a little
//! outside its edge are moved onto the edge; the largest such move is reported.
//!
//! Distances are computed at vertices only. Time and memory grow faster than for the heat
//! method (Surazhsky et al. measure `O(n^1.5)` windows for MMP); it is meant for meshes of up to
//! some ten thousand vertices.
//!
//! # References
//!
//! J. S. B. Mitchell, D. M. Mount, C. H. Papadimitriou. *The discrete geodesic problem.* SIAM
//! Journal on Computing 16(4), 647–668, 1987. <https://doi.org/10.1137/0216045>
//!
//! V. Surazhsky, T. Surazhsky, D. Kirsanov, S. J. Gortler, H. Hoppe. *Fast exact and
//! approximate geodesics on meshes.* ACM Transactions on Graphics 24(3), 553–560, 2005.
//! <https://doi.org/10.1145/1073204.1073228>
//!
//! S.-Q. Xin, G.-J. Wang. *Improving Chen and Han's algorithm on the discrete geodesic
//! problem.* ACM Transactions on Graphics 28(4), Article 104, 2009.
//! <https://doi.org/10.1145/1559755.1559761>

use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::f64::consts::TAU;
use std::fmt;

use ddg_core::geometry::triangle_area;
use ddg_core::{Geometry, GeometryError, TriMesh};

/// Options of [`exact_distance`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ExactOptions {
    /// Child windows narrower than this fraction of their edge are not created. Default `1e-12`.
    pub min_window_width: f64,
    /// The limit on the number of windows popped from the queue; reaching it is an error
    /// ([`ExactError::WindowLimit`]). `None`: no limit.
    pub max_windows: Option<usize>,
}

impl Default for ExactOptions {
    fn default() -> Self {
        ExactOptions {
            min_window_width: 1e-12,
            max_windows: None,
        }
    }
}

/// What a run of [`exact_distance`] did.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct ExactReport {
    /// Windows created (put on the queue).
    pub n_created: usize,
    /// Windows popped and propagated across a face.
    pub n_propagated: usize,
    /// Windows popped and dropped by the vertex filter.
    pub n_filtered: usize,
    /// Windows popped and dropped because they cannot improve any vertex any more.
    pub n_beyond: usize,
    /// Children not created because they were thinner than the limit.
    pub n_thin: usize,
    /// Vertices that lit their faces as (pseudo-)sources, counting repeats after an
    /// improvement.
    pub n_pseudo_sources: usize,
    /// The largest number of entries on the queue at one time.
    pub max_queue: usize,
    /// The largest distance by which an end point of a child window was moved onto its edge.
    pub max_clamp: f64,
}

/// The distances and the report. See [`exact_distance`].
#[derive(Debug, Clone, PartialEq)]
pub struct ExactDistance {
    /// The geodesic distance of every vertex from the nearest source (`+∞` for a vertex in a
    /// component without a source).
    pub distances: Vec<f64>,
    /// Counts.
    pub report: ExactReport,
}

/// Why [`exact_distance`] failed.
#[derive(Debug, Clone, PartialEq)]
pub enum ExactError {
    /// The geometry does not belong to the mesh.
    Geometry(GeometryError),
    /// No source was given.
    NoSource,
    /// A source is not a vertex of the mesh.
    SourceOutOfRange {
        /// The source.
        source: usize,
        /// The number of vertices.
        n_vertices: usize,
    },
    /// An option is out of range.
    InvalidOptions(String),
    /// The limit on popped windows was reached.
    WindowLimit {
        /// The limit.
        limit: usize,
    },
}

impl fmt::Display for ExactError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ExactError::Geometry(e) => write!(f, "{e}"),
            ExactError::NoSource => write!(f, "no source vertex was given"),
            ExactError::SourceOutOfRange { source, n_vertices } => write!(
                f,
                "source {source} is not a vertex (the mesh has {n_vertices} vertices)"
            ),
            ExactError::InvalidOptions(m) => write!(f, "invalid options: {m}"),
            ExactError::WindowLimit { limit } => write!(
                f,
                "the exact distance stopped after {limit} windows (ExactOptions::max_windows)"
            ),
        }
    }
}

impl std::error::Error for ExactError {}

impl From<GeometryError> for ExactError {
    fn from(e: GeometryError) -> Self {
        ExactError::Geometry(e)
    }
}

#[derive(Debug, Clone, Copy)]
struct Window {
    /// The halfedge: the window propagates into its face.
    h: u32,
    b0: f64,
    b1: f64,
    /// The pseudo-source in the frame of `h`, `s[1] < 0`.
    s: [f64; 2],
    sigma: f64,
}

#[derive(Debug, Clone, Copy)]
enum Event {
    Window(Window),
    Vertex(u32),
}

/// A queue entry, ordered so that the smallest key comes first out of a `BinaryHeap`.
#[derive(Debug, Clone, Copy)]
struct Entry {
    key: f64,
    /// Insertion counter: ties are broken first in, first out (determinism).
    order: u64,
    event: Event,
}

impl PartialEq for Entry {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}
impl Eq for Entry {}
impl PartialOrd for Entry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Entry {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .key
            .total_cmp(&self.key)
            .then_with(|| other.order.cmp(&self.order))
    }
}

fn cross(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}

fn norm(a: [f64; 2]) -> f64 {
    a[0].hypot(a[1])
}

/// The apex of the triangle on `[0, ℓ_ab] × {0}` with sides `ℓ_bc`, `ℓ_ca`, above the x-axis.
fn apex(l_ab: f64, l_bc: f64, l_ca: f64) -> [f64; 2] {
    let x = (l_ab * l_ab + l_ca * l_ca - l_bc * l_bc) / (2.0 * l_ab);
    // The height from the area (Kahan's formula), not from `ℓ_ca² − x²`, which cancels in
    // needles.
    let area = triangle_area(l_ab, l_bc, l_ca).unwrap_or(0.0);
    [x, 2.0 * area / l_ab]
}

/// The smallest distance `σ + |(x, 0) − s|` over `x ∈ [b0, b1]`.
fn window_key(w: &Window) -> f64 {
    let dx = if w.s[0] < w.b0 {
        w.b0 - w.s[0]
    } else if w.s[0] > w.b1 {
        w.s[0] - w.b1
    } else {
        0.0
    };
    w.sigma + dx.hypot(w.s[1])
}

struct State<'a> {
    mesh: &'a TriMesh,
    lengths: &'a [f64],
    /// Whether a vertex may be a pseudo-source (saddle, flat or boundary).
    bends: Vec<bool>,
    distance: Vec<f64>,
    queue: BinaryHeap<Entry>,
    order: u64,
    report: ExactReport,
    min_width: f64,
    n_unreached: usize,
    max_distance: f64,
}

impl State<'_> {
    fn push(&mut self, key: f64, event: Event) {
        self.order += 1;
        self.queue.push(Entry {
            key,
            order: self.order,
            event,
        });
        self.report.max_queue = self.report.max_queue.max(self.queue.len());
    }

    fn length(&self, h: usize) -> f64 {
        self.lengths[self.mesh.edge(h)]
    }

    /// Records a path of length `d` to `v`; a vertex that may bend paths is queued as a
    /// pseudo-source when its distance improves.
    fn reach(&mut self, v: usize, d: f64) {
        if d < self.distance[v] {
            let first = self.distance[v] == f64::INFINITY;
            self.distance[v] = d;
            if first {
                self.n_unreached -= 1;
                if self.n_unreached == 0 {
                    self.max_distance = self.distance.iter().fold(0.0, |m, &x| m.max(x));
                }
            }
            if self.bends[v] {
                self.push(d, Event::Vertex(v as u32));
            }
        }
    }

    fn push_window(&mut self, w: Window) {
        let width = w.b1 - w.b0;
        if width <= self.min_width * self.length(w.h as usize) {
            self.report.n_thin += 1;
            return;
        }
        self.report.n_created += 1;
        self.push(window_key(&w), Event::Window(w));
    }

    /// Lights every face around `v` from `v` at distance `d`.
    fn emit(&mut self, v: usize, d: f64) {
        self.report.n_pseudo_sources += 1;
        let outgoing: Vec<usize> = self.mesh.outgoing_halfedges(v).collect();
        for h0 in outgoing {
            let o = self.mesh.next(h0);
            let back = self.mesh.next(o);
            self.reach(self.mesh.to_vertex(h0), d + self.length(h0));
            self.reach(self.mesh.from_vertex(back), d + self.length(back));
            if let Some(t) = self.mesh.twin(o) {
                let l = self.length(o);
                let p = apex(l, self.length(back), self.length(h0));
                // In the frame of the twin (from the head of `o` to its tail) `v` lies below.
                self.push_window(Window {
                    h: t as u32,
                    b0: 0.0,
                    b1: l,
                    s: [l - p[0], -p[1]],
                    sigma: d,
                });
            }
        }
    }

    fn clamp(&mut self, mu: f64, len: f64) -> f64 {
        let c = mu.clamp(0.0, len);
        self.report.max_clamp = self.report.max_clamp.max((c - mu).abs());
        c
    }

    fn propagate(&mut self, w: Window) {
        let h = w.h as usize;
        let (a, b) = (self.mesh.from_vertex(h), self.mesh.to_vertex(h));
        let l = self.length(h);
        let tol = 1e-12 * l.max(w.sigma);
        let d0 = (w.b0 - w.s[0]).hypot(w.s[1]);
        let d1 = (w.b1 - w.s[0]).hypot(w.s[1]);
        if self.distance[a] + w.b1 < w.sigma + d1 - tol
            || self.distance[b] + (l - w.b0) < w.sigma + d0 - tol
        {
            self.report.n_filtered += 1;
            return;
        }
        self.report.n_propagated += 1;
        let eps = 1e-12 * l;
        if w.b0 <= eps {
            self.reach(a, w.sigma + d0);
        }
        if w.b1 >= l - eps {
            self.reach(b, w.sigma + d1);
        }
        let hn = self.mesh.next(h); // b -> c
        let hp = self.mesh.next(hn); // c -> a
        let c_vertex = self.mesh.to_vertex(hn);
        let (l_bc, l_ca) = (self.length(hn), self.length(hp));
        let c = apex(l, l_bc, l_ca);
        let s = w.s;
        // Where the ray from `s` through `c` meets the x-axis.
        let xc = s[0] + (c[0] - s[0]) * (-s[1]) / (c[1] - s[1]);
        if xc >= w.b0 - eps && xc <= w.b1 + eps {
            self.reach(c_vertex, w.sigma + norm([c[0] - s[0], c[1] - s[1]]));
        }
        let ray = |x: f64| [x - s[0], -s[1]];
        // Left part: onto the side from `a` to `c`.
        if w.b0 < xc {
            if let Some(t) = self.mesh.twin(hp) {
                let x1 = w.b1.min(xc);
                let u = [c[0] / l_ca, c[1] / l_ca];
                let hit = |x: f64| cross(s, ray(x)) / cross(u, ray(x));
                let mu0 = hit(w.b0);
                let mu1 = if xc <= w.b1 { l_ca } else { hit(x1) };
                let (mu0, mu1) = (self.clamp(mu0, l_ca), self.clamp(mu1, l_ca));
                // Rotate so that `a -> c` is the x-axis; the face stays below.
                let s2 = [u[0] * s[0] + u[1] * s[1], -u[1] * s[0] + u[0] * s[1]];
                self.push_window(Window {
                    h: t as u32,
                    b0: mu0,
                    b1: mu1,
                    s: s2,
                    sigma: w.sigma,
                });
            }
        }
        // Right part: onto the side from `c` to `b`.
        if w.b1 > xc {
            if let Some(t) = self.mesh.twin(hn) {
                let x0 = w.b0.max(xc);
                let u = [(l - c[0]) / l_bc, -c[1] / l_bc];
                let sc = [s[0] - c[0], s[1] - c[1]];
                let hit = |x: f64| cross(sc, ray(x)) / cross(u, ray(x));
                let mu0 = if xc >= w.b0 { 0.0 } else { hit(x0) };
                let mu1 = hit(w.b1);
                let (mu0, mu1) = (self.clamp(mu0, l_bc), self.clamp(mu1, l_bc));
                let s2 = [u[0] * sc[0] + u[1] * sc[1], -u[1] * sc[0] + u[0] * sc[1]];
                self.push_window(Window {
                    h: t as u32,
                    b0: mu0,
                    b1: mu1,
                    s: s2,
                    sigma: w.sigma,
                });
            }
        }
    }
}

/// The exact geodesic distance on the polyhedral surface `(mesh, geometry)` from the nearest of
/// `sources` (vertex ids) to every vertex. See the [module documentation](self).
///
/// # Errors
/// [`ExactError::NoSource`], [`ExactError::SourceOutOfRange`], [`ExactError::Geometry`] when
/// the geometry is not that of the mesh, [`ExactError::InvalidOptions`], and
/// [`ExactError::WindowLimit`].
pub fn exact_distance(
    mesh: &TriMesh,
    geometry: &Geometry,
    sources: &[usize],
    options: &ExactOptions,
) -> Result<ExactDistance, ExactError> {
    geometry.check_mesh(mesh)?;
    if sources.is_empty() {
        return Err(ExactError::NoSource);
    }
    if let Some(&source) = sources.iter().find(|&&s| s >= mesh.n_vertices()) {
        return Err(ExactError::SourceOutOfRange {
            source,
            n_vertices: mesh.n_vertices(),
        });
    }
    if !(options.min_window_width.is_finite() && options.min_window_width >= 0.0) {
        return Err(ExactError::InvalidOptions(format!(
            "min_window_width must be a finite number >= 0, not {}",
            options.min_window_width
        )));
    }
    let n = mesh.n_vertices();
    let mut angle_sum = vec![0.0; n];
    for (f, face) in mesh.faces().iter().enumerate() {
        for k in 0..3 {
            angle_sum[face[k] as usize] += geometry.corner_angles()[f][k];
        }
    }
    let bends: Vec<bool> = (0..n)
        .map(|v| mesh.is_boundary_vertex(v) || angle_sum[v] >= TAU - 1e-9)
        .collect();
    let mut state = State {
        mesh,
        lengths: geometry.edge_lengths(),
        bends,
        distance: vec![f64::INFINITY; n],
        queue: BinaryHeap::new(),
        order: 0,
        report: ExactReport::default(),
        min_width: options.min_window_width,
        n_unreached: n,
        max_distance: f64::INFINITY,
    };
    for &s in sources {
        if state.distance[s] > 0.0 {
            state.distance[s] = 0.0;
            state.n_unreached -= 1;
            state.push(0.0, Event::Vertex(s as u32));
        }
    }
    let mut n_popped = 0usize;
    let mut since_max = 0usize;
    while let Some(entry) = state.queue.pop() {
        match entry.event {
            Event::Vertex(v) => {
                let v = v as usize;
                // A stale entry: the vertex was reached by a shorter path since.
                if entry.key == state.distance[v] {
                    state.emit(v, entry.key);
                }
            }
            Event::Window(w) => {
                n_popped += 1;
                if let Some(limit) = options.max_windows
                    && n_popped > limit
                {
                    return Err(ExactError::WindowLimit { limit });
                }
                if state.n_unreached == 0 {
                    since_max += 1;
                    if since_max >= 4096 {
                        since_max = 0;
                        state.max_distance = state.distance.iter().fold(0.0, |m, &x| m.max(x));
                    }
                    if entry.key >= state.max_distance {
                        state.report.n_beyond += 1;
                        continue;
                    }
                }
                state.propagate(w);
            }
        }
    }
    Ok(ExactDistance {
        distances: state.distance,
        report: state.report,
    })
}
