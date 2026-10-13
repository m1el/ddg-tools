//! Tracing queries (SPEC §4): over the input mesh and over the intrinsic triangulation, from a
//! vertex, from a point of an edge and from a point of a face.
//!
//! On flat convex shapes a traced path is a straight segment in the plane, which gives every
//! output an exact answer to compare with. On curved shapes the checks are the length of the
//! path in space and that tracing back from the end arrives at the start.

mod common;

use std::f64::consts::{PI, TAU};

use common::{Shape, dist, jittered_grid, long_torus, position, sheared_parallelogram, tube};
use ddg_core::generators::{self, SplitMix64};
use ddg_intrinsic::{
    DelaunayOptions, IntrinsicError, IntrinsicTriangulation, SurfacePoint, Trace, TraceOptions,
    TraceStop,
};

const RECORD: TraceOptions = TraceOptions {
    record_crossings: true,
    max_crossings: None,
};

fn wrap(angle: f64) -> f64 {
    let a = angle.rem_euclid(TAU);
    if a >= TAU { 0.0 } else { a }
}

fn angle_distance(a: f64, b: f64) -> f64 {
    let d = (a - b).abs() % TAU;
    d.min(TAU - d)
}

/// A triangulation of a flat shape, seen in the plane: where its points are and where the
/// frame of each point points to.
trait Planar {
    fn trace(
        &self,
        start: SurfacePoint,
        direction: f64,
        length: f64,
    ) -> Result<Trace, IntrinsicError>;
    fn n_vertices(&self) -> usize;
    fn n_edges(&self) -> usize;
    fn n_faces(&self) -> usize;
    fn face(&self, f: usize) -> [usize; 3];
    /// Tail and head of the edge, in the direction its parameter runs.
    fn edge(&self, e: usize) -> [usize; 2];
    fn vertex_position(&self, v: usize) -> [f64; 2];
    /// The world angle of direction zero at a vertex, its angle sum, and whether it is on the
    /// boundary.
    fn vertex_frame(&self, v: usize) -> (f64, f64, bool);

    fn point(&self, p: SurfacePoint) -> [f64; 2] {
        match p {
            SurfacePoint::Vertex(v) => self.vertex_position(v as usize),
            SurfacePoint::Edge { edge, t } => {
                let [a, b] = self.edge(edge as usize).map(|v| self.vertex_position(v));
                [a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])]
            }
            SurfacePoint::Face { face, bary } => {
                let q = self.face(face as usize).map(|v| self.vertex_position(v));
                [0, 1].map(|k| bary[0] * q[0][k] + bary[1] * q[1][k] + bary[2] * q[2][k])
            }
        }
    }

    /// The direction, in the frame of point `p`, that has the given angle in the plane.
    fn to_frame(&self, p: SurfacePoint, world: f64) -> f64 {
        let towards = |a: [f64; 2], b: [f64; 2]| (b[1] - a[1]).atan2(b[0] - a[0]);
        match p {
            SurfacePoint::Vertex(v) => {
                let (zero, theta, boundary) = self.vertex_frame(v as usize);
                let full = if boundary { PI } else { TAU };
                let phi = wrap(world - zero) * full / theta;
                if phi >= TAU { phi - TAU } else { phi }
            }
            SurfacePoint::Edge { edge, .. } => {
                let [a, b] = self.edge(edge as usize).map(|v| self.vertex_position(v));
                wrap(world - towards(a, b))
            }
            SurfacePoint::Face { face, .. } => {
                let q = self.face(face as usize).map(|v| self.vertex_position(v));
                wrap(world - towards(q[0], q[1]))
            }
        }
    }
}

struct Input<'a>(&'a Shape, &'a IntrinsicTriangulation);
struct Intrinsic<'a>(&'a Shape, &'a IntrinsicTriangulation);

fn xy(p: [f64; 3]) -> [f64; 2] {
    [p[0], p[1]]
}

/// The frame of an input vertex: direction zero is its first outgoing halfedge.
fn input_vertex_frame(shape: &Shape, tri: &IntrinsicTriangulation, v: usize) -> (f64, f64, bool) {
    let h = shape.mesh.outgoing_halfedge(v);
    let (a, b) = (shape.positions[v], shape.positions[shape.mesh.to_vertex(h)]);
    (
        (b[1] - a[1]).atan2(b[0] - a[0]),
        tri.input().vertex_angle_sum(v),
        shape.mesh.is_boundary_vertex(v),
    )
}

impl Planar for Input<'_> {
    fn trace(
        &self,
        start: SurfacePoint,
        direction: f64,
        length: f64,
    ) -> Result<Trace, IntrinsicError> {
        self.1.input().trace(start, direction, length, &RECORD)
    }
    fn n_vertices(&self) -> usize {
        self.0.mesh.n_vertices()
    }
    fn n_edges(&self) -> usize {
        self.0.mesh.n_edges()
    }
    fn n_faces(&self) -> usize {
        self.0.mesh.n_faces()
    }
    fn face(&self, f: usize) -> [usize; 3] {
        self.0.mesh.faces()[f].map(|v| v as usize)
    }
    fn edge(&self, e: usize) -> [usize; 2] {
        self.0.mesh.edge_vertices(e)
    }
    fn vertex_position(&self, v: usize) -> [f64; 2] {
        xy(self.0.positions[v])
    }
    fn vertex_frame(&self, v: usize) -> (f64, f64, bool) {
        input_vertex_frame(self.0, self.1, v)
    }
}

impl Planar for Intrinsic<'_> {
    fn trace(
        &self,
        start: SurfacePoint,
        direction: f64,
        length: f64,
    ) -> Result<Trace, IntrinsicError> {
        self.1.trace(start, direction, length, &RECORD)
    }
    fn n_vertices(&self) -> usize {
        self.1.n_vertices()
    }
    fn n_edges(&self) -> usize {
        self.1.n_edges()
    }
    fn n_faces(&self) -> usize {
        self.1.n_faces()
    }
    fn face(&self, f: usize) -> [usize; 3] {
        self.1.faces()[f].map(|v| v as usize)
    }
    fn edge(&self, e: usize) -> [usize; 2] {
        self.1.edge_vertices(e)
    }
    fn vertex_position(&self, v: usize) -> [f64; 2] {
        xy(self.0.positions[v])
    }
    fn vertex_frame(&self, v: usize) -> (f64, f64, bool) {
        // The frame of a vertex of the intrinsic triangulation is that of the input vertex.
        input_vertex_frame(self.0, self.1, v)
    }
}

fn random_bary(rng: &mut SplitMix64) -> [f64; 3] {
    let (mut a, mut b) = (rng.next_f64(), rng.next_f64());
    if a + b > 1.0 {
        (a, b) = (1.0 - a, 1.0 - b);
    }
    // Stay off the boundary of the face.
    let c = 1.0 - a - b;
    let shrink = |x: f64| 0.02 + 0.94 * x;
    [shrink(a), shrink(b), shrink(c)]
}

fn random_point(surface: &impl Planar, kind: u64, rng: &mut SplitMix64) -> SurfacePoint {
    match kind % 3 {
        0 => SurfacePoint::Vertex((rng.next_u64() % surface.n_vertices() as u64) as u32),
        1 => SurfacePoint::Edge {
            edge: (rng.next_u64() % surface.n_edges() as u64) as u32,
            t: rng.uniform(0.05, 0.95),
        },
        _ => SurfacePoint::Face {
            face: (rng.next_u64() % surface.n_faces() as u64) as u32,
            bary: random_bary(rng),
        },
    }
}

/// Traces from random points towards random points of a flat convex shape and compares
/// everything with the straight segment: positions to `tolerance` (the shapes have size 1),
/// barycentric coordinates and directions to `1000 · tolerance`. Returns the number of
/// crossings seen.
fn check_straight_lines(
    surface: &impl Planar,
    name: &str,
    seed: u64,
    n: usize,
    tolerance: f64,
) -> usize {
    let mut rng = SplitMix64::new(seed);
    let mut crossings = 0;
    for query in 0..n {
        let start = random_point(surface, query as u64, &mut rng);
        // The end: a point strictly inside a face, so that the segment is inside the shape.
        let target = random_point(surface, 2, &mut rng);
        let (a, b) = (surface.point(start), surface.point(target));
        let length = (b[0] - a[0]).hypot(b[1] - a[1]);
        if length < 1e-3 {
            continue;
        }
        let world = (b[1] - a[1]).atan2(b[0] - a[0]);
        let direction = surface.to_frame(start, world);
        let trace = surface
            .trace(start, direction, length)
            .unwrap_or_else(|e| panic!("{name}: query {query} from {start}: {e}"));
        let context = format!("{name}: query {query} from {start} to {target}: {trace:?}");
        assert_eq!(trace.stop, TraceStop::Completed, "{context}");
        assert_eq!(trace.length, length, "{context}");
        assert_eq!(trace.n_crossings, trace.crossings.len(), "{context}");
        assert!(trace.clamp < 1000.0 * tolerance, "{context}");
        let end = surface.point(trace.end);
        assert!(
            (end[0] - b[0]).hypot(end[1] - b[1]) < tolerance,
            "{context}: ended at {end:?}, expected {b:?}"
        );
        // The end is the target, in the same face (the target is interior to its face).
        match (trace.end, target) {
            (
                SurfacePoint::Face { face, bary },
                SurfacePoint::Face {
                    face: f,
                    bary: expected,
                },
            ) => {
                assert_eq!(face, f, "{context}");
                for k in 0..3 {
                    assert!(
                        (bary[k] - expected[k]).abs() < 1000.0 * tolerance,
                        "{context}"
                    );
                }
            }
            _ => panic!("{context}: the end is not a point of a face"),
        }
        // Every crossing lies on the segment, at its distance from the start.
        let mut previous = 0.0;
        for c in &trace.crossings {
            let p = surface.point(SurfacePoint::Edge {
                edge: c.edge,
                t: c.t,
            });
            let s = c.distance / length;
            let on_segment = [a[0] + s * (b[0] - a[0]), a[1] + s * (b[1] - a[1])];
            assert!(
                (p[0] - on_segment[0]).hypot(p[1] - on_segment[1]) < tolerance,
                "{context}: {c:?}"
            );
            assert!(
                c.distance >= previous - tolerance && c.distance <= length + tolerance,
                "{context}: {c:?}"
            );
            previous = c.distance;
        }
        crossings += trace.n_crossings;
        // The way back, in the frame of the end face.
        let back = surface.to_frame(trace.end, world + PI);
        assert!(
            angle_distance(trace.back_direction, back) < 1000.0 * tolerance,
            "{context}: expected the back direction {back}"
        );
    }
    crossings
}

#[test]
fn on_flat_convex_shapes_paths_are_straight_segments() {
    let shapes = [
        Shape::generated("grid 5x5", generators::grid_disk(5, 5).unwrap(), true),
        jittered_grid(8, 0.3, 21),
        sheared_parallelogram(6, 1.0),
        sheared_parallelogram(7, 4.0),
    ];
    for (s, shape) in shapes.iter().enumerate() {
        let mut tri = shape.triangulation();
        let n = check_straight_lines(&Input(shape, &tri), &shape.name, 40 + s as u64, 600, 1e-12);
        assert!(n > 1000, "{}: only {n} crossings", shape.name);
        // The same over the intrinsic triangulation: fresh, Delaunay, flipped at random.
        let name = format!("{} (intrinsic)", shape.name);
        check_straight_lines(&Intrinsic(shape, &tri), &name, 50 + s as u64, 300, 1e-12);
        tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        let name = format!("{} (intrinsic, Delaunay)", shape.name);
        let n = check_straight_lines(&Intrinsic(shape, &tri), &name, 80 + s as u64, 600, 1e-12);
        assert!(n > 500, "{name}: only {n} crossings");
        // Random flips make triangles as thin as a flip allows; each such flip moves the
        // surface that the lengths describe by up to 1e-10 (SPEC N5).
        let mut rng = SplitMix64::new(60 + s as u64);
        for _ in 0..4 * tri.n_edges() {
            let e = (rng.next_u64() % tri.n_edges() as u64) as usize;
            let _ = tri.flip_edge(e);
        }
        let name = format!("{} (intrinsic, random flips)", shape.name);
        let n = check_straight_lines(&Intrinsic(shape, &tri), &name, 70 + s as u64, 600, 1e-8);
        assert!(n > 500, "{name}: only {n} crossings");
    }
}

/// A path that reaches the boundary stops there.
#[test]
fn a_path_that_leaves_the_surface_stops_at_the_boundary() {
    let shape = jittered_grid(6, 0.3, 4);
    let tri = shape.triangulation();
    let surface = Input(&shape, &tri);
    let mut rng = SplitMix64::new(8);
    let mut seen = 0;
    for query in 0..300 {
        let start = SurfacePoint::Face {
            face: (rng.next_u64() % shape.mesh.n_faces() as u64) as u32,
            bary: random_bary(&mut rng),
        };
        let a = surface.point(start);
        let world = rng.uniform(0.0, TAU);
        let (dx, dy) = (world.cos(), world.sin());
        // Distance to the border of the unit square along the direction.
        let hit = [
            if dx > 0.0 {
                (1.0 - a[0]) / dx
            } else {
                -a[0] / dx
            },
            if dy > 0.0 {
                (1.0 - a[1]) / dy
            } else {
                -a[1] / dy
            },
        ];
        let to_border = hit[0].min(hit[1]);
        let trace = surface
            .trace(start, surface.to_frame(start, world), to_border + 0.5)
            .unwrap();
        let context = format!("query {query}: {trace:?}");
        assert_eq!(trace.stop, TraceStop::Boundary, "{context}");
        assert!(
            (trace.length - to_border).abs() < 1e-13,
            "{context}: {to_border}"
        );
        let SurfacePoint::Edge { edge, .. } = trace.end else {
            panic!("{context}");
        };
        assert!(shape.mesh.is_boundary_edge(edge as usize), "{context}");
        let end = surface.point(trace.end);
        let expected = [a[0] + to_border * dx, a[1] + to_border * dy];
        assert!(
            (end[0] - expected[0]).hypot(end[1] - expected[1]) < 1e-13,
            "{context}"
        );
        let back = surface.to_frame(trace.end, world + PI);
        assert!(
            angle_distance(trace.back_direction, back) < 1e-10,
            "{context}"
        );
        seen += 1;
        // The same path, shorter: it ends inside.
        let inside = surface
            .trace(start, surface.to_frame(start, world), 0.999 * to_border)
            .unwrap();
        assert_eq!(inside.stop, TraceStop::Completed);
        assert!(inside.n_crossings <= trace.n_crossings, "{context}");
    }
    assert_eq!(seen, 300);
    // At a boundary vertex a direction beyond π points out of the surface.
    let corner = shape.mesh.boundary_vertices()[3];
    let out = surface
        .trace(SurfacePoint::Vertex(corner), 1.25 * PI, 0.3)
        .unwrap();
    assert_eq!(out.stop, TraceStop::Boundary);
    assert_eq!((out.length, out.end), (0.0, SurfacePoint::Vertex(corner)));
    // A point of a boundary edge, heading out.
    let e = (0..shape.mesh.n_edges())
        .find(|&e| shape.mesh.is_boundary_edge(e))
        .unwrap();
    let on_edge = SurfacePoint::Edge {
        edge: e as u32,
        t: 0.3,
    };
    let h = shape.mesh.edge_halfedge(e);
    let forward = shape.mesh.from_vertex(h) == shape.mesh.edge_vertices(e)[0];
    // The face lies on the left of the halfedge.
    let (inwards, outwards) = if forward {
        (0.5 * PI, 1.5 * PI)
    } else {
        (1.5 * PI, 0.5 * PI)
    };
    let out = surface.trace(on_edge, outwards, 0.1).unwrap();
    assert_eq!(
        (out.stop, out.length, out.end),
        (TraceStop::Boundary, 0.0, on_edge)
    );
    let inside = surface.trace(on_edge, inwards, 0.05).unwrap();
    assert_eq!(inside.stop, TraceStop::Completed);
    assert!(matches!(inside.end, SurfacePoint::Face { .. }));
}

/// Paths along edges: from a vertex along a signpost, and from a point of an edge along it.
#[test]
fn paths_along_edges_end_on_the_edge_or_at_its_vertex() {
    for shape in [jittered_grid(5, 0.3, 2), long_torus(8, 3, 4.0, 0.2, 0.37)] {
        let mut tri = shape.triangulation();
        tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        for h in 0..tri.n_halfedges() {
            let (a, b, e) = (tri.from_vertex(h), tri.to_vertex(h), tri.edge(h));
            let length = tri.edge_length(e);
            let forward = tri.edge_halfedge(e) == h;
            let start = SurfacePoint::Vertex(a as u32);
            // Part of the way: a point of the edge.
            let part = tri
                .trace(start, tri.signpost(h), 0.25 * length, &RECORD)
                .unwrap();
            let context = format!("{}: halfedge {h}: {part:?}", shape.name);
            assert_eq!(part.stop, TraceStop::Completed, "{context}");
            assert_eq!(part.along_edge, Some(e as u32), "{context}");
            assert_eq!((part.n_crossings, part.first_face), (0, None), "{context}");
            let expected_t = if forward { 0.25 } else { 0.75 };
            match part.end {
                SurfacePoint::Edge { edge, t } => {
                    assert_eq!(edge as usize, e, "{context}");
                    assert!((t - expected_t).abs() < 1e-15, "{context}");
                }
                _ => panic!("{context}"),
            }
            // In the frame of the edge the way back points against the direction of travel.
            assert_eq!(
                part.back_direction,
                if forward { PI } else { 0.0 },
                "{context}"
            );
            // Too far: the path stops at the vertex.
            let beyond = tri
                .trace(start, tri.signpost(h), 1.5 * length, &RECORD)
                .unwrap();
            assert_eq!(beyond.stop, TraceStop::Vertex, "{context}");
            assert_eq!(beyond.end, SurfacePoint::Vertex(b as u32), "{context}");
            assert_eq!(beyond.length, length, "{context}");
            // From the point of the edge onwards, along the edge and back.
            let on_edge = part.end;
            let onwards = if forward { 0.0 } else { PI };
            let rest = tri.trace(on_edge, onwards, 0.75 * length, &RECORD).unwrap();
            assert_eq!(rest.along_edge, Some(e as u32), "{context}");
            match rest.end {
                // 0.25 + 0.75 need not be exactly the length: either the vertex or the edge's
                // very end.
                SurfacePoint::Vertex(v) => assert_eq!(v as usize, b, "{context}"),
                SurfacePoint::Edge { t, .. } => {
                    assert!(
                        (t - if forward { 1.0 } else { 0.0 }).abs() < 1e-15,
                        "{context}"
                    )
                }
                _ => panic!("{context}: {rest:?}"),
            }
            let back = tri
                .trace(on_edge, wrap(onwards + PI), length, &RECORD)
                .unwrap();
            assert_eq!(back.stop, TraceStop::Vertex, "{context}");
            assert_eq!(back.end, SurfacePoint::Vertex(a as u32), "{context}");
            assert!(
                (back.length - 0.25 * length).abs() < 1e-15 * length,
                "{context}"
            );
            assert_eq!(back.back_direction, tri.signpost(h), "{context}");
        }
    }
}

/// On curved shapes: the traced path has the requested length in space, and tracing back
/// from its end along the reported direction arrives at the start.
#[test]
fn tracing_back_arrives_at_the_start() {
    let shapes = [
        Shape::generated("icosphere 2", generators::icosphere(2).unwrap(), false),
        common::jittered_icosphere(2, 0.35, 5),
        Shape::generated(
            "torus 12x6",
            generators::torus(12, 6, 1.0, 0.4).unwrap(),
            false,
        ),
        long_torus(8, 3, 4.0, 0.2, 0.37),
        Shape::generated(
            "bumpy grid",
            generators::bumpy_grid_disk(10, 10, 0.25).unwrap(),
            false,
        ),
        tube(3, 6, 0.02, 10.0, 0.37, 1.0),
    ];
    for (s, shape) in shapes.iter().enumerate() {
        let mut tri = shape.triangulation();
        let h = shape.mean_edge_length();
        let mut rng = SplitMix64::new(300 + s as u64);
        // Over the input mesh.
        let mut completed = 0;
        for query in 0..400 {
            let start = SurfacePoint::Face {
                face: (rng.next_u64() % shape.mesh.n_faces() as u64) as u32,
                bary: random_bary(&mut rng),
            };
            let (direction, length) = (rng.uniform(0.0, TAU), rng.uniform(0.2, 6.0) * h);
            let trace = tri
                .input()
                .trace(start, direction, length, &RECORD)
                .unwrap();
            let context = format!("{}: query {query} from {start}: {trace:?}", shape.name);
            // The length of the path in space.
            let in_space = common::path_length(shape, start, &trace.crossings, trace.end);
            assert!(
                (in_space - trace.length).abs() < 1e-12 * h.max(length),
                "{context}: {in_space}"
            );
            if trace.stop != TraceStop::Completed {
                assert_eq!(trace.stop, TraceStop::Boundary, "{context}");
                continue;
            }
            completed += 1;
            let back = tri
                .input()
                .trace(trace.end, trace.back_direction, length, &RECORD)
                .unwrap();
            assert_eq!(back.n_crossings, trace.n_crossings, "{context}\n{back:?}");
            let gap = dist(position(shape, back.end), position(shape, start));
            assert!(
                gap < 1e-11 * h.max(length),
                "{context}\n{back:?}: gap {gap}"
            );
            // And the direction back at the start is the opposite of the one it left in.
            assert!(
                angle_distance(back.back_direction, direction) < 1e-9,
                "{context}\n{back:?}"
            );
        }
        assert!(completed > 100, "{}", shape.name);
        // Over the intrinsic Delaunay triangulation, which has self-edges on some of these.
        tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        let mut completed = 0;
        for query in 0..400 {
            let face = (rng.next_u64() % tri.n_faces() as u64) as u32;
            let start_bary = random_bary(&mut rng);
            let start = SurfacePoint::Face {
                face,
                bary: start_bary,
            };
            let (direction, length) = (rng.uniform(0.0, TAU), rng.uniform(0.2, 6.0) * h);
            let trace = tri.trace(start, direction, length, &RECORD).unwrap();
            let context = format!(
                "{} (intrinsic): query {query} from {start}: {trace:?}",
                shape.name
            );
            if trace.stop != TraceStop::Completed {
                continue;
            }
            completed += 1;
            let back = tri
                .trace(trace.end, trace.back_direction, length, &RECORD)
                .unwrap();
            assert_eq!(back.n_crossings, trace.n_crossings, "{context}\n{back:?}");
            match back.end {
                SurfacePoint::Face { face: f, bary } => {
                    assert_eq!(f, face, "{context}\n{back:?}");
                    // Barycentric coordinates are lengths over the size of the face.
                    let size = tri
                        .face_lengths(face as usize)
                        .iter()
                        .fold(f64::INFINITY, |a, &b| a.min(b));
                    for k in 0..3 {
                        assert!(
                            (bary[k] - start_bary[k]).abs() < 1e-10 * (length / size).max(1.0),
                            "{context}\n{back:?}"
                        );
                    }
                }
                _ => panic!("{context}\n{back:?}"),
            }
            assert!(
                angle_distance(back.back_direction, direction) < 1e-9,
                "{context}"
            );
        }
        assert!(completed > 100, "{} (intrinsic)", shape.name);
    }
}

/// From a vertex over the input mesh: a signpost of the intrinsic triangulation is a direction
/// of the input mesh, and the frame of an interior vertex turns once in 2π.
#[test]
fn directions_at_a_cone_vertex_are_normalized() {
    let shape = Shape::generated("icosphere 1", generators::icosphere(1).unwrap(), false);
    let tri = shape.triangulation();
    let input = tri.input();
    for v in 0..shape.mesh.n_vertices() {
        let theta = input.vertex_angle_sum(v);
        assert!(
            theta < TAU - 0.01,
            "a vertex of the icosphere is a cone point"
        );
        let outgoing: Vec<usize> = shape.mesh.outgoing_halfedges(v).collect();
        for (n, &h) in outgoing.iter().enumerate() {
            let length = shape.geometry.edge_lengths()[shape.mesh.edge(h)];
            let next = outgoing[(n + 1) % outgoing.len()];
            let corner = shape.geometry.corner_angles()[h / 3][h % 3];
            // The signposts of consecutive halfedges differ by the normalized corner angle.
            let step = wrap(input.signpost(next) - input.signpost(h));
            assert!((step - corner * TAU / theta).abs() < 1e-14);
            // Halfway into the corner, in normalized angle, is halfway in the triangle: the
            // path of length 0.3 ℓ ends at that point of the face.
            let direction = wrap(input.signpost(h) + 0.5 * corner * TAU / theta);
            let trace = input
                .trace(
                    SurfacePoint::Vertex(v as u32),
                    direction,
                    0.3 * length,
                    &RECORD,
                )
                .unwrap();
            assert_eq!(trace.first_face, Some((h / 3) as u32));
            assert_eq!((trace.stop, trace.n_crossings), (TraceStop::Completed, 0));
            let end = position(&shape, trace.end);
            let (a, b) = (shape.positions[v], shape.positions[shape.mesh.to_vertex(h)]);
            let c = shape.positions[shape.mesh.to_vertex(shape.mesh.next(h))];
            // In space: the point at distance 0.3 ℓ from a, on the bisector of the corner.
            assert!((dist(a, end) - 0.3 * length).abs() < 1e-14);
            let unit = |p: [f64; 3]| {
                let d = dist(a, p);
                [0, 1, 2].map(|k| (p[k] - a[k]) / d)
            };
            let (ub, uc, ue) = (unit(b), unit(c), unit(end));
            let cos = |p: [f64; 3], q: [f64; 3]| p[0] * q[0] + p[1] * q[1] + p[2] * q[2];
            assert!((cos(ub, ue) - (0.5 * corner).cos()).abs() < 1e-14);
            assert!((cos(uc, ue) - (0.5 * corner).cos()).abs() < 1e-14);
            // The way back from the end arrives at the vertex.
            let back = input
                .trace(trace.end, trace.back_direction, 0.3 * length, &RECORD)
                .unwrap();
            assert!(dist(position(&shape, back.end), a) < 1e-14, "{back:?}");
        }
    }
}

#[test]
fn invalid_queries_are_errors() {
    let shape = jittered_grid(4, 0.2, 1);
    let tri = shape.triangulation();
    let input = tri.input();
    let options = TraceOptions::default();
    let inside = SurfacePoint::Face {
        face: 3,
        bary: [0.2, 0.3, 0.5],
    };
    // Points.
    let bad_points = [
        (SurfacePoint::Vertex(1000), "the vertex does not exist"),
        (
            SurfacePoint::Edge { edge: 1000, t: 0.5 },
            "the edge does not exist",
        ),
        (
            SurfacePoint::Edge { edge: 2, t: 1.5 },
            "the edge parameter is not in [0, 1]",
        ),
        (
            SurfacePoint::Edge {
                edge: 2,
                t: f64::NAN,
            },
            "the edge parameter is not in [0, 1]",
        ),
        (
            SurfacePoint::Face {
                face: 1000,
                bary: [0.2, 0.3, 0.5],
            },
            "the face does not exist",
        ),
        (
            SurfacePoint::Face {
                face: 3,
                bary: [-0.1, 0.6, 0.5],
            },
            "a barycentric coordinate is negative or not finite",
        ),
        (
            SurfacePoint::Face {
                face: 3,
                bary: [0.2, 0.3, 0.6],
            },
            "the barycentric coordinates do not sum to one",
        ),
        (
            SurfacePoint::Edge { edge: 2, t: 0.0 },
            "a tracing query cannot start at the end of an edge; start at the vertex",
        ),
        (
            SurfacePoint::Face {
                face: 3,
                bary: [0.0, 1.0, 0.0],
            },
            "a tracing query cannot start at the corner of a face; start at the vertex",
        ),
    ];
    for (point, reason) in bad_points {
        for result in [
            input.trace(point, 1.0, 0.1, &options),
            tri.trace(point, 1.0, 0.1, &options),
        ] {
            match result.unwrap_err() {
                IntrinsicError::InvalidPoint {
                    point: p,
                    reason: r,
                } => {
                    assert_eq!(r, reason);
                    // NaN is not equal to itself; compare what is printed.
                    assert_eq!(p.to_string(), point.to_string());
                }
                other => panic!("{other}"),
            }
        }
    }
    let message = input
        .trace(SurfacePoint::Edge { edge: 2, t: 1.5 }, 1.0, 0.1, &options)
        .unwrap_err()
        .to_string();
    assert_eq!(
        message,
        "the point t = 1.5 of edge 2 is not a valid point: the edge parameter is not in [0, 1]"
    );
    // Directions and lengths.
    for direction in [-0.1, TAU, 7.0, f64::NAN, f64::INFINITY] {
        let error = input.trace(inside, direction, 0.1, &options).unwrap_err();
        assert!(
            matches!(error, IntrinsicError::InvalidDirection { .. }),
            "{error}"
        );
        assert!(error.to_string().contains("[0, 2π)"), "{error}");
    }
    for length in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let error = input.trace(inside, 1.0, length, &options).unwrap_err();
        assert!(
            matches!(error, IntrinsicError::InvalidLength { .. }),
            "{error}"
        );
    }
    // A limit on the crossings.
    let limited = TraceOptions {
        record_crossings: false,
        max_crossings: Some(2),
    };
    let far = tri.input().trace(inside, 0.3, 0.6, &RECORD).unwrap();
    assert!(far.n_crossings > 2, "{far:?}");
    match input.trace(inside, 0.3, 0.6, &limited).unwrap_err() {
        IntrinsicError::CrossingLimit {
            start,
            limit,
            travelled,
        } => {
            assert_eq!((start, limit), (inside, 2));
            assert_eq!(travelled, far.crossings[2].distance);
        }
        other => panic!("{other}"),
    }
    let enough = TraceOptions {
        record_crossings: false,
        max_crossings: Some(far.n_crossings),
    };
    let same = input.trace(inside, 0.3, 0.6, &enough).unwrap();
    assert_eq!(same.end, far.end);
    assert!(same.crossings.is_empty() && same.n_crossings == far.n_crossings);
}
