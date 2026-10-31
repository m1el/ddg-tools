//! Pins the exact results of the signpost representation: the state after flips, insertions,
//! edge splits, Delaunay flipping, refinement, the optimal Delaunay iteration, vertex moves and
//! the simplicial export, as hashes of every stored bit (connectivity, lengths, signposts,
//! angle sums, places on the input mesh, the fingerprint of SPEC N15), and a few edge lengths
//! as bit patterns.
//!
//! The values were recorded on the code before the connectivity was factored out into
//! `complex::Complex` (SPEC N50); the refactoring had to leave every one of them unchanged.
//! To print the table again: `DDG_PIN_PRINT=1 cargo test --release -p ddg-intrinsic --test
//! bitwise_pins -- --nocapture`.

mod common;

use common::{
    Shape, inner_bary, jittered_grid, jittered_icosphere, long_torus, needle_tetrahedron,
};
use common::{load_data, needles_and_caps, sheared_parallelogram};
use ddg_core::generators::{self, SplitMix64};
use ddg_intrinsic::{
    DelaunayOptions, EdgeTraceOptions, ExportOptions, IntrinsicTriangulation,
    OptimalDelaunayOptions, RefineOptions, SurfacePoint,
};

struct Fnv(u64);

impl Fnv {
    fn new() -> Fnv {
        Fnv(0xcbf2_9ce4_8422_2325)
    }
    fn w(&mut self, x: u64) {
        for byte in x.to_le_bytes() {
            self.0 ^= u64::from(byte);
            self.0 = self.0.wrapping_mul(0x0000_0100_0000_01b3);
        }
    }
    fn point(&mut self, p: SurfacePoint) {
        match p {
            SurfacePoint::Vertex(v) => {
                self.w(0);
                self.w(u64::from(v));
            }
            SurfacePoint::Edge { edge, t } => {
                self.w(1);
                self.w(u64::from(edge));
                self.w(t.to_bits());
            }
            SurfacePoint::Face { face, bary } => {
                self.w(2);
                self.w(u64::from(face));
                for b in bary {
                    self.w(b.to_bits());
                }
            }
        }
    }
}

/// Every stored bit of a triangulation, through its public accessors.
fn state_hash(tri: &IntrinsicTriangulation) -> u64 {
    let mut h = Fnv::new();
    h.w(tri.fingerprint());
    h.w(tri.n_vertices() as u64);
    h.w(tri.n_edges() as u64);
    h.w(tri.n_faces() as u64);
    for face in tri.faces() {
        for &v in face {
            h.w(u64::from(v));
        }
    }
    for g in 0..tri.n_halfedges() {
        h.w(tri.twin(g).map_or(u64::MAX, |t| t as u64));
        h.w(tri.edge(g) as u64);
        h.w(tri.signpost(g).to_bits());
    }
    for e in 0..tri.n_edges() {
        h.w(tri.edge_halfedge(e) as u64);
        h.w(tri.edge_length(e).to_bits());
        h.w(tri.edge_on_input(e).map_or(u64::MAX, |i| i as u64));
    }
    for v in 0..tri.n_vertices() {
        h.w(tri.outgoing_halfedge(v) as u64);
        h.w(tri.vertex_angle_sum(v).to_bits());
        h.point(tri.vertex_location(v));
    }
    h.0
}

fn shapes() -> Vec<Shape> {
    vec![
        jittered_icosphere(2, 0.35, 5),
        Shape::generated(
            "torus 5x3",
            generators::torus(5, 3, 1.0, 0.25).unwrap(),
            false,
        ),
        long_torus(8, 3, 4.0, 0.2, 0.37),
        needle_tetrahedron(8.0, 4.0),
        jittered_grid(10, 0.3, 3),
        Shape::generated(
            "annulus 16x3",
            generators::annulus(16, 3, 0.3, 1.0).unwrap(),
            true,
        ),
        sheared_parallelogram(8, 3.0),
        needles_and_caps(6),
    ]
}

fn record(shape: &Shape, out: &mut Vec<(String, u64)>) {
    let mut put = |what: &str, value: u64| out.push((format!("{}/{what}", shape.name), value));
    let mut rng = SplitMix64::new(0x5eed);
    let mut tri = shape.triangulation();
    put("new", state_hash(&tri));

    // Random flips, some refused.
    let mut n_flips = 0u64;
    for _ in 0..3 * tri.n_edges() {
        let e = (rng.next_u64() % tri.n_edges() as u64) as usize;
        n_flips += u64::from(tri.flip_edge(e).is_ok());
    }
    put("flips/count", n_flips);
    put("flips", state_hash(&tri));

    // Insertions in faces and splits of edges (boundary edges included).
    let mut n_inserted = 0u64;
    for _ in 0..12 {
        let f = (rng.next_u64() % tri.n_faces() as u64) as usize;
        n_inserted += u64::from(tri.insert_vertex(f, inner_bary(&mut rng)).is_ok());
    }
    put("insertions/count", n_inserted);
    put("insertions", state_hash(&tri));
    let mut n_split = 0u64;
    for _ in 0..12 {
        let e = (rng.next_u64() % tri.n_edges() as u64) as usize;
        n_split += u64::from(tri.split_edge(e, rng.uniform(0.2, 0.8)).is_ok());
    }
    put("splits/count", n_split);
    put("splits", state_hash(&tri));

    // Delaunay flipping.
    let report = tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
    put("delaunay/n_flips", report.n_flips as u64);
    put("delaunay", state_hash(&tri));
    for e in 0..3 {
        put(&format!("delaunay/length{e}"), tri.edge_length(e).to_bits());
    }
    let last = tri.n_edges() - 1;
    put("delaunay/length_last", tri.edge_length(last).to_bits());

    // Edge traces and the overlay.
    match tri.trace_edges(&EdgeTraceOptions::default()) {
        Ok(traces) => {
            let mut h = Fnv::new();
            for e in 0..traces.n_edges() {
                for c in traces.crossings(e) {
                    h.w(c.halfedge as u64);
                    h.w(c.edge as u64);
                    h.w(c.t.to_bits());
                    h.w(c.distance.to_bits());
                }
                h.w(traces.mismatch(e).to_bits());
            }
            put("traces", h.0);
            match tri.overlay(&traces) {
                Ok(overlay) => {
                    let mut h = Fnv::new();
                    for v in &overlay.vertices {
                        h.point(v.on_input);
                        h.point(v.on_triangulation);
                    }
                    for &c in &overlay.corners {
                        h.w(u64::from(c));
                    }
                    put("overlay", h.0);
                }
                Err(_) => put("overlay/error", 1),
            }
        }
        Err(_) => put("traces/error", 1),
    }

    // Moves of the inserted vertices towards their targets.
    let mut n_moved = 0u64;
    for v in tri.n_input_vertices()..tri.n_vertices() {
        if tri.is_boundary_vertex(v) {
            continue;
        }
        if let Ok(goal) = tri.odt_target(v, true) {
            n_moved += u64::from(tri.move_vertex(v, goal).is_ok());
        }
    }
    put("moves/count", n_moved);
    put("moves", state_hash(&tri));

    // Refinement, the optimal Delaunay iteration and the export, from the input.
    let mut refined = shape.triangulation();
    let report = refined.refine(&RefineOptions::default()).unwrap();
    put("refine/n_inserted", report.n_inserted as u64);
    put("refine", state_hash(&refined));
    let export = refined
        .export_simplicial(&ExportOptions::default())
        .unwrap();
    put("export/n_splits", export.n_splits as u64);
    put("export", export.fingerprint);
    put("export/state", state_hash(&refined));
    let mut odt = shape.triangulation();
    let report = odt
        .optimal_delaunay(&OptimalDelaunayOptions::default())
        .unwrap();
    put("odt/n_inserted", report.n_inserted as u64);
    put("odt", state_hash(&odt));
}

/// Recorded on the code before the refactoring (see the module documentation).
const PINS: &[(&str, u64)] = &[
    ("jittered icosphere 2 amp 0.35/new", 0xb436ceb6f37ef9e2),
    (
        "jittered icosphere 2 amp 0.35/flips/count",
        0x0000000000000348,
    ),
    ("jittered icosphere 2 amp 0.35/flips", 0x563e911bc5457de4),
    (
        "jittered icosphere 2 amp 0.35/insertions/count",
        0x000000000000000c,
    ),
    (
        "jittered icosphere 2 amp 0.35/insertions",
        0x109f3cb763bcb120,
    ),
    (
        "jittered icosphere 2 amp 0.35/splits/count",
        0x000000000000000c,
    ),
    ("jittered icosphere 2 amp 0.35/splits", 0xa79da18157c95a17),
    (
        "jittered icosphere 2 amp 0.35/delaunay/n_flips",
        0x000000000000011d,
    ),
    ("jittered icosphere 2 amp 0.35/delaunay", 0x1ed79e76e78a224b),
    (
        "jittered icosphere 2 amp 0.35/delaunay/length0",
        0x3fd3afb25661d00b,
    ),
    (
        "jittered icosphere 2 amp 0.35/delaunay/length1",
        0x3fd0462246ecd8df,
    ),
    (
        "jittered icosphere 2 amp 0.35/delaunay/length2",
        0x3fd7a29d53113a72,
    ),
    (
        "jittered icosphere 2 amp 0.35/delaunay/length_last",
        0x3fc75a90f096bb46,
    ),
    ("jittered icosphere 2 amp 0.35/traces", 0x627c8637b742532b),
    ("jittered icosphere 2 amp 0.35/overlay", 0x56910f662413903b),
    (
        "jittered icosphere 2 amp 0.35/moves/count",
        0x0000000000000018,
    ),
    ("jittered icosphere 2 amp 0.35/moves", 0x6771a2c8b4d4bebb),
    (
        "jittered icosphere 2 amp 0.35/refine/n_inserted",
        0x0000000000000020,
    ),
    ("jittered icosphere 2 amp 0.35/refine", 0xd4cfb8fce0afc0ea),
    (
        "jittered icosphere 2 amp 0.35/export/n_splits",
        0x0000000000000000,
    ),
    ("jittered icosphere 2 amp 0.35/export", 0x8be0f1c009aa84eb),
    (
        "jittered icosphere 2 amp 0.35/export/state",
        0xd4cfb8fce0afc0ea,
    ),
    (
        "jittered icosphere 2 amp 0.35/odt/n_inserted",
        0x00000000000000f6,
    ),
    ("jittered icosphere 2 amp 0.35/odt", 0x6f6082ea29d82095),
    ("torus 5x3/new", 0x2c505c0853207684),
    ("torus 5x3/flips/count", 0x000000000000004e),
    ("torus 5x3/flips", 0xd5793ff7a2e076f3),
    ("torus 5x3/insertions/count", 0x000000000000000c),
    ("torus 5x3/insertions", 0xbc465f7d7f94450f),
    ("torus 5x3/splits/count", 0x000000000000000c),
    ("torus 5x3/splits", 0x76561f99e737f5a4),
    ("torus 5x3/delaunay/n_flips", 0x000000000000005a),
    ("torus 5x3/delaunay", 0x64001a87b5f5c1f7),
    ("torus 5x3/delaunay/length0", 0x3fe7da892f84fb92),
    ("torus 5x3/delaunay/length1", 0x3fdbb67ae8584cac),
    ("torus 5x3/delaunay/length2", 0x3fe7da892f84fb7f),
    ("torus 5x3/delaunay/length_last", 0x3fbac075a11f42d3),
    ("torus 5x3/traces", 0x196840e75eea76bc),
    ("torus 5x3/overlay", 0x853e7f230293cda1),
    ("torus 5x3/moves/count", 0x0000000000000015),
    ("torus 5x3/moves", 0x15075fb902e23344),
    ("torus 5x3/refine/n_inserted", 0x0000000000000005),
    ("torus 5x3/refine", 0x35306e4e833dbd7f),
    ("torus 5x3/export/n_splits", 0x000000000000001e),
    ("torus 5x3/export", 0xaad20438859b1632),
    ("torus 5x3/export/state", 0x92e98afb5b06572c),
    ("torus 5x3/odt/n_inserted", 0x0000000000000014),
    ("torus 5x3/odt", 0x0d096be29eb93eae),
    (
        "long torus 8x3 R 4 r 0.2 twist 0.37/new",
        0x8f687500c68d8922,
    ),
    (
        "long torus 8x3 R 4 r 0.2 twist 0.37/flips/count",
        0x0000000000000062,
    ),
    (
        "long torus 8x3 R 4 r 0.2 twist 0.37/flips",
        0x08fd64ee892db3d6,
    ),
    (
        "long torus 8x3 R 4 r 0.2 twist 0.37/insertions/count",
        0x000000000000000c,
    ),
    (
        "long torus 8x3 R 4 r 0.2 twist 0.37/insertions",
        0xdb4facf39cdeb163,
    ),
    (
        "long torus 8x3 R 4 r 0.2 twist 0.37/splits/count",
        0x000000000000000c,
    ),
    (
        "long torus 8x3 R 4 r 0.2 twist 0.37/splits",
        0xa5b8460a43b5b042,
    ),
    (
        "long torus 8x3 R 4 r 0.2 twist 0.37/delaunay/n_flips",
        0x0000000000000089,
    ),
    (
        "long torus 8x3 R 4 r 0.2 twist 0.37/delaunay",
        0xb3f35811e94c0d2c,
    ),
    (
        "long torus 8x3 R 4 r 0.2 twist 0.37/delaunay/length0",
        0x3fd62b9586ad0a27,
    ),
    (
        "long torus 8x3 R 4 r 0.2 twist 0.37/delaunay/length1",
        0x3fd62b9586ad0a26,
    ),
    (
        "long torus 8x3 R 4 r 0.2 twist 0.37/delaunay/length2",
        0x3fe5834a5ca5104e,
    ),
    (
        "long torus 8x3 R 4 r 0.2 twist 0.37/delaunay/length_last",
        0x3fee06b76d68ed7a,
    ),
    (
        "long torus 8x3 R 4 r 0.2 twist 0.37/traces",
        0x8dd5e6a84aa8e7ad,
    ),
    (
        "long torus 8x3 R 4 r 0.2 twist 0.37/overlay",
        0xf3989eb10d26a788,
    ),
    (
        "long torus 8x3 R 4 r 0.2 twist 0.37/moves/count",
        0x000000000000000a,
    ),
    (
        "long torus 8x3 R 4 r 0.2 twist 0.37/moves",
        0x1fd3db69604e2d43,
    ),
    (
        "long torus 8x3 R 4 r 0.2 twist 0.37/refine/n_inserted",
        0x0000000000000036,
    ),
    (
        "long torus 8x3 R 4 r 0.2 twist 0.37/refine",
        0xb28c1109e3c951df,
    ),
    (
        "long torus 8x3 R 4 r 0.2 twist 0.37/export/n_splits",
        0x000000000000005f,
    ),
    (
        "long torus 8x3 R 4 r 0.2 twist 0.37/export",
        0x729d1f9cc881a8e4,
    ),
    (
        "long torus 8x3 R 4 r 0.2 twist 0.37/export/state",
        0xea38a9d10a08bf7b,
    ),
    (
        "long torus 8x3 R 4 r 0.2 twist 0.37/odt/n_inserted",
        0x000000000000001a,
    ),
    (
        "long torus 8x3 R 4 r 0.2 twist 0.37/odt",
        0x6658867081b75510,
    ),
    ("needle tetrahedron h 8 lift 4/new", 0x0a5f6ea343cf371b),
    (
        "needle tetrahedron h 8 lift 4/flips/count",
        0x000000000000000d,
    ),
    ("needle tetrahedron h 8 lift 4/flips", 0x6b60fba71565b64a),
    (
        "needle tetrahedron h 8 lift 4/insertions/count",
        0x000000000000000c,
    ),
    (
        "needle tetrahedron h 8 lift 4/insertions",
        0xcd14c6e93e997689,
    ),
    (
        "needle tetrahedron h 8 lift 4/splits/count",
        0x000000000000000c,
    ),
    ("needle tetrahedron h 8 lift 4/splits", 0xf8f8d1e25d681f94),
    (
        "needle tetrahedron h 8 lift 4/delaunay/n_flips",
        0x000000000000005c,
    ),
    ("needle tetrahedron h 8 lift 4/delaunay", 0xbf76231b4053fbee),
    (
        "needle tetrahedron h 8 lift 4/delaunay/length0",
        0x40084bdc3efe8486,
    ),
    (
        "needle tetrahedron h 8 lift 4/delaunay/length1",
        0x3ffd611cab6c82f8,
    ),
    (
        "needle tetrahedron h 8 lift 4/delaunay/length2",
        0x4001e7392602a608,
    ),
    (
        "needle tetrahedron h 8 lift 4/delaunay/length_last",
        0x3fdf49c0353baa4c,
    ),
    ("needle tetrahedron h 8 lift 4/traces", 0xf83452366cfcf3e2),
    ("needle tetrahedron h 8 lift 4/overlay", 0xcc00e512f28f6c45),
    (
        "needle tetrahedron h 8 lift 4/moves/count",
        0x0000000000000008,
    ),
    ("needle tetrahedron h 8 lift 4/moves", 0xb99c3300c40866c7),
    (
        "needle tetrahedron h 8 lift 4/refine/n_inserted",
        0x0000000000000001,
    ),
    ("needle tetrahedron h 8 lift 4/refine", 0xd3f6a03b5b187339),
    (
        "needle tetrahedron h 8 lift 4/export/n_splits",
        0x0000000000000008,
    ),
    ("needle tetrahedron h 8 lift 4/export", 0xed03bfde87c85eeb),
    (
        "needle tetrahedron h 8 lift 4/export/state",
        0x67131cbed6c021fb,
    ),
    (
        "needle tetrahedron h 8 lift 4/odt/n_inserted",
        0x0000000000000000,
    ),
    ("needle tetrahedron h 8 lift 4/odt", 0x4b45c4b04406b126),
    ("jittered grid 10 amp 0.3/new", 0xfef04434764a4b53),
    ("jittered grid 10 amp 0.3/flips/count", 0x00000000000001ed),
    ("jittered grid 10 amp 0.3/flips", 0xf7512815bc764765),
    (
        "jittered grid 10 amp 0.3/insertions/count",
        0x000000000000000c,
    ),
    ("jittered grid 10 amp 0.3/insertions", 0x21ce75b52702d76f),
    ("jittered grid 10 amp 0.3/splits/count", 0x000000000000000c),
    ("jittered grid 10 amp 0.3/splits", 0xff43c8c750d0b055),
    (
        "jittered grid 10 amp 0.3/delaunay/n_flips",
        0x00000000000000cc,
    ),
    ("jittered grid 10 amp 0.3/delaunay", 0x6afb35c9bbde247d),
    (
        "jittered grid 10 amp 0.3/delaunay/length0",
        0x3fb999999999999a,
    ),
    (
        "jittered grid 10 amp 0.3/delaunay/length1",
        0x3fb999999999999a,
    ),
    (
        "jittered grid 10 amp 0.3/delaunay/length2",
        0x3fc16292a41ae2cf,
    ),
    (
        "jittered grid 10 amp 0.3/delaunay/length_last",
        0x3fb773218554e314,
    ),
    ("jittered grid 10 amp 0.3/traces", 0x05ea1ba0aca64489),
    ("jittered grid 10 amp 0.3/overlay", 0x1409f84ee31b3a39),
    ("jittered grid 10 amp 0.3/moves/count", 0x0000000000000014),
    ("jittered grid 10 amp 0.3/moves", 0xe48984ab3110ad74),
    (
        "jittered grid 10 amp 0.3/refine/n_inserted",
        0x0000000000000005,
    ),
    ("jittered grid 10 amp 0.3/refine", 0xf55969c76bda94c4),
    (
        "jittered grid 10 amp 0.3/export/n_splits",
        0x0000000000000000,
    ),
    ("jittered grid 10 amp 0.3/export", 0x5fbe812ccf8634c4),
    ("jittered grid 10 amp 0.3/export/state", 0xf55969c76bda94c4),
    (
        "jittered grid 10 amp 0.3/odt/n_inserted",
        0x0000000000000088,
    ),
    ("jittered grid 10 amp 0.3/odt", 0xb55f999f2b85a5aa),
    ("annulus 16x3/new", 0x5168c7f68ee91cb0),
    ("annulus 16x3/flips/count", 0x00000000000000e6),
    ("annulus 16x3/flips", 0x002f10f2426f263e),
    ("annulus 16x3/insertions/count", 0x000000000000000c),
    ("annulus 16x3/insertions", 0xc4e92690b717600e),
    ("annulus 16x3/splits/count", 0x000000000000000c),
    ("annulus 16x3/splits", 0x5e00a1ca941b7ebf),
    ("annulus 16x3/delaunay/n_flips", 0x0000000000000060),
    ("annulus 16x3/delaunay", 0xdf03e81d58cc88fa),
    ("annulus 16x3/delaunay/length0", 0x3fcdddddddddddde),
    ("annulus 16x3/delaunay/length1", 0x3fd8f8b83c69a60b),
    ("annulus 16x3/delaunay/length2", 0x3fda7a5a9c8a0c86),
    ("annulus 16x3/delaunay/length_last", 0x3fc8f91a6e06ce3d),
    ("annulus 16x3/traces", 0x75374cf32bfb8ca8),
    ("annulus 16x3/overlay", 0xbd3d3c5cb859b886),
    ("annulus 16x3/moves/count", 0x0000000000000012),
    ("annulus 16x3/moves", 0x5067894e5f3a7054),
    ("annulus 16x3/refine/n_inserted", 0x0000000000000010),
    ("annulus 16x3/refine", 0x2426ef3f869e74b0),
    ("annulus 16x3/export/n_splits", 0x0000000000000000),
    ("annulus 16x3/export", 0x4ae4b589a7ce67cf),
    ("annulus 16x3/export/state", 0x2426ef3f869e74b0),
    ("annulus 16x3/odt/n_inserted", 0x0000000000000050),
    ("annulus 16x3/odt", 0x31ace29ea0ba65b2),
    ("sheared parallelogram 8 sigma 3/new", 0x1abb321da2567885),
    (
        "sheared parallelogram 8 sigma 3/flips/count",
        0x00000000000000e7,
    ),
    ("sheared parallelogram 8 sigma 3/flips", 0xbb98cd480948ae97),
    (
        "sheared parallelogram 8 sigma 3/insertions/count",
        0x000000000000000c,
    ),
    (
        "sheared parallelogram 8 sigma 3/insertions",
        0x129c48492e80115a,
    ),
    (
        "sheared parallelogram 8 sigma 3/splits/count",
        0x000000000000000c,
    ),
    ("sheared parallelogram 8 sigma 3/splits", 0xba13fdb9229d3e4c),
    (
        "sheared parallelogram 8 sigma 3/delaunay/n_flips",
        0x00000000000000ea,
    ),
    (
        "sheared parallelogram 8 sigma 3/delaunay",
        0x310746a643a181d3,
    ),
    (
        "sheared parallelogram 8 sigma 3/delaunay/length0",
        0x3fc0000000000000,
    ),
    (
        "sheared parallelogram 8 sigma 3/delaunay/length1",
        0x3fd94c583ada5b53,
    ),
    (
        "sheared parallelogram 8 sigma 3/delaunay/length2",
        0x3faf00a02e5ad4c9,
    ),
    (
        "sheared parallelogram 8 sigma 3/delaunay/length_last",
        0x3fb5bb15584a0123,
    ),
    ("sheared parallelogram 8 sigma 3/traces", 0x6f4578c406a4ab95),
    (
        "sheared parallelogram 8 sigma 3/overlay",
        0xacb2616b11e80bad,
    ),
    (
        "sheared parallelogram 8 sigma 3/moves/count",
        0x0000000000000015,
    ),
    ("sheared parallelogram 8 sigma 3/moves", 0xf1c383bd5a427898),
    (
        "sheared parallelogram 8 sigma 3/refine/n_inserted",
        0x0000000000000040,
    ),
    ("sheared parallelogram 8 sigma 3/refine", 0xbeed95c25c136279),
    (
        "sheared parallelogram 8 sigma 3/export/n_splits",
        0x0000000000000000,
    ),
    ("sheared parallelogram 8 sigma 3/export", 0xa988691bbdd5e49c),
    (
        "sheared parallelogram 8 sigma 3/export/state",
        0xbeed95c25c136279,
    ),
    (
        "sheared parallelogram 8 sigma 3/odt/n_inserted",
        0x0000000000000010,
    ),
    ("sheared parallelogram 8 sigma 3/odt", 0x909116ef74db9edc),
    ("needles and caps 6/new", 0xc99fa6f2e9c4c76b),
    ("needles and caps 6/flips/count", 0x0000000000000073),
    ("needles and caps 6/flips", 0x9284d0b302652c56),
    ("needles and caps 6/insertions/count", 0x000000000000000c),
    ("needles and caps 6/insertions", 0x1ac4d69f1e7b1d66),
    ("needles and caps 6/splits/count", 0x000000000000000c),
    ("needles and caps 6/splits", 0x7929a55980a4099f),
    ("needles and caps 6/delaunay/n_flips", 0x0000000000000063),
    ("needles and caps 6/delaunay", 0xa27a3921c1bfe5f4),
    ("needles and caps 6/delaunay/length0", 0x3fc08e2ddfc9110d),
    ("needles and caps 6/delaunay/length1", 0x3f80c9714fbcda3b),
    ("needles and caps 6/delaunay/length2", 0x3fb2c121b577d392),
    (
        "needles and caps 6/delaunay/length_last",
        0x3faee1d295c0f17a,
    ),
    ("needles and caps 6/traces", 0x9decfb1d96f92447),
    ("needles and caps 6/overlay", 0x6467cf1d8bde302e),
    ("needles and caps 6/moves/count", 0x0000000000000012),
    ("needles and caps 6/moves", 0xb50f48f008dbc876),
    ("needles and caps 6/refine/n_inserted", 0x00000000000000bb),
    ("needles and caps 6/refine", 0xb6df282bba22d440),
    ("needles and caps 6/export/n_splits", 0x0000000000000000),
    ("needles and caps 6/export", 0x2bb6c0bd6ebc50fc),
    ("needles and caps 6/export/state", 0xb6df282bba22d440),
    ("needles and caps 6/odt/n_inserted", 0x0000000000000036),
    ("needles and caps 6/odt", 0xce7f87d165597069),
];

/// Recorded on the code before the refactoring, on the meshes of `data/` (spot, camelbad).
const MODEL_PINS: &[(&str, u64)] = &[
    ("spot_triangulated/new", 0xdda05b2f9dd2f78a),
    ("spot_triangulated/flips/count", 0x0000000000003a44),
    ("spot_triangulated/flips", 0x644458500eef45fd),
    ("spot_triangulated/insertions/count", 0x000000000000000c),
    ("spot_triangulated/insertions", 0x2d0f7c12e9cf2c4f),
    ("spot_triangulated/splits/count", 0x000000000000000c),
    ("spot_triangulated/splits", 0x6adcc6ce6418bf09),
    ("spot_triangulated/delaunay/n_flips", 0x0000000000001035),
    ("spot_triangulated/delaunay", 0x550c6034a07407fc),
    ("spot_triangulated/delaunay/length0", 0x3fac38d40ce74dbf),
    ("spot_triangulated/delaunay/length1", 0x3fb412f00dda0502),
    ("spot_triangulated/delaunay/length2", 0x3fad124b8940a1d6),
    ("spot_triangulated/delaunay/length_last", 0x3fb094a44e1d30e1),
    ("spot_triangulated/traces", 0xfbe93c0e9d48eebd),
    ("spot_triangulated/overlay", 0xd59c37644144b108),
    ("spot_triangulated/moves/count", 0x0000000000000018),
    ("spot_triangulated/moves", 0xd4ad0692ac57dc9e),
    ("spot_triangulated/refine/n_inserted", 0x00000000000000ff),
    ("spot_triangulated/refine", 0x26ad58e112702073),
    ("spot_triangulated/export/n_splits", 0x0000000000000000),
    ("spot_triangulated/export", 0x4b6cdc42d132a557),
    ("spot_triangulated/export/state", 0x26ad58e112702073),
    ("spot_triangulated/odt/n_inserted", 0x0000000000001296),
    ("spot_triangulated/odt", 0x22377ab0aa27a21a),
    ("camelbad/new", 0x0522449d73f8c9b9),
    ("camelbad/flips/count", 0x000000000000dedc),
    ("camelbad/flips", 0xab5792aaec3c0037),
    ("camelbad/insertions/count", 0x000000000000000c),
    ("camelbad/insertions", 0x2dd0471d99fd2970),
    ("camelbad/splits/count", 0x000000000000000c),
    ("camelbad/splits", 0x2a1a8b56d655892a),
    ("camelbad/delaunay/n_flips", 0x0000000000004b2a),
    ("camelbad/delaunay", 0x819b0deac0a49f37),
    ("camelbad/delaunay/length0", 0x3f8957128b9a3c33),
    ("camelbad/delaunay/length1", 0x3f891a582953ab0c),
    ("camelbad/delaunay/length2", 0x3f87cef41f0595fd),
    ("camelbad/delaunay/length_last", 0x3f973c3e161fbfc1),
    ("camelbad/traces", 0xc43fde2235a156df),
    ("camelbad/overlay", 0xebd3537e5269006d),
    ("camelbad/moves/count", 0x0000000000000018),
    ("camelbad/moves", 0xc87bf909e7e09250),
    ("camelbad/refine/n_inserted", 0x0000000000001049),
    ("camelbad/refine", 0x445520dc235ab9ba),
    ("camelbad/export/n_splits", 0x0000000000000000),
    ("camelbad/export", 0xb09fa521e88c08f3),
    ("camelbad/export/state", 0x445520dc235ab9ba),
    ("camelbad/odt/n_inserted", 0x0000000000003bd8),
    ("camelbad/odt", 0xfdf2aa7a66984ad2),
];

fn compare(got: &[(String, u64)], pins: &[(&str, u64)], table: &str) {
    if std::env::var_os("DDG_PIN_PRINT").is_some() || pins.is_empty() {
        println!("const {table}: &[(&str, u64)] = &[");
        for (name, value) in got {
            println!("    ({name:?}, {value:#018x}),");
        }
        println!("];");
    }
    assert!(!pins.is_empty(), "no pins recorded in {table}");
    assert_eq!(got.len(), pins.len(), "number of pinned values in {table}");
    let mut wrong = Vec::new();
    for ((name, value), (pin_name, pin)) in got.iter().zip(pins) {
        assert_eq!(name, pin_name);
        if value != pin {
            wrong.push(format!("{name}: {value:#018x}, pinned {pin:#018x}"));
        }
    }
    assert!(wrong.is_empty(), "changed:\n{}", wrong.join("\n"));
}

#[test]
fn results_are_bitwise_those_recorded_before_the_refactoring() {
    let mut got = Vec::new();
    for shape in shapes() {
        record(&shape, &mut got);
    }
    compare(&got, PINS, "PINS");
}

#[test]
fn results_on_models_are_bitwise_those_recorded_before_the_refactoring() {
    let shapes: Vec<Shape> = [
        "models/spot/spot/spot_triangulated.obj",
        "optimal-cones/optimal-cones-data/data/camelbad.obj",
    ]
    .iter()
    .filter_map(|p| load_data(p))
    .collect();
    if shapes.len() < 2 {
        eprintln!("SKIPPED: the pins on models need spot and camelbad in data/");
        return;
    }
    let mut got = Vec::new();
    for shape in &shapes {
        record(shape, &mut got);
    }
    compare(&got, MODEL_PINS, "MODEL_PINS");
}
