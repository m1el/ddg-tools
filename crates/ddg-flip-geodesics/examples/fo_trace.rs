//! Prints every iteration of FlipOut on one vertex pair (vertex, side, angle, flips, refused
//! flips, what happened, running length), then the report and the certificate's failures.
//!
//! Usage: `cargo run --release --example fo_trace -- MESH.obj A B [signposts|integers]`.

use ddg_core::io::read_obj;
use ddg_core::{Geometry, TriMesh};
use ddg_flip_geodesics::{CertifyOptions, GeodesicNetwork, ShortenOptions, dijkstra_path};
use ddg_integer_coords::IntegerTriangulation;
use ddg_intrinsic::{IntrinsicTriangulation, Triangulation};

fn run<T: Triangulation + Clone>(tri: &T, a: usize, b: usize) {
    let path = dijkstra_path(tri, a, b).expect("dijkstra");
    let mut net = GeodesicNetwork::new(tri.clone(), std::slice::from_ref(&path), &[], None)
        .expect("a Dijkstra path is a valid network");
    println!(
        "# segments {} length {:.17e}",
        path.segments.len(),
        net.length()
    );
    let report = net
        .shorten_with(&ShortenOptions::default(), |s, _| {
            println!(
                "it {} v {} {:?} seg ({}, {}) angle {:.17e} flips {} refused {} {:?} length {:.17e}",
                s.iteration,
                s.vertex,
                s.side,
                s.joint.incoming.0,
                s.joint.outgoing.0,
                s.angle,
                s.flips,
                s.refused_flips,
                s.shortening,
                s.length
            );
        })
        .expect("shorten");
    println!("# {report:?}");
    let cert = net.certify(&CertifyOptions::default()).expect("certify");
    println!(
        "# geodesic {} min_angle {:.17e} failures {:?}",
        cert.geodesic, cert.min_angle, cert.failures
    );
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let obj = read_obj(&args[1]).unwrap_or_else(|e| panic!("{e}"));
    let mut faces = Vec::new();
    for f in 0..obj.n_faces() {
        let corners = &obj.corner_verts[obj.face_corners(f)];
        for k in 1..corners.len().saturating_sub(1) {
            faces.push([corners[0], corners[k], corners[k + 1]]);
        }
    }
    let mesh = TriMesh::new(obj.positions.len(), faces).unwrap_or_else(|e| panic!("{e}"));
    let geometry =
        Geometry::from_positions(&mesh, &obj.positions).unwrap_or_else(|e| panic!("{e}"));
    let a: usize = args[2].parse().expect("vertex A");
    let b: usize = args[3].parse().expect("vertex B");
    match args.get(4).map(String::as_str).unwrap_or("signposts") {
        "integers" => run(
            &IntegerTriangulation::new(&mesh, &geometry).expect("integer triangulation"),
            a,
            b,
        ),
        _ => run(
            &IntrinsicTriangulation::new(&mesh, &geometry).expect("signpost triangulation"),
            a,
            b,
        ),
    }
}
