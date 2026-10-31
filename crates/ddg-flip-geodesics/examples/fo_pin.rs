//! Fingerprints FlipOut's results on given vertex pairs, to check that a change of the
//! shortening code keeps every result bitwise (not a timing harness: it observes every
//! iteration and hashes the final state).
//!
//! Usage: `cargo run --release --example fo_pin -- MESH.obj PAIRS.txt [signposts|integers]`.
//! Prints one line per pair:
//! `a b length_bits flips iterations status refused stalled sweeps stale blocked loop_moves
//! contractions twice_incident degree_one min_angle_bits geodesic failures steps_hash state_hash`,
//! `steps_hash` over every iteration's `ShortenStep` (joint, side, vertex, angle, flips, refused,
//! shortening, running length; `Debug` prints `f64` exactly) and `state_hash` over the final
//! connectivity, edge-length bits, paths and layering.

use std::hash::{DefaultHasher, Hash, Hasher};

use ddg_core::io::read_obj;
use ddg_core::{Geometry, TriMesh};
use ddg_flip_geodesics::{CertifyOptions, GeodesicNetwork, ShortenOptions, dijkstra_path};
use ddg_integer_coords::IntegerTriangulation;
use ddg_intrinsic::{IntrinsicTriangulation, Triangulation};

fn run<T: Triangulation + Clone>(tri: &T, pairs: &[(usize, usize)]) {
    for &(a, b) in pairs {
        let path = match dijkstra_path(tri, a, b) {
            Ok(p) => p,
            Err(e) => {
                println!("{a} {b} skip {e}");
                continue;
            }
        };
        let mut net = GeodesicNetwork::new(tri.clone(), std::slice::from_ref(&path), &[], None)
            .expect("a Dijkstra path is a valid network");
        let mut steps = DefaultHasher::new();
        let report = net
            .shorten_with(&ShortenOptions::default(), |step, _| {
                format!("{step:?}").hash(&mut steps);
            })
            .expect("shorten");
        let cert = net.certify(&CertifyOptions::default()).expect("certify");
        let mut state = DefaultHasher::new();
        let t = net.triangulation();
        let c = t.complex();
        for h in 0..c.n_halfedges() {
            (c.from_vertex(h), c.edge(h)).hash(&mut state);
        }
        for e in 0..c.n_edges() {
            t.edge_length(e).to_bits().hash(&mut state);
        }
        format!("{:?}", net.paths()).hash(&mut state);
        format!("{:?}", net.network().layering()).hash(&mut state);
        let r = &report;
        println!(
            "{a} {b} {:016x} {} {} {:?} {} {} {} {} {} {} {} {} {} {:016x} {} {} {:016x} {:016x}",
            net.length().to_bits(),
            r.flips,
            r.iterations,
            r.status,
            r.refused_flips,
            r.stalled,
            r.sweeps,
            r.stale_entries,
            r.blocked_entries,
            r.loop_perturbations,
            r.contractions,
            r.twice_incident_wedges,
            r.degree_one_wedges,
            cert.min_angle.to_bits(),
            u8::from(cert.geodesic),
            cert.failures.len(),
            steps.finish(),
            state.finish(),
        );
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("usage: fo_pin MESH.obj PAIRS.txt [signposts|integers]");
        std::process::exit(2);
    }
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
    let pairs: Vec<(usize, usize)> = std::fs::read_to_string(&args[2])
        .unwrap_or_else(|e| panic!("{}: {e}", args[2]))
        .lines()
        .filter_map(|l| {
            let mut it = l.split_whitespace().map(|x| x.parse::<usize>());
            Some((it.next()?.ok()?, it.next()?.ok()?))
        })
        .collect();
    match args.get(3).map(String::as_str).unwrap_or("signposts") {
        "integers" => run(
            &IntegerTriangulation::new(&mesh, &geometry).expect("integer triangulation"),
            &pairs,
        ),
        _ => run(
            &IntrinsicTriangulation::new(&mesh, &geometry).expect("signpost triangulation"),
            &pairs,
        ),
    }
}
