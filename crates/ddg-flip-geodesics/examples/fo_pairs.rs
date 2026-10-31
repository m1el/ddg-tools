//! Times FlipOut on given vertex pairs of a mesh, for comparison with the reference driver
//! (geometry-central's `FlipEdgeNetwork`, SPEC §7 layer 5) on the same pairs.
//!
//! Usage: `cargo run --release --example fo_pairs -- MESH.obj PAIRS.txt [signposts|integers]`
//! (`PAIRS.txt`: one `a b` vertex pair per line). Prints one line per pair:
//! `a b dijkstra_ms construct_ms shorten_ms extract_ms length flips iterations min_angle geodesic points`,
//! the same columns as the reference driver. `construct` is the clone of the triangulation plus
//! building the network (the reference builds its signpost triangulation there).

use std::time::Instant;

use ddg_core::io::read_obj;
use ddg_core::{Geometry, TriMesh};
use ddg_flip_geodesics::{
    CertifyOptions, ExtractOptions, GeodesicNetwork, ShortenOptions, dijkstra_path,
};
use ddg_integer_coords::IntegerTriangulation;
use ddg_intrinsic::{IntrinsicTriangulation, Triangulation};

fn ms(t: Instant) -> f64 {
    t.elapsed().as_secs_f64() * 1e3
}

fn run<T: Triangulation + Clone + Sync>(tri: &T, pairs: &[(usize, usize)]) {
    let extract = ExtractOptions {
        certify_input: false,
        ..ExtractOptions::default()
    };
    for &(a, b) in pairs {
        let t0 = Instant::now();
        let path = match dijkstra_path(tri, a, b) {
            Ok(p) => p,
            Err(e) => {
                println!("{a} {b} skip {e}");
                continue;
            }
        };
        let t_dij = ms(t0);
        let t0 = Instant::now();
        let mut net = GeodesicNetwork::new(tri.clone(), std::slice::from_ref(&path), &[], None)
            .expect("a Dijkstra path is a valid network");
        let t_con = ms(t0);
        let t0 = Instant::now();
        let report = net.shorten(&ShortenOptions::default()).expect("shorten");
        let t_sh = ms(t0);
        let t0 = Instant::now();
        let poly = net.extract(&extract).expect("extract");
        let t_ex = ms(t0);
        let cert = net.certify(&CertifyOptions::default()).expect("certify");
        println!(
            "{a} {b} {t_dij:.6} {t_con:.6} {t_sh:.6} {t_ex:.6} {:.17e} {} {} {:.6e} {} {}",
            net.length(),
            report.flips,
            report.iterations,
            cert.min_angle,
            u8::from(cert.geodesic),
            poly[0].points.len()
        );
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("usage: fo_pairs MESH.obj PAIRS.txt [signposts|integers]");
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
    println!("# V {} F {}", mesh.n_vertices(), mesh.n_faces());
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
