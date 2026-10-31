//! Profiles `GeodesicNetwork::shorten` alone: for every pair the network is built first,
//! then `shorten` runs `REPS` times on clones of it, with a `perf record --control` FIFO
//! (environment `PERF_CTL`, `PERF_ACK`) enabled only around the shortening.
//!
//! Usage: `perf record -D -1 --control fifo:CTL,ACK … fo_prof MESH.obj PAIRS.txt
//! [signposts|integers] [REPS]`. Prints the total shortening time.

use std::io::{Read, Write};
use std::time::Instant;

use ddg_core::io::read_obj;
use ddg_core::{Geometry, TriMesh};
use ddg_flip_geodesics::{GeodesicNetwork, ShortenOptions, dijkstra_path};
use ddg_integer_coords::IntegerTriangulation;
use ddg_intrinsic::{IntrinsicTriangulation, Triangulation};

struct Ctl(Option<(std::fs::File, std::fs::File)>);

impl Ctl {
    fn send(&mut self, what: &[u8]) {
        if let Some((ctl, ack)) = &mut self.0 {
            ctl.write_all(what).expect("perf control");
            let mut buf = [0u8; 5];
            let _ = ack.read(&mut buf);
        }
    }
}

fn run<T: Triangulation + Clone>(tri: &T, pairs: &[(usize, usize)], reps: usize, ctl: &mut Ctl) {
    let mut total = 0.0;
    for &(a, b) in pairs {
        let Ok(path) = dijkstra_path(tri, a, b) else {
            continue;
        };
        let net = GeodesicNetwork::new(tri.clone(), std::slice::from_ref(&path), &[], None)
            .expect("a Dijkstra path is a valid network");
        for _ in 0..reps {
            let mut n = net.clone();
            ctl.send(b"enable\n");
            let t0 = Instant::now();
            n.shorten(&ShortenOptions::default()).expect("shorten");
            total += t0.elapsed().as_secs_f64();
            ctl.send(b"disable\n");
        }
    }
    println!("shorten total {:.3} ms", total * 1e3);
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
    let pairs: Vec<(usize, usize)> = std::fs::read_to_string(&args[2])
        .unwrap_or_else(|e| panic!("{}: {e}", args[2]))
        .lines()
        .filter_map(|l| {
            let mut it = l.split_whitespace().map(|x| x.parse::<usize>());
            Some((it.next()?.ok()?, it.next()?.ok()?))
        })
        .collect();
    let reps: usize = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(1);
    let mut ctl = Ctl(
        match (std::env::var("PERF_CTL"), std::env::var("PERF_ACK")) {
            (Ok(c), Ok(a)) => Some((
                std::fs::OpenOptions::new()
                    .write(true)
                    .open(c)
                    .expect("ctl fifo"),
                std::fs::File::open(a).expect("ack fifo"),
            )),
            _ => None,
        },
    );
    match args.get(3).map(String::as_str).unwrap_or("signposts") {
        "integers" => run(
            &IntegerTriangulation::new(&mesh, &geometry).expect("integer triangulation"),
            &pairs,
            reps,
            &mut ctl,
        ),
        _ => run(
            &IntrinsicTriangulation::new(&mesh, &geometry).expect("signpost triangulation"),
            &pairs,
            reps,
            &mut ctl,
        ),
    }
}
