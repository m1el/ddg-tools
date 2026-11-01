//! Bitwise pins of `exact_distance` on `TriMesh`es (N37): the generalization to Δ-complexes
//! (N37b) shares the propagation code and must not change a bit here. The hashes were
//! recorded before the refactoring.

mod common;

use common::{Shape, cube, jittered_grid};
use ddg_core::generators::{annulus, bumpy_grid_disk, icosphere, torus};
use ddg_intrinsic_apps::exact::{ExactOptions, exact_distance};

fn fnv(h: &mut u64, x: u64) {
    for k in 0..8 {
        *h ^= (x >> (8 * k)) & 0xff;
        *h = h.wrapping_mul(0x100000001b3);
    }
}

/// A hash of the bits of the distances and of every field of the report.
fn digest(shape: &Shape, sources: &[usize]) -> u64 {
    let r = exact_distance(&shape.mesh, &shape.geometry, sources, &ExactOptions::default())
        .unwrap_or_else(|e| panic!("{}: {e}", shape.name));
    let mut h = 0xcbf29ce484222325u64;
    for d in &r.distances {
        fnv(&mut h, d.to_bits());
    }
    let p = r.report;
    for x in [
        p.n_created,
        p.n_propagated,
        p.n_filtered,
        p.n_beyond,
        p.n_thin,
        p.n_pseudo_sources,
        p.max_queue,
    ] {
        fnv(&mut h, x as u64);
    }
    fnv(&mut h, p.max_clamp.to_bits());
    h
}

fn pinned_cases() -> Vec<(Shape, Vec<usize>)> {
    vec![
        (Shape::generated("icosphere2", icosphere(2).unwrap()), vec![0]),
        (Shape::generated("icosphere3", icosphere(3).unwrap()), vec![3, 100]),
        (Shape::generated("torus", torus(12, 7, 1.0, 0.4).unwrap()), vec![5]),
        (
            Shape::generated("bumpy", bumpy_grid_disk(14, 9, 0.3).unwrap()),
            vec![0, 40],
        ),
        (
            Shape::generated("annulus", annulus(16, 4, 0.5, 1.0).unwrap()),
            vec![2],
        ),
        (cube(3), vec![0]),
        (jittered_grid(10, 0.3, 5), vec![17]),
    ]
}

const PINS: &[(&str, u64)] = &[
    ("icosphere2", 0x0feaa2bcb516a688),
    ("icosphere3", 0x66e2370e8292dd30),
    ("torus", 0x3151707a468b9cbd),
    ("bumpy", 0xb79039628b035a54),
    ("annulus", 0xceea75405b75dac9),
    ("cube 3", 0x503e9db324e45dc3),
    ("jittered grid 10", 0xa1336b42fbcc30d8),
];

#[test]
fn trimesh_entry_point_is_bitwise_pinned() {
    let mut bad = Vec::new();
    for (shape, sources) in pinned_cases() {
        let d = digest(&shape, &sources);
        println!("PIN {} {:#018x}", shape.name, d);
        match PINS.iter().find(|(n, _)| *n == shape.name) {
            Some((_, p)) if *p == d => {}
            _ => bad.push(shape.name.clone()),
        }
    }
    assert!(bad.is_empty(), "pin mismatch: {bad:?}");
}
