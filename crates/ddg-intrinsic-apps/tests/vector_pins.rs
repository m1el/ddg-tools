//! Bitwise pins of the tangent fields on signposts (Integer Coordinates I26): every bit of the
//! connection matrices, the smoothest fields, their singularities, the harmonic fields, the
//! transport and the halfedge directions on a few triangulations, hashed. They were recorded
//! with the signpost-only code before the tangent fields became generic over the
//! representation; the generic code must give the signposts the same bits.

mod common;

use common::*;
use ddg_intrinsic::{DelaunayOptions, IntrinsicTriangulation, RefineOptions};
use ddg_intrinsic_apps::vector::*;

/// FNV-1a over 64-bit words.
struct Hash(u64);

impl Hash {
    fn new() -> Hash {
        Hash(0xcbf2_9ce4_8422_2325)
    }
    fn word(&mut self, x: u64) {
        for b in x.to_le_bytes() {
            self.0 ^= b as u64;
            self.0 = self.0.wrapping_mul(0x0100_0000_01b3);
        }
    }
    fn real(&mut self, x: f64) {
        self.word(x.to_bits());
    }
    fn reals(&mut self, xs: &[f64]) {
        self.word(xs.len() as u64);
        for &x in xs {
            self.real(x);
        }
    }
    fn complex(&mut self, xs: &[C64]) {
        self.word(xs.len() as u64);
        for x in xs {
            self.real(x.re);
            self.real(x.im);
        }
    }
}

fn triangulations() -> Vec<(String, IntrinsicTriangulation, Vec<[f64; 3]>)> {
    let mut out = Vec::new();
    for shape in [
        jittered_icosphere(2, 0.12, 5),
        bumpy_sheared_grid(8, 0.2, 1.5),
        sheared_torus(12, 6, 1.0, 0.4, 2),
    ] {
        let mut tri = shape.triangulation();
        out.push((
            format!("{} input", shape.name),
            tri.clone(),
            shape.positions.clone(),
        ));
        tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        out.push((
            format!("{} iDT", shape.name),
            tri.clone(),
            shape.positions.clone(),
        ));
        tri.refine(&RefineOptions::default()).unwrap();
        out.push((
            format!("{} refined", shape.name),
            tri,
            shape.positions.clone(),
        ));
    }
    out
}

/// The hash of everything the tangent fields compute on `tri`.
fn fields_hash(tri: &IntrinsicTriangulation, positions: &[[f64; 3]]) -> u64 {
    let mut h = Hash::new();
    h.reals(&halfedge_directions(tri));
    h.reals(&tangent_scales(tri));
    let transport: Vec<f64> = (0..tri.n_halfedges())
        .map(|e| transport_angle(tri, e))
        .collect();
    h.reals(&transport);
    let holonomy: Vec<f64> = (0..tri.n_faces()).map(|f| face_holonomy(tri, f)).collect();
    h.reals(&holonomy);
    let boundary = (0..tri.n_vertices()).any(|v| tri.is_boundary_vertex(v));
    for n_sym in [1u32, 2, 4] {
        for kind in [
            ConnectionLaplacian::FiniteDifference,
            ConnectionLaplacian::FiniteElement,
        ] {
            let m = connection_matrices(tri, n_sym, kind).unwrap();
            h.reals(m.stiffness.values());
            h.reals(m.mass.values());
            let options = DirectionFieldOptions {
                n_sym,
                laplacian: kind,
                max_iterations: 50,
                ..Default::default()
            };
            match smoothest_field(tri, &options) {
                Ok((field, report)) => {
                    h.complex(&field.values);
                    h.real(report.eigenvalue);
                    h.real(report.residual);
                    h.word(report.iterations as u64);
                    let s = singularities(tri, &field).unwrap();
                    h.word(s.index_sum as u64);
                    for x in &s.faces {
                        h.word(((x.face as u64) << 32) | (x.index as u32 as u64));
                    }
                    h.real(s.max_rounding);
                    let v3 = field.input_vectors_3d(tri, positions).unwrap();
                    for p in v3 {
                        h.reals(&p);
                    }
                }
                Err(e) => h.word(e.to_string().len() as u64),
            }
            if boundary {
                for bv in [BoundaryValues::InwardNormal, BoundaryValues::Tangent] {
                    let options = HarmonicFieldOptions {
                        n_sym,
                        laplacian: kind,
                        boundary: bv,
                    };
                    let (field, report) = harmonic_field(tri, &options).unwrap();
                    h.complex(&field.values);
                    h.real(report.residual);
                    let r = check_maximum_principle(tri, &field, 1e-9);
                    h.word(r.n_checked as u64);
                    h.word(r.n_outside_cone as u64);
                    h.real(r.min_margin);
                }
            }
        }
    }
    h.0
}

/// Recorded before the generic code (commit a492320, signpost-only `vector`).
const PINS: [(&str, u64); 9] = [
    ("jittered_icosphere_2 input", 0x44c44ba3cbad7f4a),
    ("jittered_icosphere_2 iDT", 0x7cc66799eb62ac0d),
    ("jittered_icosphere_2 refined", 0x3bff3b29c2447265),
    ("bumpy_sheared_grid_8 input", 0x920f35245e350f87),
    ("bumpy_sheared_grid_8 iDT", 0x438291217b4fdae0),
    ("bumpy_sheared_grid_8 refined", 0xe8f22cfc8a95b6ee),
    ("sheared torus 12x6 shift 2 input", 0x1b1c4e0044ee0da9),
    ("sheared torus 12x6 shift 2 iDT", 0x1789c7cd21c9985a),
    ("sheared torus 12x6 shift 2 refined", 0x1789c7cd21c9985a),
];

#[test]
fn signpost_tangent_fields_are_bitwise_pinned() {
    let tris = triangulations();
    let mut mismatches = Vec::new();
    for ((name, tri, positions), (pin_name, pin)) in tris.iter().zip(PINS) {
        let hash = fields_hash(tri, positions);
        println!("    (\"{name}\", 0x{hash:016x}),");
        if name != pin_name || hash != pin {
            mismatches.push(name.clone());
        }
    }
    assert!(mismatches.is_empty(), "changed: {mismatches:?}");
}
