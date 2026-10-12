//! Differential test against the reference C++ implementation.
//!
//! Reads the output of the driver `bff_ref` (`reference/bff_ref.cpp`: a small program around
//! the reference library that dumps, with full precision, the positions it used, the UV of
//! every corner, the cut it chose and where it started the boundary) and flattens the same
//! mesh with the same cones and the same cut. `reference/run.sh` runs everything.
//!
//! ```text
//! reference_diff cones out.txt cones.txt stored.obj|- [more triples ...]
//! reference_diff auto  out.txt [more ...]
//! reference_diff disk  out.txt [more ...]
//! ```
//!
//! One line per case. All residuals are after the best similarity per chart, relative to the
//! diagonal of the chart's bounding box in the reference layout.
//!
//! * `V`, `F`, `cones`, `cutE` — sizes; `cutE` is the number of edges of the driver's cut.
//! * `mollify` — the amount the reference adds to every edge length (zero unless a triangle
//!   is degenerate to within `1e-6` of the mean edge length). The geometry used here is
//!   mollified in the same way.
//! * `exact: rms, max` — this crate with its default exact solves against the driver.
//! * `shifted: rms, max` — this crate with the reference's diagonal shift of `1e-8` against
//!   the driver: the two programs then compute the same thing, and this is rounding error.
//! * `scale` — scale of the fitted similarity in the shifted run (driver size / our size).
//!
//! In mode `cones`, `stored.obj` is a file with a UV map of the same mesh to compare the
//! *driver's* output with (`-` for none); this says how well today's reference reproduces
//! that file, and nothing about this crate:
//!
//! * `stored cut` — `same` if the driver chose exactly the cut that is visible in the stored
//!   UV map, otherwise `common/stored` edge counts; a line below the table then describes
//!   the difference.
//! * `ref-file: rms, max` — the driver's layout against the stored one (same cut only).
//! * `metric p1, p99` — the ratio of the lengths of every triangle edge in the driver's
//!   layout and in the stored one, divided by its median: 1st and 99th percentile. The flat
//!   metric does not depend on the cut, so this compares maps with different cuts too.

#[path = "../tests/common/mod.rs"]
mod common;

use common::{compare_layouts, length_ratio_quantiles};
use ddg_bff::{
    ConeFlattenInput, ConeFlattenOptions, CutSource, DiskBff, DiskOptions, Extension,
    GaussBonnetPolicy, REFERENCE_DIAGONAL_SHIFT, REFERENCE_MOLLIFY_DELTA, ScaleFactorOptions,
    flatten_with_cones, mollify_edge_lengths,
};
use ddg_core::geometry::edge_length;
use ddg_core::io::read_cones;
use ddg_core::{Geometry, TriMesh};

struct Reference {
    positions: Vec<[f64; 3]>,
    faces: Vec<[u32; 3]>,
    corner_uvs: Vec<[f64; 2]>,
    cut: Vec<[u32; 2]>,
    /// The wedge at which the reference starts its boundary walk: (vertex, face).
    start: Option<(u32, u32)>,
}

fn parse(path: &str) -> Reference {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let mut r = Reference {
        positions: Vec::new(),
        faces: Vec::new(),
        corner_uvs: Vec::new(),
        cut: Vec::new(),
        start: None,
    };
    for line in text.lines() {
        let mut tokens = line.split_ascii_whitespace();
        let numbers = |tokens: std::str::SplitAsciiWhitespace<'_>| -> Vec<f64> {
            tokens
                .map(|t| {
                    t.parse()
                        .unwrap_or_else(|_| panic!("{path}: bad number {t}"))
                })
                .collect()
        };
        match tokens.next() {
            Some("f") => {
                let x = numbers(tokens);
                r.faces.push([x[0] as u32, x[1] as u32, x[2] as u32]);
                r.corner_uvs
                    .extend([[x[3], x[4]], [x[5], x[6]], [x[7], x[8]]]);
            }
            Some("c") => {
                let x = numbers(tokens);
                r.cut.push([x[0] as u32, x[1] as u32]);
            }
            Some("start") => {
                let x = numbers(tokens);
                r.start = Some((x[0] as u32, x[1] as u32));
            }
            Some("p") => {
                let x = numbers(tokens);
                r.positions.push([x[0], x[1], x[2]]);
            }
            _ => {}
        }
    }
    r
}

/// The geometry the reference builds: edge lengths from positions, mollified.
fn reference_geometry(mesh: &TriMesh, positions: &[[f64; 3]]) -> (Geometry, f64) {
    let lengths: Vec<f64> = mesh
        .edges()
        .iter()
        .map(|&[a, b]| edge_length(positions[a as usize], positions[b as usize]))
        .collect();
    let mollified = mollify_edge_lengths(mesh, &lengths, REFERENCE_MOLLIFY_DELTA).unwrap();
    let epsilon = mollified.epsilon;
    (
        Geometry::from_edge_lengths(mesh, mollified.edge_lengths).unwrap_or_else(|e| panic!("{e}")),
        epsilon,
    )
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mode = args.first().map(String::as_str).unwrap_or("");
    let step = if mode == "cones" { 3 } else { 1 };
    if !matches!(mode, "cones" | "auto" | "disk") || args.len() < 1 + step {
        eprintln!(
            "usage: reference_diff cones out.txt cones.txt stored.obj|- ... | auto out.txt ... | disk out.txt ..."
        );
        std::process::exit(2);
    }
    println!(
        "{:<20} {:>7} {:>7} {:>6} {:>6} {:>8} | {:>8} {:>8} | {:>8} {:>8} | {:>8}{}",
        "case",
        "V",
        "F",
        "cones",
        "cutE",
        "mollify",
        "exact",
        "max",
        "shifted",
        "max",
        "scale",
        if mode == "cones" {
            format!(
                " | {:>11} {:>8} {:>8} | {:>7} {:>7}",
                "stored cut", "ref-file", "max", "met p1", "p99"
            )
        } else {
            String::new()
        }
    );
    let mut notes: Vec<String> = Vec::new();
    for case in args[1..].chunks(step) {
        let reference = parse(&case[0]);
        let name = std::path::Path::new(&case[0])
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("?")
            .trim_end_matches(".out.txt")
            .to_string();
        let mesh = TriMesh::new(reference.positions.len(), reference.faces.clone())
            .unwrap_or_else(|e| panic!("{name}: {e}"));
        let (geometry, epsilon) = reference_geometry(&mesh, &reference.positions);
        let mut line = format!(
            "{:<20} {:>7} {:>7}",
            name,
            mesh.n_vertices(),
            mesh.n_faces()
        );
        let mut results = Vec::new();
        let mut n_cones = 0;
        for shift in [0.0, REFERENCE_DIAGONAL_SHIFT] {
            let (uvs, face_chart, n_charts) = if mode == "cones" {
                let cones =
                    read_cones(&case[1], mesh.n_vertices()).unwrap_or_else(|e| panic!("{e}"));
                n_cones = cones.len();
                let cut: Vec<u32> = reference
                    .cut
                    .iter()
                    .map(|&[a, b]| mesh.find_edge(a as usize, b as usize).expect("cut edge") as u32)
                    .collect();
                // Start the boundary curve where the reference starts it. This matters only
                // with the diagonal shift, which makes the exterior angles miss 2π slightly.
                let anchors: Vec<u32> = reference
                    .start
                    .iter()
                    .map(|&(v, f)| {
                        let k = reference.faces[f as usize]
                            .iter()
                            .position(|&x| x == v)
                            .expect("start vertex is a corner of the start face");
                        3 * f + k as u32
                    })
                    .collect();
                let flat = flatten_with_cones(
                    &ConeFlattenInput {
                        mesh: &mesh,
                        geometry: &geometry,
                        cones: &cones,
                        scale_factors: None,
                        cut: CutSource::Edges(&cut),
                        anchors: &anchors,
                    },
                    &ConeFlattenOptions {
                        scale_factors: ScaleFactorOptions {
                            gauss_bonnet: GaussBonnetPolicy::SpreadOverVertices,
                            diagonal_shift: shift,
                            ..ScaleFactorOptions::default()
                        },
                        max_length_adjustment: f64::INFINITY,
                        ..ConeFlattenOptions::default()
                    },
                )
                .unwrap_or_else(|e| panic!("{name}: {e}"));
                (
                    flat.corner_uvs.clone(),
                    flat.face_chart.clone(),
                    flat.n_charts(),
                )
            } else {
                let mut bff = DiskBff::new(
                    &mesh,
                    &geometry,
                    DiskOptions {
                        diagonal_shift: shift,
                        boundary_start: reference.start.map(|(v, _)| v),
                        // The reference closes its curves whatever it takes.
                        max_length_adjustment: f64::INFINITY,
                        ..DiskOptions::default()
                    },
                )
                .unwrap_or_else(|e| panic!("{name}: {e}"));
                let flat = if mode == "auto" {
                    let zero = vec![0.0; bff.boundary().len()];
                    bff.flatten_with_scale_factors(&zero, Extension::Conjugate)
                } else {
                    bff.flatten_to_disk(10)
                }
                .unwrap_or_else(|e| panic!("{name}: {e}"));
                let corner: Vec<[f64; 2]> = mesh
                    .faces()
                    .iter()
                    .flat_map(|f| f.map(|v| flat.uv[v as usize]))
                    .collect();
                (corner, vec![0u32; mesh.n_faces()], 1)
            };
            results.push(compare_layouts(
                &face_chart,
                n_charts,
                &uvs,
                &reference.corner_uvs,
            ));
        }
        line += &format!(
            " {:>6} {:>6} {:>8.1e} | {:8.1e} {:8.1e} | {:8.1e} {:8.1e} | {:8.2e}",
            n_cones,
            reference.cut.len(),
            epsilon,
            results[0].rms,
            results[0].max,
            results[1].rms,
            results[1].max,
            results[1].scale,
        );
        if mode == "cones" && case[2] != "-" {
            // The driver's output against a stored UV map of the same mesh.
            let obj = ddg_core::io::read_obj(&case[2]).unwrap_or_else(|e| panic!("{e}"));
            let stored = obj.corner_uvs().expect("stored UVs");
            let mut stored_cut = ddg_bff::discontinuous_edges(&mesh, &stored).unwrap();
            let mut driver_cut: Vec<u32> = reference
                .cut
                .iter()
                .map(|&[a, b]| mesh.find_edge(a as usize, b as usize).expect("cut edge") as u32)
                .collect();
            stored_cut.sort_unstable();
            driver_cut.sort_unstable();
            let q = length_ratio_quantiles(&reference.corner_uvs, &stored);
            if stored_cut == driver_cut {
                let cut = ddg_bff::cut_mesh(&mesh, &driver_cut).unwrap();
                let a = compare_layouts(
                    cut.mesh.face_components(),
                    cut.mesh.n_components(),
                    &reference.corner_uvs,
                    &stored,
                );
                line += &format!(" | {:>11} {:8.1e} {:8.1e}", "same", a.rms, a.max);
            } else {
                let only_stored: Vec<u32> = stored_cut
                    .iter()
                    .copied()
                    .filter(|e| driver_cut.binary_search(e).is_err())
                    .collect();
                let only_driver: Vec<u32> = driver_cut
                    .iter()
                    .copied()
                    .filter(|e| stored_cut.binary_search(e).is_err())
                    .collect();
                line += &format!(
                    " | {:>11} {:>8} {:>8}",
                    format!(
                        "{}/{}",
                        stored_cut.len() - only_stored.len(),
                        stored_cut.len()
                    ),
                    "-",
                    "-"
                );
                let length = |edges: &[u32]| {
                    edges
                        .iter()
                        .map(|&e| geometry.edge_lengths()[e as usize])
                        .sum::<f64>()
                };
                let list = |edges: &[u32]| {
                    if edges.len() > 6 {
                        String::new()
                    } else {
                        format!(
                            " {:?}",
                            edges
                                .iter()
                                .map(|&e| mesh.edges()[e as usize])
                                .collect::<Vec<_>>()
                        )
                    }
                };
                notes.push(format!(
                    "{name}: {} edges only in the stored cut (total length {:.6}){}, {} only in \
                     the driver's (total length {:.6}){}",
                    only_stored.len(),
                    length(&only_stored),
                    list(&only_stored),
                    only_driver.len(),
                    length(&only_driver),
                    list(&only_driver),
                ));
            }
            line += &format!(" | {:7.4} {:7.4}", q[1], q[3]);
        }
        println!("{line}");
    }
    for note in &notes {
        println!("{note}");
    }
}
