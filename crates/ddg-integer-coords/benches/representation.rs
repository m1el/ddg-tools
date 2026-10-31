//! WP10, SPEC 10.2: benchmarks of the integer representation next to the signposts' on the same
//! inputs: a run of flips (1000 edges, each flipped and back), Delaunay flipping (IDT), Delaunay refinement to 25 degrees, the common
//! subdivision, and attribute transfer (integers only: `Transfer`). The inputs are mollified
//! (`1e-5`) as in the paper's evaluation. Meshes: `data/models` (spot, bob, blub), the optimal
//! cones dataset's camelbad, and the converted Thingi10k meshes 44395 (Fig. 15) and 54674
//! (Fig. 16) when `data/thingi10k/converted/` exists; missing files are skipped with a notice.
//!
//! ```text
//! cargo bench -p ddg-integer-coords --bench representation                # everything
//! cargo bench -p ddg-integer-coords --bench representation -- refine      # one operation
//! ```
//!
//! Single thread; run nothing else at the same time when the numbers matter. The signposts'
//! own retriangulation benches are `cargo bench -p ddg-intrinsic --bench retriangulation`.

use std::path::PathBuf;
use std::time::Duration;

use criterion::{BatchSize, Criterion, criterion_group, criterion_main};
use ddg_core::io::read_obj;
use ddg_core::mollify::mollify_edge_lengths;
use ddg_core::{Geometry, TriMesh};
use ddg_integer_coords::{
    InitOptions, IntegerRefineOptions, IntegerTriangulation, PAPER_MOLLIFY_EPSILON, Transfer,
};
use ddg_intrinsic::{
    DelaunayOptions, EdgeTraceOptions, IntrinsicTriangulation, RefineOptions, Triangulation,
};

struct Input {
    name: String,
    mesh: TriMesh,
    lengths: Vec<f64>,
}

fn from_file(name: &str, relative: &str) -> Option<Input> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../data")
        .join(relative);
    if !path.is_file() {
        eprintln!("SKIPPED: {} not found", path.display());
        return None;
    }
    let obj = read_obj(&path).ok()?;
    let mut faces = Vec::new();
    for f in 0..obj.n_faces() {
        let c = &obj.corner_verts[obj.face_corners(f)];
        for k in 1..c.len().saturating_sub(1) {
            faces.push([c[0], c[k], c[k + 1]]);
        }
    }
    let mesh = TriMesh::new(obj.positions.len(), faces).ok()?;
    let geometry = Geometry::from_positions(&mesh, &obj.positions).ok()?;
    let lengths = geometry.edge_lengths().to_vec();
    Some(Input {
        name: name.to_string(),
        mesh,
        lengths,
    })
}

fn integer(i: &Input) -> IntegerTriangulation {
    IntegerTriangulation::with_options(
        &i.mesh,
        &i.lengths,
        &InitOptions {
            mollify: Some(PAPER_MOLLIFY_EPSILON),
        },
    )
    .expect("init")
    .0
}

fn signpost(i: &Input) -> IntrinsicTriangulation {
    let m = mollify_edge_lengths(&i.mesh, &i.lengths, PAPER_MOLLIFY_EPSILON).expect("mollify");
    let g = Geometry::from_edge_lengths(&i.mesh, m.edge_lengths).expect("geometry");
    IntrinsicTriangulation::new(&i.mesh, &g).expect("init")
}

fn bench(c: &mut Criterion) {
    let inputs: Vec<Input> = [
        from_file("spot", "models/spot/spot/spot_triangulated.obj"),
        from_file("bob", "models/bob/bob/bob_isotropic.obj"),
        from_file("blub", "models/blub/blub/blub_triangulated.obj"),
        from_file(
            "camelbad",
            "optimal-cones/optimal-cones-data/data/camelbad.obj",
        ),
        from_file("thingi 44395 (Fig. 15)", "thingi10k/converted/44395.obj"),
        from_file("thingi 54674 (Fig. 16)", "thingi10k/converted/54674.obj"),
    ]
    .into_iter()
    .flatten()
    .collect();
    let delaunay = DelaunayOptions::default();
    for inp in inputs {
        let ints = integer(&inp);
        let sign = signpost(&inp);
        let mut g = c.benchmark_group(format!("representation/{}", inp.name));
        g.sample_size(10).measurement_time(Duration::from_secs(8));

        // Each of 1000 flippable edges flipped and flipped back (2000 flips), from the Delaunay state.
        let mut flat_i = ints.clone();
        flat_i.flip_to_delaunay(&delaunay).unwrap();
        let mut flat_s = sign.clone();
        flat_s.flip_to_delaunay(&delaunay).unwrap();
        let flippable: Vec<usize> = (0..flat_i.complex().n_edges())
            .filter(|&e| matches!(Triangulation::can_flip(&flat_i, e), Ok(Ok(()))))
            .take(1000)
            .collect();
        g.bench_function("2000 flips and flips back, integers", |b| {
            b.iter_batched(
                || flat_i.clone(),
                |mut t| {
                    for &e in &flippable {
                        t.flip_edge(e).unwrap();
                        t.flip_edge(e).unwrap();
                    }
                    t
                },
                BatchSize::LargeInput,
            )
        });
        g.bench_function("2000 flips and flips back, signposts", |b| {
            b.iter_batched(
                || flat_s.clone(),
                |mut t| {
                    for &e in &flippable {
                        t.flip_edge(e).unwrap();
                        t.flip_edge(e).unwrap();
                    }
                    t
                },
                BatchSize::LargeInput,
            )
        });
        g.bench_function("idt integers (with certificate)", |b| {
            b.iter_batched(
                || ints.clone(),
                |mut t| t.flip_to_delaunay_certified(&delaunay).unwrap(),
                BatchSize::LargeInput,
            )
        });
        g.bench_function("idt integers (no certificate)", |b| {
            b.iter_batched(
                || ints.clone(),
                |mut t| t.flip_to_delaunay(&delaunay).unwrap(),
                BatchSize::LargeInput,
            )
        });
        g.bench_function("idt signposts", |b| {
            b.iter_batched(
                || sign.clone(),
                |mut t| t.flip_to_delaunay(&delaunay).unwrap(),
                BatchSize::LargeInput,
            )
        });
        g.bench_function("refine25 integers (Chew, certified)", |b| {
            b.iter_batched(
                || ints.clone(),
                |mut t| t.refine(&IntegerRefineOptions::default()).unwrap(),
                BatchSize::LargeInput,
            )
        });
        g.bench_function("refine25 signposts (Ruppert)", |b| {
            b.iter_batched(
                || sign.clone(),
                |mut t| t.refine(&RefineOptions::default()).unwrap(),
                BatchSize::LargeInput,
            )
        });
        // Common subdivision and transfer of the refined states.
        let mut ri = ints.clone();
        ri.refine(&IntegerRefineOptions::default()).unwrap();
        let mut rs = sign.clone();
        rs.refine(&RefineOptions::default()).unwrap();
        g.bench_function("common subdivision integers (refined)", |b| {
            b.iter(|| ri.common_subdivision().unwrap())
        });
        g.bench_function(
            "common subdivision signposts (traces + overlay, refined)",
            |b| {
                b.iter(|| {
                    let traces = rs.trace_edges(&EdgeTraceOptions::default()).unwrap();
                    rs.overlay(&traces).unwrap()
                })
            },
        );
        g.bench_function("transfer operators integers (refined)", |b| {
            b.iter(|| Transfer::from_state(&ri).unwrap())
        });
        let tr = Transfer::from_state(&ri).unwrap();
        let f: Vec<f64> = (0..ri.complex().n_vertices())
            .map(|v| (v as f64 * 0.37).sin())
            .collect();
        g.bench_function("transfer solve T1 to T0 integers", |b| {
            b.iter(|| tr.to_input(&f).unwrap())
        });
        g.finish();
    }
}

criterion_group!(benches, bench);
criterion_main!(benches);
