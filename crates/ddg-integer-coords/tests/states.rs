//! Integer coordinates read off signpost triangulations after random flips, Delaunay flipping,
//! insertions and edge splits (WP1 item 5): every conversion either fails with an error or
//! passes the certificate; the combinatorial tracing agrees with the signposts' traces; and the
//! certificate rejects every single mutation of the integers.
//!
//! Run with `--nocapture` to see the tables.

mod common;

use std::collections::BTreeMap;
use std::sync::OnceLock;

use common::{SCENARIOS, SEEDS, Scenario, all_shapes, convert, state};
use ddg_core::generators::SplitMix64;
use ddg_integer_coords::{
    ConversionReport, Crossing, IntegerError, IntegerTriangulation, NO_ROUNDABOUT, Parts,
    PathPiece, ValidityItem,
};
use ddg_intrinsic::{EdgeTraceOptions, IntrinsicTriangulation};

struct Job {
    shape: String,
    scenario: Scenario,
    seed: u64,
    signposts: IntrinsicTriangulation,
    result: Result<(IntegerTriangulation, ConversionReport), IntegerError>,
}

/// Every (shape, scenario, seed), computed once per test binary, in parallel.
fn jobs() -> &'static [Job] {
    static JOBS: OnceLock<Vec<Job>> = OnceLock::new();
    JOBS.get_or_init(|| {
        let shapes = all_shapes();
        let mut specs = Vec::new();
        for (i, _) in shapes.iter().enumerate() {
            for scenario in SCENARIOS {
                for seed in SEEDS {
                    specs.push((i, scenario, seed));
                }
            }
        }
        let threads = std::thread::available_parallelism()
            .map_or(4, |n| n.get())
            .min(48);
        let next = std::sync::atomic::AtomicUsize::new(0);
        let mut out: Vec<Option<Job>> = (0..specs.len()).map(|_| None).collect();
        let results = std::sync::Mutex::new(&mut out);
        std::thread::scope(|scope| {
            for _ in 0..threads {
                scope.spawn(|| {
                    loop {
                        let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                        if i >= specs.len() {
                            break;
                        }
                        let (s, scenario, seed) = specs[i];
                        let shape = &shapes[s];
                        let signposts = state(shape, scenario, seed);
                        let result = convert(&signposts);
                        let job = Job {
                            shape: shape.name.clone(),
                            scenario,
                            seed,
                            signposts,
                            result,
                        };
                        results.lock().unwrap()[i] = Some(job);
                    }
                });
            }
        });
        out.into_iter().map(|j| j.unwrap()).collect()
    })
}

#[test]
fn every_conversion_certifies_or_fails_with_an_error() {
    let jobs = jobs();
    let mut failures = Vec::new();
    let mut n_ok = 0;
    eprintln!(
        "{:<40} {:<26} {:>5} {:>7} {:>7} {:>9} {:>5} {:>6} {:>6}",
        "shape", "scenario", "seed", "|V|", "|E|", "Σn⁺", "on E", "snap", "drop"
    );
    for job in jobs {
        match &job.result {
            Ok((t, report)) => {
                n_ok += 1;
                // The certificate is part of the conversion; run it once more, standalone.
                let again = t.certify().expect("a returned state certifies");
                assert_eq!(again, report.certificate);
                assert_eq!(report.certificate.n_crossings as u128, t.total_crossings());
                eprintln!(
                    "{:<40} {:<26} {:>5} {:>7} {:>7} {:>9} {:>5} {:>6} {:>6}",
                    job.shape,
                    format!("{:?}", job.scenario),
                    job.seed,
                    t.complex().n_vertices(),
                    t.complex().n_edges(),
                    report.certificate.n_crossings,
                    report.certificate.n_inserted_on_edges,
                    report.n_along_by_endpoints,
                    report.n_arrival_crossings_dropped,
                );
            }
            Err(e) => failures.push(format!(
                "{} / {:?} / seed {}: {e}",
                job.shape, job.scenario, job.seed
            )),
        }
    }
    eprintln!(
        "conversions: {} of {} certified, {} failed",
        n_ok,
        jobs.len(),
        failures.len()
    );
    for f in &failures {
        eprintln!("  FAILED {f}");
    }
    assert!(
        failures.is_empty(),
        "{} conversions failed:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

#[test]
fn the_initial_state_is_the_input_and_converts_to_itself() {
    for shape in all_shapes() {
        let t = IntegerTriangulation::new(&shape.mesh, &shape.geometry).unwrap();
        let report = t
            .certify()
            .unwrap_or_else(|e| panic!("{}: {e}", shape.name));
        assert_eq!(report.n_crossings, 0);
        assert_eq!(report.n_along, shape.mesh.n_edges());
        assert!(t.normal_coordinates().iter().all(|&n| n == -1));
        let (c, _) = convert(&shape.triangulation()).unwrap();
        assert_eq!(
            c.normal_coordinates(),
            t.normal_coordinates(),
            "{}",
            shape.name
        );
        assert_eq!(c.roundabouts(), t.roundabouts(), "{}", shape.name);
        assert_eq!(c.vertex_locations(), t.vertex_locations());
    }
}

/// SPEC §3.1 / Alg. 2 + 3: every crossing of every edge, traced in both directions and
/// identified by the roundabouts (walking back past inserted vertices), names the input edge
/// that the signposts' trace of that edge crosses at that position.
#[test]
fn identified_crossings_match_the_signpost_traces() {
    let mut n_checked = 0u64;
    let mut n_walked_back = 0u64;
    for job in jobs() {
        let Ok((t, _)) = &job.result else { continue };
        let traces = job
            .signposts
            .trace_edges(&EdgeTraceOptions::default())
            .unwrap();
        let c = t.complex();
        for e in 0..c.n_edges() {
            let n = t.normal_coordinate(e);
            if n <= 0 {
                continue;
            }
            let h = c.edge_halfedge(e);
            for p in 0..n {
                let segment = t
                    .extract_curve(Crossing::new(h, p))
                    .unwrap_or_else(|err| panic!("{} {:?}: {err}", job.shape, job.scenario));
                if !t.is_input_vertex(segment.start_vertex as usize) {
                    n_walked_back += 1;
                }
                let id = t
                    .identify_traced_edge(&segment)
                    .unwrap_or_else(|err| panic!("{} {:?}: {err}", job.shape, job.scenario));
                let expected = traces.crossings(e)[p as usize].edge;
                assert_eq!(
                    id.input_edge, expected,
                    "{} {:?} seed {}: crossing {p} of edge {e}",
                    job.shape, job.scenario, job.seed
                );
                n_checked += 1;
            }
        }
    }
    eprintln!(
        "identified {n_checked} crossings, {n_walked_back} by walking back past inserted vertices"
    );
    assert!(n_checked > 10_000);
}

/// Alg. 4: the paths of all input edges use every crossing once, and their pieces join up.
#[test]
fn extracted_input_edges_partition_the_crossings() {
    for job in jobs() {
        let Ok((t, _)) = &job.result else { continue };
        let mesh = t.input().mesh();
        let mut total = 0u128;
        let mut pieces_through_inserted = 0usize;
        for ie in 0..mesh.n_edges() {
            for h in [
                Some(mesh.edge_halfedge(ie)),
                mesh.twin(mesh.edge_halfedge(ie)),
            ]
            .into_iter()
            .flatten()
            {
                let path = t
                    .extract_edge(h)
                    .unwrap_or_else(|err| panic!("{} {:?}: {err}", job.shape, job.scenario));
                assert_eq!(path.vertices[0] as usize, mesh.from_vertex(h));
                assert_eq!(*path.vertices.last().unwrap() as usize, mesh.to_vertex(h));
                assert_eq!(path.pieces.len() + 1, path.vertices.len());
                let crossings: usize = path
                    .pieces
                    .iter()
                    .map(|p| match p {
                        PathPiece::Along { .. } => 0,
                        PathPiece::Crossings(s) => s.crossings.len(),
                    })
                    .sum();
                if h == mesh.edge_halfedge(ie) {
                    total += crossings as u128;
                    pieces_through_inserted += path.vertices.len() - 2;
                }
            }
        }
        assert_eq!(
            total,
            t.total_crossings(),
            "{} {:?}",
            job.shape,
            job.scenario
        );
        let _ = pieces_through_inserted;
    }
}

/// The kinds of mutation of WP1 item 5.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Mutation {
    NormalPlus,
    NormalMinus,
    RoundaboutPlus,
    RoundaboutMinus,
    RoundaboutSwap,
}

/// Every single mutation of `n` (±1) or of a roundabout (±1 modulo `deg₀` at an interior vertex,
/// ±1 at a boundary vertex; a swap of two different roundabouts at one vertex) is rejected, and
/// the error names a place. Tallies which item caught them.
#[test]
fn the_certificate_rejects_every_mutation() {
    let mut tally: BTreeMap<(Mutation, String), usize> = BTreeMap::new();
    let mut n_mutations = 0usize;
    let mut escaped = Vec::new();
    for (j, job) in jobs().iter().enumerate() {
        let Ok((t, _)) = &job.result else { continue };
        let mut rng = SplitMix64::new(0x5eed + j as u64);
        let parts = t.clone().into_parts();
        let ne = parts.normal.len();
        let nh = parts.roundabouts.len();
        let per_kind = 12;
        let mut tries: Vec<(Mutation, Parts)> = Vec::new();
        for _ in 0..per_kind {
            let e = (rng.next_u64() % ne as u64) as usize;
            for (kind, d) in [(Mutation::NormalPlus, 1), (Mutation::NormalMinus, -1)] {
                let mut p = parts.clone();
                p.normal[e] += d;
                tries.push((kind, p));
            }
            // A halfedge leaving an input vertex.
            let h = loop {
                let h = (rng.next_u64() % nh as u64) as usize;
                if parts.roundabouts[h] != NO_ROUNDABOUT {
                    break h;
                }
            };
            let a = parts.complex.from_vertex(h);
            let local = t.local_indices();
            let deg = local.degree(a);
            for (kind, d) in [
                (Mutation::RoundaboutPlus, 1),
                (Mutation::RoundaboutMinus, -1),
            ] {
                let mut p = parts.clone();
                let r = p.roundabouts[h] as i64 + d;
                p.roundabouts[h] = if local.is_boundary(a) {
                    r.max(0) as u32 // 0 − 1 stays out of the way: tested as +1 instead
                } else {
                    r.rem_euclid(deg as i64) as u32
                };
                if p.roundabouts[h] == parts.roundabouts[h] {
                    p.roundabouts[h] += 1;
                }
                tries.push((kind, p));
            }
            // Two halfedges at the same vertex with different roundabouts.
            let outs: Vec<usize> = parts.complex.outgoing_halfedges(a).collect();
            if let Some(&g) = outs
                .iter()
                .find(|&&g| parts.roundabouts[g] != parts.roundabouts[h])
            {
                let mut p = parts.clone();
                p.roundabouts.swap(h, g);
                tries.push((Mutation::RoundaboutSwap, p));
            }
        }
        for (kind, p) in tries {
            n_mutations += 1;
            let caught = match IntegerTriangulation::from_parts(p) {
                Err(e) => Some(e),
                Ok(m) => m.certify().err(),
            };
            match caught {
                Some(e) => {
                    let item = match e.validity_item() {
                        Some(item) => format!("{item:?}"),
                        None => format!("other: {e}"),
                    };
                    assert!(!e.to_string().is_empty());
                    *tally.entry((kind, item)).or_default() += 1;
                }
                None => escaped.push(format!("{} {:?} {:?}", job.shape, job.scenario, kind)),
            }
        }
    }
    eprintln!("{n_mutations} mutations:");
    for ((kind, item), count) in &tally {
        eprintln!("  {kind:?} caught by {item}: {count}");
    }
    assert!(escaped.is_empty(), "not caught: {escaped:#?}");
    let _ = ValidityItem::Ranges;
}
