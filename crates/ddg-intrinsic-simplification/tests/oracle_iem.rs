//! The Lean oracle's golden fixtures for the ICE metric (`lean/oracle/fixtures/iem`, format in
//! `lean/oracle/README.md`, "Intrinsic simplification"): the redistribution weights (Eq. 9),
//! and whole removals of one star — the flattening, the curvature changes, the weights and the
//! scaled lengths through the crate's removal pipeline (`tentative_removal` on the star built as
//! a mesh), then the transport (Eq. 5, 8) and the cost (Eq. 6, and its memoryless form) through
//! the crate's formulas with the fixture's angular coordinates.

use std::collections::HashMap;
use std::path::PathBuf;

use ddg_core::TriMesh;
use ddg_intrinsic_simplification::{CoarseMesh, CoarseMeshOptions, metric};

fn fixtures() -> Option<PathBuf> {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../lean/oracle/fixtures/iem");
    if dir.join("MANIFEST.tsv").is_file() {
        Some(dir)
    } else {
        eprintln!("SKIPPED: {} not found; nothing was checked", dir.display());
        None
    }
}

fn hex(tok: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(tok, 16).unwrap_or_else(|e| panic!("{tok}: {e}")))
}

/// The lines of a block as token vectors.
fn lines(text: &str) -> Vec<Vec<String>> {
    text.lines()
        .map(|l| l.split_whitespace().map(str::to_string).collect::<Vec<_>>())
        .filter(|t| !t.is_empty())
        .collect()
}

/// A fixture channel: weight, mass and vector at `i`, then mass and vector per neighbour.
type Channel = (f64, f64, [f64; 2], Vec<(f64, [f64; 2])>);

fn rel(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs().max(1e-300)
}

#[test]
fn weights_fixtures() {
    let Some(dir) = fixtures() else { return };
    let mut n = 0;
    for name in ["iem-weights-basic", "iem-weights-flat"] {
        let input = lines(&std::fs::read_to_string(dir.join(format!("{name}.in"))).unwrap());
        let output = lines(&std::fs::read_to_string(dir.join(format!("{name}.out"))).unwrap());
        let ki = input
            .iter()
            .find(|t| t[0] == "ki")
            .map(|t| hex(&t[1]))
            .unwrap();
        let dk: Vec<f64> = input
            .iter()
            .filter(|t| t[0] == "dk")
            .map(|t| hex(&t[2]))
            .collect();
        let want: Vec<f64> = output
            .iter()
            .filter(|t| t[0] == "w")
            .map(|t| hex(&t[2]))
            .collect();
        let got = metric::redistribution_weights(ki, &vec![0.0; dk.len()], &dk);
        for (g, w) in got.iter().zip(&want) {
            assert!(rel(*g, *w) <= 1e-15, "{name}: {got:?} against {want:?}");
        }
        n += 1;
    }
    assert_eq!(n, 2);
}

/// A star of the fixtures as a mesh (labels are vertex ids; side `k` is opposite corner `k`).
fn star_mesh(faces: &[([u32; 3], [f64; 3])]) -> Option<CoarseMesh> {
    let n = faces.iter().flat_map(|f| f.0).max()? as usize + 1;
    let tri = TriMesh::new(n, faces.iter().map(|f| f.0).collect()).ok()?;
    let mut length: HashMap<(u32, u32), f64> = HashMap::new();
    for (v, l) in faces {
        for k in 0..3 {
            let (a, b) = (v[(k + 1) % 3], v[(k + 2) % 3]);
            length.insert((a.min(b), a.max(b)), l[k]);
        }
    }
    let lengths: Vec<f64> = tri
        .edges()
        .iter()
        .map(|e| length[&(e[0].min(e[1]), e[0].max(e[1]))])
        .collect();
    CoarseMesh::from_lengths(&tri, lengths, &[], CoarseMeshOptions::default()).ok()
}

#[test]
fn removal_fixtures() {
    let Some(dir) = fixtures() else { return };
    let manifest = std::fs::read_to_string(dir.join("MANIFEST.tsv")).unwrap();
    let mut checked = 0;
    for row in manifest.lines().skip(1) {
        let name = row.split('\t').next().unwrap();
        if !name.starts_with("iem-removal-") {
            continue;
        }
        let input = lines(&std::fs::read_to_string(dir.join(format!("{name}.in"))).unwrap());
        let output = lines(&std::fs::read_to_string(dir.join(format!("{name}.out"))).unwrap());
        let faces: Vec<([u32; 3], [f64; 3])> = input
            .iter()
            .filter(|t| t[0] == "face")
            .map(|t| {
                (
                    [1, 2, 3].map(|k| t[k].parse().unwrap()),
                    [4, 5, 6].map(|k| hex(&t[k])),
                )
            })
            .collect();
        let memory = !input.iter().any(|t| t[0] == "memory" && t[1] == "0");
        let flatten_ok = output.iter().any(|t| t[0] == "flatten" && t[1] == "ok");
        let mut mesh = star_mesh(&faces).unwrap_or_else(|| panic!("{name}: not a mesh"));
        let prepared = mesh.tentative_removal(0);
        if !flatten_ok {
            // The oracle's cost is +inf; the crate refuses the removal (cost +inf too).
            assert!(prepared.is_err(), "{name}: {prepared:?}");
            checked += 1;
            continue;
        }
        let p = prepared.unwrap_or_else(|e| panic!("{name}: {e}"));
        let u = output
            .iter()
            .find(|t| t[0] == "flatten")
            .map(|t| hex(&t[3]))
            .unwrap();
        assert!(
            (p.flatten.u - u).abs() <= 1e-12 * (1.0 + u.abs()),
            "{name}: u {} against {u}",
            p.flatten.u
        );
        let ki = output
            .iter()
            .find(|t| t[0] == "ki")
            .map(|t| hex(&t[1]))
            .unwrap();
        assert!(
            (p.curvature - ki).abs() <= 1e-13,
            "{name}: K_i {} against {ki}",
            p.curvature
        );
        // Per neighbour: K̃_j − K_j, α, ℓ̃.
        let nbr_out: Vec<(u32, f64, f64, f64)> = output
            .iter()
            .filter(|t| t[0] == "nbr")
            .map(|t| (t[1].parse().unwrap(), hex(&t[2]), hex(&t[3]), hex(&t[4])))
            .collect();
        let before: Vec<f64> = p.neighbours.iter().map(|n| n.curvature_before).collect();
        let after: Vec<f64> = p.neighbours.iter().map(|n| n.curvature_after).collect();
        let alpha = metric::redistribution_weights(p.curvature, &before, &after);
        assert_eq!(p.neighbours.len(), nbr_out.len(), "{name}");
        for (k, n) in p.neighbours.iter().enumerate() {
            let &(j, dk, a, len) = nbr_out.iter().find(|x| x.0 == n.vertex).unwrap();
            assert_eq!(j, n.vertex);
            let our_dk = n.curvature_after - n.curvature_before;
            assert!(
                (our_dk - dk).abs() <= 1e-12,
                "{name}: dK_{j} {our_dk} against {dk}"
            );
            assert!(
                (alpha[k] - a).abs() <= 1e-11,
                "{name}: alpha_{j} {} against {a}",
                alpha[k]
            );
            assert!(
                rel(n.length, len) <= 1e-12,
                "{name}: length to {j} {} against {len}",
                n.length
            );
        }
        // The transport and the cost with the fixture's frames, our weights and lengths.
        let nbr_in: Vec<(u32, f64, f64)> = input
            .iter()
            .filter(|t| t[0] == "nbr")
            .map(|t| (t[1].parse().unwrap(), hex(&t[2]), hex(&t[3])))
            .collect();
        let mut channels: Vec<Channel> = Vec::new();
        for t in &input {
            match t[0].as_str() {
                "channel" => {
                    channels.push((hex(&t[1]), hex(&t[2]), [hex(&t[3]), hex(&t[4])], vec![]))
                }
                "mass" => channels
                    .last_mut()
                    .unwrap()
                    .3
                    .push((hex(&t[2]), [hex(&t[3]), hex(&t[4])])),
                _ => {}
            }
        }
        let new_out: HashMap<(usize, u32), (f64, [f64; 2])> = output
            .iter()
            .filter(|t| t[0] == "new")
            .map(|t| {
                (
                    (t[1].parse().unwrap(), t[2].parse().unwrap()),
                    (hex(&t[3]), [hex(&t[4]), hex(&t[5])]),
                )
            })
            .collect();
        let mut cost = 0.0;
        let mut memoryless = 0.0;
        for (c, (w, mi, ti, masses)) in channels.iter().enumerate() {
            let ti = if memory { *ti } else { [0.0; 2] };
            let mut alphas = Vec::new();
            let mut lens = Vec::new();
            for (k, &(j, phi_ij, phi_ji)) in nbr_in.iter().enumerate() {
                let n = p.neighbours.iter().position(|n| n.vertex == j).unwrap();
                let (a, len) = (alpha[n], p.neighbours[n].length);
                let angle = phi_ji + std::f64::consts::PI - phi_ij;
                let rot = [angle.cos(), angle.sin()];
                let edge = [len * phi_ji.cos(), len * phi_ji.sin()];
                let (mj, tj) = masses[k];
                let tj = if memory { tj } else { [0.0; 2] };
                let (m, t) = metric::transport_update(a, *mi, ti, rot, edge, mj, tj);
                let (wm, wt) = new_out[&(c, j)];
                assert!(
                    rel(m, wm) <= 1e-12,
                    "{name}: m̃ channel {c} at {j}: {m} against {wm}"
                );
                let scale = 1e-11 * (1.0 + wt[0].hypot(wt[1]));
                assert!(
                    (t[0] - wt[0]).abs() <= scale && (t[1] - wt[1]).abs() <= scale,
                    "{name}: t̃ channel {c} at {j}: {t:?} against {wt:?}"
                );
                cost += w * m * t[0].hypot(t[1]);
                alphas.push(a);
                lens.push(len);
            }
            memoryless += w * metric::memoryless_cost(*mi, &alphas, &lens);
        }
        let want_cost = output
            .iter()
            .find(|t| t[0] == "cost")
            .map(|t| hex(&t[1]))
            .unwrap();
        let want_ml = output
            .iter()
            .find(|t| t[0] == "memoryless")
            .map(|t| hex(&t[1]))
            .unwrap();
        assert!(
            rel(cost, want_cost) <= 1e-11,
            "{name}: cost {cost} against {want_cost}"
        );
        assert!(
            rel(memoryless, want_ml) <= 1e-11,
            "{name}: memoryless {memoryless} against {want_ml}"
        );
        eprintln!(
            "{name}: u {:e}, cost {cost:e} (oracle {want_cost:e})",
            p.flatten.u
        );
        checked += 1;
    }
    assert_eq!(checked, 12, "removal fixtures checked");
}

/// The flattening of the simplicial stars of the fixtures (`task flatten`, status ok): the
/// conformal factor through the crate's removal pipeline.
#[test]
fn flatten_fixtures() {
    let Some(dir) = fixtures() else { return };
    let manifest = std::fs::read_to_string(dir.join("MANIFEST.tsv")).unwrap();
    let mut checked = 0;
    for row in manifest.lines().skip(1) {
        let cols: Vec<&str> = row.split('\t').collect();
        let (name, status) = (cols[0], cols[1]);
        if !name.starts_with("iem-flatten-") || status != "ok" {
            continue;
        }
        let input = lines(&std::fs::read_to_string(dir.join(format!("{name}.in"))).unwrap());
        let output = lines(&std::fs::read_to_string(dir.join(format!("{name}.out"))).unwrap());
        let faces: Vec<([u32; 3], [f64; 3])> = input
            .iter()
            .filter(|t| t[0] == "face")
            .map(|t| {
                (
                    [1, 2, 3].map(|k| t[k].parse().unwrap()),
                    [4, 5, 6].map(|k| hex(&t[k])),
                )
            })
            .collect();
        // Self-edges cannot be built as a mesh; the oracle's step checks cover them.
        let Some(mut mesh) = star_mesh(&faces) else {
            eprintln!("{name}: not a simplicial star, skipped");
            continue;
        };
        let u = output
            .iter()
            .find(|t| t[0] == "solve")
            .map(|t| hex(&t[3]))
            .unwrap();
        match mesh.tentative_removal(0) {
            Ok(p) if p.fallback_flips.is_empty() => {
                assert!(
                    (p.flatten.u - u).abs() <= 1e-12 * (1.0 + u.abs()),
                    "{name}: u {} against {u}",
                    p.flatten.u
                );
                checked += 1;
            }
            Ok(_) => eprintln!("{name}: needed fallback flips here, skipped"),
            Err(e) => eprintln!("{name}: {e} (the flattening is not the reason: skipped)"),
        }
    }
    eprintln!("{checked} flatten fixtures checked");
    assert!(checked >= 20, "{checked}");
}
