//! Compares `triangulate_polygons` with triangulations recorded from a real Blender.
//!
//! * `golden_fixture_matches_blender` checks the committed fixture
//!   `tests/data/blender_tessellation_5.0.1.txt` (a small dump made by
//!   `scripts/blender_tessellation_dump.py` with Blender 5.0.1) and needs no Blender.
//! * `external_dump_agreement` reads a larger dump named by the environment variable
//!   `DDG_BLENDER_DUMP` and prints the agreement per polygon family; it is a no-op when the
//!   variable is unset. `scripts/blender-tessellation-diff.sh` produces the dump and runs it.

#![cfg(feature = "blender-tessellation")]

use std::collections::BTreeMap;

use ddg_core::triangulate::triangulate_polygons;

struct Group {
    name: String,
    /// "exact" or "degenerate".
    expect: String,
    positions: Vec<[f64; 3]>,
    corner_verts: Vec<u32>,
    face_starts: Vec<u32>,
    /// Blender's loop triangles.
    triangles: Vec<[u32; 3]>,
}

fn parse_dump(text: &str) -> Vec<Group> {
    let mut lines = text
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty());
    let mut groups = Vec::new();
    while let Some(header) = lines.next() {
        let mut h = header.split_whitespace();
        assert_eq!(h.next(), Some("group"), "bad header {header:?}");
        let name = h.next().unwrap().to_string();
        let expect = h.next().unwrap().to_string();
        let mut count = |key: &str| -> usize {
            let line = lines.next().unwrap();
            let (k, v) = line.split_once(' ').unwrap();
            assert_eq!(k, key);
            v.parse().unwrap()
        };
        let nv = count("nv");
        let positions = (0..nv)
            .map(|_| {
                let mut it = lines
                    .next()
                    .unwrap()
                    .split_whitespace()
                    .map(|t| t.parse::<f64>().unwrap());
                [it.next().unwrap(), it.next().unwrap(), it.next().unwrap()]
            })
            .collect();
        let mut count = |key: &str| -> usize {
            let line = lines.next().unwrap();
            let (k, v) = line.split_once(' ').unwrap();
            assert_eq!(k, key);
            v.parse().unwrap()
        };
        let np = count("np");
        let mut corner_verts = Vec::new();
        let mut face_starts = Vec::new();
        for _ in 0..np {
            let mut it = lines
                .next()
                .unwrap()
                .split_whitespace()
                .map(|t| t.parse::<u32>().unwrap());
            let k = it.next().unwrap();
            face_starts.push(corner_verts.len() as u32);
            corner_verts.extend(it);
            assert_eq!(corner_verts.len() as u32 - face_starts.last().unwrap(), k);
        }
        let mut count = |key: &str| -> usize {
            let line = lines.next().unwrap();
            let (k, v) = line.split_once(' ').unwrap();
            assert_eq!(k, key);
            v.parse().unwrap()
        };
        let nt = count("nt");
        let triangles = (0..nt)
            .map(|_| {
                let mut it = lines
                    .next()
                    .unwrap()
                    .split_whitespace()
                    .map(|t| t.parse::<u32>().unwrap());
                [it.next().unwrap(), it.next().unwrap(), it.next().unwrap()]
            })
            .collect();
        assert_eq!(lines.next(), Some("end"));
        groups.push(Group {
            name,
            expect,
            positions,
            corner_verts,
            face_starts,
            triangles,
        });
    }
    groups
}

#[derive(Default, Clone, Copy)]
struct Tally {
    polygons: usize,
    /// Same triangles, same order, same corner order.
    identical: usize,
    /// Same set of triangles (as corner sets), but listed or rotated differently.
    same_triangles: usize,
    /// A different triangulation.
    different: usize,
}

impl Tally {
    fn add(&mut self, other: Tally) {
        self.polygons += other.polygons;
        self.identical += other.identical;
        self.same_triangles += other.same_triangles;
        self.different += other.different;
    }
}

fn compare(group: &Group) -> Tally {
    let ours = triangulate_polygons(&group.positions, &group.corner_verts, &group.face_starts)
        .unwrap_or_else(|e| panic!("group {}: {e}", group.name));
    assert_eq!(
        ours.triangle_corners.len(),
        group.triangles.len(),
        "group {}",
        group.name
    );
    let mut tally = Tally::default();
    for (p, &start) in group.face_starts.iter().enumerate() {
        let end = group
            .face_starts
            .get(p + 1)
            .map_or(group.corner_verts.len(), |&e| e as usize);
        let first = start as usize - 2 * p;
        let count = end - start as usize - 2;
        let a = &ours.triangle_corners[first..first + count];
        let b = &group.triangles[first..first + count];
        tally.polygons += 1;
        if a == b {
            tally.identical += 1;
            continue;
        }
        let canonical = |tris: &[[u32; 3]]| {
            let mut v: Vec<[u32; 3]> = tris
                .iter()
                .map(|t| {
                    let mut s = *t;
                    s.sort_unstable();
                    s
                })
                .collect();
            v.sort_unstable();
            v
        };
        if canonical(a) == canonical(b) {
            tally.same_triangles += 1;
        } else {
            tally.different += 1;
        }
    }
    tally
}

/// Compares every group and returns the totals for the "exact" and the "degenerate" families.
fn report(groups: &[Group], print: bool) -> (Tally, Tally) {
    // Merge "<name>_alt" (same polygons, other normal-cache state in Blender) into "<name>".
    let mut by_family: BTreeMap<(String, String), Tally> = BTreeMap::new();
    for g in groups {
        let family = g.name.strip_suffix("_alt").unwrap_or(&g.name).to_string();
        by_family
            .entry((g.expect.clone(), family))
            .or_default()
            .add(compare(g));
    }
    let mut exact = Tally::default();
    let mut degenerate = Tally::default();
    if print {
        println!(
            "{:<11} {:<26} {:>9} {:>10} {:>12} {:>10} {:>11}",
            "expectation",
            "family",
            "polygons",
            "identical",
            "same tris",
            "different",
            "identical %"
        );
    }
    for ((expect, family), tally) in &by_family {
        if print {
            println!(
                "{:<11} {:<26} {:>9} {:>10} {:>12} {:>10} {:>10.3}%",
                expect,
                family,
                tally.polygons,
                tally.identical,
                tally.same_triangles,
                tally.different,
                100.0 * tally.identical as f64 / tally.polygons as f64
            );
        }
        if expect == "exact" {
            exact.add(*tally)
        } else {
            degenerate.add(*tally)
        }
    }
    if print {
        for (label, t) in [("exact", exact), ("degenerate", degenerate)] {
            println!(
                "{:<11} {:<26} {:>9} {:>10} {:>12} {:>10} {:>10.3}%",
                label,
                "TOTAL",
                t.polygons,
                t.identical,
                t.same_triangles,
                t.different,
                100.0 * t.identical as f64 / t.polygons.max(1) as f64
            );
        }
    }
    (exact, degenerate)
}

#[test]
fn golden_fixture_matches_blender() {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/data/blender_tessellation_5.0.1.txt"
    );
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let groups = parse_dump(&text);
    let (exact, degenerate) = report(&groups, false);
    assert!(
        exact.polygons > 300,
        "fixture is unexpectedly small: {} polygons",
        exact.polygons
    );
    assert_eq!(
        exact.identical,
        exact.polygons,
        "non-degenerate polygons must be triangulated exactly like Blender 5.0.1 \
         ({} of {} differ)",
        exact.polygons - exact.identical,
        exact.polygons
    );
    assert_eq!(
        degenerate.identical,
        degenerate.polygons,
        "degenerate polygons of the fixture matched Blender 5.0.1 when it was recorded \
         ({} of {} differ now)",
        degenerate.polygons - degenerate.identical,
        degenerate.polygons
    );
}

#[test]
fn external_dump_agreement() {
    let Ok(path) = std::env::var("DDG_BLENDER_DUMP") else {
        eprintln!(
            "SKIPPED: DDG_BLENDER_DUMP is not set (see scripts/blender-tessellation-diff.sh)"
        );
        return;
    };
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
    for line in text.lines().take_while(|l| l.starts_with('#')) {
        println!("{line}");
    }
    let groups = parse_dump(&text);
    let (exact, _degenerate) = report(&groups, true);
    assert_eq!(
        exact.identical,
        exact.polygons,
        "{} of {} non-degenerate polygons are triangulated differently from Blender",
        exact.polygons - exact.identical,
        exact.polygons
    );
}
