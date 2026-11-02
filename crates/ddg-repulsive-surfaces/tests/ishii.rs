//! Fig. 16: the genus-two handlebody-knots of Ishii et al. 2012 (`generate::ishii`). Run with `--release`.

use ddg_repulsive_surfaces::ccd;
use ddg_repulsive_surfaces::generate::ishii::{self, GraphKind};

#[test]
fn all_21_entries_are_transcribed_with_the_printed_crossing_numbers() {
    let all = ishii::all_entries();
    assert_eq!(all.len(), 21);
    for e in &all {
        let expected = e.name.chars().next().unwrap().to_digit(10).unwrap() as usize;
        assert_eq!(e.crossings, expected, "{}", e.name);
        assert_eq!(e.edges.len(), 3, "{}", e.name);
        // projection reproduces the count and the over/under sequence along every edge
        ishii::check_projection(e).unwrap();
        // every crossing has a clear vertical separation
        for c in ishii::projection_crossings(e) {
            let z = |(k, s): (usize, f64)| {
                let p = &e.edges[k].points;
                let i = s.floor() as usize;
                p[i][2] + (s - i as f64) * (p[i + 1][2] - p[i][2])
            };
            assert!((z(c.a) - z(c.b)).abs() > 0.35, "{}: crossings too close in z", e.name);
        }
        // handcuff graph: a loop at each vertex and a bridge; theta: three vertex-to-vertex edges
        let loops = e.edges.iter().filter(|x| x.from == x.to).count();
        match e.kind {
            GraphKind::Theta => assert_eq!(loops, 0, "{}", e.name),
            GraphKind::Handcuff => assert_eq!(loops, 2, "{}", e.name),
        }
    }
}

#[test]
fn mirror_image_fails_the_diagram_check() {
    // swapping over and under at every crossing must be detected (the check is not vacuous)
    for name in ["4_1", "6_12", "6_13"] {
        let mut e = ishii::entry(name).unwrap();
        ishii::check_projection(&e).unwrap();
        for edge in &mut e.edges {
            for p in &mut edge.points {
                p[2] = -p[2];
            }
        }
        assert!(ishii::check_projection(&e).is_err(), "{name}");
    }
    // 6_12 and 6_13 have the same projection but different crossing signs
    let (a, b) = (ishii::entry("6_12").unwrap(), ishii::entry("6_13").unwrap());
    assert_ne!(a.edges.iter().map(|e| e.recorded.clone()).collect::<Vec<_>>(), b.edges.iter().map(|e| e.recorded.clone()).collect::<Vec<_>>());
}

#[test]
fn every_entry_meshes_as_an_embedded_genus_two_tube() {
    // 14 000 faces is the smallest size at which every tube is resolved (checked 10 000, 12 000)
    let names: Vec<&str> = ishii::NAMES.to_vec();
    std::thread::scope(|s| {
        let hs: Vec<_> = names
            .iter()
            .map(|n| {
                s.spawn(move || {
                    let e = ishii::entry(n).unwrap();
                    let (mesh, radius) = ishii::tube(&e, 14_000).unwrap_or_else(|x| panic!("{n}: {x}"));
                    let surf = mesh.surface().unwrap();
                    assert_eq!(surf.n_components(), 1, "{n}");
                    assert_eq!(1 - surf.euler_characteristic() / 2, 2, "{n}: genus");
                    assert!(!ccd::intersects(&surf, &mesh.positions).unwrap(), "{n}: intersects");
                    assert!(radius > 0.05 && radius <= 0.08, "{n}: radius {radius}");
                })
            })
            .collect();
        for h in hs {
            h.join().unwrap();
        }
    });
}
