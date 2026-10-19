//! Bad input gives errors, never panics (SPEC §3.1, §3.6, §3.9, §7.1).

mod energy_common;
use energy_common::*;

use ddg_core::TriMesh;
use ddg_developability::{DevError, Energy, EnergyOptions, evaluate, evaluate_energy};

#[test]
fn wrong_length_and_non_finite_positions() {
    let (m, p) = interior_star(5, 1, 0.3);
    for energy in Energy::ALL {
        let o = EnergyOptions::new(energy);
        let err = evaluate(&m, &p[..4], &o).unwrap_err();
        assert_eq!(
            err,
            DevError::LengthMismatch {
                what: "positions",
                expected: 6,
                got: 4
            }
        );
        assert!(err.to_string().contains("expected 6"));
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut q = p.clone();
            q[3][1] = bad;
            q[5][0] = bad;
            assert_eq!(
                evaluate(&m, &q, &o).unwrap_err(),
                DevError::NonFinitePosition { vertex: 3 }
            );
            assert_eq!(
                evaluate_energy(&m, &q, &o).unwrap_err(),
                DevError::NonFinitePosition { vertex: 3 }
            );
        }
    }
}

#[test]
fn degenerate_faces() {
    let (m, p) = interior_star(6, 2, 0.3);
    for energy in Energy::ALL {
        let o = EnergyOptions::new(energy);
        // collinear corners: vertex 4 on the segment from 0 to 3 makes face 2 = (0, 3, 4) flat
        let mut q = p.clone();
        q[4] = scale(add(q[0], q[3]), 0.5);
        let err = evaluate(&m, &q, &o).unwrap_err();
        assert_eq!(err, DevError::DegenerateFace { face: 2 });
        assert!(err.to_string().contains("face 2"));
        // coincident vertices
        let mut q = p.clone();
        q[2] = q[1];
        assert_eq!(
            evaluate(&m, &q, &o).unwrap_err(),
            DevError::DegenerateFace { face: 0 }
        );
        // extreme scales: the normal under/overflows
        let tiny: Vec<_> = p.iter().map(|&x| scale(x, 1e-170)).collect();
        assert!(matches!(
            evaluate(&m, &tiny, &o),
            Err(DevError::DegenerateFace { .. })
        ));
        let huge: Vec<_> = p.iter().map(|&x| scale(x, 1e170)).collect();
        assert!(matches!(
            evaluate(&m, &huge, &o),
            Err(DevError::DegenerateFace { .. })
        ));
    }
}

#[test]
fn invalid_threads_option() {
    let (m, p) = interior_star(5, 1, 0.3);
    let o = EnergyOptions {
        threads: Some(0),
        ..Default::default()
    };
    assert!(matches!(
        o.validate(),
        Err(DevError::InvalidOption {
            name: "threads",
            ..
        })
    ));
    assert!(matches!(
        evaluate(&m, &p, &o),
        Err(DevError::InvalidOption {
            name: "threads",
            ..
        })
    ));
    assert!(EnergyOptions::default().validate().is_ok());
    let d = EnergyOptions::default();
    assert_eq!(
        (
            d.energy,
            d.include_boundary,
            d.omit_valence_three,
            d.threads
        ),
        (Energy::Covariance, true, true, None)
    );
}

/// Valence-3 star with zero vector area at the centre (the ring is collinear); boundary vertices
/// 1 and 3 have antiparallel face normals. The intrinsic energies are undefined at 0, 1, 3.
fn zero_vector_area_star() -> Mesh {
    let p = vec![
        [0.0, 0.0, 0.0],
        [-1.0, 1.0, -1.0],
        [0.0, 1.0, 0.0],
        [1.0, 1.0, 1.0],
    ];
    let m = TriMesh::new(4, vec![[0, 1, 2], [0, 2, 3], [0, 3, 1]]).unwrap();
    (m, p)
}

#[test]
fn undefined_intrinsic_energy() {
    let (m, p) = zero_vector_area_star();
    for energy in [Energy::Intrinsic, Energy::IntrinsicMax] {
        // all vertices uncounted: NaN entries, no error
        let o = EnergyOptions {
            energy,
            include_boundary: false,
            omit_valence_three: true,
            threads: None,
        };
        let e = evaluate(&m, &p, &o).unwrap();
        assert!(e.counted.iter().all(|&c| !c));
        // vertex 2 has two equal normals and is defined; 0, 1, 3 are not
        for v in [0, 1, 3] {
            assert!(
                e.vertex_energy[v].is_nan() && e.vertex_margin[v].is_nan(),
                "vertex {v}"
            );
        }
        assert!(e.vertex_energy[2].abs() < 1e-15);
        assert_eq!(e.energy, 0.0);
        assert!(e.gradient.iter().all(|g| *g == [0.0; 3]));
        // centre counted: zero vector area
        let o2 = EnergyOptions {
            omit_valence_three: false,
            ..o.clone()
        };
        match evaluate(&m, &p, &o2).unwrap_err() {
            DevError::UndefinedEnergy { vertex: 0, reason } => {
                assert!(reason.contains("vector area"))
            }
            other => panic!("{other:?}"),
        }
        // boundary counted: antiparallel normals at vertex 1
        let o3 = EnergyOptions {
            include_boundary: true,
            ..o.clone()
        };
        match evaluate_energy(&m, &p, &o3).unwrap_err() {
            e @ DevError::UndefinedEnergy { vertex: 1, reason } => {
                assert!(reason.contains("antipode"));
                assert!(e.to_string().contains("vertex 1"));
            }
            other => panic!("{other:?}"),
        }
    }
    // the extrinsic energies are defined there
    for energy in [
        Energy::Combinatorial,
        Energy::CombinatorialMax,
        Energy::Covariance,
        Energy::CovarianceMax,
    ] {
        let o = EnergyOptions {
            energy,
            include_boundary: true,
            omit_valence_three: false,
            threads: None,
        };
        let e = evaluate(&m, &p, &o).unwrap();
        assert!(e.vertex_energy.iter().all(|x| x.is_finite()));
    }
}

#[test]
fn random_garbage_never_panics() {
    let (m, p) = perturbed_icosphere(1, 0.1, 3);
    let mut r = ddg_core::generators::SplitMix64::new(99);
    for trial in 0..200 {
        let mut q = p.clone();
        // collapse random vertices onto others, fold, or scale wildly
        for _ in 0..(trial % 5) {
            let a = (r.next_u64() % q.len() as u64) as usize;
            let b = (r.next_u64() % q.len() as u64) as usize;
            q[a] = match trial % 3 {
                0 => q[b],
                1 => scale(q[b], -1.0),
                _ => [r.uniform(-1.0, 1.0) * 1e150, 0.0, r.uniform(-1.0, 1.0)],
            };
        }
        for energy in Energy::ALL {
            let o = EnergyOptions {
                energy,
                include_boundary: true,
                omit_valence_three: false,
                threads: Some(2),
            };
            if let Ok(e) = evaluate(&m, &q, &o) {
                assert!(
                    e.energy.is_finite(),
                    "trial {trial} {energy:?}: energy {}",
                    e.energy
                );
                assert!(
                    e.gradient.iter().all(|g| g.iter().all(|x| x.is_finite())),
                    "trial {trial} {energy:?}"
                );
            }
        }
    }
}
