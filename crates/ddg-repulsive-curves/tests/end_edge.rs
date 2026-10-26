//! End-edge collapse of the authors' `TangentInterpolation` scene and the guard against it
//! (WP3c, `SPEC.md` §17j).

use ddg_repulsive_curves::scene::{BackendChoice, Scene, read_scene};
use ddg_repulsive_curves::{FlowOptions, RcError, Session, Status};

fn scene() -> Scene {
    read_scene("scenes/TangentInterpolation/scene.txt").unwrap()
}

fn session(scene: &Scene, tweak: impl Fn(&mut FlowOptions)) -> Session {
    let mut o = scene.options(BackendChoice::Exact, Some(4)).unwrap();
    tweak(&mut o);
    Session::new(scene.network.clone(), scene.positions.clone(), o).unwrap()
}

fn min_edge(s: &Session) -> f64 {
    s.network()
        .edge_geometry(s.positions())
        .length
        .iter()
        .cloned()
        .fold(f64::MAX, f64::min)
}

/// The scene has no length constraint (penalty only) and a fixed tangent at both ends: an end
/// edge shrinks by about a factor 0.6 per step. By default the step that would take it below
/// 1e-3 of the mean edge length fails with the named error, the curve stays that of the last
/// accepted step, and the metric never became singular.
#[test]
fn the_default_guard_names_the_collapsing_edge() {
    let sc = scene();
    let mut s = session(&sc, |_| {});
    let l0 = s.mean_edge_length0();
    let mut last_ok = None;
    let mut error = None;
    for _ in 0..80 {
        match s.step() {
            Ok(_) => last_ok = Some((s.positions().to_vec(), min_edge(&s))),
            Err(e) => {
                error = Some(e);
                break;
            }
        }
    }
    let e = error.expect("the scene collapses within 80 steps");
    let RcError::EdgeCollapse {
        edge,
        length,
        length_before,
        threshold,
        step,
    } = e.clone()
    else {
        panic!("expected EdgeCollapse, got {e}");
    };
    let (kept, kept_min) = last_ok.unwrap();
    assert!((threshold - 1e-3 * l0).abs() < 1e-15 * l0);
    assert!(length < threshold && length < length_before);
    assert!(
        edge == 0 || edge == sc.network.n_edges() - 1,
        "an end edge, got {edge}"
    );
    assert!(step > 5 && step < 60, "step {step}");
    // Not taken: the curve is that of the last accepted step, and its edges are above the bound.
    assert!(s.positions() == kept.as_slice());
    assert!(kept_min >= threshold, "{kept_min} < {threshold}");
    assert_eq!(s.status(), Status::Error(e.clone()));
    let msg = e.to_string();
    assert!(
        msg.contains("edge") && msg.contains("min_edge_length_ratio"),
        "{msg}"
    );
}

/// Without the guard the same run ends later with the metric singular: the old behaviour.
#[test]
fn without_the_guard_the_metric_becomes_singular() {
    let sc = scene();
    let mut s = session(&sc, |o| o.min_edge_length_ratio = None);
    let mut error = None;
    for _ in 0..100 {
        if let Err(e) = s.step() {
            error = Some(e);
            break;
        }
    }
    assert!(
        matches!(error, Some(RcError::SingularMetric { .. })),
        "{error:?}"
    );
}

/// Fixed edge lengths never trigger the guard, and the option is validated.
#[test]
fn guard_is_quiet_with_fixed_lengths_and_validates_its_option() {
    let n = 40usize;
    let pos: Vec<[f64; 3]> = (0..n)
        .map(|k| {
            let t = std::f64::consts::TAU * k as f64 / n as f64;
            [
                t.sin() + 2.0 * (2.0 * t).sin(),
                t.cos() - 2.0 * (2.0 * t).cos(),
                -(3.0 * t).sin(),
            ]
        })
        .collect();
    let edges = (0..n as u32).map(|k| [k, (k + 1) % n as u32]).collect();
    let net = ddg_repulsive_curves::CurveNetwork::new(n, edges).unwrap();
    let mut s = Session::new(
        net,
        pos,
        FlowOptions {
            threads: Some(2),
            ..FlowOptions::default()
        },
    )
    .unwrap();
    for _ in 0..30 {
        s.step().unwrap();
    }
    assert!(min_edge(&s) > 0.3 * s.mean_edge_length0());
    for bad in [0.0, 1.0, -1.0, f64::NAN, f64::INFINITY] {
        let o = FlowOptions {
            min_edge_length_ratio: Some(bad),
            ..FlowOptions::default()
        };
        assert!(
            matches!(
                o.validate(),
                Err(RcError::InvalidOption {
                    name: "min_edge_length_ratio",
                    ..
                })
            ),
            "{bad}"
        );
    }
    assert_eq!(FlowOptions::default().min_edge_length_ratio, Some(1e-3));
}
