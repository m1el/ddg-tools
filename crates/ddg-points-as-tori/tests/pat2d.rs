//! The 2D variant (SPEC.md §8.7): classical circle fit, 2D blend, grid, marching squares, output.

use ddg_points_as_tori::analytic::{Kind, Shape2d};
use ddg_points_as_tori::generate::d2::Instance2d;
use ddg_points_as_tori::generate::rng::Rng;
use ddg_points_as_tori::pat2d::*;
use ddg_points_as_tori::{GradientMode, LambdaRule, Neighborhood, NormalPolicy, PatError};

fn shape(kind: Kind, p: &[f64]) -> Shape2d {
    Shape2d::new(kind, p).expect("shape")
}

fn cloud_of(s: &Shape2d, n: usize, seed: u64) -> Cloud2 {
    let inst = Instance2d::Analytic(s.clone());
    let c = inst.sample_uniform(n, &mut Rng::new(seed)).expect("sample");
    Cloud2::new(c.positions, c.normals, NormalPolicy::Normalize).expect("cloud").0
}

fn build(cloud: &Cloud2, fit: &ClassicalFit, threads: usize, eval: Eval2dOptions) -> (Pat2dField, BuildReport2) {
    let mut o = BuildOptions2 { eval, ..Default::default() }.with_threads(Some(threads));
    o.eval.threads = Some(threads);
    Pat2d::build(cloud, fit, &o).expect("build")
}

fn padded(b: [[f64; 2]; 2], f: f64) -> [[f64; 2]; 2] {
    let w = (b[1][0] - b[0][0]).max(b[1][1] - b[0][1]);
    [[b[0][0] - f * w, b[0][1] - f * w], [b[1][0] + f * w, b[1][1] + f * w]]
}

/// Mean and max |φ − exact| on a grid over the padded box; mean also within `band` of the surface.
fn grid_error(s: &Shape2d, f: &Pat2dField, res: usize) -> (f64, f64, f64) {
    let g = evaluate_grid(f, padded(s.bbox(), 0.25), res, res).expect("grid");
    let nodes = g.nodes();
    let (mut sum, mut mx, mut bsum, mut bn) = (0.0, 0.0f64, 0.0, 0usize);
    for (p, v) in nodes.iter().zip(&g.values) {
        let e = (v - s.distance(*p)).abs();
        sum += e;
        mx = mx.max(e);
        if s.distance(*p).abs() < 0.1 {
            bsum += e;
            bn += 1;
        }
    }
    (sum / nodes.len() as f64, mx, bsum / bn.max(1) as f64)
}

#[test]
fn fit_recovers_a_circle() {
    let s = shape(Kind::Circle, &[0.5]);
    let cloud = cloud_of(&s, 200, 1);
    let (field, rep) = build(&cloud, &ClassicalFit::default(), 4, Eval2dOptions::default());
    assert_eq!(rep.fit.lines, 0);
    assert_eq!(rep.fit.tokens, 16);
    for (i, k) in field.circles().kappa.iter().enumerate() {
        assert!((k + 2.0).abs() < 0.1, "kappa {k} at {i}");
    }
    // contact: g(o) = 0 and the gradient at o is n*
    for i in 0..cloud.len() {
        let c = field.circles();
        let g = ddg_points_as_tori::pat2d::circle_sdf(&c.primitives[i], c.contact[i]);
        assert!(g.abs() < 1e-12, "g(o) = {g}");
        let (gr, _) = ddg_points_as_tori::pat2d::circle_sdf_gradient(&c.primitives[i], c.contact[i]);
        assert!((gr[0] - c.normal[i][0]).abs() < 1e-9 && (gr[1] - c.normal[i][1]).abs() < 1e-9);
    }
}

#[test]
fn concave_curve_gets_a_positive_curvature_and_the_right_sign() {
    // the inside of a ring seen from the hole: normals point toward the center
    let n = 200;
    let pts: Vec<[f64; 2]> = (0..n).map(|i| i as f64 / n as f64 * std::f64::consts::TAU).map(|t| [0.5 * t.cos(), 0.5 * t.sin()]).collect();
    let nrm: Vec<[f64; 2]> = pts.iter().map(|p| [-p[0] / 0.5, -p[1] / 0.5]).collect();
    let (cloud, _) = Cloud2::new(pts, nrm, NormalPolicy::Normalize).unwrap();
    let (field, _) = build(&cloud, &ClassicalFit::default(), 2, Eval2dOptions::default());
    for k in &field.circles().kappa {
        assert!((k - 2.0).abs() < 0.1);
    }
    // inside the hole is "outside the shape": positive; the distance to the circle line
    let v = field.distance([0.0, 0.0]).unwrap();
    assert!((v - 0.5).abs() < 0.02, "{v}");
    let v = field.distance([0.9, 0.0]).unwrap();
    assert!((v + 0.4).abs() < 0.02, "{v}");
}

#[test]
fn straight_line_is_a_half_plane() {
    let n = 50;
    let pts: Vec<[f64; 2]> = (0..n).map(|i| [i as f64 * 0.01, 0.2]).collect();
    let nrm = vec![[0.0, 1.0]; n];
    let (cloud, _) = Cloud2::new(pts, nrm, NormalPolicy::default()).unwrap();
    let (field, rep) = build(&cloud, &ClassicalFit::default(), 2, Eval2dOptions::default());
    assert_eq!(rep.fit.lines, n);
    for x in [[0.1, 0.5], [0.3, -0.4], [0.25, 0.2]] {
        let (v, g) = field.distance_gradient(x).unwrap();
        assert!((v - (x[1] - 0.2)).abs() < 1e-12, "{v}");
        assert!((g[0]).abs() < 1e-9 && (g[1] - 1.0).abs() < 1e-9);
    }
}

#[test]
fn grid_error_on_analytic_shapes() {
    let cases: Vec<(&str, Shape2d)> = vec![
        ("circle", shape(Kind::Circle, &[0.5])),
        ("ellipse", shape(Kind::Ellipse, &[0.6, 0.35])),
        ("box", shape(Kind::Box, &[0.5, 0.3])),
        ("hexagon", shape(Kind::Hexagon, &[0.5])),
        ("heart", shape(Kind::Heart, &[])),
        ("pill", shape(Kind::Pill, &[0.5, 0.2])),
    ];
    for (name, s) in cases {
        let cloud = cloud_of(&s, 600, 3);
        let (field, _) = build(&cloud, &ClassicalFit::default(), 4, Eval2dOptions::default());
        let (mean, max, band) = grid_error(&s, &field, 120);
        let w = s.bbox()[1][0] - s.bbox()[0][0];
        eprintln!("{name}: mean {mean:.5} max {max:.4} band {band:.5} (width {w:.2})");
        assert!(mean < 0.003 * w, "{name}: mean error {mean}");
        assert!(band < 0.001 * w, "{name}: band error {band}");
    }
}

#[test]
fn gradients_match_finite_differences() {
    let s = shape(Kind::Ellipse, &[0.6, 0.35]);
    let cloud = cloud_of(&s, 120, 5);
    let mut cases = vec![
        (Neighborhood::All, LambdaRule::PerQuery { c: 8.0 }, GradientMode::Exact),
        (Neighborhood::All, LambdaRule::Constant { lambda: 12.0 }, GradientMode::Frozen),
        (Neighborhood::All, LambdaRule::Constant { lambda: 12.0 }, GradientMode::Exact),
    ];
    cases.push((Neighborhood::All, LambdaRule::PerQuery { c: 64.0 }, GradientMode::Exact));
    for (nb, lam, gm) in cases {
        let (field, _) = build(&cloud, &ClassicalFit::default(), 2, Eval2dOptions { neighborhood: nb, lambda: lam, gradient: gm, ..Default::default() });
        let mut rng = Rng::new(9);
        for _ in 0..60 {
            let x = [rng.uniform(-1.0, 1.0), rng.uniform(-1.0, 1.0)];
            let (_, g) = field.distance_gradient(x).unwrap();
            let h = 1e-6;
            let fd = [
                (field.distance([x[0] + h, x[1]]).unwrap() - field.distance([x[0] - h, x[1]]).unwrap()) / (2.0 * h),
                (field.distance([x[0], x[1] + h]).unwrap() - field.distance([x[0], x[1] - h]).unwrap()) / (2.0 * h),
            ];
            assert!((g[0] - fd[0]).abs() < 1e-5 * (1.0 + fd[0].abs()) && (g[1] - fd[1]).abs() < 1e-5 * (1.0 + fd[1].abs()), "{lam:?} {gm:?} at {x:?}: {g:?} vs {fd:?}");
        }
    }
}

#[test]
fn default_gradient_is_close_to_finite_differences_with_knn() {
    // with 32-NN the set changes discontinuously; away from those events the frozen gradient is
    // the reference's (approximate) one. Check the median agreement.
    let s = shape(Kind::Hexagon, &[0.5]);
    let cloud = cloud_of(&s, 400, 6);
    let (field, _) = build(&cloud, &ClassicalFit::default(), 2, Eval2dOptions { gradient: GradientMode::Exact, ..Default::default() });
    let mut rng = Rng::new(10);
    let mut errs = Vec::new();
    for _ in 0..200 {
        let x = [rng.uniform(-0.8, 0.8), rng.uniform(-0.8, 0.8)];
        let (_, g) = field.distance_gradient(x).unwrap();
        let h = 1e-7;
        let fd = [
            (field.distance([x[0] + h, x[1]]).unwrap() - field.distance([x[0] - h, x[1]]).unwrap()) / (2.0 * h),
            (field.distance([x[0], x[1] + h]).unwrap() - field.distance([x[0], x[1] - h]).unwrap()) / (2.0 * h),
        ];
        errs.push(((g[0] - fd[0]).powi(2) + (g[1] - fd[1]).powi(2)).sqrt());
    }
    errs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    assert!(errs[errs.len() / 2] < 1e-6, "median {}", errs[errs.len() / 2]);
    assert!(errs[errs.len() * 9 / 10] < 1e-3, "p90 {}", errs[errs.len() * 9 / 10]);
}

#[test]
fn bitwise_independent_of_the_thread_count() {
    let s = shape(Kind::Heart, &[]);
    let cloud = cloud_of(&s, 500, 7);
    let mut reference: Option<(Vec<u64>, Vec<u64>)> = None;
    for threads in [1usize, 3, 8, 32] {
        let (field, _) = build(&cloud, &ClassicalFit::default(), threads, Eval2dOptions::default());
        let g = evaluate_grid(&field, padded(s.bbox(), 0.2), 90, 90).unwrap();
        let mut grad = vec![[0.0; 2]; g.values.len()];
        field.evaluate(&g.nodes(), None, Some(&mut grad)).unwrap();
        let bits = (
            g.values.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
            field.circles().kappa.iter().chain(grad.iter().flatten()).map(|v| v.to_bits()).collect::<Vec<_>>(),
        );
        match &reference {
            None => reference = Some(bits),
            Some(r) => assert!(*r == bits, "threads = {threads} differs"),
        }
    }
}

#[test]
fn marching_squares_traces_the_circle() {
    let s = shape(Kind::Circle, &[0.5]);
    let cloud = cloud_of(&s, 300, 8);
    let (field, _) = build(&cloud, &ClassicalFit::default(), 4, Eval2dOptions::default());
    let g = evaluate_grid(&field, [[-1.0, -1.0], [1.0, 1.0]], 129, 129).unwrap();
    let segs = marching_squares(&g, 0.0);
    let polys = chain_segments(&segs);
    assert_eq!(polys.len(), 1);
    assert!(polys[0].closed);
    let mut area = 0.0;
    let p = &polys[0].points;
    for i in 0..p.len() {
        let (a, b) = (p[i], p[(i + 1) % p.len()]);
        area += 0.5 * (a[0] * b[1] - b[0] * a[1]);
        assert!((a[0].hypot(a[1]) - 0.5).abs() < 0.01);
    }
    // inside on the left: counter-clockwise, positive area
    assert!((area - std::f64::consts::PI * 0.25).abs() < 0.01, "area {area}");
    // an offset level is a bigger circle
    let polys = chain_segments(&marching_squares(&g, 0.2));
    assert_eq!(polys.len(), 1);
    assert!((polys[0].points[0][0].hypot(polys[0].points[0][1]) - 0.7).abs() < 0.02);
}

#[test]
fn marching_squares_saddle_and_border() {
    // f = x*y on [-1,1]^2 at level 0: two lines crossing; with a grid that puts 4 crossings in one cell
    let g = Grid2d { origin: [-1.0, -1.0], spacing: [1.0, 1.0], nx: 3, ny: 3, values: vec![1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0] };
    let segs = marching_squares(&g, 0.0);
    assert_eq!(segs.len(), 8); // two per cell, four cells
    let polys = chain_segments(&segs);
    assert!(!polys.is_empty());
    let empty = Grid2d { origin: [0.0; 2], spacing: [1.0; 2], nx: 2, ny: 2, values: vec![1.0; 4] };
    assert!(marching_squares(&empty, 0.0).is_empty());
    let nan = Grid2d { origin: [0.0; 2], spacing: [1.0; 2], nx: 2, ny: 2, values: vec![f64::NAN, 1.0, -1.0, 1.0] };
    assert!(marching_squares(&nan, 0.0).is_empty());
}

/// The seam for WP9b: any predictor of `(a0, a1, a2)` works. This one knows the circle exactly.
struct ExactCirclePredictor {
    radius: f64,
}

impl QuadraticPredictor for ExactCirclePredictor {
    fn neighborhood_size(&self) -> usize {
        4
    }
    fn predict(&self, cloud: &Cloud2, _n: &ddg_points_as_tori::Neighborhoods, _t: Option<usize>) -> Result<Predicted2, PatError> {
        // a circle seen from the outside with outward normals: Q(s) = -r + sqrt(r^2 - s^2) ≈ -s^2/(2r), a2 = -1/(2r)
        Ok(Predicted2 { coefficients: vec![Quad2 { a0: 0.0, a1: 0.0, a2: -0.5 / self.radius }; cloud.len()], sigma: vec![f64::NAN; cloud.len()], fallbacks: 0 })
    }
}

#[test]
fn a_custom_predictor_plugs_in() {
    let s = shape(Kind::Circle, &[0.5]);
    let cloud = cloud_of(&s, 100, 2);
    let (field, rep) = Pat2d::build(&cloud, &ExactCirclePredictor { radius: 0.5 }, &BuildOptions2::default().with_threads(Some(2))).unwrap();
    assert_eq!(rep.fit.tokens, 4);
    for x in [[0.0f64, 0.0], [0.9, 0.1], [0.3, -0.2]] {
        let want = (x[0] * x[0] + x[1] * x[1]).sqrt() - 0.5;
        assert!((field.distance(x).unwrap() - want).abs() < 1e-12);
    }
    // a predictor that answers with the wrong number of rows is caught
    struct Short;
    impl QuadraticPredictor for Short {
        fn neighborhood_size(&self) -> usize {
            4
        }
        fn predict(&self, _c: &Cloud2, _n: &ddg_points_as_tori::Neighborhoods, _t: Option<usize>) -> Result<Predicted2, PatError> {
            Ok(Predicted2 { coefficients: vec![], sigma: vec![], fallbacks: 0 })
        }
    }
    assert!(matches!(Pat2d::build(&cloud, &Short, &BuildOptions2::default()), Err(PatError::SizeMismatch { .. })));
}

#[test]
fn bad_input_is_an_error() {
    assert!(matches!(Cloud2::new(vec![[0.0; 2]], vec![[1.0, 0.0]], NormalPolicy::default()), Err(PatError::TooFewPoints { .. })));
    assert!(matches!(Cloud2::new(vec![[0.0; 2], [1.0, 0.0]], vec![[1.0, 0.0]], NormalPolicy::default()), Err(PatError::LengthMismatch { .. })));
    assert!(matches!(Cloud2::new(vec![[0.0; 2], [f64::NAN, 0.0]], vec![[1.0, 0.0]; 2], NormalPolicy::default()), Err(PatError::NonFinite { index: 1, .. })));
    assert!(matches!(Cloud2::new(vec![[0.0; 2], [1.0, 0.0]], vec![[0.0, 0.0], [1.0, 0.0]], NormalPolicy::default()), Err(PatError::ZeroNormal { index: 0 })));
    assert!(matches!(Cloud2::new(vec![[0.0; 2], [1.0, 0.0]], vec![[2.0, 0.0], [1.0, 0.0]], NormalPolicy::default()), Err(PatError::NonUnitNormal { index: 0, .. })));
    let (cloud, _) = Cloud2::new(vec![[0.0; 2], [1.0, 0.0], [0.5, 0.5]], vec![[0.0, 1.0]; 3], NormalPolicy::default()).unwrap();
    assert!(Pat2d::build(&cloud, &ClassicalFit { k: 2, bandwidth: 1.0 }, &BuildOptions2::default()).is_err());
    let (field, _) = Pat2d::build(&cloud, &ClassicalFit::default(), &BuildOptions2::default().with_threads(Some(1))).unwrap();
    assert!(field.distance([f64::NAN, 0.0]).is_err());
    // a query exactly at every neighbor: all points equal
    let (dup, _) = Cloud2::new(vec![[0.0; 2]; 3], vec![[0.0, 1.0]; 3], NormalPolicy::default()).unwrap();
    let (f2, _) = Pat2d::build(&dup, &ClassicalFit::default(), &BuildOptions2::default().with_threads(Some(1))).unwrap();
    assert!(matches!(f2.distance([0.0, 0.0]), Err(PatError::DegenerateQuery { .. })));
    assert!(evaluate_grid(&field, [[0.0, 0.0], [0.0, 1.0]], 4, 4).is_err());
    assert!(Raster::new(0, 3, [[0.0; 2], [1.0; 2]], [0; 3]).is_err());
}

#[test]
fn png_and_svg_are_written() {
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../runs_tmp/pat-wp9a/test-out");
    let mut r = Raster::new(40, 30, [[0.0, 0.0], [4.0, 3.0]], [255, 255, 255]).unwrap();
    r.dot([2.0, 1.5], 3.0, [200, 0, 0]);
    r.line([0.0, 0.0], [4.0, 3.0], 2.0, [0, 0, 200]);
    r.circle_arc([2.0, 1.5], 1.0, [1.0, 0.0], 3.0, 1.5, [0, 120, 0]);
    let p = dir.join("t.png");
    write_png(&p, &r).unwrap();
    let b = std::fs::read(&p).unwrap();
    assert_eq!(&b[..8], &[0x89, b'P', b'N', b'G', 0x0d, 0x0a, 0x1a, 0x0a]);
    assert!(b.len() > 3 * 40 * 30);
    let poly = Polyline { points: vec![[0.0, 0.0], [1.0, 1.0], [2.0, 0.0]], closed: true };
    write_svg_contours(&dir.join("t.svg"), [[0.0, 0.0], [4.0, 3.0]], 200, &[poly], &[[1.0, 1.0]]).unwrap();
    assert!(std::fs::read_to_string(dir.join("t.svg")).unwrap().contains("<polygon"));
}

/// Measurements for NOTES-WP9a: error against k, bandwidth and sample count.
/// `cargo test --release --test pat2d -- --ignored --nocapture measure`
#[test]
#[ignore]
fn measure() {
    let shapes: Vec<(&str, Shape2d)> = vec![
        ("circle", shape(Kind::Circle, &[0.5])),
        ("ellipse", shape(Kind::Ellipse, &[0.6, 0.35])),
        ("box", shape(Kind::Box, &[0.5, 0.3])),
        ("hexagon", shape(Kind::Hexagon, &[0.5])),
        ("heart", shape(Kind::Heart, &[])),
        ("pill", shape(Kind::Pill, &[0.5, 0.2])),
        ("moon", shape(Kind::Moon, &[0.5, 0.6, 0.35])),
        ("pie", shape(Kind::Pie, &[0.6, 1.0])),
    ];
    println!("shape\tn\tk\tbw\tmean\tmax\tband");
    for (name, s) in &shapes {
        for n in [200usize, 800] {
            for (k, bw) in [(8usize, 0.6), (16, 0.6), (16, 1e9), (32, 0.6), (32, 0.4), (64, 0.5)] {
                let cloud = cloud_of(s, n, 3);
                let bwv = if bw > 1e8 { f64::INFINITY } else { bw };
                let (field, _) = build(&cloud, &ClassicalFit { k, bandwidth: bwv }, 16, Eval2dOptions::default());
                let (mean, max, band) = grid_error(s, &field, 160);
                println!("{name}\t{n}\t{k}\t{bw}\t{mean:.5}\t{max:.4}\t{band:.5}");
            }
        }
    }
}
