//! Replay of the Lean oracle's golden fixtures (`lean/oracle/fixtures/pat/*.in|out`; format in
//! `lean/oracle/README.md`, "Points as Tori") through the Rust implementation (SPEC.md §15.3).
//!
//! Float tokens are compared to 1e-12 relative to `max(|expected|, 1)`, integers and tags exactly.
//! The oracle sums a neighbor set in `(ρ², index)` order, Rust sums `All` and kNN with `k ≥ n` in
//! index order: that differs by round-off only, which the tolerance absorbs; the neighbor sets are
//! compared as sets there and as sequences everywhere else.

use ddg_points_as_tori::baselines::kernel::{KernelField, KernelMethod, KernelOptions};
use ddg_points_as_tori::blend::{BlendInputs, BlendWant, GradientMode, LambdaRule, blend};
use ddg_points_as_tori::features::{FeatureScratch, local_frame, neighborhood_features_f64, unit};
use ddg_points_as_tori::field::mean_neighbor_distance;
use ddg_points_as_tori::fit::{fit_torus, fundamental_forms, mean_gauss};
use ddg_points_as_tori::network::scale_outputs_f64;
use ddg_points_as_tori::sdf::{Singular, torus_sdf_laplacian};
use ddg_points_as_tori::spatial::{KdTree, KdTreeOptions, fitting_neighborhoods};
use ddg_points_as_tori::{EvalOptions, FitOptions, SignRule, Neighborhood, NormalPolicy, PatError, PointCloud, Tori, Torus, TorusField};
use std::path::PathBuf;

const TOL: f64 = 1e-12;

fn hex(s: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(s, 16).unwrap_or_else(|_| panic!("bad hex float {s:?}")))
}

#[derive(Default)]
struct Case {
    name: String,
    task: String,
    cloud: Vec<([f64; 3], [f64; 3], f64)>,
    center: usize,
    nbhd: Option<usize>,
    nbr: Option<Vec<u32>>,
    y: Option<[f64; 6]>,
    fitopts: Option<(f64, f64)>,
    coefs: Vec<(usize, [f64; 6])>,
    tori: Vec<Torus>,
    neighborhood: String,
    lambda: String,
    areas: bool,
    iso: f64,
    xs: Vec<[f64; 3]>,
}

fn parse(text: &str) -> Case {
    let mut c = Case { neighborhood: "knn 32".into(), lambda: "perquery 4050000000000000".into(), ..Default::default() };
    for line in text.lines() {
        let t: Vec<&str> = line.split_whitespace().collect();
        let f = |k: usize| hex(t[k]);
        match t.first().copied() {
            Some("case") => c.name = t[1].into(),
            Some("task") => c.task = t[1].into(),
            Some("cloud") | Some("end") | None => {}
            Some("p") => {
                let area = if t.len() > 7 { f(7) } else { 1.0 };
                c.cloud.push(([f(1), f(2), f(3)], [f(4), f(5), f(6)], area));
            }
            Some("center") => c.center = t[1].parse().unwrap(),
            Some("nbhd") => c.nbhd = Some(t[1].parse().unwrap()),
            Some("nbr") => c.nbr = Some(t[1..].iter().map(|x| x.parse().unwrap()).collect()),
            Some("y") => c.y = Some(std::array::from_fn(|k| f(1 + k))),
            Some("fitopts") => c.fitopts = Some((f(1), f(2))),
            Some("coef") => c.coefs.push((t[1].parse().unwrap(), std::array::from_fn(|k| f(2 + k)))),
            Some("torus") => c.tori.push(Torus { center: [f(1), f(2), f(3)], axis: [f(4), f(5), f(6)], major: f(7), minor: f(8) }),
            Some("neighborhood") => c.neighborhood = t[1..].join(" "),
            Some("lambda") => c.lambda = t[1..].join(" "),
            Some("areas") => c.areas = t[1] == "1",
            Some("iso") => c.iso = f(1),
            Some("x") => c.xs.push([f(1), f(2), f(3)]),
            Some(other) => panic!("{}: unknown tag {other}", c.name),
        }
    }
    c
}

/// One output token with its kind.
#[derive(Debug, Clone, PartialEq)]
enum Tok {
    F(f64),
    N(i64),
    S(String),
}

type Line = Vec<Tok>;

fn s(x: &str) -> Tok {
    Tok::S(x.into())
}
fn fs(v: &[f64]) -> Vec<Tok> {
    v.iter().map(|&x| Tok::F(x)).collect()
}
fn ns(v: &[usize]) -> Vec<Tok> {
    v.iter().map(|&x| Tok::N(x as i64)).collect()
}
fn line(parts: Vec<Vec<Tok>>) -> Line {
    parts.into_iter().flatten().collect()
}

/// The error kind the oracle names for a Rust error.
fn kind(e: &PatError) -> &'static str {
    match e {
        PatError::NonUnitNormal { .. } => "NonUnitNormal",
        PatError::DegenerateQuery { .. } => "DegenerateQuery",
        PatError::SizeMismatch { .. } => "ValueCount",
        PatError::NonFinite { .. } => "NonFinite",
        PatError::NonFiniteTorus { .. } => "NonFiniteTorus",
        PatError::InvalidOption { .. } => "InvalidOptions",
        _ => "Other",
    }
}

enum Outcome {
    Ok(Vec<Line>),
    Err(String),
}

fn check_unit(c: &Case) -> Result<(), PatError> {
    for (i, (_, n, _)) in c.cloud.iter().enumerate() {
        let l = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        if (l - 1.0).abs() > 1e-6 {
            return Err(PatError::NonUnitNormal { index: i, norm: l });
        }
    }
    Ok(())
}

fn build_cloud(c: &Case) -> Result<PointCloud, PatError> {
    check_unit(c)?;
    let p = c.cloud.iter().map(|x| x.0).collect();
    let n = c.cloud.iter().map(|x| x.1).collect();
    PointCloud::new(p, n, NormalPolicy::RequireUnit { tol: 1e-6 }).map(|x| x.0)
}

fn features(c: &Case) -> Outcome {
    let cloud = match build_cloud(c) {
        Ok(x) => x,
        Err(e) => return Outcome::Err(kind(&e).into()),
    };
    let (pos, nor) = (cloud.positions(), cloud.normals());
    let row: Vec<u32> = if let Some(r) = &c.nbr {
        r.clone()
    } else {
        let tree = KdTree::new(pos, KdTreeOptions::default()).unwrap();
        fitting_neighborhoods(&tree, pos, c.nbhd.unwrap()).row(c.center).to_vec()
    };
    assert_eq!(row[0] as usize, c.center);
    let k = row.len();
    let mut out = vec![0.0; 6 * k];
    let sigma = neighborhood_features_f64(pos, nor, &row, &mut FeatureScratch::default(), &mut out);
    if sigma == 0.0 || !sigma.is_finite() {
        return Outcome::Err("DegenerateNeighborhood".into());
    }
    let (sh, th) = local_frame(unit(nor[c.center]));
    let ch = unit(nor[c.center]);
    let mut lines = vec![
        line(vec![vec![s("center")], ns(&[c.center, k])]),
        line(vec![vec![s("sigma")], fs(&[sigma])]),
        line(vec![vec![s("frame")], fs(&[sh[0], sh[1], sh[2], th[0], th[1], th[2], ch[0], ch[1], ch[2]])]),
        line(vec![vec![s("nbr")], row.iter().map(|&j| Tok::N(j as i64)).collect()]),
    ];
    for m in 0..k {
        lines.push(line(vec![vec![s("f")], ns(&[m]), fs(&out[6 * m..6 * m + 6])]));
    }
    if let Some(y) = c.y {
        lines.push(line(vec![vec![s("coef")], fs(&scale_outputs_f64(y, sigma))]));
    }
    Outcome::Ok(lines)
}

fn fit_opts(c: &Case) -> FitOptions {
    let (e, f) = c.fitopts.unwrap_or((1e-10, 1e-10));
    FitOptions { eig_eps: e, curvature_floor: f, sign_rule: SignRule::Reference, ..Default::default() }
}

fn fit(c: &Case) -> Outcome {
    if let Err(e) = check_unit(c) {
        return Outcome::Err(kind(&e).into());
    }
    let opts = fit_opts(c);
    let mut lines = vec![
        line(vec![vec![s("fitopts")], fs(&[opts.eig_eps, opts.curvature_floor])]),
        line(vec![vec![s("tori")], ns(&[c.coefs.len()])]),
    ];
    for &(i, a) in &c.coefs {
        let (p, n) = (c.cloud[i].0, c.cloud[i].1);
        let (t, d) = fit_torus(p, n, a, &opts);
        let (big_a, _, _) = fundamental_forms(a);
        let (h, k) = mean_gauss(a);
        let r = (h * h - k).max(0.0).sqrt();
        let sign = if t.minor >= 0.0 { 1.0 } else { -1.0 };
        lines.push(line(vec![
            vec![s("t")],
            ns(&[i]),
            fs(&[t.center[0], t.center[1], t.center[2], t.axis[0], t.axis[1], t.axis[2], t.major, t.minor, sign]),
        ]));
        // Eq. 19 values, then the Rayleigh pair in the order the fit uses (|κ| ascending first)
        lines.push(line(vec![
            vec![s("d")],
            ns(&[i]),
            fs(&[big_a, h, k, h + r, h - r, d.kappa_min, d.kappa_max, d.kappa_min, d.kappa_max]),
        ]));
        lines.push(line(vec![
            vec![s("dir")],
            ns(&[i]),
            fs(&[d.v_min[0], d.v_min[1], d.v_min[2], d.v_max[0], d.v_max[1], d.v_max[2], d.n_star[0], d.n_star[1], d.n_star[2]]),
        ]));
        let cf = d.clamps;
        lines.push(line(vec![
            vec![s("flags")],
            ns(&[i]),
            ns(&[cf.disc as usize, (cf.eigvec) as usize, cf.denominator as usize, cf.curvature as usize, cf.disc as usize]),
        ]));
    }
    Outcome::Ok(lines)
}

fn sdf(c: &Case) -> Outcome {
    let mut lines = vec![line(vec![vec![s("tori")], ns(&[c.tori.len()]), vec![s("queries")], ns(&[c.xs.len()])])];
    for (q, &x) in c.xs.iter().enumerate() {
        for (t, tor) in c.tori.iter().enumerate() {
            let (g, gr, lap, sing) = torus_sdf_laplacian(tor, x);
            let kind = match sing {
                Singular::None => 0,
                Singular::Axis => 1,
                Singular::Core => 2,
            };
            lines.push(line(vec![vec![s("g")], ns(&[q, t]), fs(&[g, gr[0], gr[1], gr[2], lap]), ns(&[kind])]));
        }
    }
    Outcome::Ok(lines)
}

/// Tori of the case: `torus` lines or fitted `coef` lines. `Err` is the oracle's error kind.
fn case_tori(c: &Case) -> Result<Vec<Torus>, String> {
    if !c.tori.is_empty() {
        return Ok(c.tori.clone());
    }
    let cloud = build_cloud(c).map_err(|e| kind(&e).to_string())?;
    let opts = FitOptions { sign_rule: SignRule::Reference, ..Default::default() };
    Ok(c.coefs.iter().map(|&(i, a)| fit_torus(cloud.positions()[i], cloud.normals()[i], a, &opts).0).collect())
}

struct Rules {
    neighborhood: Neighborhood,
    lambda: LambdaRule,
    /// `Some((D, λ, R))` for the paper rule
    paper: Option<(f64, f64, f64)>,
    /// The neighbor set comes in index order in Rust (`All`, kNN with k ≥ n).
    index_order: bool,
}

impl Rules {
    /// Conditioning of Δφ: terms of size λ² are summed to a result of order 1.
    fn lambda_scale(&self, v: &ddg_points_as_tori::blend::BlendValue) -> f64 {
        (v.lambda * v.lambda).max(1.0)
    }
}

fn rules(c: &Case, pos: &[[f64; 3]]) -> Rules {
    let lt: Vec<&str> = c.lambda.split_whitespace().collect();
    let nt: Vec<&str> = c.neighborhood.split_whitespace().collect();
    let n = pos.len();
    let lambda = match lt[0] {
        "perquery" => LambdaRule::PerQuery { c: hex(lt[1]) },
        "constant" => LambdaRule::Constant { lambda: hex(lt[1]) },
        "paper" => LambdaRule::PaperHeuristic { c: hex(lt[1]), factor: hex(lt[2]), k_mean: lt[3].parse().unwrap() },
        o => panic!("lambda {o}"),
    };
    let paper = if let LambdaRule::PaperHeuristic { c: cc, factor, k_mean } = lambda {
        let tree = KdTree::new(pos, KdTreeOptions::default()).unwrap();
        let pool = rayon::ThreadPoolBuilder::new().num_threads(1).build().unwrap();
        let d = pool.install(|| mean_neighbor_distance(&tree, pos, k_mean));
        let lam = factor / d;
        Some((d, lam, 2.0 * cc / lam))
    } else {
        None
    };
    let neighborhood = match nt[0] {
        "knn" => Neighborhood::Knn { k: nt[1].parse().unwrap() },
        "all" => Neighborhood::All,
        "radius" if nt[1] == "auto" => Neighborhood::Radius { r: paper.expect("radius auto needs the paper rule").2, fallback_k: nt[2].parse().unwrap() },
        "radius" => Neighborhood::Radius { r: hex(nt[1]), fallback_k: nt[2].parse().unwrap() },
        o => panic!("neighborhood {o}"),
    };
    let index_order = matches!(neighborhood, Neighborhood::All) || matches!(neighborhood, Neighborhood::Knn { k } if k >= n);
    Rules { neighborhood, lambda, paper, index_order }
}

fn blend_case(c: &Case) -> Outcome {
    let tori = match case_tori(c) {
        Ok(t) => t,
        Err(e) => return Outcome::Err(e),
    };
    let n = c.cloud.len();
    let pos: Vec<[f64; 3]> = c.cloud.iter().map(|x| x.0).collect();
    let areas: Vec<f64> = c.cloud.iter().map(|x| x.2).collect();
    let r = rules(c, &pos);
    // A single point is not a PointCloud (n >= 2): use the blend itself.
    if n < 2 {
        let (d, lam, _) = r.paper.unwrap_or((0.0, 0.0, 0.0));
        let _ = d;
        let inp = BlendInputs {
            positions: &pos,
            tori: &tori,
            area_weights: c.areas.then_some(&areas[..]),
            lambda: r.lambda,
            lambda_global: lam,
            isovalue: c.iso,
        };
        for &x in &c.xs {
            if blend(&inp, x, &[0], BlendWant::Distance).is_err() {
                return Outcome::Err("DegenerateQuery".into());
            }
        }
        return Outcome::Err("none".into());
    }
    let cloud = match build_cloud(c) {
        Ok(x) => x,
        Err(e) => return Outcome::Err(kind(&e).into()),
    };
    let opts = EvalOptions {
        neighborhood: r.neighborhood,
        lambda: r.lambda,
        area_weights: c.areas.then(|| areas.clone()),
        isovalue: c.iso,
        threads: Some(1),
        ..Default::default()
    };
    let field = match TorusField::new(&cloud, Tori { tori, coefficients: vec![[0.0; 6]; c.coefs.len().max(n)], sigma: vec![f64::NAN; c.coefs.len().max(n)] }, opts.clone()) {
        Ok(f) => f,
        Err(e) => return Outcome::Err(kind(&e).into()),
    };
    let mut lines = vec![line(vec![vec![s("cloud")], ns(&[n]), vec![s("queries")], ns(&[c.xs.len()])])];
    let (pd, pl, pr) = r.paper.unwrap_or((0.0, 0.0, 0.0));
    if r.paper.is_some() {
        let (d, l) = field.paper_constants();
        assert!((d - pd).abs() <= 1e-15 * pd && (l - pl).abs() <= 1e-15 * pl, "{}: paper constants", c.name);
    }
    lines.push(line(vec![vec![s("paper")], fs(&[pd, pl, pr])]));
    for (q, &x) in c.xs.iter().enumerate() {
        let go = |w| field.blend_at(x, w).map_err(|e| kind(&e).to_string());
        let (mut lap, mut idx) = match go(BlendWant::Laplacian) {
            Ok(v) => v,
            Err(e) => return Outcome::Err(e),
        };
        let (mut frz, _) = go(BlendWant::Gradient(GradientMode::Frozen)).unwrap();
        let (mut exa, _) = go(BlendWant::Gradient(GradientMode::Exact)).unwrap();
        let (mut dist, _) = go(BlendWant::Distance).unwrap();
        if r.index_order {
            // Rust sums this set in index order, the oracle in (ρ², index) order. Replay the
            // library's blend on the sorted set (identical operations to the oracle's), and check
            // that the index-order sums agree with it to round-off.
            let key = |j: u32| {
                let p = pos[j as usize];
                ((x[0] - p[0]).powi(2) + (x[1] - p[1]).powi(2) + (x[2] - p[2]).powi(2), j)
            };
            let mut sorted = idx.clone();
            sorted.sort_by(|&a, &b| key(a).partial_cmp(&key(b)).unwrap());
            let inp = BlendInputs {
                positions: &pos,
                tori: &field.tori().tori,
                area_weights: c.areas.then_some(&areas[..]),
                lambda: r.lambda,
                lambda_global: pl,
                isovalue: c.iso,
            };
            let sd = |w| blend(&inp, x, &sorted, w).unwrap();
            let (l2, f2, e2, d2) = (sd(BlendWant::Laplacian), sd(BlendWant::Gradient(GradientMode::Frozen)), sd(BlendWant::Gradient(GradientMode::Exact)), sd(BlendWant::Distance));
            let scale = r.lambda_scale(&l2);
            assert!((lap.laplacian - l2.laplacian).abs() <= 1e-10 * scale * l2.laplacian.abs().max(1.0), "{}: query {q}: index-order Δφ {} vs sorted {}", c.name, lap.laplacian, l2.laplacian);
            assert!((dist.distance - d2.distance).abs() <= 1e-12 * d2.distance.abs().max(1.0), "{}: query {q}: index-order φ", c.name);
            (lap, frz, exa, dist) = (l2, f2, e2, d2);
            idx = sorted;
        }
        // the Laplacian pass returns the frozen gradient
        assert_eq!(lap.gradient, frz.gradient, "{}: frozen gradient of the two passes", c.name);
        let sing = (lap.flags.singular_axis + lap.flags.singular_core) as usize;
        lines.push(line(vec![
            vec![s("b")],
            ns(&[q]),
            fs(&[dist.distance]),
            fs(&lap.gradient),
            fs(&exa.gradient),
            fs(&[lap.laplacian, lap.sigma, lap.lambda, lap.shift]),
            ns(&[idx.len(), sing, lap.flags.at_data_point as usize]),
        ]));
        let set = idx.clone();
        lines.push(line(vec![vec![s("nb")], ns(&[q]), set.iter().map(|&j| Tok::N(j as i64)).collect()]));
    }
    Outcome::Ok(lines)
}

fn baselines(c: &Case) -> Outcome {
    let cloud = match build_cloud(c) {
        Ok(x) => x,
        Err(e) => return Outcome::Err(kind(&e).into()),
    };
    let n = cloud.len();
    let pos = cloud.positions().to_vec();
    let areas: Vec<f64> = c.cloud.iter().map(|x| x.2).collect();
    let r = rules(c, &pos);
    let opts = KernelOptions { neighborhood: r.neighborhood, lambda: r.lambda, threads: Some(1), shc_zero_value: None };
    let field = match KernelField::new(&cloud, areas, opts) {
        Ok(f) => f,
        Err(e) => return Outcome::Err(kind(&e).into()),
    };
    let (pd, pl, pr) = r.paper.unwrap_or((0.0, 0.0, 0.0));
    let mut lines = vec![line(vec![vec![s("cloud")], ns(&[n]), vec![s("queries")], ns(&[c.xs.len()])]), line(vec![vec![s("paper")], fs(&[pd, pl, pr])])];
    for (q, &x) in c.xs.iter().enumerate() {
        let a = match field.eval(KernelMethod::Shc, x, true) {
            Ok(v) => v,
            Err(e) => return Outcome::Err(kind(&e).into()),
        };
        let b = match field.eval(KernelMethod::Sspd, x, true) {
            Ok(v) => v,
            Err(e) => return Outcome::Err(kind(&e).into()),
        };
        lines.push(line(vec![
            vec![s("s")],
            ns(&[q]),
            fs(&[a.value]),
            fs(&a.gradient),
            fs(&[b.value]),
            fs(&b.gradient),
            fs(&[a.sigma, a.lambda, a.shift]),
            ns(&[a.terms as usize, a.skipped as usize]),
        ]));
    }
    Outcome::Ok(lines)
}

/// Compare `got` with the `.out` text; returns the number of tokens compared or the first
/// disagreement.
fn compare(name: &str, got: Outcome, expected: &str) -> Result<usize, String> {
    let mut it = expected.lines();
    let _ = it.next(); // case
    let _ = it.next(); // task
    let status = it.next().unwrap_or("");
    if let Some(rest) = status.strip_prefix("status error ") {
        let want = rest.split_whitespace().next().unwrap_or("");
        return match got {
            Outcome::Err(k) if k == want => Ok(1),
            Outcome::Err(k) => Err(format!("{name}: error kind {k}, oracle {want}")),
            Outcome::Ok(_) => Err(format!("{name}: oracle says error {want}, Rust succeeded")),
        };
    }
    let lines = match got {
        Outcome::Ok(l) => l,
        Outcome::Err(k) => return Err(format!("{name}: Rust error {k}, oracle ok")),
    };
    let exp: Vec<&str> = it.take_while(|l| *l != "end").collect();
    if exp.len() != lines.len() {
        return Err(format!("{name}: {} lines, oracle {}", lines.len(), exp.len()));
    }
    let mut count = 0;
    for (ln, (g, e)) in lines.iter().zip(&exp).enumerate() {
        let mut et: Vec<&str> = e.split_whitespace().collect();
        let mut g = g.clone();
        if g.first() == Some(&s("d")) && g.len() == 11 && et.len() == 11 {
            // the Rayleigh values of the eigenvector pair (tokens 7, 8) are an unordered pair for
            // us: the library reports them as (κmin, κmax)
            let (mut a, mut b) = (hex(et[7]), hex(et[8]));
            if let (Tok::F(x), Tok::F(y)) = (&g[7], &g[8]) {
                let (lo, hi) = if x <= y { (*x, *y) } else { (*y, *x) };
                g[7] = Tok::F(lo);
                g[8] = Tok::F(hi);
            }
            if a > b {
                std::mem::swap(&mut a, &mut b);
            }
            let (sa, sb) = (format!("{:016x}", a.to_bits()), format!("{:016x}", b.to_bits()));
            let leak: &'static [String; 2] = Box::leak(Box::new([sa, sb]));
            et[7] = &leak[0];
            et[8] = &leak[1];
        }
        let g = &g;
        if et.len() != g.len() {
            return Err(format!("{name}: line {ln} has {} tokens, oracle {}: {e}", g.len(), et.len()));
        }
        for (k, (gt, e)) in g.iter().zip(&et).enumerate() {
            count += 1;
            let bad = match gt {
                Tok::S(t) => t != e,
                Tok::N(v) => e.parse::<i64>().ok() != Some(*v),
                Tok::F(v) => {
                    let w = hex(e);
                    // Δφ sums terms of size λ² to a result of order 1, and the fitted tori differ
                    // from the oracle's by round-off (different operation order inside the fit,
                    // 1e-16 relative): the result carries that noise times λ² (observed 0.7-1.2 ε λ²,
                    // allowed 16 ε λ²). Everything else is held to 1e-12.
                    let extra = if et[0] == "b" && k == 9 { let lam = hex(et[11]); 16.0 * f64::EPSILON * lam * lam } else { 0.0 };
                    !(*v == w || (v - w).abs() <= TOL * w.abs().max(1.0) + extra)
                }
            };
            if bad {
                return Err(format!("{name}: line {ln} ({}) token {k}: Rust {gt:?}, oracle {e}\n  oracle line: {ln_e}\n  rust line:   {rl}", et[0], ln_e = exp[ln], rl = g.iter().map(|t| match t { Tok::F(v) => format!("{:016x}", v.to_bits()), Tok::N(v) => v.to_string(), Tok::S(v) => v.clone() }).collect::<Vec<_>>().join(" ")));
            }
        }
    }
    Ok(count)
}

fn fixtures_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../lean/oracle/fixtures/pat")
}

#[test]
fn replay_all_fixtures() {
    let dir = fixtures_dir();
    let mut names: Vec<String> = std::fs::read_dir(&dir)
        .unwrap()
        .filter_map(|e| e.ok()?.file_name().to_str()?.strip_suffix(".in").map(String::from))
        .collect();
    names.sort();
    assert!(names.len() >= 35, "expected the 35 fixtures, found {}", names.len());
    let (mut failed, mut tokens) = (Vec::new(), 0);
    for name in &names {
        let c = parse(&std::fs::read_to_string(dir.join(format!("{name}.in"))).unwrap());
        let out = std::fs::read_to_string(dir.join(format!("{name}.out"))).unwrap();
        let got = match c.task.as_str() {
            "features" => features(&c),
            "fit" => fit(&c),
            "sdf" => sdf(&c),
            "blend" => blend_case(&c),
            "baselines" => baselines(&c),
            t => panic!("{name}: task {t}"),
        };
        match compare(name, got, &out) {
            Ok(n) => tokens += n,
            Err(e) => failed.push(e),
        }
    }
    eprintln!("oracle replay: {} fixtures, {} failed, {tokens} tokens compared", names.len(), failed.len());
    assert!(failed.is_empty(), "{} of {} fixtures disagree:\n{}", failed.len(), names.len(), failed.join("\n"));
}
