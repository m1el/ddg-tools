//! Study of estimators of the normal derivative `∂u/∂n` on the cut (SPEC §3.5, §3.7) where the
//! exact values are known, on cuts chosen so that crossings come arbitrarily close to mesh
//! vertices.
//!
//! * Unit sphere cut along circles (SPEC §6 item 2): `∂u/∂n = tan(ψ₀/2)/R` on the side of the
//!   cap of polar angle `ψ₀`, `cot(ψ₀/2)/R` on the other side, `R = 1/√(4π)`. Constant along the
//!   cut, so this measures robustness but cannot see the bias of a smoother.
//! * Flat annulus `0.4 ≤ r ≤ 1` cut along the circle `r = ρ`, with the source
//!   `K = K₀ + κ cos(mθ)` as a curvature override: both regions are annuli and
//!   `u = K₀ r²/4 + A ln r + B + (κ r²/(4 − m²) + C r^m + D r^−m) cos(mθ)` in each. Here `∂u/∂n`
//!   varies along the cut with `m` periods.
//!
//! In both, `s = (∂u₀/∂n)² − (∂u₁/∂n)²` (Dirichlet energy). Errors are absolute errors divided
//! by the largest exact value on the cut, in percent.
//!
//! ```text
//! cargo run --release -p ddg-variational-cuts --example estimator_study -- 3 6
//! ```
//!
//! The estimators (A nodal, with sliver smoothing or interpolation over slivers; B segment
//! gradient; C clusters; H smoothed, `κ = c² h^(1+γ)/|seg|^γ` with `γ = 1`, H2 with `γ = 2`) are
//! computed here from the public fields `flux`, `dual_length`, `log_scale`, `triangles`, `cut`
//! of an `Evaluation`. The last row is what the library returns with its default options; it
//! must agree with `H  c=0.50`.

use std::f64::consts::PI;

use ddg_core::generators::{SplitMix64, annulus, icosphere};
use ddg_core::sparse::{Method, Parallelism, Symbolic};
use ddg_core::{Geometry, SparseMatrix, TriMesh};
use ddg_variational_cuts::{
    CurvatureSource, Cut, CutProblem, EdgeKind, Energy, EvalOptions, Evaluation, Evaluator,
    smooth_slivers,
};

struct Sphere {
    mesh: TriMesh,
    geometry: Geometry,
    positions: Vec<[f64; 3]>,
}

fn sphere(level: u32) -> Sphere {
    let g = icosphere(level).unwrap();
    let mesh = g.topology().unwrap();
    let geometry = Geometry::from_positions(&mesh, &g.positions).unwrap();
    Sphere {
        mesh,
        geometry,
        positions: g.positions,
    }
}

/// The cap `{x · axis ≥ c}` is region 0.
#[derive(Clone, Copy)]
struct Case {
    axis: [f64; 3],
    c: f64,
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn random_axis(rng: &mut SplitMix64) -> [f64; 3] {
    loop {
        let v = [
            rng.uniform(-1.0, 1.0),
            rng.uniform(-1.0, 1.0),
            rng.uniform(-1.0, 1.0),
        ];
        let n = dot(v, v).sqrt();
        if n > 0.2 && n <= 1.0 {
            return v.map(|x| x / n);
        }
    }
}

/// For a segment: the fractions of its two points measured from the lone vertex of its face.
fn fractions(mesh: &TriMesh, cut: &Cut, segment: usize) -> [f64; 2] {
    let seg = cut.segments[segment];
    let face = mesh.faces()[seg.face as usize];
    let r = face.map(|v| cut.vertex_region[v as usize]);
    let k = if r[1] == r[2] {
        0
    } else if r[0] == r[2] {
        1
    } else {
        2
    };
    let a = face[k];
    seg.points.map(|p| {
        let point = cut.points[p as usize];
        if mesh.edges()[point.edge as usize][0] == a {
            point.t
        } else {
            1.0 - point.t
        }
    })
}

fn find(parent: &mut [u32], mut x: u32) -> u32 {
    while parent[x as usize] != x {
        parent[x as usize] = parent[parent[x as usize] as usize];
        x = parent[x as usize];
    }
    x
}

/// Per cut point: `[value on the side of the lower region, on the side of the higher region]`.
type PerSide = Vec<[f64; 2]>;

fn sides(e: &Evaluation, f: impl Fn(usize) -> f64) -> PerSide {
    (0..e.cut.points.len())
        .map(|k| e.triangles.point_vertices_by_region(k).map(&f))
        .collect()
}

/// (A) `v/ℓ` per cut vertex.
fn nodal(e: &Evaluation) -> PerSide {
    sides(e, |v| e.flux[v] / e.dual_length[v])
}

/// (C) pooled over runs of crossings within `theta` of the same mesh vertex.
fn clustered(s: &Sphere, e: &Evaluation, theta: f64) -> (PerSide, usize) {
    let n = e.cut.points.len();
    let mut parent: Vec<u32> = (0..n as u32).collect();
    for (i, seg) in e.cut.segments.iter().enumerate() {
        let [sp, sq] = fractions(&s.mesh, &e.cut, i);
        if sp < theta && sq < theta {
            let (a, b) = (
                find(&mut parent, seg.points[0]),
                find(&mut parent, seg.points[1]),
            );
            parent[a as usize] = b;
        }
    }
    let v = sides(e, |v| e.flux[v]);
    let l = sides(e, |v| e.dual_length[v]);
    let mut sum = vec![[0.0f64; 4]; n];
    for k in 0..n {
        let root = find(&mut parent, k as u32) as usize;
        sum[root][0] += v[k][0];
        sum[root][1] += v[k][1];
        sum[root][2] += l[k][0];
        sum[root][3] += l[k][1];
    }
    let mut merged = 0;
    let out = (0..n)
        .map(|k| {
            let root = find(&mut parent, k as u32) as usize;
            if root != k {
                merged += 1;
            }
            [sum[root][0] / sum[root][2], sum[root][1] / sum[root][3]]
        })
        .collect();
    (out, merged)
}

/// (B) gradient of the piecewise linear `u` in the triangle on each cut segment, brought to the
/// cut vertices by the length-weighted mean of the adjacent segments.
fn segment_gradient(e: &Evaluation) -> PerSide {
    let tri = &e.triangles;
    let n = tri.n_vertices();
    let (mut num, mut den) = (vec![0.0; n], vec![0.0; n]);
    for t in 0..tri.n_faces() {
        for k in 0..3 {
            if tri.edge_kind[t][k] == EdgeKind::Cut {
                let (i, j, apex) = (
                    tri.faces[t][k] as usize,
                    tri.faces[t][(k + 1) % 3] as usize,
                    tri.faces[t][(k + 2) % 3] as usize,
                );
                let length = tri.len2[t][k].sqrt();
                // Outward derivative: (0 − u_apex) / height, height = 2 area / length.
                let d = -e.log_scale[apex] * length / (2.0 * tri.area[t]);
                for v in [i, j] {
                    num[v] += 0.5 * length * d;
                    den[v] += 0.5 * length;
                }
            }
        }
    }
    sides(e, |v| num[v] / den[v])
}

/// (H) screened projection: `(diag(ℓ) + Σ_seg κ_seg [[1,−1],[−1,1]]) d = v` with
/// `κ_seg = c² h_f^(1+γ) / |seg|^γ`, `h_f²` the mean squared edge length of the parent face.
fn screened(s: &Sphere, e: &Evaluation, c: f64, gamma: f64) -> PerSide {
    let n = e.cut.points.len();
    if n == 0 {
        return Vec::new();
    }
    let v = sides(e, |v| e.flux[v]);
    let l = sides(e, |v| e.dual_length[v]);
    let lengths = s.geometry.edge_lengths();
    let scale = e.scale.length_scale;
    let mut out = vec![[0.0; 2]; n];
    for side in 0..2 {
        let mut triplets: Vec<(u32, u32, f64)> = (0..n).map(|k| (k as u32, k as u32, l[k][side])).collect();
        for (seg, length) in e.cut.segments.iter().zip(&e.triangles.segment_length) {
            let f = seg.face as usize;
            let h2: f64 = (0..3)
                .map(|k| (lengths[s.mesh.edge(3 * f + k)] * scale).powi(2))
                .sum::<f64>()
                / 3.0;
            let w = c * c * h2.powf(0.5 * (1.0 + gamma)) / length.powf(gamma);
            let [p, q] = seg.points;
            triplets.push((p, p, w));
            triplets.push((q, q, w));
            triplets.push((p.max(q), p.min(q), -w));
        }
        let m = SparseMatrix::from_triplets(n, n, &triplets).unwrap();
        let f = Symbolic::analyze(&m, Method::Llt)
            .unwrap()
            .factor(&m, Parallelism::Sequential)
            .unwrap();
        let mut rhs: Vec<f64> = v.iter().map(|x| x[side]).collect();
        f.solve_in_place(&mut rhs, Parallelism::Sequential).unwrap();
        for k in 0..n {
            out[k][side] = rhs[k];
        }
    }
    out
}

/// Neighbours of every cut point along the cut, with segment lengths.
fn neighbours(e: &Evaluation) -> Vec<Vec<(usize, f64)>> {
    let mut out = vec![Vec::new(); e.cut.points.len()];
    for (seg, &length) in e.cut.segments.iter().zip(&e.triangles.segment_length) {
        let [p, q] = seg.points.map(|x| x as usize);
        out[p].push((q, length));
        out[q].push((p, length));
    }
    out
}

/// Sliver values replaced by the arclength-linear interpolation between the nearest non-sliver
/// points on each side.
fn interpolate_slivers(e: &Evaluation, values: &[f64]) -> Vec<f64> {
    let nb = neighbours(e);
    let sliver = &e.triangles.sliver_point;
    let walk = |start: usize, first: (usize, f64)| -> Option<(f64, f64)> {
        let (mut prev, mut cur, mut dist) = (start, first.0, first.1);
        loop {
            if cur == start {
                return None;
            }
            if !sliver[cur] {
                return Some((dist, values[cur]));
            }
            let next = nb[cur].iter().find(|&&(x, _)| x != prev)?;
            prev = cur;
            cur = next.0;
            dist += next.1;
        }
    };
    (0..values.len())
        .map(|k| {
            if !sliver[k] {
                return values[k];
            }
            let found: Vec<(f64, f64)> = nb[k].iter().filter_map(|&f| walk(k, f)).collect();
            match found.as_slice() {
                [(d1, v1), (d2, v2)] => (d2 * v1 + d1 * v2) / (d1 + d2),
                [(_, v1)] => *v1,
                _ => values[k],
            }
        })
        .collect()
}

#[derive(Clone, Copy)]
enum Post {
    None,
    Smooth(usize),
    Interpolate,
}

fn post(e: &Evaluation, values: &[f64], post: Post) -> Vec<f64> {
    match post {
        Post::None => values.to_vec(),
        Post::Smooth(passes) => smooth_slivers(
            &e.cut,
            &e.triangles.segment_length,
            &e.triangles.sliver_point,
            values,
            passes,
        ),
        Post::Interpolate => interpolate_slivers(e, values),
    }
}

/// max and sum of relative errors, count.
#[derive(Clone, Copy, Default)]
struct Stat {
    max: f64,
    sum: f64,
    n: usize,
}

impl Stat {
    fn add(&mut self, value: f64, exact: f64, largest: f64) {
        let err = (value - exact).abs() / largest;
        let err = if err.is_finite() { err } else { f64::INFINITY };
        self.max = self.max.max(err);
        self.sum += err;
        self.n += 1;
    }
    fn mean(&self) -> f64 {
        self.sum / self.n.max(1) as f64
    }
}

const METHODS: [&str; 15] = [
    "A0 nodal v/l", "A1 +1 pass", "A3 +3 passes", "Ai +interp", "B  segment", "C  th=0.05",
    "C  th=0.10", "C  th=0.25", "C  th=0.40", "H  c=0.25", "H  c=0.35", "H  c=0.50",
    "H  c=1.00", "H2 c=0.25", "library default",
];

/// For each method: stats of d (lower side), d (higher side), s.
type Row = [[Stat; 3]; METHODS.len()];

/// `exact[k]` = exact `[d on the lower side, d on the higher side, s]` at cut point `k`.
fn accumulate(s: &Sphere, e: &Evaluation, exact: &[[f64; 3]], row: &mut Row) {
    let largest = [0, 1, 2].map(|j| exact.iter().map(|x| x[j].abs()).fold(0.0, f64::max));
    let a = nodal(e);
    let estimates: [(PerSide, Post); METHODS.len()] = [
        (a.clone(), Post::None),
        (a.clone(), Post::Smooth(1)),
        (a.clone(), Post::Smooth(3)),
        (a, Post::Interpolate),
        (segment_gradient(e), Post::None),
        (clustered(s, e, 0.05).0, Post::None),
        (clustered(s, e, 0.10).0, Post::None),
        (clustered(s, e, 0.25).0, Post::None),
        (clustered(s, e, 0.40).0, Post::None),
        (screened(s, e, 0.25, 1.0), Post::None),
        (screened(s, e, 0.35, 1.0), Post::None),
        (screened(s, e, 0.50, 1.0), Post::None),
        (screened(s, e, 1.00, 1.0), Post::None),
        (screened(s, e, 0.25, 2.0), Post::None),
        // What the library computes with its default options (FluxRecovery::Smoothed, c = 0.5).
        (sides(e, |v| e.normal_derivative[v]), Post::None),
    ];
    for (m, (d, p)) in estimates.iter().enumerate() {
        for side in 0..2 {
            let values: Vec<f64> = d.iter().map(|x| x[side]).collect();
            for (k, v) in post(e, &values, *p).into_iter().enumerate() {
                row[m][side].add(v, exact[k][side], largest[side]);
            }
        }
        let speed: Vec<f64> = d.iter().map(|x| x[0] * x[0] - x[1] * x[1]).collect();
        for (k, v) in post(e, &speed, *p).into_iter().enumerate() {
            row[m][2].add(v, exact[k][2], largest[2]);
        }
    }
}

#[derive(Default)]
struct Counts {
    cuts: usize,
    points: usize,
    clamped: usize,
    slivers: usize,
    merged: usize,
    smallest: f64,
}

fn print_table(title: &str, counts: &Counts, row: &Row) {
    println!();
    println!(
        "{title}: {} cuts, {} cut vertices, {} clamped, smallest offset {:.1e}, {} sliver vertices, \
         {} merged at th=0.25",
        counts.cuts, counts.points, counts.clamped, counts.smallest, counts.slivers, counts.merged
    );
    println!(
        "{:<16} {:>21} {:>21} {:>21}",
        "estimator", "d0 max / mean", "d1 max / mean", "s max / mean"
    );
    for (m, name) in METHODS.iter().enumerate() {
        let f = |st: Stat| format!("{:>10.3e} / {:>8.2e}", 100.0 * st.max, 100.0 * st.mean());
        println!("{:<16} {} {} {}", name, f(row[m][0]), f(row[m][1]), f(row[m][2]));
    }
}

fn count(s: &Sphere, e: &Evaluation, counts: &mut Counts) {
    counts.cuts += 1;
    counts.points += e.report.n_cut_points;
    counts.clamped += e.report.n_clamped_crossings;
    counts.slivers += e.report.n_sliver_points;
    counts.smallest = if counts.cuts == 1 {
        e.report.smallest_crossing_offset
    } else {
        counts.smallest.min(e.report.smallest_crossing_offset)
    };
    counts.merged += clustered(s, e, 0.25).1;
}

fn sphere_study(lo: u32, hi: u32, options: &EvalOptions<'_>) {
    let problem = CutProblem::new(Energy::Dirichlet, 0.0);
    for level in lo..=hi {
        let s = sphere(level);
        let mut evaluator = Evaluator::new(&s.mesh, &s.geometry).unwrap();
        let h = s.geometry.edge_lengths().iter().sum::<f64>() / s.mesh.n_edges() as f64;
        let mut rng = SplitMix64::new(1000 + u64::from(level));
        // A vertex whose z is closest to 0.55, for the latitude circles.
        let ring = (0..s.positions.len())
            .min_by(|&i, &j| {
                (s.positions[i][2] - 0.55)
                    .abs()
                    .total_cmp(&(s.positions[j][2] - 0.55).abs())
            })
            .unwrap();
        let z = [0.0, 0.0, 1.0];
        let mut families: Vec<(&str, Vec<Case>)> = Vec::new();
        families.push((
            "latitude through a vertex ring (phi = 0 at vertices)",
            vec![Case {
                axis: z,
                c: s.positions[ring][2],
            }],
        ));
        families.push((
            "latitude next to a vertex ring",
            [1e-6, -1e-6, 1e-4, -1e-4, 1e-3, -1e-3, 1e-2, -1e-2, 0.1, -0.1]
                .iter()
                .map(|e| Case {
                    axis: z,
                    c: s.positions[ring][2] + e * h,
                })
                .collect(),
        ));
        families.push((
            "random tilted planes",
            (0..20)
                .map(|_| Case {
                    axis: random_axis(&mut rng),
                    c: rng.uniform(0.3, 0.8),
                })
                .collect(),
        ));
        families.push((
            "tilted planes through a vertex (phi = 0 at a vertex)",
            (0..20)
                .map(|_| {
                    let axis = random_axis(&mut rng);
                    let target = rng.uniform(0.3, 0.8);
                    let v = (0..s.positions.len())
                        .min_by(|&i, &j| {
                            (dot(s.positions[i], axis) - target)
                                .abs()
                                .total_cmp(&(dot(s.positions[j], axis) - target).abs())
                        })
                        .unwrap();
                    Case {
                        axis,
                        c: dot(s.positions[v], axis),
                    }
                })
                .collect(),
        ));
        for (name, cases) in &families {
            let mut row: Row = [[Stat::default(); 3]; METHODS.len()];
            let mut counts = Counts::default();
            for case in cases {
                let phi: Vec<f64> = s
                    .positions
                    .iter()
                    .map(|p| case.c - dot(*p, case.axis))
                    .collect();
                let e = evaluator.evaluate(&phi, &problem, options).unwrap();
                let r = 1.0 / (4.0 * PI).sqrt();
                let psi = case.c.acos();
                let d = [(0.5 * psi).tan() / r, 1.0 / ((0.5 * psi).tan() * r)];
                let exact = vec![[d[0], d[1], d[0] * d[0] - d[1] * d[1]]; e.cut.points.len()];
                count(&s, &e, &mut counts);
                accumulate(&s, &e, &exact, &mut row);
            }
            print_table(
                &format!("SPHERE level {level} ({} vertices), {name}", s.positions.len()),
                &counts,
                &row,
            );
        }
    }
}

/// Radial derivative at `rho` of the solution of `div grad u = k0 + kappa cos(m θ)` on the
/// annulus `a ≤ r ≤ b` with `u = 0` on both circles: `(mode 0, mode m)`.
fn annulus_derivative(k0: f64, kappa: f64, m: f64, a: f64, b: f64, rho: f64) -> (f64, f64) {
    // Mode 0: u = k0 r²/4 + A ln r + B.
    let a0 = -k0 * (b * b - a * a) / (4.0 * (b / a).ln());
    let mode0 = 0.5 * k0 * rho + a0 / rho;
    // Mode m: f = c r² + A r^m + B r^−m, c = kappa/(4 − m²).
    let c = kappa / (4.0 - m * m);
    let (am, bm, an, bn) = (a.powf(m), b.powf(m), a.powf(-m), b.powf(-m));
    let det = am * bn - an * bm;
    let (ra, rb) = (-c * a * a, -c * b * b);
    let ca = (ra * bn - an * rb) / det;
    let cb = (am * rb - ra * bm) / det;
    let modem = 2.0 * c * rho + m * ca * rho.powf(m - 1.0) - m * cb * rho.powf(-m - 1.0);
    (mode0, modem)
}

fn annulus_study(lo: u32, hi: u32, options: &EvalOptions<'_>) {
    let (r_in, r_out, k0, kappa) = (0.4, 1.0, 4.0, 6.0);
    for level in lo..=hi {
        let scale_up = 1usize << (level - 3);
        let (n_around, n_radial) = (60 * scale_up, 8 * scale_up);
        let dr = (r_out - r_in) / n_radial as f64;
        let ring = r_out - dr * (5 * n_radial / 8) as f64; // = 0.625, a ring of vertices
        // Meshes: the regular polar grid and perturbed ones.
        let regular = annulus(n_around, n_radial, r_in, r_out).unwrap();
        let build = |positions: Vec<[f64; 3]>| -> Sphere {
            let mesh = regular.topology().unwrap();
            let geometry = Geometry::from_positions(&mesh, &positions).unwrap();
            Sphere {
                mesh,
                geometry,
                positions,
            }
        };
        let perturbed = |seed: u64| -> Vec<[f64; 3]> {
            let mut rng = SplitMix64::new(seed);
            let mut out = Vec::new();
            for i in 0..n_around {
                for j in 0..=n_radial {
                    let interior = j > 0 && j < n_radial;
                    let amp = if interior { 0.3 } else { 0.0 };
                    let r = r_out - dr * j as f64 + amp * dr * rng.uniform(-1.0, 1.0);
                    let theta =
                        2.0 * PI * (i as f64 + amp * rng.uniform(-1.0, 1.0)) / n_around as f64;
                    out.push([r * theta.cos(), r * theta.sin(), 0.0]);
                }
            }
            out
        };
        let families: Vec<(&str, Vec<(Sphere, f64)>)> = vec![
            (
                "regular polar grid, cut between two rings",
                vec![(build(regular.positions.clone()), ring + 0.37 * dr)],
            ),
            (
                "regular polar grid, cut on / next to a ring of vertices",
                [0.0, 1e-6, -1e-6, 1e-4, -1e-4, 1e-3, -1e-3, 1e-2, -1e-2, 0.1, -0.1]
                    .iter()
                    .map(|e| (build(regular.positions.clone()), ring + e * dr))
                    .collect(),
            ),
            (
                "perturbed grid",
                (0..6).map(|seed| (build(perturbed(seed)), ring)).collect(),
            ),
        ];
        for m in [4.0, 12.0] {
            for (name, cases) in &families {
                let mut row: Row = [[Stat::default(); 3]; METHODS.len()];
                let mut counts = Counts::default();
                let mut vertices = 0;
                for (s, rho) in cases {
                    vertices = s.positions.len();
                    let mut evaluator = Evaluator::new(&s.mesh, &s.geometry).unwrap();
                    let unit = evaluator.scale().length_scale;
                    let density: Vec<f64> = s
                        .mesh
                        .faces()
                        .iter()
                        .map(|face| {
                            let c = [0, 1].map(|k| {
                                face.iter().map(|&v| s.positions[v as usize][k]).sum::<f64>() / 3.0
                            });
                            (k0 + kappa * (m * c[1].atan2(c[0])).cos()) / (unit * unit)
                        })
                        .collect();
                    let problem = CutProblem {
                        curvature: CurvatureSource::PerFace(density),
                        ..CutProblem::new(Energy::Dirichlet, 0.0)
                    };
                    let phi: Vec<f64> = s
                        .positions
                        .iter()
                        .map(|p| (p[0] * p[0] + p[1] * p[1]).sqrt() - rho)
                        .collect();
                    let e = evaluator.evaluate(&phi, &problem, options).unwrap();
                    let inner = annulus_derivative(k0, kappa, m, r_in, *rho, *rho);
                    let outer = annulus_derivative(k0, kappa, m, *rho, r_out, *rho);
                    let exact: Vec<[f64; 3]> = e
                        .cut
                        .point_positions(&s.mesh, &s.positions)
                        .unwrap()
                        .iter()
                        .map(|p| {
                            let cos = (m * p[1].atan2(p[0])).cos();
                            let d0 = (inner.0 + inner.1 * cos) / unit;
                            let d1 = -(outer.0 + outer.1 * cos) / unit;
                            [d0, d1, d0 * d0 - d1 * d1]
                        })
                        .collect();
                    count(s, &e, &mut counts);
                    accumulate(s, &e, &exact, &mut row);
                }
                print_table(
                    &format!(
                        "ANNULUS {n_around} x {n_radial} ({vertices} vertices), m = {m}, {name}"
                    ),
                    &counts,
                    &row,
                );
            }
        }
    }
}

fn main() {
    let args: Vec<u32> = std::env::args()
        .skip(1)
        .map(|a| a.parse().expect("levels are integers"))
        .collect();
    let (lo, hi) = match args.as_slice() {
        [a, b] => (*a, *b),
        _ => (3, 5),
    };
    let options = EvalOptions::default();
    println!("errors in percent of the largest exact value on the cut: max / mean over all cut");
    println!("vertices of the family. d0 = du/dn on the side of region 0, d1 on the side of");
    println!("region 1, s = d0^2 - d1^2.");
    sphere_study(lo, hi, &options);
    annulus_study(lo, hi, &options);
}
