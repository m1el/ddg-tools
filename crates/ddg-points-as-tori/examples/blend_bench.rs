//! Single-threaded cost of the blend (WP4b): `blend_bench [n_points=20000] [n_queries=200000]`.
//! Synthetic ellipsoid cloud with osculating-sphere-like coefficients; prints queries/s for
//! distance only and for distance + gradient + Laplacian. Run alone.

use ddg_points_as_tori::{EvalOptions, EvalOut, NormalPolicy, PointCloud, Tori, TorusField, fit_tori, FitOptions};

fn lcg(s: &mut u64) -> f64 {
    *s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
    (*s >> 11) as f64 / (1u64 << 53) as f64
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let n: usize = a.get(1).map(|s| s.parse().unwrap()).unwrap_or(20000);
    let nq: usize = a.get(2).map(|s| s.parse().unwrap()).unwrap_or(200000);
    let mut s = 7;
    let (mut p, mut nr) = (Vec::new(), Vec::new());
    let ax = [0.6, 0.45, 0.3];
    while p.len() < n {
        let v = [lcg(&mut s) * 2.0 - 1.0, lcg(&mut s) * 2.0 - 1.0, lcg(&mut s) * 2.0 - 1.0];
        let l = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        if !(0.1..=1.0).contains(&l) {
            continue;
        }
        let u = [v[0] / l, v[1] / l, v[2] / l];
        p.push([u[0] * ax[0], u[1] * ax[1], u[2] * ax[2]]);
        let g = [u[0] / ax[0], u[1] / ax[1], u[2] / ax[2]];
        let gl = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
        nr.push([g[0] / gl, g[1] / gl, g[2] / gl]);
    }
    let (cloud, _) = PointCloud::new(p, nr, NormalPolicy::default()).unwrap();
    let co: Vec<[f32; 6]> = (0..n).map(|i| [0.0, 0.0, 0.0, (i % 7) as f32 * 0.03 - 0.5, 0.05, -1.5 + (i % 5) as f32 * 0.1]).collect();
    let (tori, _): (Tori, _) = fit_tori(&cloud, &co, &FitOptions::default()).unwrap();
    let field = TorusField::new(&cloud, tori, EvalOptions { threads: Some(1), ..Default::default() }).unwrap();
    let qs: Vec<[f64; 3]> = (0..nq).map(|_| [lcg(&mut s) - 0.5, lcg(&mut s) - 0.5, lcg(&mut s) * 0.6 - 0.3]).map(|q| [q[0] * 1.4, q[1] * 1.1, q[2]]).collect();
    let (mut d, mut g, mut l) = (vec![0.0; nq], vec![[0.0; 3]; nq], vec![0.0; nq]);
    for rep in 0..4 {
        let t0 = std::time::Instant::now();
        field.evaluate(&qs, EvalOut { distance: Some(&mut d), gradient: None, laplacian: None }).unwrap();
        let t1 = t0.elapsed().as_secs_f64();
        let t0 = std::time::Instant::now();
        field.evaluate(&qs, EvalOut { distance: Some(&mut d), gradient: Some(&mut g), laplacian: Some(&mut l) }).unwrap();
        let t2 = t0.elapsed().as_secs_f64();
        println!("rep {rep}: distance {:.0} q/s ({:.2} us/q), distance+gradient+Laplacian {:.0} q/s ({:.2} us/q); checksum {:.12e}", nq as f64 / t1, t1 / nq as f64 * 1e6, nq as f64 / t2, t2 / nq as f64 * 1e6, d.iter().sum::<f64>() + l.iter().sum::<f64>());
    }
}
