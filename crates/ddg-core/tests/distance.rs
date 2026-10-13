//! Tests of the upwind (Hopf–Lax) distance solver: fixed point, plane waves, point sources,
//! classes, relabeling, bounds and bad input.

// Index loops mirror the subscripts of the formulas being checked.
#![allow(clippy::needless_range_loop)]

use ddg_core::distance::{DistanceError, triangle_update, upwind_distance};
use ddg_core::generators::{self, GeneratedMesh, SplitMix64};
use ddg_core::{Geometry, TriMesh};
use proptest::prelude::*;

fn build(positions: &[[f64; 3]], faces: &[[u32; 3]]) -> (TriMesh, Geometry) {
    let mesh = TriMesh::new(positions.len(), faces.to_vec()).unwrap();
    let geometry = Geometry::from_positions(&mesh, positions).unwrap();
    (mesh, geometry)
}

fn grid(n: usize) -> GeneratedMesh {
    generators::grid_disk(n, n).unwrap()
}

/// Flat grid sheared in x: `x ← x + s y`. Largest angle: 117° for `s = 0.5`, 135° for 1, 153° for 2, 162° for 3.
fn sheared_grid(n: usize, s: f64) -> GeneratedMesh {
    let mut m = grid(n);
    for p in &mut m.positions {
        p[0] += s * p[1];
    }
    m
}

/// Flat grid squeezed in y: needle-like triangles.
fn needle_grid(n: usize, squeeze: f64) -> GeneratedMesh {
    let mut m = grid(n);
    for p in &mut m.positions {
        p[1] *= squeeze;
    }
    m
}

fn irregular_grid(n: usize, seed: u64) -> GeneratedMesh {
    let mut m = grid(n);
    let a = 0.2 / n as f64;
    generators::perturb(&mut m.positions, [a, a, 0.0], seed);
    m
}

fn bumpy(n: usize, seed: u64) -> GeneratedMesh {
    let mut m = generators::bumpy_grid_disk(n, n, 0.3).unwrap();
    let a = 0.15 / n as f64;
    generators::perturb(&mut m.positions, [a, a, a], seed);
    m
}

fn close(a: f64, b: f64, tol: f64) -> bool {
    if a.is_infinite() || b.is_infinite() {
        return a == b;
    }
    (a - b).abs() <= tol * b.abs().max(a.abs()).max(1e-300) + 1e-15
}

/// The update of vertex `i` recomputed from scratch with the public local update.
fn reference_update(
    mesh: &TriMesh,
    geometry: &Geometry,
    values: &[f64],
    class: Option<&[u32]>,
    i: usize,
) -> f64 {
    let usable =
        |j: u32| values[j as usize].is_finite() && class.is_none_or(|c| c[j as usize] == c[i]);
    let len =
        |a: u32, b: u32| geometry.edge_lengths()[mesh.find_edge(a as usize, b as usize).unwrap()];
    let mut best = f64::INFINITY;
    for face in mesh.faces() {
        let Some(c) = face.iter().position(|&v| v as usize == i) else {
            continue;
        };
        let (j, k) = (face[(c + 1) % 3], face[(c + 2) % 3]);
        let (lij, lik, ljk) = (len(i as u32, j), len(i as u32, k), len(j, k));
        let cand = match (usable(j), usable(k)) {
            (true, true) => {
                triangle_update(values[j as usize], values[k as usize], lij, lik, ljk).unwrap()
            }
            (true, false) => values[j as usize] + lij,
            (false, true) => values[k as usize] + lik,
            (false, false) => f64::INFINITY,
        };
        best = best.min(cand);
    }
    best
}

fn assert_fixed_point(
    mesh: &TriMesh,
    geometry: &Geometry,
    seeds: &[(u32, f64)],
    class: Option<&[u32]>,
    what: &str,
) -> usize {
    let r = upwind_distance(mesh, geometry, seeds, class).unwrap();
    assert!(
        r.converged,
        "{what}: not converged after {} re-openings",
        r.n_reopened
    );
    let is_seed = |v: usize| seeds.iter().any(|s| s.0 as usize == v);
    for &(v, x) in seeds {
        assert_eq!(r.values[v as usize], x, "{what}: seed value changed");
    }
    for i in 0..mesh.n_vertices() {
        if is_seed(i) {
            continue;
        }
        let want = reference_update(mesh, geometry, &r.values, class, i);
        assert!(
            close(r.values[i], want, 1e-12),
            "{what}: vertex {i}: value {} but update {want}",
            r.values[i]
        );
    }
    r.n_reopened
}

// ---------------------------------------------------------------------------------------------
// 1. Fixed point

#[test]
fn fixed_point_on_regular_sheared_needle_and_random_meshes() {
    let n = 10;
    let meshes: Vec<(&str, GeneratedMesh)> = vec![
        ("regular", grid(n)),
        ("sheared 2", sheared_grid(n, 2.0)),
        ("sheared 2.5", sheared_grid(n, 2.5)),
        ("needle", needle_grid(n, 0.02)),
        ("irregular", irregular_grid(n, 7)),
        ("bumpy", bumpy(n, 3)),
        ("icosphere", generators::icosphere(3).unwrap()),
    ];
    for (name, m) in meshes {
        let (mesh, geometry) = build(&m.positions, &m.faces);
        let nv = mesh.n_vertices() as u32;
        let seeds = [(0, 0.0), (nv / 2, 0.0), (nv - 1, 0.1)];
        let reopened = assert_fixed_point(&mesh, &geometry, &seeds, None, name);
        println!("{name}: {nv} vertices, {reopened} re-openings");
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]

    #[test]
    fn fixed_point_random(
        n in 3usize..12,
        shear in 0.0..2.0f64,
        mseed in any::<u64>(),
        picks in proptest::collection::vec((any::<u32>(), 0.0..2.0f64), 1..5),
    ) {
        let mut m = sheared_grid(n, shear);
        let a = 0.15 / n as f64;
        generators::perturb(&mut m.positions, [a, a, a], mseed);
        let (mesh, geometry) = build(&m.positions, &m.faces);
        let nv = mesh.n_vertices() as u32;
        let mut seeds: Vec<(u32, f64)> = Vec::new();
        for (v, x) in picks {
            let v = v % nv;
            if seeds.iter().all(|s| s.0 != v) {
                seeds.push((v, x));
            }
        }
        assert_fixed_point(&mesh, &geometry, &seeds, None, "random");
    }

    #[test]
    fn local_update_properties(
        a in 0.2..3.0f64, b in 0.2..3.0f64, c in 0.2..3.0f64,
        tj in 0.0..5.0f64, tk in 0.0..5.0f64, d in 0.0..1.0f64,
    ) {
        prop_assume!(a + b > c * 1.001 && a + c > b * 1.001 && b + c > a * 1.001);
        let u = triangle_update(tj, tk, a, b, c).unwrap();
        let ends = (tj + a).min(tk + b);
        prop_assert!(u <= ends);
        prop_assert!(u >= tj.min(tk));
        // Monotone in both arguments.
        prop_assert!(triangle_update(tj + d, tk, a, b, c).unwrap() >= u);
        prop_assert!(triangle_update(tj, tk + d, a, b, c).unwrap() >= u);
        // Symmetric under swapping j and k.
        let sw = triangle_update(tk, tj, b, a, c).unwrap();
        prop_assert!((sw - u).abs() <= 1e-12 * u.max(1.0));
    }
}

// ---------------------------------------------------------------------------------------------
// 2. Plane wave

/// Exact values `n·x` seeded on the first two vertex layers (`i = 0, 1` of the grid); returns
/// the maximum error over vertices with a verified upwind chain, and how many there are.
fn plane_wave(m: &GeneratedMesh, n: usize, angle: f64) -> (f64, usize, usize, f64) {
    let (mesh, geometry) = build(&m.positions, &m.faces);
    let dir = [angle.cos(), angle.sin()];
    let exact: Vec<f64> = m
        .positions
        .iter()
        .map(|p| dir[0] * p[0] + dir[1] * p[1] + 1.0)
        .collect();
    // Seed the first two layers on the side the wave comes from: grid columns i = 0, 1 for
    // waves moving mostly along x, rows j = 0, 1 for waves moving mostly along y.
    let along_x = dir[0] >= dir[1];
    let seeds: Vec<(u32, f64)> = (0..exact.len() as u32)
        .filter(|&v| {
            let (i, j) = (v as usize / (n + 1), v as usize % (n + 1));
            if along_x { i < 2 } else { j < 2 }
        })
        .map(|v| (v, exact[v as usize]))
        .collect();
    let r = upwind_distance(&mesh, &geometry, &seeds, None).unwrap();
    // "Good" vertices: some triangle (i, j, k) with good j, k whose edge jk contains the foot
    // of the characteristic through i. Their value must be exact. Computed by iteration.
    let nv = mesh.n_vertices();
    let mut good = vec![false; nv];
    for s in &seeds {
        good[s.0 as usize] = true;
    }
    let mut changed = true;
    while changed {
        changed = false;
        for face in mesh.faces() {
            for c in 0..3 {
                let (i, j, k) = (
                    face[c] as usize,
                    face[(c + 1) % 3] as usize,
                    face[(c + 2) % 3] as usize,
                );
                if good[i] || !good[j] || !good[k] {
                    continue;
                }
                // Solve p_i - s dir = q_j + t (q_k - q_j).
                let (p, qj, qk) = (m.positions[i], m.positions[j], m.positions[k]);
                let e = [qk[0] - qj[0], qk[1] - qj[1]];
                let rhs = [p[0] - qj[0], p[1] - qj[1]];
                // [dir, e] [s, -t]^T = rhs  -> 2x2 solve
                let det = dir[0] * e[1] - dir[1] * e[0];
                if det.abs() < 1e-14 {
                    continue;
                }
                let s = (rhs[0] * e[1] - rhs[1] * e[0]) / det;
                let t = (dir[0] * rhs[1] - dir[1] * rhs[0]) / det;
                if s > 0.0 && (0.0..=1.0).contains(&t) {
                    good[i] = true;
                    changed = true;
                }
            }
        }
    }
    let mut max_err: f64 = 0.0;
    let mut worst_below: f64 = 0.0;
    let mut n_good = 0;
    for i in 0..nv {
        let err = r.values[i] - exact[i];
        worst_below = worst_below.min(err);
        if good[i] {
            n_good += 1;
            max_err = max_err.max(err.abs());
        }
    }
    (max_err, n_good, nv, worst_below)
}

#[test]
fn plane_wave_is_exact_where_an_upwind_triangle_exists() {
    let n = 24;
    let cases: Vec<(&str, GeneratedMesh)> = vec![
        ("regular", grid(n)),
        ("irregular 1", irregular_grid(n, 11)),
        ("irregular 2", irregular_grid(n, 12)),
        ("sheared 3", {
            // Shear after perturbing; the plane wave is still n·x.
            let mut m = irregular_grid(n, 5);
            for p in &mut m.positions {
                p[0] += 3.0 * p[1];
            }
            m
        }),
    ];
    for (name, m) in &cases {
        // The seed layers of the sheared mesh are slanted lines: only shallow waves leave them.
        let angles: &[f64] = if name.starts_with("sheared") {
            &[4.0, 12.0]
        } else {
            &[10.0, 33.0, 57.0]
        };
        for &angle_deg in angles {
            let (err, n_good, nv, below) = plane_wave(m, n, angle_deg.to_radians());
            println!(
                "plane wave {name} {angle_deg} deg: max error {err:.2e} on {n_good}/{nv} vertices, lowest signed error {below:.2e}"
            );
            assert!(err < 1e-13, "{name} {angle_deg}: error {err}");
            // Never below the exact value (T is 1-Lipschitz, the update cannot undercut it).
            assert!(below > -1e-13);
            // For a sheared mesh the wave direction is skewed relative to the layers, so the
            // fraction of verified vertices is smaller.
            assert!(
                n_good as f64 > 0.4 * nv as f64,
                "{name}: only {n_good}/{nv} verified"
            );
        }
    }
}

// ---------------------------------------------------------------------------------------------
// 3. Point sources

#[test]
fn point_source_on_a_flat_grid() {
    println!("point source on a flat grid (centre vertex), max and rms error");
    println!(
        "{:>6} {:>10} {:>12} {:>7} {:>12} {:>7}",
        "n", "h", "max err", "order", "rms err", "order"
    );
    let mut prev: Option<(f64, f64)> = None;
    let mut first_max = 0.0;
    let mut last_max = 0.0;
    for (idx, n) in [8usize, 16, 32, 64, 128].into_iter().enumerate() {
        let m = grid(n);
        let (mesh, geometry) = build(&m.positions, &m.faces);
        let centre = (n / 2) * (n + 1) + n / 2;
        let r = upwind_distance(&mesh, &geometry, &[(centre as u32, 0.0)], None).unwrap();
        assert_eq!(r.n_unreached, 0);
        let c = m.positions[centre];
        let (mut max_err, mut sq) = (0.0f64, 0.0f64);
        for (p, &t) in m.positions.iter().zip(&r.values) {
            let d = ((p[0] - c[0]).powi(2) + (p[1] - c[1]).powi(2)).sqrt();
            let e = (t - d).abs();
            max_err = max_err.max(e);
            sq += e * e;
        }
        let rms = (sq / m.positions.len() as f64).sqrt();
        let order = |now: f64, before: Option<f64>| match before {
            Some(b) => format!("{:7.2}", (b / now).log2()),
            None => "      -".to_string(),
        };
        println!(
            "{:>6} {:>10.5} {:>12.4e} {} {:>12.4e} {}",
            n,
            1.0 / n as f64,
            max_err,
            order(max_err, prev.map(|p| p.0)),
            rms,
            order(rms, prev.map(|p| p.1))
        );
        prev = Some((max_err, rms));
        if idx == 0 {
            first_max = max_err;
        }
        last_max = max_err;
    }
    // A point source is the worst case for a first-order scheme: the error is concentrated in
    // the cone of directions at the source and decays slower than h (observed order 0.6 rising
    // to 0.75, cf. the smooth-start tests below, which show order 1).
    assert!(last_max < first_max / 4.0, "{first_max} -> {last_max}");
}

fn great_circle(p: [f64; 3], q: [f64; 3]) -> f64 {
    let dot = p[0] * q[0] + p[1] * q[1] + p[2] * q[2];
    let cx = p[1] * q[2] - p[2] * q[1];
    let cy = p[2] * q[0] - p[0] * q[2];
    let cz = p[0] * q[1] - p[1] * q[0];
    (cx * cx + cy * cy + cz * cz).sqrt().atan2(dot)
}

#[test]
fn geodesic_distance_on_icospheres() {
    println!("geodesic distance from vertex 0 on icosphere(level) vs. great circle");
    println!(
        "{:>5} {:>8} {:>10} {:>12} {:>7} {:>12} {:>7} {:>9}",
        "level", "verts", "h", "max err", "order", "rms err", "order", "reopened"
    );
    let mut prev: Option<(f64, f64)> = None;
    let mut maxes = Vec::new();
    for level in 2..=6u32 {
        let m = generators::icosphere(level).unwrap();
        let (mesh, geometry) = build(&m.positions, &m.faces);
        let r = upwind_distance(&mesh, &geometry, &[(0, 0.0)], None).unwrap();
        assert_eq!(r.n_unreached, 0);
        let p0 = m.positions[0];
        let (mut max_err, mut sq) = (0.0f64, 0.0f64);
        for (p, &t) in m.positions.iter().zip(&r.values) {
            let e = (t - great_circle(p0, *p)).abs();
            max_err = max_err.max(e);
            sq += e * e;
        }
        let rms = (sq / m.positions.len() as f64).sqrt();
        let h = geometry.edge_lengths().iter().cloned().fold(0.0, f64::max);
        let order = |now: f64, before: Option<f64>| match before {
            Some(b) => format!("{:7.2}", (b / now).log2()),
            None => "      -".to_string(),
        };
        println!(
            "{:>5} {:>8} {:>10.5} {:>12.4e} {} {:>12.4e} {} {:>9}",
            level,
            m.positions.len(),
            h,
            max_err,
            order(max_err, prev.map(|p| p.0)),
            rms,
            order(rms, prev.map(|p| p.1)),
            r.n_reopened
        );
        prev = Some((max_err, rms));
        maxes.push(max_err);
    }
    // Four halvings of h; a point source converges with observed order 0.4 rising to 0.7.
    assert!(maxes[0] / maxes[4] > 3.0, "{maxes:?}");
    for w in maxes.windows(2) {
        assert!(w[1] < w[0], "{maxes:?}");
    }
}

/// Seeds on a disc around `centre` carry the exact distance; the front then expands without the
/// singularity of a point source, and the error of the first-order scheme is O(h).
#[test]
fn smooth_start_converges_with_first_order() {
    // Flat grid: disc of radius 0.15 around the centre.
    println!("flat grid, exact seeds within 0.15 of the centre");
    println!(
        "{:>6} {:>12} {:>7} {:>12} {:>7}",
        "n", "max err", "order", "rms err", "order"
    );
    let mut prev: Option<(f64, f64)> = None;
    let mut orders = Vec::new();
    for n in [16usize, 32, 64, 128] {
        let m = grid(n);
        let (mesh, geometry) = build(&m.positions, &m.faces);
        let c = m.positions[(n / 2) * (n + 1) + n / 2];
        let dist = |p: [f64; 3]| ((p[0] - c[0]).powi(2) + (p[1] - c[1]).powi(2)).sqrt();
        let seeds: Vec<(u32, f64)> = (0..m.positions.len())
            .filter(|&v| dist(m.positions[v]) < 0.15)
            .map(|v| (v as u32, dist(m.positions[v])))
            .collect();
        let r = upwind_distance(&mesh, &geometry, &seeds, None).unwrap();
        let (mut mx, mut sq, mut cnt) = (0.0f64, 0.0f64, 0usize);
        for (v, p) in m.positions.iter().enumerate() {
            if dist(*p) < 0.15 {
                continue;
            }
            let e = (r.values[v] - dist(*p)).abs();
            mx = mx.max(e);
            sq += e * e;
            cnt += 1;
        }
        let rms = (sq / cnt as f64).sqrt();
        let (o1, o2) = prev.map_or((f64::NAN, f64::NAN), |p| {
            ((p.0 / mx).log2(), (p.1 / rms).log2())
        });
        println!("{n:>6} {mx:>12.4e} {o1:>7.2} {rms:>12.4e} {o2:>7.2}");
        if prev.is_some() {
            orders.push(o1.min(o2));
        }
        prev = Some((mx, rms));
    }
    assert!(orders.iter().all(|&o| o > 0.85), "{orders:?}");

    // Icosphere: exact great-circle values within 0.4 radians of vertex 0.
    println!("icosphere, exact seeds within 0.4 rad of vertex 0");
    println!(
        "{:>5} {:>12} {:>7} {:>12} {:>7}",
        "level", "max err", "order", "rms err", "order"
    );
    let mut prev: Option<(f64, f64)> = None;
    let mut orders = Vec::new();
    for level in 2..=6u32 {
        let m = generators::icosphere(level).unwrap();
        let (mesh, geometry) = build(&m.positions, &m.faces);
        let p0 = m.positions[0];
        let seeds: Vec<(u32, f64)> = (0..m.positions.len())
            .filter(|&v| great_circle(p0, m.positions[v]) < 0.4)
            .map(|v| (v as u32, great_circle(p0, m.positions[v])))
            .collect();
        let r = upwind_distance(&mesh, &geometry, &seeds, None).unwrap();
        let (mut mx, mut sq, mut cnt) = (0.0f64, 0.0f64, 0usize);
        for (v, p) in m.positions.iter().enumerate() {
            let d = great_circle(p0, *p);
            if d < 0.4 {
                continue;
            }
            let e = (r.values[v] - d).abs();
            mx = mx.max(e);
            sq += e * e;
            cnt += 1;
        }
        let rms = (sq / cnt as f64).sqrt();
        let (o1, o2) = prev.map_or((f64::NAN, f64::NAN), |p| {
            ((p.0 / mx).log2(), (p.1 / rms).log2())
        });
        println!("{level:>5} {mx:>12.4e} {o1:>7.2} {rms:>12.4e} {o2:>7.2}");
        if prev.is_some() {
            orders.push(o1.min(o2));
        }
        prev = Some((mx, rms));
    }
    // The rms error is first order; the maximum error sits at the antipode, where the fronts
    // meet (a singularity of the distance), and converges more slowly (observed 0.6 to 0.7).
    assert!(orders.iter().all(|&o| o > 0.5), "{orders:?}");
}

// ---------------------------------------------------------------------------------------------
// 4. Classes

#[test]
fn class_restriction_isolates_the_two_sides() {
    let n = 20;
    let m = irregular_grid(n, 21);
    let (mesh, geometry) = build(&m.positions, &m.faces);
    let nv = mesh.n_vertices();
    // Two classes separated by a slanted line; seeds in both.
    let class: Vec<u32> = m
        .positions
        .iter()
        .map(|p| u32::from(p[1] > 0.3 + 0.4 * p[0]))
        .collect();
    let pick = |c: u32, k: usize| (0..nv).filter(|&v| class[v] == c).nth(k).unwrap() as u32;
    let seeds_a = [(pick(0, 3), 0.0), (pick(0, 90), 0.2)];
    let seeds_b1 = [(pick(1, 5), 0.0)];
    let seeds_b2 = [(pick(1, 40), 0.5), (pick(1, 77), 0.0)];
    let run = |a: &[(u32, f64)], b: &[(u32, f64)]| {
        let mut s = a.to_vec();
        s.extend_from_slice(b);
        upwind_distance(&mesh, &geometry, &s, Some(&class))
            .unwrap()
            .values
    };
    let t1 = run(&seeds_a, &seeds_b1);
    let t2 = run(&seeds_a, &seeds_b2);
    let t3 = run(&seeds_a, &[]);
    for v in 0..nv {
        if class[v] == 0 {
            assert_eq!(t1[v], t2[v], "class 0 changed at {v}");
            assert_eq!(t1[v], t3[v]);
        } else {
            assert!(t3[v].is_infinite(), "class 1 reached without a seed");
            assert!(t1[v].is_finite() && t2[v].is_finite());
        }
    }
    // And the other way round.
    let u1 = run(&[seeds_a[0]], &seeds_b2);
    let u2 = run(&[seeds_a[1]], &seeds_b2);
    for v in 0..nv {
        if class[v] == 1 {
            assert_eq!(u1[v], u2[v]);
        }
    }
    // The classes are fixed points of their own restricted scheme.
    let mut all = seeds_a.to_vec();
    all.extend_from_slice(&seeds_b2);
    assert_fixed_point(&mesh, &geometry, &all, Some(&class), "classes");
}

// ---------------------------------------------------------------------------------------------
// 5. Relabeling

#[test]
fn relabeling_does_not_change_the_result() {
    let n = 12;
    let meshes: Vec<(&str, GeneratedMesh)> = vec![
        ("sheared 2.5", sheared_grid(n, 2.5)),
        ("needle", needle_grid(n, 0.03)),
        ("bumpy", bumpy(n, 9)),
        ("icosphere", generators::icosphere(3).unwrap()),
    ];
    for (name, m) in meshes {
        let (mesh, geometry) = build(&m.positions, &m.faces);
        let nv = mesh.n_vertices();
        let mut rng = SplitMix64::new(99);
        let mut perm: Vec<u32> = (0..nv as u32).collect();
        for i in (1..nv).rev() {
            perm.swap(i, (rng.next_u64() % (i as u64 + 1)) as usize);
        }
        // Permute vertices, shuffle faces, rotate each face's corner list.
        let mut positions = vec![[0.0; 3]; nv];
        for v in 0..nv {
            positions[perm[v] as usize] = m.positions[v];
        }
        let mut faces: Vec<[u32; 3]> = m
            .faces
            .iter()
            .map(|f| {
                let r = (rng.next_u64() % 3) as usize;
                [
                    perm[f[r] as usize],
                    perm[f[(r + 1) % 3] as usize],
                    perm[f[(r + 2) % 3] as usize],
                ]
            })
            .collect();
        for i in (1..faces.len()).rev() {
            faces.swap(i, (rng.next_u64() % (i as u64 + 1)) as usize);
        }
        let (mesh2, geometry2) = build(&positions, &faces);
        let seeds = [(1u32, 0.0), (nv as u32 / 3, 0.05), (nv as u32 - 2, 0.0)];
        let seeds2: Vec<(u32, f64)> = seeds.iter().map(|&(v, x)| (perm[v as usize], x)).collect();
        let r1 = upwind_distance(&mesh, &geometry, &seeds, None).unwrap();
        let r2 = upwind_distance(&mesh2, &geometry2, &seeds2, None).unwrap();
        for v in 0..nv {
            assert!(
                close(r1.values[v], r2.values[perm[v] as usize], 1e-12),
                "{name}: vertex {v}: {} vs {}",
                r1.values[v],
                r2.values[perm[v] as usize]
            );
        }
        println!("{name}: re-openings {} / {}", r1.n_reopened, r2.n_reopened);
    }
}

// ---------------------------------------------------------------------------------------------
// 6. Bounds, re-openings, bad input

#[test]
fn neighbour_bound_and_reopening_counts() {
    let n = 24;
    let meshes: Vec<(&str, GeneratedMesh)> = vec![
        ("regular", grid(n)),
        ("sheared 2", sheared_grid(n, 2.0)),
        ("sheared 3", sheared_grid(n, 3.0)),
        ("needle", needle_grid(n, 0.02)),
        ("bumpy", bumpy(n, 4)),
    ];
    for (name, m) in meshes {
        let (mesh, geometry) = build(&m.positions, &m.faces);
        let centre = (mesh.n_vertices() / 2) as u32;
        let seeds = [(centre, 0.0)];
        let r = upwind_distance(&mesh, &geometry, &seeds, None).unwrap();
        for i in 0..mesh.n_vertices() {
            if i as u32 == centre {
                continue;
            }
            for j in mesh.vertex_neighbors(i) {
                let l = geometry.edge_lengths()[mesh.find_edge(i, j).unwrap()];
                assert!(
                    r.values[i] <= r.values[j] + l + 1e-14,
                    "{name}: T[{i}] = {} > T[{j}] + l = {}",
                    r.values[i],
                    r.values[j] + l
                );
            }
        }
        println!(
            "{name}: {} vertices, {} re-openings",
            mesh.n_vertices(),
            r.n_reopened
        );
    }
}

#[test]
fn bad_input_is_an_error() {
    let m = grid(3);
    let (mesh, geometry) = build(&m.positions, &m.faces);
    let n = mesh.n_vertices();
    let run = |seeds: &[(u32, f64)], class: Option<&[u32]>| {
        upwind_distance(&mesh, &geometry, seeds, class)
    };
    assert!(matches!(
        run(&[(n as u32, 0.0)], None),
        Err(DistanceError::SeedOutOfRange { seed: 0, vertex, n_vertices }) if vertex as usize == n && n_vertices == n
    ));
    assert!(matches!(
        run(&[(1, 0.0), (2, 0.0), (1, 3.0)], None),
        Err(DistanceError::DuplicateSeed {
            vertex: 1,
            first: 0,
            second: 2
        })
    ));
    for bad in [-1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            run(&[(0, 0.0), (4, bad)], None),
            Err(DistanceError::InvalidSeedValue {
                seed: 1,
                vertex: 4,
                ..
            })
        ));
    }
    assert!(matches!(
        run(&[(0, 0.0)], Some(&[0, 1])),
        Err(DistanceError::ClassCountMismatch { expected, got: 2 }) if expected == n
    ));
    let other = grid(5);
    let (_, other_geometry) = build(&other.positions, &other.faces);
    assert!(matches!(
        upwind_distance(&mesh, &other_geometry, &[(0, 0.0)], None),
        Err(DistanceError::Geometry(_))
    ));
    // Messages name the element.
    let msg = run(&[(n as u32, 0.0)], None).unwrap_err().to_string();
    assert!(msg.contains(&format!("vertex {n}")), "{msg}");
    // No seeds: everything is unreached, no error.
    let r = run(&[], None).unwrap();
    assert_eq!(r.n_unreached, n);
    assert_eq!(r.n_reopened, 0);
}

#[test]
fn disconnected_components_and_zero_values() {
    // Two triangles sharing no vertex: the second one is unreachable.
    let positions = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [5.0, 0.0, 0.0],
        [6.0, 0.0, 0.0],
        [5.0, 1.0, 0.0],
    ];
    let faces = vec![[0, 1, 2], [3, 4, 5]];
    let (mesh, geometry) = build(&positions, &faces);
    let r = upwind_distance(&mesh, &geometry, &[(0, 0.0)], None).unwrap();
    assert_eq!(r.n_unreached, 3);
    assert_eq!(&r.values[..3], &[0.0, 1.0, 1.0]);
    assert!(r.values[3..].iter().all(|x| x.is_infinite()));
}
