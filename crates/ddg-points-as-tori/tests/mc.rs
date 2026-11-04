//! Marching cubes (SPEC.md §8.6) and CSG (§8.2) on analytic SDFs. Run with `--release`.

use ddg_points_as_tori::PatError;
use ddg_points_as_tori::csg::{Csg, ImplicitField};
use ddg_points_as_tori::mc::{Band, McMesh, McOptions, MeshCheck, Resolution, marching_cubes};

struct Sphere(f64);
impl ImplicitField for Sphere {
    fn distance(&self, x: [f64; 3]) -> Result<f64, PatError> {
        Ok((x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt() - self.0)
    }
    fn distance_gradient(&self, x: [f64; 3]) -> Result<(f64, [f64; 3]), PatError> {
        let l = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
        Ok((l - self.0, [x[0] / l, x[1] / l, x[2] / l]))
    }
    fn bbox(&self) -> [[f64; 3]; 2] {
        [[-1.0; 3], [1.0; 3]]
    }
}

struct Torus(f64, f64);
impl ImplicitField for Torus {
    fn distance(&self, x: [f64; 3]) -> Result<f64, PatError> {
        Ok(self.distance_gradient(x)?.0)
    }
    fn distance_gradient(&self, x: [f64; 3]) -> Result<(f64, [f64; 3]), PatError> {
        let q = (x[0] * x[0] + x[1] * x[1]).sqrt();
        let (dq, dz) = (q - self.0, x[2]);
        let l = (dq * dq + dz * dz).sqrt();
        let g = [dq / l * x[0] / q, dq / l * x[1] / q, dz / l];
        Ok((l - self.1, g))
    }
    fn bbox(&self) -> [[f64; 3]; 2] {
        [[-1.0; 3], [1.0; 3]]
    }
}

struct Cuboid([f64; 3], [f64; 3]); // center, half extents
impl ImplicitField for Cuboid {
    fn distance(&self, x: [f64; 3]) -> Result<f64, PatError> {
        Ok(self.distance_gradient(x)?.0)
    }
    fn distance_gradient(&self, x: [f64; 3]) -> Result<(f64, [f64; 3]), PatError> {
        let q = [0, 1, 2].map(|i| (x[i] - self.0[i]).abs() - self.1[i]);
        let o = q.map(|v| v.max(0.0));
        let lo = (o[0] * o[0] + o[1] * o[1] + o[2] * o[2]).sqrt();
        let m = q[0].max(q[1]).max(q[2]);
        let d = lo + m.min(0.0);
        let s = |i: usize| if x[i] >= self.0[i] { 1.0 } else { -1.0 };
        let g = if lo > 0.0 {
            [0, 1, 2].map(|i| o[i] / lo * s(i))
        } else {
            let a = (0..3).max_by(|a, b| q[*a].partial_cmp(&q[*b]).unwrap()).unwrap();
            let mut g = [0.0; 3];
            g[a] = s(a);
            g
        };
        Ok((d, g))
    }
    fn bbox(&self) -> [[f64; 3]; 2] {
        [[0, 1, 2].map(|i| self.0[i] - self.1[i] - 0.1), [0, 1, 2].map(|i| self.0[i] + self.1[i] + 0.1)]
    }
}

/// Values prescribed on the grid points of a unit box; the field is a lookup by rounding.
pub struct Lattice {
    pub n: [usize; 3],
    pub vals: Vec<f64>,
}
impl ImplicitField for Lattice {
    fn distance(&self, x: [f64; 3]) -> Result<f64, PatError> {
        let i = [0, 1, 2].map(|a| (x[a] * self.n[a] as f64).round() as usize);
        Ok(self.vals[i[0] + (self.n[0] + 1) * (i[1] + (self.n[1] + 1) * i[2])])
    }
    fn distance_gradient(&self, x: [f64; 3]) -> Result<(f64, [f64; 3]), PatError> {
        Ok((self.distance(x)?, [0.0; 3]))
    }
    fn bbox(&self) -> [[f64; 3]; 2] {
        [[0.0; 3], [1.0; 3]]
    }
}

fn opts(res: usize) -> McOptions {
    McOptions { resolution: Resolution::Longest(res), threads: Some(4), ..Default::default() }
}

fn manifold(rep: &ddg_points_as_tori::mc::McReport) -> (usize, i64, usize) {
    match &rep.mesh {
        MeshCheck::Manifold { boundary_edges, euler_characteristic, components } => (*boundary_edges, *euler_characteristic, *components),
        other => panic!("mesh is not manifold: {other:?}"),
    }
}

#[test]
fn sphere_is_watertight_with_euler_2() {
    let s = Sphere(0.6);
    for res in [16usize, 33, 64] {
        let (m, rep) = marching_cubes(&s, &opts(res)).unwrap();
        let (b, chi, comp) = manifold(&rep);
        assert_eq!((b, chi, comp), (0, 2, 1), "res {res}");
        let h = rep.spacing[0];
        let err = m.positions.iter().map(|p| s.distance(*p).unwrap().abs()).fold(0.0, f64::max);
        assert!(err <= 0.5 * h, "res {res}: vertex error {err} > h/2 = {}", 0.5 * h);
        // outward orientation: positive volume close to 4/3 pi r^3
        let vol = m.signed_volume();
        let exact = 4.0 / 3.0 * std::f64::consts::PI * 0.6f64.powi(3);
        assert!(vol > 0.0 && (vol - exact).abs() < 0.05 * exact, "volume {vol} vs {exact}");
        assert!(!rep.bound_violated);
    }
}

#[test]
fn torus_has_euler_0_and_refinement_reduces_error() {
    let t = Torus(0.6, 0.2);
    let (m0, rep0) = marching_cubes(&t, &opts(64)).unwrap();
    assert_eq!(manifold(&rep0), (0, 0, 1));
    let (m1, rep1) = marching_cubes(&t, &McOptions { refine_steps: 3, normals: true, ..opts(64) }).unwrap();
    assert_eq!(manifold(&rep1), (0, 0, 1));
    assert_eq!(m0.faces, m1.faces);
    let e = |m: &McMesh| m.positions.iter().map(|p| t.distance(*p).unwrap().abs()).fold(0.0, f64::max);
    assert!(e(&m1) < 1e-4 && e(&m1) < 0.05 * e(&m0), "{} vs {}", e(&m1), e(&m0));
    assert_eq!(rep1.zero_normals, 0);
    // normals point outward: dot with (p - nearest core point) positive; check against the gradient
    let nrm = m1.normals.as_ref().unwrap();
    for (p, n) in m1.positions.iter().zip(nrm).step_by(97) {
        let (_, g) = t.distance_gradient(*p).unwrap();
        assert!(g[0] * n[0] + g[1] * n[1] + g[2] * n[2] > 0.999999);
    }
}

#[test]
fn csg_of_boxes_union_and_difference() {
    let a = Csg::leaf(Cuboid([-0.15, 0.0, 0.0], [0.4, 0.3, 0.3]));
    let b = Csg::leaf(Cuboid([0.15, 0.0, 0.0], [0.4, 0.3, 0.3]));
    let u = a.clone().union(b.clone());
    let (m, rep) = marching_cubes(&u, &McOptions { bbox: Some([[-0.7; 3], [0.7; 3]]), ..opts(60) }).unwrap();
    assert_eq!(manifold(&rep), (0, 2, 1));
    let exact = 1.1 * 0.36; // [-0.55, 0.55] x [-0.3, 0.3]^2
    assert!((m.signed_volume() - exact).abs() < 0.03 * exact, "{} vs {exact}", m.signed_volume());
    // difference: a box with a bite, still a ball topologically; subtracting a thin through-box gives genus 1
    let hole = Csg::leaf(Cuboid([0.0, 0.0, 0.0], [0.1, 0.1, 1.0]));
    let d = u.subtract(hole);
    let (_, rep) = marching_cubes(&d, &McOptions { bbox: Some([[-0.7; 3], [0.7; 3]]), ..opts(70) }).unwrap();
    assert_eq!(manifold(&rep), (0, 0, 1));
}

#[test]
fn narrow_band_equals_dense_and_is_sparse() {
    let t = Torus(0.6, 0.2);
    let (md, rd) = marching_cubes(&t, &McOptions { band: Band::Dense, ..opts(64) }).unwrap();
    let (mb, rb) = marching_cubes(&t, &opts(64)).unwrap();
    assert_eq!(md.positions, mb.positions);
    let sorted = |m: &McMesh| {
        let mut f: Vec<[u32; 3]> = m.faces.clone();
        f.sort();
        f
    };
    assert_eq!(sorted(&md), sorted(&mb));
    assert!(rb.grid_values * 4 < rd.grid_values * 3, "{} vs {}", rb.grid_values, rd.grid_values);
    assert_eq!(rb.lipschitz_bound, Some(2.0));
    assert!(rb.max_edge_slope <= 1.0 + 1e-9);
}

#[test]
fn deterministic_across_thread_counts() {
    let u = Csg::leaf(Sphere(0.5)).union(Csg::leaf(Torus(0.5, 0.15))).offset(0.02);
    let base = marching_cubes(&u, &McOptions { threads: Some(1), refine_steps: 2, normals: true, ..opts(48) }).unwrap().0;
    for th in [3, 8, 16] {
        let m = marching_cubes(&u, &McOptions { threads: Some(th), refine_steps: 2, normals: true, ..opts(48) }).unwrap().0;
        assert_eq!(m.faces, base.faces);
        for (a, b) in m.positions.iter().zip(&base.positions) {
            assert!(a.iter().zip(b).all(|(x, y)| x.to_bits() == y.to_bits()));
        }
        for (a, b) in m.normals.as_ref().unwrap().iter().zip(base.normals.as_ref().unwrap()) {
            assert!(a.iter().zip(b).all(|(x, y)| x.to_bits() == y.to_bits()));
        }
    }
}

/// splitmix64 for the random lattices.
pub fn rnd(state: &mut u64) -> f64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^= z >> 31;
    (z >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0
}

#[test]
fn random_lattices_are_manifold_with_boundary_only_on_the_box() {
    let mut st = 12345u64;
    let mut ambiguous = 0usize;
    let mut cases = 0usize;
    for trial in 0..30000 {
        let n = [2 + trial % 3, 2 + (trial / 3) % 3, 2 + (trial / 9) % 3];
        let count = (n[0] + 1) * (n[1] + 1) * (n[2] + 1);
        let vals: Vec<f64> = (0..count).map(|_| rnd(&mut st)).collect();
        let f = Lattice { n, vals };
        let o = McOptions { resolution: Resolution::Cells(n), band: Band::Dense, threads: Some(2), ..Default::default() };
        let (m, rep) = marching_cubes(&f, &o).unwrap();
        ambiguous += rep.ambiguous_cell_faces;
        cases += 1;
        if m.faces.is_empty() {
            continue;
        }
        match &rep.mesh {
            MeshCheck::Manifold { .. } => {}
            other => panic!("trial {trial}: {other:?}"),
        }
        // boundary edges must lie in a face of the box: both endpoints share a box face coordinate
        let tm = m.to_trimesh().unwrap();
        for e in 0..tm.n_edges() {
            if tm.is_boundary_edge(e) {
                let [a, b] = tm.edge_vertices(e);
                let (pa, pb) = (m.positions[a], m.positions[b]);
                let on_face = (0..3).any(|ax| (pa[ax] == 0.0 && pb[ax] == 0.0) || (pa[ax] == 1.0 && pb[ax] == 1.0) || (pa[ax].abs() < 1e-12 && pb[ax].abs() < 1e-12) || ((pa[ax] - 1.0).abs() < 1e-12 && (pb[ax] - 1.0).abs() < 1e-12));
                assert!(on_face, "trial {trial}: boundary edge {pa:?} {pb:?} is inside the box");
            }
        }
    }
    assert!(ambiguous > 1000, "only {ambiguous} ambiguous faces in {cases} lattices");
}

#[test]
fn ambiguous_face_decider_connects_the_stronger_diagonal() {
    // One face z=0 of a single cell with the checkerboard (-,+;+,-): inside corners (0,0) and (1,1).
    // vals index: i + 2 j + 4 k. Top face (k=1) all outside.
    let make = |a: f64, b: f64| {
        let mut v = vec![1.0; 8];
        v[0] = -a; // (0,0,0)
        v[3] = -a; // (1,1,0)
        v[1] = b; // (1,0,0)
        v[2] = b; // (0,1,0)
        Lattice { n: [1, 1, 1], vals: v }
    };
    let o = McOptions { resolution: Resolution::Cells([1, 1, 1]), band: Band::Dense, threads: Some(1), ..Default::default() };
    // strongly negative inside values: the inside diagonal is connected -> one loop (a tunnel-like band): 1 component
    let (_, r1) = marching_cubes(&make(10.0, 0.1), &o).unwrap();
    // weakly negative: the inside corners are separated -> two caps
    let (_, r2) = marching_cubes(&make(0.1, 10.0), &o).unwrap();
    assert!(r1.ambiguous_cell_faces >= 1 && r2.ambiguous_cell_faces >= 1);
    assert_eq!(manifold(&r1).2, 1);
    assert_eq!(manifold(&r2).2, 2);
}

#[test]
fn csg_gradients_match_finite_differences() {
    let a = Csg::leaf(Sphere(0.5));
    let b = Csg::leaf(Torus(0.5, 0.2));
    let c = Csg::leaf(Cuboid([0.1, 0.0, 0.0], [0.3, 0.3, 0.3]));
    let trees = [
        a.clone().union(b.clone()),
        a.clone().intersect(c.clone()),
        b.clone().subtract(a.clone()).offset(0.05),
        Csg::SmoothUnion(Box::new(a.clone()), Box::new(b.clone()), 0.3),
        Csg::SmoothIntersect(Box::new(a.clone()), Box::new(c.clone()), 0.3),
        Csg::SmoothSubtract(Box::new(c.clone()), Box::new(a.clone()), 0.3),
        a.clone().union(b.clone()).subtract(c.clone()).offset(-0.02),
    ];
    let mut st = 99u64;
    let h = 1e-6;
    for (ti, t) in trees.iter().enumerate() {
        let mut checked = 0;
        for _ in 0..400 {
            let x = [rnd(&mut st) * 0.9, rnd(&mut st) * 0.9, rnd(&mut st) * 0.9];
            let (d, g) = t.distance_gradient(x).unwrap();
            assert_eq!(d, t.distance(x).unwrap());
            let mut fd = [0.0; 3];
            for i in 0..3 {
                let (mut p, mut m) = (x, x);
                p[i] += h;
                m[i] -= h;
                fd[i] = (t.distance(p).unwrap() - t.distance(m).unwrap()) / (2.0 * h);
            }
            // skip points within the FD stencil of a crease (value of the two arguments nearly tied) or the cuboid's medial surface
            let err = (0..3).map(|i| (g[i] - fd[i]).abs()).fold(0.0, f64::max);
            if err < 1e-5 {
                checked += 1;
            } else {
                // a mismatch is allowed only next to a crease: perturb and compare the two one-sided results
                let near = (0..3).any(|i| {
                    let mut p = x;
                    p[i] += 1e-4;
                    let (_, g2) = t.distance_gradient(p).unwrap();
                    (g2[i] - g[i]).abs() > 1e-3
                }) || (0..3).any(|i| {
                    let mut p = x;
                    p[i] -= 1e-4;
                    let (_, g2) = t.distance_gradient(p).unwrap();
                    (g2[i] - g[i]).abs() > 1e-3
                });
                assert!(near, "tree {ti} at {x:?}: g = {g:?}, fd = {fd:?}");
            }
        }
        assert!(checked > 300, "tree {ti}: only {checked} smooth points");
    }
}

#[test]
fn csg_ties_go_to_a_and_batches_agree() {
    let s = Csg::leaf(Sphere(0.5));
    let t = Csg::leaf(Sphere(0.5));
    // identical arguments: gradient of A (equal anyway); use different leaves with equal value instead
    let c = Csg::leaf(Cuboid([0.0; 3], [0.5, 10.0, 10.0]));
    let x = [0.5, 0.0, 0.0]; // sphere: 0; cuboid: 0 -> tie
    let u = s.clone().union(c.clone());
    let (d, g) = u.distance_gradient(x).unwrap();
    assert_eq!(d, 0.0);
    assert_eq!(g, [1.0, 0.0, 0.0]);
    let _ = t;
    let tree = s.clone().union(c).subtract(Csg::leaf(Sphere(0.3))).offset(0.01);
    let xs: Vec<[f64; 3]> = (0..1000).map(|i| [(i as f64 * 0.37).sin(), (i as f64 * 0.11).cos(), (i as f64 * 0.7).sin() * 0.5]).collect();
    let mut d = vec![0.0; xs.len()];
    let mut g = vec![[0.0; 3]; xs.len()];
    tree.distance_gradients(&xs, &mut d, &mut g).unwrap();
    let mut d2 = vec![0.0; xs.len()];
    tree.distances(&xs, &mut d2).unwrap();
    for i in 0..xs.len() {
        let (a, b) = tree.distance_gradient(xs[i]).unwrap();
        assert_eq!(a, d[i]);
        assert_eq!(a, d2[i]);
        assert_eq!(b, g[i]);
    }
    assert!(Csg::SmoothUnion(Box::new(s.clone()), Box::new(s), 0.0).distance([0.0; 3]).is_err());
}

#[test]
fn bad_input_is_an_error_not_a_panic() {
    let s = Sphere(0.5);
    for o in [
        McOptions { resolution: Resolution::Longest(0), ..opts(8) },
        McOptions { resolution: Resolution::Cells([4, 0, 4]), ..opts(8) },
        McOptions { bbox: Some([[0.0; 3], [0.0, 1.0, 1.0]]), ..opts(8) },
        McOptions { bbox: Some([[f64::NAN; 3], [1.0; 3]]), ..opts(8) },
        McOptions { isovalue: f64::INFINITY, ..opts(8) },
        McOptions { band: Band::Lipschitz { bound: -1.0, block: 8 }, ..opts(8) },
        McOptions { resolution: Resolution::Longest(100_000), ..opts(8) },
    ] {
        assert!(marching_cubes(&s, &o).is_err());
    }
    // an empty level set is an empty mesh, not an error
    let (m, _) = marching_cubes(&s, &McOptions { isovalue: 5.0, ..opts(8) }).unwrap();
    assert!(m.faces.is_empty());
}

#[test]
fn offset_isovalue_gives_the_offset_sphere() {
    let s = Sphere(0.4);
    let (m, rep) = marching_cubes(&s, &McOptions { isovalue: 0.15, refine_steps: 2, ..opts(48) }).unwrap();
    assert_eq!(manifold(&rep), (0, 2, 1));
    assert!(m.positions.iter().all(|p| (s.distance(*p).unwrap() - 0.15).abs() < 1e-6));
}

#[test]
fn growth_repairs_a_too_small_bound() {
    // With L = 0.7 the first pass skips blocks that contain the surface; the growth rule must
    // still deliver the whole connected surface (same vertices and faces as the dense grid).
    let t = Torus(0.6, 0.2);
    let (md, _) = marching_cubes(&t, &McOptions { band: Band::Dense, ..opts(64) }).unwrap();
    let (mb, rb) = marching_cubes(&t, &McOptions { band: Band::Lipschitz { bound: 0.7, block: 8 }, ..opts(64) }).unwrap();
    assert_eq!(manifold(&rb), (0, 0, 1));
    assert_eq!(md.positions, mb.positions);
    assert_eq!(md.faces.len(), mb.faces.len());
    assert!(rb.band_rounds >= 2 && rb.blocks_grown > 0, "{} rounds, {} grown", rb.band_rounds, rb.blocks_grown);
}
