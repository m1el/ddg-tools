//! Shared helpers of the integration tests: golden files of the authors' program
//! (`runs_tmp/pat-ref/`, gitignored, written by the reference-build agent; see its README.md).

#![allow(dead_code)]

use ddg_points_as_tori::npy::{NpyArray, read_npy};
use ddg_points_as_tori::{NormalPolicy, PointCloud, Torus};
use std::path::PathBuf;

/// The golden directory, or `None` (with a message) when it is absent.
pub fn golden_root() -> Option<PathBuf> {
    let root = std::env::var_os("PAT_REF_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../runs_tmp/pat-ref"));
    if root.join("cpu_run1").is_dir() && root.join("inputs").is_dir() {
        Some(root)
    } else {
        eprintln!(
            "SKIP: golden files of the reference not found at {} (set PAT_REF_DIR; they are written by the reference-build agent and are gitignored)",
            root.display()
        );
        None
    }
}

/// The 12 golden inputs.
pub const INPUTS: [&str; 12] = [
    "synth_sphere",
    "synth_torus",
    "ABC0_1",
    "ABC0_6",
    "ABC17_1",
    "ABC17_3",
    "block",
    "knot",
    "tooth",
    "femur",
    "snail",
    "holey-bunny",
];

pub fn load(path: PathBuf) -> NpyArray {
    read_npy(&path).unwrap_or_else(|e| panic!("{e}"))
}

pub fn to3(v: &[f64]) -> Vec<[f64; 3]> {
    v.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect()
}

/// One golden input with its CPU outputs.
pub struct Golden {
    pub name: String,
    pub dir: PathBuf,
    pub cloud: PointCloud,
    pub root: PathBuf,
}

impl Golden {
    pub fn open(root: &std::path::Path, name: &str) -> Golden {
        let inp = root.join("inputs").join(name);
        let p = to3(&load(inp.join("points.npy")).data);
        let n = to3(&load(inp.join("normals.npy")).data);
        let (cloud, _) = PointCloud::new(p, n, NormalPolicy::RequireUnit { tol: 1e-6 }).expect("golden cloud validates");
        Golden { name: name.to_string(), dir: root.join("cpu_run1").join(name), cloud, root: root.to_path_buf() }
    }

    pub fn arr(&self, f: &str) -> NpyArray {
        load(self.dir.join(format!("{f}.npy")))
    }

    /// The reference's tori with the axis from its stored angles (as its evaluator does).
    pub fn tori(&self) -> Vec<Torus> {
        let c = self.arr("tori_centers").data;
        let a = self.arr("tori_axes_angles").data;
        let big = self.arr("tori_major_radius").data;
        let small = self.arr("tori_minor_radius").data;
        (0..big.len())
            .map(|i| Torus {
                center: [c[3 * i], c[3 * i + 1], c[3 * i + 2]],
                axis: Torus::axis_from_angles(a[2 * i], a[2 * i + 1]),
                major: big[i],
                minor: small[i],
            })
            .collect()
    }

    /// The golden coefficients (f32).
    pub fn coeffs(&self) -> Vec<[f32; 6]> {
        let c = self.arr("net_coeffs").f32_data.expect("net_coeffs is f32");
        c.chunks_exact(6).map(|r| [r[0], r[1], r[2], r[3], r[4], r[5]]).collect()
    }

    pub fn queries(&self) -> Vec<[f64; 3]> {
        to3(&load(self.root.join("inputs/queries_random.npy")).data)
    }

    /// The 64³ grid, C order of `(ix, iy, iz)`.
    pub fn grid(&self) -> Vec<[f64; 3]> {
        let ax = load(self.root.join("inputs/grid_axis.npy")).data;
        let mut g = Vec::with_capacity(ax.len().pow(3));
        for &x in &ax {
            for &y in &ax {
                for &z in &ax {
                    g.push([x, y, z]);
                }
            }
        }
        g
    }
}

/// Max, median and 99.9th percentile of |a − b|.
pub fn diff_stats(a: &[f64], b: &[f64]) -> (f64, f64, f64) {
    assert_eq!(a.len(), b.len());
    let mut d: Vec<f64> = a.iter().zip(b).map(|(x, y)| (x - y).abs()).collect();
    d.sort_by(|x, y| x.partial_cmp(y).expect("finite differences"));
    let q = |f: f64| d[((d.len() - 1) as f64 * f).round() as usize];
    (d[d.len() - 1], q(0.5), q(0.999))
}

/// A deterministic LCG in [0, 1).
pub fn lcg(seed: &mut u64) -> f64 {
    *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
    (*seed >> 11) as f64 / (1u64 << 53) as f64
}

/// `n` points on the ellipsoid `(x/a)² + (y/b)² + (z/c)² = 1` with exact outward unit normals,
/// pseudo-random (seeded).
pub fn ellipsoid_cloud(n: usize, seed: u64, axes: [f64; 3]) -> PointCloud {
    let mut s = seed;
    let mut p = Vec::with_capacity(n);
    let mut nr = Vec::with_capacity(n);
    while p.len() < n {
        let v = [2.0 * lcg(&mut s) - 1.0, 2.0 * lcg(&mut s) - 1.0, 2.0 * lcg(&mut s) - 1.0];
        let l = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        if !(0.1..=1.0).contains(&l) {
            continue;
        }
        let u = [v[0] / l, v[1] / l, v[2] / l];
        let x = [axes[0] * u[0], axes[1] * u[1], axes[2] * u[2]];
        let g = [x[0] / (axes[0] * axes[0]), x[1] / (axes[1] * axes[1]), x[2] / (axes[2] * axes[2])];
        let gl = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
        p.push(x);
        nr.push([g[0] / gl, g[1] / gl, g[2] / gl]);
    }
    PointCloud::new(p, nr, NormalPolicy::default()).expect("valid cloud").0
}

/// Bit patterns of a slice of f64.
pub fn bits(v: &[f64]) -> Vec<u64> {
    v.iter().map(|x| x.to_bits()).collect()
}

/// `FitOptions` pinned to the authors' sign rule (`SignRule::Reference`), for parity, golden and
/// oracle tests. The crate default is `SignRule::Geometric`.
pub fn ref_fit() -> ddg_points_as_tori::FitOptions {
    ddg_points_as_tori::FitOptions { sign_rule: ddg_points_as_tori::SignRule::Reference, ..Default::default() }
}

/// `BuildOptions` pinned to the authors' sign rule (`SignRule::Reference`).
pub fn ref_build() -> ddg_points_as_tori::BuildOptions {
    ddg_points_as_tori::BuildOptions { fit: ref_fit(), ..Default::default() }
}
