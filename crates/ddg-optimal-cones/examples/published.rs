//! Investigation of the authors' published cones
//! (`data/optimal-cones/optimal-cones-data/data`: `<name>.obj`, `<name>-cones.txt`).
//!
//! ```text
//! cargo run --release -p ddg-optimal-cones --example published -- COMMAND [--data DIR] [mesh names...]
//!
//!   facts    topology, curvature sums, the cone files, and what the published UV maps say
//!            about the scale factor on the boundary (no solve)
//!   refit    the angles that minimize the distortion at the file's cone POSITIONS, with no
//!            cone penalty, under three boundary models; compared with the file's angles
//!   adjoint  the adjoint state of the file's cones in the canonical problem, at the cones
//!            and elsewhere
//!   iterate  penalized problem, merge, make the cones free sites, repeat  [--lambdas a,b,c]
//!   sweep    solves the penalized problem for a range of λ, merges clusters, re-optimizes the
//!            angles, and compares positions and angles with the file
//!            [--lambdas a,b,c] [--boundary fixed|free]
//! ```
//!
//! What these establish is summarized in the crate README ("The published dataset").
//!
//! Boundary models of `refit` (they differ only on meshes with boundary):
//!
//! * `fixed`: `u = 0` on the boundary (the canonical problem of this crate; paper Fig. 10 b).
//! * `constant`: `u = c` on the boundary with the constant `c` optimized.
//! * `free`: every boundary value of `u` optimized (paper Sec. 8.3, Fig. 10 c).
//!
//! All three are instances of one linear-algebra problem. Let `E` be the vertices that carry a
//! value of `u` (interior vertices; plus one vertex standing for the whole boundary; or all
//! vertices), `L` and `M` the Laplacian and mass on `E`, and `N ⊂ E` the vertices whose
//! curvature after flattening is prescribed (`ν_i`, zero at an interior vertex without a
//! cone). Minimizing `½ uᵀMu` subject to `(L u)_i = ν_i − κ_i` for `i ∈ N` gives
//!
//! ```text
//! (L_NE M⁻¹ L_EN) φ_N = κ_N − ν_N,      u = −M⁻¹ L_EN φ_N,      ν_i = κ_i + (L u)_i  (i ∉ N)
//! ```
//!
//! where `φ` (zero outside `N`) is the adjoint state: `L φ = −M u`. A vertex outside `N` is a
//! free cone site: its angle is whatever minimizes the distortion, and `φ = 0` there.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use ddg_core::io::{Cone, read_cones, read_obj};
use ddg_core::operators::{
    angle_defects, boundary_geodesic_curvatures, cotan_laplacian, lumped_mass_normalized,
};
use ddg_core::sparse::{Method, Parallelism, RefineOptions, Symbolic};
use ddg_core::{Geometry, SparseMatrix, TriMesh};
use ddg_optimal_cones::{ConeProblem, SolveOptions, Status, evaluate_cones, solve};

struct Dataset {
    mesh: TriMesh,
    geometry: Geometry,
    cones: Vec<Cone>,
    /// Per-corner UV coordinates of the published flattening.
    uvs: Vec<[f64; 2]>,
    positions: Vec<[f64; 3]>,
}

fn load(dir: &Path, name: &str) -> Dataset {
    let obj = read_obj(dir.join(format!("{name}.obj"))).expect("readable OBJ");
    let faces = obj.triangles().expect("triangle mesh");
    let mesh = TriMesh::new(obj.positions.len(), faces).expect("valid mesh");
    let geometry = Geometry::from_positions(&mesh, &obj.positions).expect("valid geometry");
    let cones = read_cones(dir.join(format!("{name}-cones.txt")), mesh.n_vertices())
        .expect("readable cone file");
    Dataset {
        uvs: obj.corner_uvs().expect("every corner has a UV"),
        positions: obj.positions.clone(),
        mesh,
        geometry,
        cones,
    }
}

fn main() {
    let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../data/optimal-cones/optimal-cones-data/data");
    let mut names: Vec<String> = Vec::new();
    let mut lambdas: Vec<f64> = Vec::new();
    let mut boundary = Boundary::Fixed;
    let mut args = std::env::args().skip(1);
    let command = args.next().unwrap_or_else(|| {
        eprintln!(
            "usage: published facts|refit|adjoint|iterate|sweep [--data DIR] [--lambdas a,b,c] \
             [--boundary fixed|free] [mesh names...]"
        );
        std::process::exit(2);
    });
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--data" => dir = PathBuf::from(args.next().expect("--data needs a value")),
            "--boundary" => {
                boundary = match args.next().expect("--boundary needs a value").as_str() {
                    "fixed" => Boundary::Fixed,
                    "free" => Boundary::Free,
                    other => panic!("--boundary: expected fixed or free, got {other}"),
                }
            }
            "--lambdas" => {
                lambdas = args
                    .next()
                    .expect("--lambdas needs a value")
                    .split(',')
                    .map(|s| s.parse().expect("--lambdas: comma separated numbers"))
                    .collect()
            }
            name => names.push(name.to_string()),
        }
    }
    if names.is_empty() {
        names = std::fs::read_dir(&dir)
            .expect("readable dataset directory")
            .filter_map(|e| e.ok())
            .filter_map(|e| {
                e.file_name()
                    .to_str()
                    .and_then(|n| n.strip_suffix(".obj"))
                    .map(str::to_string)
            })
            .collect();
        names.sort();
    }
    match command.as_str() {
        "facts" => facts(&dir, &names),
        "refit" => refit(&dir, &names),
        "sweep" => sweep(&dir, &names, &lambdas, boundary),
        "adjoint" => adjoint(&dir, &names),
        "iterate" => iterate(&dir, &names, &lambdas),
        other => {
            eprintln!("unknown command {other}");
            std::process::exit(2);
        }
    }
}

// ----- facts -----

fn facts(dir: &Path, names: &[String]) {
    println!(
        "{:<20} {:>6} {:>6} {:>3} {:>5} {:>3} {:>8} {:>9} {:>9} {:>5} {:>4} {:>4} {:>9} {:>9} | {:>8} {:>8} {:>8} {:>8} {:>8}",
        "mesh", "V", "F", "bl", "bE", "chi", "area", "sum Omega", "sum k", "cones", "pos", "neg",
        "sum mu", "sum |mu|", "uB min", "uB max", "u mean", "E", "E var"
    );
    for name in names {
        let data = load(dir, name);
        let (mesh, geometry) = (&data.mesh, &data.geometry);
        let omega = angle_defects(mesh, geometry).unwrap();
        let kappa = boundary_geodesic_curvatures(mesh, geometry).unwrap();
        let area = geometry.total_area();
        let positive = data.cones.iter().filter(|c| c.angle > 0.0).count();
        let sum: f64 = data.cones.iter().map(|c| c.angle).sum();
        let abs: f64 = data.cones.iter().map(|c| c.angle.abs()).sum();

        // The published UV map: per-face log scale factor u_f = ½ log(area_uv / area_3d), and
        // on every boundary edge the log of the length ratio. The scale of a UV map is
        // arbitrary; on a mesh with boundary it is fixed here by the mean boundary value.
        let nf = mesh.n_faces();
        let areas = geometry.face_areas();
        let mut u_face = vec![0.0; nf];
        let mut flipped = 0;
        for f in 0..nf {
            let [a, b, c] = [data.uvs[3 * f], data.uvs[3 * f + 1], data.uvs[3 * f + 2]];
            let uv_area = 0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]));
            if uv_area <= 0.0 {
                flipped += 1;
            }
            u_face[f] = 0.5 * (uv_area.abs() / areas[f]).ln();
        }
        let mut boundary_u = Vec::new();
        for h in 0..3 * nf {
            if mesh.twin(h).is_none() {
                let (p, q) = (data.uvs[h], data.uvs[3 * (h / 3) + (h % 3 + 1) % 3]);
                let uv_length = ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2)).sqrt();
                boundary_u.push((uv_length / geometry.edge_lengths()[mesh.edge(h)]).ln());
            }
        }
        let finite: Vec<usize> = (0..nf).filter(|&f| u_face[f].is_finite()).collect();
        let weighted = |g: &dyn Fn(usize) -> f64| -> f64 {
            finite.iter().map(|&f| areas[f] * g(f)).sum::<f64>()
                / finite.iter().map(|&f| areas[f]).sum::<f64>()
        };
        let mean_all = weighted(&|f| u_face[f]);
        let (offset, ub_min, ub_max) = if boundary_u.is_empty() {
            (mean_all, f64::NAN, f64::NAN)
        } else {
            let mean = boundary_u.iter().sum::<f64>() / boundary_u.len() as f64;
            (
                mean,
                boundary_u.iter().cloned().fold(f64::INFINITY, f64::min) - mean,
                boundary_u.iter().cloned().fold(f64::NEG_INFINITY, f64::max) - mean,
            )
        };
        let mean = weighted(&|f| u_face[f] - offset);
        let energy = 0.5 * weighted(&|f| (u_face[f] - offset).powi(2));
        println!(
            "{:<20} {:>6} {:>6} {:>3} {:>5} {:>3} {:>8.4} {:>9.4} {:>9.4} {:>5} {:>4} {:>4} {:>9.4} {:>9.4} | {:>8.4} {:>8.4} {:>8.4} {:>8.4} {:>8.4}{}",
            name,
            mesh.n_vertices(),
            nf,
            mesh.boundary_loops().len(),
            mesh.n_boundary_edges(),
            mesh.euler_characteristic(),
            area,
            omega.iter().sum::<f64>(),
            kappa.iter().sum::<f64>(),
            data.cones.len(),
            positive,
            data.cones.len() - positive,
            sum,
            abs,
            ub_min,
            ub_max,
            mean,
            energy,
            energy - 0.5 * mean * mean,
            if flipped > 0 || finite.len() < nf {
                format!("   ({flipped} flipped, {} degenerate UV triangles)", nf - finite.len())
            } else {
                String::new()
            }
        );
    }
    println!(
        "\nuB min/max: log(UV length / length) of the boundary edges minus its mean (constant = \
         the map has one scale factor on the whole boundary).\nu mean, E = ½ mean(u²), E var = E − \
         ½ mean(u)²: area-weighted over the faces, of u_f = ½ log(UV area / area) with the scale \
         fixed by u = 0 on the boundary\n(closed meshes: by mean u = 0)."
    );
}

// ----- iterate -----

/// A scheme that would produce cones with zero adjoint: solve the penalized problem, merge
/// the clusters, make the cones found so far free sites (Sec. 8.1 of the paper: no cost
/// there), solve again with the same `λ` elsewhere, and repeat until no further cone appears.
/// Canonical problem (`u = 0` on the boundary), through the public interface.
fn iterate(dir: &Path, names: &[String], lambdas: &[f64]) {
    for name in names {
        let data = load(dir, name);
        let nv = data.mesh.n_vertices();
        let area = data.geometry.total_area();
        println!("{name}: file has {} cones; area {area:.4}", data.cones.len());
        for &lambda in lambdas {
            let mut sites: Vec<u32> = Vec::new();
            let mut history = Vec::new();
            let mut last: Vec<Cone> = Vec::new();
            for _round in 0..12 {
                let mut costs = vec![1.0; nv];
                for &v in &sites {
                    costs[v as usize] = 0.0;
                }
                let problem = ConeProblem {
                    cone_costs: Some(costs),
                    ..ConeProblem::new(&data.mesh, &data.geometry, lambda)
                };
                let solution = solve(&problem, &SolveOptions::default()).expect("valid problem");
                assert_eq!(solution.report.status, Status::Converged);
                let merged = merge(&data, &solution.measure, false);
                let new: Vec<u32> = merged
                    .iter()
                    .map(|c| c.vertex)
                    .filter(|v| !sites.contains(v))
                    .collect();
                history.push(format!("{}+{}", sites.len(), new.len()));
                last = merged;
                if new.is_empty() {
                    break;
                }
                sites.extend(new);
            }
            // The cones of the last round: the free sites that kept an angle.
            let distances = ring_distances(&data.mesh, &data.cones, &last, 1);
            let ours: BTreeMap<u32, f64> = last.iter().map(|c| (c.vertex, c.angle)).collect();
            let worst = data
                .cones
                .iter()
                .filter_map(|c| ours.get(&c.vertex).map(|a| (a - c.angle).abs()))
                .fold(f64::NAN, f64::max);
            println!(
                "  lambda {:>10.4e}  lambda*A {:>10.4e}  rounds (free sites + new cones): {:<34} final {:>3} cones, sum |mu| {:>8.4}; file cones at the same vertex {:>3}, within one edge {:>3}, max |d| {:>9.2e}",
                lambda,
                lambda * area,
                history.join(" "),
                last.len(),
                last.iter().map(|c| c.angle.abs()).sum::<f64>(),
                distances.iter().filter(|&&d| d == 0).count(),
                distances.iter().filter(|&&d| d <= 1).count(),
                worst
            );
        }
    }
}

// ----- adjoint -----

/// The adjoint state of the file's cones in the canonical problem (`evaluate_cones`), at each
/// cone and its range elsewhere. An optimum of the penalized problem has `φ = λ·sign(µ)` at
/// its cones and `|φ| ≤ λ` elsewhere (on a closed mesh: after adding one constant).
fn adjoint(dir: &Path, names: &[String]) {
    for name in names {
        let data = load(dir, name);
        if data.cones.len() > 60 {
            continue;
        }
        let problem = ConeProblem::new(&data.mesh, &data.geometry, 0.0);
        let evaluation = evaluate_cones(&problem, &data.cones, Parallelism::Sequential)
            .expect("the file's cones are valid");
        let mut is_cone = vec![false; data.mesh.n_vertices()];
        for c in &data.cones {
            is_cone[c.vertex as usize] = true;
        }
        let elsewhere: Vec<f64> = (0..data.mesh.n_vertices())
            .filter(|&v| !is_cone[v] && !data.mesh.is_boundary_vertex(v))
            .map(|v| evaluation.adjoint[v])
            .collect();
        println!(
            "{name}: {} cones, Gauss-Bonnet defect {:.2e}; adjoint elsewhere in [{:.6}, {:.6}]",
            data.cones.len(),
            evaluation.gauss_bonnet_defect,
            elsewhere.iter().cloned().fold(f64::INFINITY, f64::min),
            elsewhere.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
        );
        for c in &data.cones {
            println!(
                "    vertex {:>6}  angle {:>9.5}  adjoint {:>10.6}",
                c.vertex + 1,
                c.angle,
                evaluation.adjoint[c.vertex as usize]
            );
        }
    }
}

// ----- refit -----

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Boundary {
    Fixed,
    Constant,
    Free,
}

/// The Laplacian, mass and curvature on the vertices that carry a value of `u`.
struct Model {
    /// For each mesh vertex its index in `E`, or `usize::MAX` if it carries no value.
    index: Vec<usize>,
    n: usize,
    laplacian: SparseMatrix,
    /// Area-normalized masses, summed over the vertices that share an index.
    mass: Vec<f64>,
    /// Angle defect at interior vertices, geodesic curvature at boundary vertices, summed
    /// likewise.
    kappa: Vec<f64>,
    /// Whether the index belongs to an interior vertex.
    interior: Vec<bool>,
}

impl Model {
    fn new(data: &Dataset, boundary: Boundary) -> Model {
        let (mesh, geometry) = (&data.mesh, &data.geometry);
        let nv = mesh.n_vertices();
        let interior_vertices = mesh.interior_vertices();
        let mut index = vec![usize::MAX; nv];
        let n = match boundary {
            Boundary::Fixed => {
                for (i, &v) in interior_vertices.iter().enumerate() {
                    index[v as usize] = i;
                }
                interior_vertices.len()
            }
            Boundary::Constant => {
                for (i, &v) in interior_vertices.iter().enumerate() {
                    index[v as usize] = i;
                }
                let has_boundary = interior_vertices.len() < nv;
                for slot in index.iter_mut().filter(|slot| **slot == usize::MAX) {
                    *slot = interior_vertices.len();
                }
                interior_vertices.len() + usize::from(has_boundary)
            }
            Boundary::Free => {
                for (v, slot) in index.iter_mut().enumerate() {
                    *slot = v;
                }
                nv
            }
        };
        let full = cotan_laplacian(mesh, geometry).unwrap();
        let triplets: Vec<(u32, u32, f64)> = full
            .triplets()
            .filter(|&(r, c, _)| index[r as usize] != usize::MAX && index[c as usize] != usize::MAX)
            .map(|(r, c, value)| (index[r as usize] as u32, index[c as usize] as u32, value))
            .collect();
        let laplacian = SparseMatrix::from_triplets(n, n, &triplets).unwrap();
        let mass_v = lumped_mass_normalized(mesh, geometry).unwrap();
        let omega = angle_defects(mesh, geometry).unwrap();
        let k = boundary_geodesic_curvatures(mesh, geometry).unwrap();
        let mut mass = vec![0.0; n];
        let mut kappa = vec![0.0; n];
        let mut interior = vec![false; n];
        for v in 0..nv {
            if index[v] == usize::MAX {
                continue;
            }
            mass[index[v]] += mass_v[v];
            kappa[index[v]] += omega[v] + k[v];
            interior[index[v]] = !mesh.is_boundary_vertex(v);
        }
        Model {
            index,
            n,
            laplacian,
            mass,
            kappa,
            interior,
        }
    }

    /// Solves the optimality system of: minimize `½ uᵀMu + Σ_pinned p_i ν_i` subject to
    /// `(L u)_i = ν_i − κ_i` at the `prescribed` vertices (index in `E`, `ν_i`), where
    /// `ν_i = κ_i + (L u)_i` is free at every other vertex:
    ///
    /// ```text
    /// M u + L φ = 0,     (L u)_i = ν_i − κ_i  (i prescribed),     φ_i = p_i  (i pinned),
    /// φ_i = 0  at the vertices in neither list (free sites)
    /// ```
    ///
    /// With `homogeneous` the curvature `κ` is taken as zero: the response to the pinned
    /// values alone. Returns `(u, φ, L u)` on `E`.
    ///
    /// The system is solved as it stands by a sparse LU factorization with refinement.
    /// (Eliminating `u` gives the system with `L_NE M⁻¹ L_EN` of the module documentation,
    /// whose condition number is that of the Laplacian squared: too inaccurate on the camels.)
    fn solve(
        &self,
        prescribed: &[(usize, f64)],
        pinned: &[(usize, f64)],
        homogeneous: bool,
    ) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let n = self.n;
        let m = prescribed.len();
        let mut local = vec![u32::MAX; n];
        for (a, &(i, _)) in prescribed.iter().enumerate() {
            local[i] = a as u32;
        }
        let mut triplets: Vec<(u32, u32, f64)> = Vec::new();
        for (e, &mass) in self.mass.iter().enumerate() {
            triplets.push((e as u32, e as u32, mass));
        }
        for (r, c, value) in self.laplacian.triplets() {
            let a = local[c as usize];
            if a != u32::MAX {
                triplets.push((r, n as u32 + a, value));
                triplets.push((n as u32 + a, r, value));
            }
        }
        let saddle = SparseMatrix::from_triplets(n + m, n + m, &triplets).unwrap();
        let mut phi = vec![0.0; n];
        for &(i, value) in pinned {
            assert!(local[i] == u32::MAX, "a vertex is both prescribed and pinned");
            phi[i] = value;
        }
        let mut rhs = vec![0.0; n + m];
        self.laplacian.mul_vec(&phi, &mut rhs[..n]).unwrap();
        rhs[..n].iter_mut().for_each(|x| *x = -*x);
        let kappa = |i: usize| if homogeneous { 0.0 } else { self.kappa[i] };
        for (a, &(i, target)) in prescribed.iter().enumerate() {
            rhs[n + a] = target - kappa(i);
        }
        let mut x = vec![0.0; n + m];
        Symbolic::analyze(&saddle, Method::Lu)
            .unwrap()
            .factor(&saddle, Parallelism::Sequential)
            .unwrap()
            .solve_refined(
                &saddle,
                &rhs,
                &mut x,
                RefineOptions {
                    max_steps: 5,
                    tolerance: 1e-15,
                },
                Parallelism::Sequential,
            )
            .unwrap();
        let u = x[..n].to_vec();
        for (&(i, _), &value) in prescribed.iter().zip(&x[n..]) {
            phi[i] = value;
        }
        let mut l_u = vec![0.0; n];
        self.laplacian.mul_vec(&u, &mut l_u).unwrap();
        // Both equations must hold: the prescribed curvatures come out as prescribed, and
        // L φ = −M u.
        let state = prescribed
            .iter()
            .map(|&(i, target)| (kappa(i) + l_u[i] - target).abs())
            .fold(0.0, f64::max);
        let mut l_phi = vec![0.0; n];
        self.laplacian.mul_vec(&phi, &mut l_phi).unwrap();
        let adjoint = (0..n)
            .map(|e| (l_phi[e] + self.mass[e] * u[e]).abs())
            .fold(0.0, f64::max);
        // The angles compared with the files have six digits.
        assert!(
            state < 1e-8 && adjoint < 1e-8,
            "inaccurate solve: state residual {state:.2e}, adjoint residual {adjoint:.2e}"
        );
        (u, phi, l_u)
    }

    fn distortion(&self, u: &[f64]) -> f64 {
        0.5 * u.iter().zip(&self.mass).map(|(u, m)| m * u * u).sum::<f64>()
    }
}

struct Refit {
    /// Least-distortion angles at the file's positions, in the order of the file.
    angles: Vec<f64>,
    distortion: f64,
    /// Largest `|φ|` over the interior vertices without a cone, for the refitted angles.
    adjoint_elsewhere: f64,
    /// For the file's own angles: the adjoint at each cone, and the distortion. `None` on a
    /// closed mesh, where the file's rounded angles do not sum to `2πχ`.
    file: Option<(Vec<f64>, f64)>,
    /// The values of `λ` (area-normalized masses) for which the file's cone vertices, with the
    /// signs of the file's angles, are exactly the support of the optimum of the penalized
    /// problem `½ uᵀMu + λ Σ |µ_i|` under this boundary model; `lo > hi` if there is none.
    support_interval: (f64, f64),
    /// The mesh vertex that determines the lower end (it would become a cone below it) and
    /// its distance in edges to the nearest cone of the file (capped at 3).
    lower_vertex: Option<(u32, usize)>,
}

fn refit_one(data: &Dataset, boundary: Boundary) -> Refit {
    let model = Model::new(data, boundary);
    let cone_index: Vec<usize> = data
        .cones
        .iter()
        .map(|c| model.index[c.vertex as usize])
        .collect();
    let mut is_cone = vec![false; model.n];
    for &i in &cone_index {
        is_cone[i] = true;
    }
    let plain: Vec<(usize, f64)> = (0..model.n)
        .filter(|&i| model.interior[i] && !is_cone[i])
        .map(|i| (i, 0.0))
        .collect();
    // λ = 0 on the support: the cones are free sites.
    let (u, phi, l_u) = model.solve(&plain, &[], false);
    let angles: Vec<f64> = cone_index.iter().map(|&i| model.kappa[i] + l_u[i]).collect();
    let closed = data.mesh.boundary_loops().is_empty();
    let file = if closed {
        None
    } else {
        let mut with_cones = plain.clone();
        with_cones.extend(cone_index.iter().zip(&data.cones).map(|(&i, c)| (i, c.angle)));
        let (u_file, phi_file, _) = model.solve(&with_cones, &[], false);
        Some((
            cone_index.iter().map(|&i| phi_file[i]).collect(),
            model.distortion(&u_file),
        ))
    };

    // The penalized optimum on the same support is affine in λ: the solution above plus λ
    // times the response to φ = sign(µ) at the cones. The support is optimal for the λ with
    // sign(µ_s)·µ_s(λ) ≥ 0 at every cone and |φ_i(λ)| ≤ λ at every other interior vertex.
    let signs: Vec<(usize, f64)> = cone_index
        .iter()
        .zip(&data.cones)
        .map(|(&i, c)| (i, c.angle.signum()))
        .collect();
    let (_, phi_1, l_u_1) = model.solve(&plain, &signs, true);
    let (mut lo, mut hi) = (0.0f64, f64::INFINITY);
    let mut lower_at = None;
    // a·λ ≥ b.
    let mut require = |a: f64, b: f64, at: Option<usize>| {
        if a.abs() <= 1e-12 {
            if b > 1e-12 {
                // Cannot hold for any λ.
                lo = f64::INFINITY;
                lower_at = at;
            }
        } else if a > 0.0 {
            if b / a > lo {
                lo = b / a;
                lower_at = at;
            }
        } else {
            hi = hi.min(b / a);
        }
    };
    for (k, &(i, sign)) in signs.iter().enumerate() {
        // sign·(µ0 + λ µ1) ≥ 0.
        require(sign * l_u_1[i], -sign * angles[k], None);
    }
    for &(i, _) in &plain {
        require(1.0 - phi_1[i], phi[i], Some(i));
        require(1.0 + phi_1[i], -phi[i], Some(i));
    }
    let lower_vertex = lower_at.and_then(|i| {
        let vertex = (0..data.mesh.n_vertices()).find(|&v| model.index[v] == i)? as u32;
        let probe = [Cone { vertex, angle: 0.0 }];
        Some((vertex, ring_distances(&data.mesh, &probe, &data.cones, 3)[0]))
    });
    Refit {
        angles,
        distortion: model.distortion(&u),
        adjoint_elsewhere: plain.iter().fold(0.0f64, |m, &(i, _)| m.max(phi[i].abs())),
        file,
        support_interval: (lo, hi),
        lower_vertex,
    }
}

fn refit(dir: &Path, names: &[String]) {
    println!(
        "{:<20} {:>5} {:<9} {:>9} {:>9} {:>10} {:>10} {:>10} {:>10} | {:>10} {:>10} {:>10} {:>10}  lower end set by",
        "mesh",
        "cones",
        "boundary",
        "max |d|",
        "rms |d|",
        "E refit",
        "E file",
        "max|phi|S",
        "max|phi|N",
        "lam lo",
        "lam hi",
        "lam lo*A",
        "lam hi*A",
    );
    for name in names {
        let data = load(dir, name);
        if data.cones.len() > 500 {
            println!("{name:<20} {:>5} skipped (the cones cover the mesh)", data.cones.len());
            continue;
        }
        let closed = data.mesh.boundary_loops().is_empty();
        let models: &[Boundary] = if closed {
            &[Boundary::Free]
        } else {
            &[Boundary::Fixed, Boundary::Constant, Boundary::Free]
        };
        let area = data.geometry.total_area();
        let mut best: Option<(f64, Boundary, Refit)> = None;
        for &boundary in models {
            let r = refit_one(&data, boundary);
            let differences: Vec<f64> = r
                .angles
                .iter()
                .zip(&data.cones)
                .map(|(a, c)| a - c.angle)
                .collect();
            let worst = differences.iter().fold(0.0f64, |m, d| m.max(d.abs()));
            let rms = (differences.iter().map(|d| d * d).sum::<f64>()
                / differences.len() as f64)
                .sqrt();
            let (phi_s, e_file) = match &r.file {
                Some((phi, e)) => (
                    format!("{:>10.3e}", phi.iter().fold(0.0f64, |m, p| m.max(p.abs()))),
                    format!("{e:>10.6}"),
                ),
                None => (format!("{:>10}", "-"), format!("{:>10}", "-")),
            };
            let (lo, hi) = r.support_interval;
            let interval = if lo <= hi {
                format!(
                    "{:>10.3e} {:>10.3e} {:>10.5} {:>10.5}  {}",
                    lo,
                    hi,
                    lo * area,
                    hi * area,
                    match r.lower_vertex {
                        Some((vertex, distance)) => format!(
                            "vertex {} ({} from a cone)",
                            vertex + 1,
                            match distance {
                                1 => "1 edge".to_string(),
                                d if d <= 3 => format!("{d} edges"),
                                _ => "more than 3 edges".to_string(),
                            }
                        ),
                        None => "-".to_string(),
                    }
                )
            } else {
                // No λ: the two requirements exclude each other. Shown as "needs λ ≥ … but ≤ …".
                format!(
                    "none: needs lam*A >= {:.5} ({}) and <= {:.5}",
                    lo * area,
                    match r.lower_vertex {
                        Some((vertex, distance)) =>
                            format!("vertex {}, {} edge(s) from a cone", vertex + 1, distance),
                        None => "-".to_string(),
                    },
                    hi * area
                )
            };
            println!(
                "{:<20} {:>5} {:<9} {:>9.2e} {:>9.2e} {:>10.6} {} {} {:>10.3e} | {}",
                name,
                data.cones.len(),
                if closed {
                    "closed".to_string()
                } else {
                    format!("{boundary:?}").to_lowercase()
                },
                worst,
                rms,
                r.distortion,
                e_file,
                phi_s,
                r.adjoint_elsewhere,
                interval
            );
            if best.as_ref().is_none_or(|(w, _, _)| worst < *w) {
                best = Some((worst, boundary, r));
            }
        }
        let (_, boundary, r) = best.expect("at least one model");
        if data.cones.len() <= 14 {
            for (cone, angle) in data.cones.iter().zip(&r.angles) {
                println!(
                    "      vertex {:>6}   file {:>9.5}   refit ({:?}) {:>9.5}   difference {:>9.2e}",
                    cone.vertex + 1,
                    cone.angle,
                    boundary,
                    angle,
                    angle - cone.angle
                );
            }
        }
    }
    println!(
        "\nmax |d|, rms |d|: refitted minus file angle over the file's cones (radians).\n\
         E: distortion ½ Σ m u² with area-normalized masses, of the refitted and of the file's angles.\n\
         max|phi|S: largest |adjoint| at the file's cones for the file's own angles (zero = the \
         angles minimize the distortion).\nmax|phi|N: largest |adjoint| at the other interior \
         vertices for the refitted angles.\nlam lo, lam hi: the file's cone vertices (with the \
         file's signs) are exactly the support of the optimum of the penalized problem\n  for \
         lambda in this interval (area-normalized masses; times the area: raw masses)."
    );
}

// ----- sweep -----

/// For every file cone, the graph distance (in edges, capped) to the nearest of `ours`.
fn ring_distances(mesh: &TriMesh, file: &[Cone], ours: &[Cone], cap: usize) -> Vec<usize> {
    let mut distance = vec![usize::MAX; mesh.n_vertices()];
    let mut frontier: Vec<usize> = ours.iter().map(|c| c.vertex as usize).collect();
    for &v in &frontier {
        distance[v] = 0;
    }
    for d in 1..=cap {
        let mut next = Vec::new();
        for &v in &frontier {
            for w in mesh.vertex_neighbors(v) {
                if distance[w] == usize::MAX {
                    distance[w] = d;
                    next.push(w);
                }
            }
        }
        frontier = next;
    }
    file.iter().map(|c| distance[c.vertex as usize].min(cap + 1)).collect()
}

/// The mesh glued to a mirror copy of itself along its boundary: a closed mesh. Vertex `v` of
/// the original keeps its index; the copy of interior vertex `v` gets a new one.
///
/// The canonical problem on this double, with the seam vertices as free cone sites
/// (`wR = 0`) and `λ/2`, is the problem on the original mesh with every boundary value of `u`
/// optimized ("free" boundary) and `λ`: the solution is mirror symmetric, each seam row of the
/// Laplacian, seam mass and seam curvature is twice the one of the original boundary vertex,
/// and the total area doubles. This lets the crate's certified solver handle that boundary
/// model through its public interface.
///
/// `None` if some face has three boundary vertices: it would coincide with its own copy.
fn double(data: &Dataset) -> Option<(TriMesh, Geometry)> {
    let mesh = &data.mesh;
    let nv = mesh.n_vertices();
    if mesh
        .faces()
        .iter()
        .any(|f| f.iter().all(|&v| mesh.is_boundary_vertex(v as usize)))
    {
        return None;
    }
    let mut copy: Vec<u32> = (0..nv as u32).collect();
    let mut next = nv as u32;
    for (v, slot) in copy.iter_mut().enumerate() {
        if !mesh.is_boundary_vertex(v) {
            *slot = next;
            next += 1;
        }
    }
    let mut faces = mesh.faces().to_vec();
    faces.extend(
        mesh.faces()
            .iter()
            .map(|f| [copy[f[0] as usize], copy[f[2] as usize], copy[f[1] as usize]]),
    );
    let doubled = TriMesh::new(next as usize, faces).expect("the double is a valid mesh");
    // Intrinsic geometry: every edge of the double has the length of the edge it copies.
    let mut original = vec![0u32; next as usize];
    for v in 0..nv {
        original[v] = v as u32;
        original[copy[v] as usize] = v as u32;
    }
    let lengths: Vec<f64> = (0..doubled.n_edges())
        .map(|e| {
            let [a, b] = doubled.edge_vertices(e);
            let edge = mesh
                .find_edge(original[a] as usize, original[b] as usize)
                .expect("the copy of an edge");
            data.geometry.edge_lengths()[edge]
        })
        .collect();
    let geometry = Geometry::from_edge_lengths(&doubled, lengths).expect("valid lengths");
    Some((doubled, geometry))
}

/// Merges clusters of cones of a per-vertex measure: a cluster is a set of cone vertices that
/// is connected through mesh edges between its members (`by_sign`: between members of the
/// same sign only). Each cluster becomes one cone with the summed angle at the member that
/// minimizes `Σ_k |µ_k| dist(v, k)²`, distances along edges inside the cluster.
fn merge(data: &Dataset, mu: &[f64], by_sign: bool) -> Vec<Cone> {
    let mesh = &data.mesh;
    let nv = mesh.n_vertices();
    let is_cone = |v: usize| mu[v].abs() > 1e-12;
    let mut visited = vec![false; nv];
    let mut cones = Vec::new();
    for seed in 0..nv {
        if visited[seed] || !is_cone(seed) {
            continue;
        }
        visited[seed] = true;
        let mut members = vec![seed];
        let mut head = 0;
        while head < members.len() {
            let v = members[head];
            head += 1;
            for w in mesh.vertex_neighbors(v) {
                if !visited[w] && is_cone(w) && (!by_sign || (mu[w] > 0.0) == (mu[seed] > 0.0)) {
                    visited[w] = true;
                    members.push(w);
                }
            }
        }
        members.sort_unstable();
        let m = members.len();
        let local = |v: usize| members.binary_search(&v).ok();
        let mut variance = vec![0.0; m];
        for source in 0..m {
            // Dijkstra inside the cluster (clusters are small; a linear scan for the minimum).
            let mut dist = vec![f64::INFINITY; m];
            let mut done = vec![false; m];
            dist[source] = 0.0;
            for _ in 0..m {
                let Some(a) = (0..m)
                    .filter(|&a| !done[a] && dist[a].is_finite())
                    .min_by(|&a, &b| dist[a].total_cmp(&dist[b]))
                else {
                    break;
                };
                done[a] = true;
                for w in mesh.vertex_neighbors(members[a]) {
                    if let Some(b) = local(w) {
                        let e = mesh.find_edge(members[a], w).expect("an edge");
                        let candidate = dist[a] + data.geometry.edge_lengths()[e];
                        if candidate < dist[b] {
                            dist[b] = candidate;
                        }
                    }
                }
            }
            for v in 0..m {
                variance[v] += mu[members[source]].abs() * dist[v] * dist[v];
            }
        }
        let best = (0..m)
            .min_by(|&a, &b| variance[a].total_cmp(&variance[b]))
            .expect("a member");
        let angle: f64 = members.iter().map(|&v| mu[v]).sum();
        if angle.abs() > 1e-12 {
            cones.push(Cone {
                vertex: members[best] as u32,
                angle,
            });
        }
    }
    cones
}

/// The least-distortion angles at the given positions under a boundary model, and the
/// distortion.
fn refit_positions(data: &Dataset, boundary: Boundary, positions: &[Cone]) -> (Vec<Cone>, f64) {
    let model = Model::new(data, boundary);
    let mut is_cone = vec![false; model.n];
    for c in positions {
        is_cone[model.index[c.vertex as usize]] = true;
    }
    let plain: Vec<(usize, f64)> = (0..model.n)
        .filter(|&i| model.interior[i] && !is_cone[i])
        .map(|i| (i, 0.0))
        .collect();
    let (u, _, l_u) = model.solve(&plain, &[], false);
    let cones = positions
        .iter()
        .map(|c| {
            let i = model.index[c.vertex as usize];
            Cone {
                vertex: c.vertex,
                angle: model.kappa[i] + l_u[i],
            }
        })
        .collect();
    (cones, model.distortion(&u))
}

fn sweep(dir: &Path, names: &[String], lambdas: &[f64], boundary: Boundary) {
    let default: Vec<f64> = (0..=12).map(|k| 10f64.powf(-3.0 + 0.25 * k as f64)).collect();
    let lambdas = if lambdas.is_empty() { &default[..] } else { lambdas };
    for name in names {
        let data = load(dir, name);
        let nv = data.mesh.n_vertices();
        let closed = data.mesh.boundary_loops().is_empty();
        let area = data.geometry.total_area();
        let file_total: f64 = data.cones.iter().map(|c| c.angle.abs()).sum();
        let free = boundary == Boundary::Free && !closed;
        println!(
            "{name}: file has {} cones, total |angle| {:.4}; area {:.4}; penalized problem with {}",
            data.cones.len(),
            file_total,
            area,
            if closed {
                "no boundary"
            } else if free {
                "free boundary values (solved on the double of the mesh)"
            } else {
                "u = 0 on the boundary"
            }
        );
        println!(
            "  {:>10} {:>10} {:>5} {:>9} {:>28} | {:>6} {:>4} {:>5} {:>9} {:>9} {:>9} | {:>6} {:>4} {:>5} {:>9} {:>9} {:>9}",
            "lambda", "lambda*A", "raw", "sum |mu|", "E   | common   max |d|",
            "merged", "same", "ring1", "sum refit", "E refit", "max |d|",
            "signed", "same", "ring1", "sum refit", "E refit", "max |d|"
        );
        println!(
            "  {:>238}",
            "| raw angles summed over the nearest file cone: max |d|, farthest | refit on the raw support, then merged: cones, same, sum, max |d|"
        );
        let doubled = if free {
            let Some(doubled) = double(&data) else {
                println!(
                    "  skipped: a face has three boundary vertices, so the mesh cannot be doubled"
                );
                continue;
            };
            Some(doubled)
        } else {
            None
        };
        for &lambda in lambdas {
            let solution = match &doubled {
                Some((mesh, geometry)) => {
                    let mut costs = vec![1.0; mesh.n_vertices()];
                    for (v, cost) in costs.iter_mut().enumerate().take(nv) {
                        if data.mesh.is_boundary_vertex(v) {
                            *cost = 0.0;
                        }
                    }
                    let problem = ConeProblem {
                        cone_costs: Some(costs),
                        ..ConeProblem::new(mesh, geometry, 0.5 * lambda)
                    };
                    solve(&problem, &SolveOptions::default())
                }
                None => solve(
                    &ConeProblem::new(&data.mesh, &data.geometry, lambda),
                    &SolveOptions::default(),
                ),
            }
            .expect("valid problem");
            if solution.report.status != Status::Converged {
                println!("  {lambda:>10.4e} not converged: {:?}", solution.report.message);
                continue;
            }
            // The measure on the interior vertices of the original mesh.
            let mu: Vec<f64> = (0..nv)
                .map(|v| {
                    if data.mesh.is_boundary_vertex(v) {
                        0.0
                    } else {
                        solution.measure[v]
                    }
                })
                .collect();
            let raw = mu.iter().filter(|m| m.abs() > 1e-12).count();
            let total: f64 = mu.iter().map(|m| m.abs()).sum();
            // On the double the distortion is per unit of the doubled area, i.e. the same.
            let mut line = format!(
                "  {:>10.4e} {:>10.4e} {:>5} {:>9.4} {:>9.5}",
                lambda,
                lambda * area,
                raw,
                total,
                solution.distortion_optimal
            );
            // The raw optimum against the file, without merging or re-optimizing.
            let file: BTreeMap<u32, f64> = data.cones.iter().map(|c| (c.vertex, c.angle)).collect();
            let common: Vec<f64> = (0..nv)
                .filter(|&v| mu[v].abs() > 1e-12)
                .filter_map(|v| file.get(&(v as u32)).map(|a| (a - mu[v]).abs()))
                .collect();
            line += &format!(
                " | {:>6} {:>9.2e}",
                common.len(),
                common.iter().cloned().fold(f64::NAN, f64::max)
            );
            if raw > 1500 {
                println!("{line} | (too many cones to merge)");
                continue;
            }
            for by_sign in [false, true] {
                let merged = merge(&data, &mu, by_sign);
                let distances = ring_distances(&data.mesh, &data.cones, &merged, 1);
                let count = |d: usize| distances.iter().filter(|&&x| x <= d).count();
                let refitted = if merged.len() <= 3 * data.cones.len().max(4) && !merged.is_empty() {
                    let (cones, distortion) = refit_positions(&data, boundary, &merged);
                    let ours: BTreeMap<u32, f64> =
                        cones.iter().map(|c| (c.vertex, c.angle)).collect();
                    let worst = data
                        .cones
                        .iter()
                        .filter_map(|c| ours.get(&c.vertex).map(|a| (a - c.angle).abs()))
                        .fold(f64::NAN, f64::max);
                    format!(
                        "{:>9.4} {:>9.5} {:>9.2e}",
                        cones.iter().map(|c| c.angle.abs()).sum::<f64>(),
                        distortion,
                        worst
                    )
                } else {
                    format!("{:>9} {:>9} {:>9}", "-", "-", "-")
                };
                line += &format!(
                    " | {:>6} {:>4} {:>5} {}",
                    merged.len(),
                    count(0),
                    count(1),
                    refitted
                );
            }
            // The raw optimum summed over the Voronoi cells of the file's cones (straight-line
            // distance): would a coarser grouping than edge-connected clusters give the file?
            {
                let mut sums = vec![0.0; data.cones.len()];
                let mut farthest = 0.0f64;
                for v in (0..nv).filter(|&v| mu[v].abs() > 1e-12) {
                    let distance = |c: &Cone| -> f64 {
                        let (p, q) = (data.positions[v], data.positions[c.vertex as usize]);
                        (0..3).map(|k| (p[k] - q[k]).powi(2)).sum::<f64>().sqrt()
                    };
                    let nearest = (0..data.cones.len())
                        .min_by(|&a, &b| distance(&data.cones[a]).total_cmp(&distance(&data.cones[b])))
                        .expect("the file has cones");
                    sums[nearest] += mu[v];
                    farthest = farthest.max(distance(&data.cones[nearest]));
                }
                let worst = sums
                    .iter()
                    .zip(&data.cones)
                    .map(|(s, c)| (s - c.angle).abs())
                    .fold(0.0, f64::max);
                line += &format!(" | {:>9.2e} {:>7.3}", worst, farthest);
            }
            // The other order: re-optimize the angles on the raw support, then merge.
            if raw > 0 && raw <= 3 * data.cones.len().max(20) {
                let support: Vec<Cone> = (0..nv)
                    .filter(|&v| mu[v].abs() > 1e-12)
                    .map(|v| Cone {
                        vertex: v as u32,
                        angle: mu[v],
                    })
                    .collect();
                let (cones, _) = refit_positions(&data, boundary, &support);
                let mut refitted = vec![0.0; nv];
                for c in &cones {
                    refitted[c.vertex as usize] = c.angle;
                }
                let merged = merge(&data, &refitted, false);
                let distances = ring_distances(&data.mesh, &data.cones, &merged, 1);
                let ours: BTreeMap<u32, f64> =
                    merged.iter().map(|c| (c.vertex, c.angle)).collect();
                let worst = data
                    .cones
                    .iter()
                    .filter_map(|c| ours.get(&c.vertex).map(|a| (a - c.angle).abs()))
                    .fold(f64::NAN, f64::max);
                line += &format!(
                    " | {:>6} {:>4} {:>9.4} {:>9.2e}",
                    merged.len(),
                    distances.iter().filter(|&&x| x == 0).count(),
                    merged.iter().map(|c| c.angle.abs()).sum::<f64>(),
                    worst
                );
            }
            println!("{line}");
        }
    }
    println!(
        "\nlambda: for area-normalized masses; lambda*A: for the masses of the OBJ as it is.\n\
         raw, sum |mu|, E: cones, total angle and distortion of the optimum of the penalized problem.\n\
         common, max |d|: cones of that optimum at a cone vertex of the file, and the largest angle \
         difference there.\n\
         merged: cones after merging every edge-connected cluster; signed: after merging the \
         positive and the negative cones separately.\nsame / ring1: number of the file's cones \
         with one of these cones at the same vertex / within one edge.\nsum refit, E refit: total \
         angle and distortion after replacing the angles at these positions by the \
         least-distortion angles (same boundary model).\nmax |d|: largest difference between a \
         refitted angle and the file's angle at the same vertex."
    );
}
