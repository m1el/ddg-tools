//! Shared helpers of the integration tests: random meshes of every topology, random problem
//! features, and a reference solution of the canonical problem by a generic convex solver.
#![allow(dead_code)]

use clarabel::algebra::CscMatrix;
use clarabel::solver::{DefaultSettings, DefaultSolver, IPSolver, SolverStatus, SupportedConeT};
use ddg_core::generators::{
    GeneratedMesh, SplitMix64, annulus, bumpy_grid_disk, icosphere, open_cylinder, perturb, torus,
};
use ddg_core::operators::{
    angle_defects, boundary_geodesic_curvatures, cotan_laplacian, lumped_mass_normalized,
};
use ddg_core::{Geometry, TriMesh};
use ddg_optimal_cones::{
    BoundaryTreatment, Cone, ConeError, ConeProblem, MethodUsed, Solution, SolveMethod,
    SolveOptions, factor_once_support, solve_with_boundary,
};

/// [`solve_with_boundary`] with the method chosen by the environment variable
/// `DDG_CONES_METHOD`: unset or `auto` leaves `options` alone; `reference` forces the reference
/// method; `fast` forces the factor-once method wherever it supports the problem (with the
/// limits for giving up lifted, so that the test measures correctness, not thresholds) and
/// leaves the problems it does not support to the reference method.
pub fn solve_t(
    problem: &ConeProblem<'_>,
    boundary: BoundaryTreatment,
    options: &SolveOptions<'_>,
) -> Result<Solution, ConeError> {
    let mut options = options.clone();
    let mode = std::env::var("DDG_CONES_METHOD").unwrap_or_default();
    let mut forced_fast = false;
    match mode.as_str() {
        "reference" => options.method = SolveMethod::Reference,
        "fast" => {
            if factor_once_support(problem, boundary)?.is_none() {
                options.method = SolveMethod::FactorOnce;
                options.factor_once.max_working_set_fraction = 1e9;
                options.factor_once.min_working_set = usize::MAX / 2;
                options.factor_once.memory_limit_bytes = usize::MAX;
                forced_fast = true;
            } else {
                options.method = SolveMethod::Reference;
            }
        }
        _ => {}
    }
    let result = solve_with_boundary(problem, boundary, &options);
    if let (true, Ok(s)) = (forced_fast, &result) {
        assert_eq!(s.report.method.used, MethodUsed::FactorOnce);
    }
    result
}

/// A mesh with its geometry.
pub struct TestMesh {
    pub name: String,
    pub positions: Vec<[f64; 3]>,
    pub mesh: TriMesh,
    pub geometry: Geometry,
}

impl TestMesh {
    pub fn from_positions(name: &str, positions: Vec<[f64; 3]>, faces: Vec<[u32; 3]>) -> TestMesh {
        let mesh = TriMesh::new(positions.len(), faces).expect("valid test mesh");
        let geometry = Geometry::from_positions(&mesh, &positions).expect("valid test geometry");
        TestMesh {
            name: name.to_string(),
            positions,
            mesh,
            geometry,
        }
    }

    pub fn interior(&self) -> Vec<usize> {
        (0..self.mesh.n_vertices())
            .filter(|&v| !self.mesh.is_boundary_vertex(v))
            .collect()
    }

    pub fn is_closed(&self) -> bool {
        self.mesh.boundary_loops().is_empty()
    }
}

pub const N_KINDS: usize = 5;

/// A small random mesh. `kind`: 0 disk, 1 annulus, 2 open cylinder, 3 sphere, 4 torus.
/// Deterministic in its arguments.
pub fn random_mesh(kind: usize, seed: u64) -> TestMesh {
    let mut rng = SplitMix64::new(seed ^ 0x9e37_79b9_7f4a_7c15);
    let pick =
        |rng: &mut SplitMix64, lo: usize, hi: usize| lo + (rng.next_u64() as usize) % (hi - lo + 1);
    let (name, base): (&str, GeneratedMesh) = match kind % N_KINDS {
        0 => {
            let (nx, ny) = (pick(&mut rng, 2, 5), pick(&mut rng, 2, 5));
            let amplitude = rng.uniform(0.0, 0.4);
            ("disk", bumpy_grid_disk(nx, ny, amplitude).unwrap())
        }
        1 => {
            let (na, nr) = (pick(&mut rng, 5, 9), pick(&mut rng, 2, 3));
            (
                "annulus",
                annulus(na, nr, rng.uniform(0.3, 0.6), 1.0).unwrap(),
            )
        }
        2 => {
            let (na, nh) = (pick(&mut rng, 4, 8), pick(&mut rng, 2, 4));
            (
                "cylinder",
                open_cylinder(na, nh, rng.uniform(0.4, 1.0), rng.uniform(0.8, 2.0)).unwrap(),
            )
        }
        3 => ("sphere", icosphere(pick(&mut rng, 0, 1) as u32).unwrap()),
        _ => {
            let (n_major, n_minor) = (pick(&mut rng, 4, 7), pick(&mut rng, 3, 5));
            (
                "torus",
                torus(n_major, n_minor, 1.0, rng.uniform(0.3, 0.6)).unwrap(),
            )
        }
    };
    // Anisotropic scaling and noise make the geometry generic; retry with less noise if a
    // triangle degenerates.
    let scale = [
        rng.uniform(0.6, 1.6),
        rng.uniform(0.6, 1.6),
        rng.uniform(0.6, 1.6),
    ];
    let mut amplitude = rng.uniform(0.0, 0.08);
    let noise_seed = rng.next_u64();
    loop {
        let mut positions = base.positions.clone();
        perturb(&mut positions, [amplitude; 3], noise_seed);
        for p in &mut positions {
            for k in 0..3 {
                p[k] *= scale[k];
            }
        }
        let mesh = TriMesh::new(positions.len(), base.faces.clone()).expect("generator mesh");
        if let Ok(geometry) = Geometry::from_positions(&mesh, &positions) {
            return TestMesh {
                name: format!("{name}-{seed}"),
                positions,
                mesh,
                geometry,
            };
        }
        amplitude *= 0.5;
    }
}

/// The optional parts of a problem.
#[derive(Debug, Clone, Default)]
pub struct Features {
    pub lambda: f64,
    pub distortion_weights: Option<Vec<f64>>,
    pub cone_costs: Option<Vec<f64>>,
    pub lower_bounds: Option<Vec<f64>>,
    pub upper_bounds: Option<Vec<f64>>,
    pub prescribed_cones: Vec<Cone>,
}

impl Features {
    pub fn plain(lambda: f64) -> Features {
        Features {
            lambda,
            ..Features::default()
        }
    }

    pub fn problem<'a>(&self, mesh: &'a TestMesh) -> ConeProblem<'a> {
        ConeProblem {
            mesh: &mesh.mesh,
            geometry: &mesh.geometry,
            lambda: self.lambda,
            distortion_weights: self.distortion_weights.clone(),
            cone_costs: self.cone_costs.clone(),
            lower_bounds: self.lower_bounds.clone(),
            upper_bounds: self.upper_bounds.clone(),
            prescribed_cones: self.prescribed_cones.clone(),
        }
    }

    /// Whether every interior vertex has a positive distortion weight (then the minimizer is
    /// unique).
    pub fn strictly_weighted(&self, mesh: &TestMesh) -> bool {
        match &self.distortion_weights {
            None => true,
            Some(w) => mesh.interior().iter().all(|&v| w[v] > 0.0),
        }
    }
}

/// Which features [`random_features`] may switch on.
#[derive(Debug, Clone, Copy)]
pub struct FeatureMask {
    pub weights: bool,
    pub zero_distortion_weights: bool,
    pub free_sites: bool,
    pub forbidden: bool,
    pub bounds: bool,
    pub prescribed: bool,
}

impl FeatureMask {
    pub const ALL: FeatureMask = FeatureMask {
        weights: true,
        zero_distortion_weights: true,
        free_sites: true,
        forbidden: true,
        bounds: true,
        prescribed: true,
    };
    /// Everything except zero distortion weights, so that the minimizer is unique.
    pub const UNIQUE: FeatureMask = FeatureMask {
        zero_distortion_weights: false,
        ..FeatureMask::ALL
    };
}

/// Random problem features, each switched on with probability about one half (if allowed).
pub fn random_features(mesh: &TestMesh, seed: u64, mask: FeatureMask) -> Features {
    let mut rng = SplitMix64::new(seed ^ 0x5851_f42d_4c95_7f2d);
    let nv = mesh.mesh.n_vertices();
    let interior = mesh.interior();
    let coin = |rng: &mut SplitMix64, p: f64| rng.next_f64() < p;
    // λ log-uniform in [1e-3, 1]: from many cones to none.
    let lambda = 10f64.powf(rng.uniform(-3.0, 0.0));
    let mut f = Features::plain(lambda);
    if mask.weights && coin(&mut rng, 0.5) {
        let mut w: Vec<f64> = (0..nv).map(|_| 2f64.powf(rng.uniform(-2.0, 2.0))).collect();
        if mask.zero_distortion_weights && coin(&mut rng, 0.3) {
            for x in &mut w {
                if rng.next_f64() < 0.3 {
                    *x = 0.0;
                }
            }
            // Keep at least one positive weight at an unknown vertex.
            if let Some(&v) = interior.first() {
                w[v] = 1.0;
            }
        }
        f.distortion_weights = Some(w);
    }
    let mut costs: Option<Vec<f64>> = None;
    if mask.weights && coin(&mut rng, 0.5) {
        costs = Some((0..nv).map(|_| 2f64.powf(rng.uniform(-2.0, 2.0))).collect());
    }
    if mask.free_sites && coin(&mut rng, 0.4) {
        let c = costs.get_or_insert_with(|| vec![1.0; nv]);
        for _ in 0..1 + (rng.next_u64() % 2) {
            c[(rng.next_u64() as usize) % nv] = 0.0;
        }
    }
    if mask.forbidden && coin(&mut rng, 0.4) {
        let c = costs.get_or_insert_with(|| vec![1.0; nv]);
        let share = rng.uniform(0.05, 0.5);
        for x in c.iter_mut() {
            if rng.next_f64() < share {
                *x = f64::INFINITY;
            }
        }
    }
    f.cone_costs = costs;
    if mask.bounds && coin(&mut rng, 0.5) {
        // A mixture of common settings: positive only, symmetric caps, per-vertex bounds.
        match rng.next_u64() % 4 {
            0 => f.lower_bounds = Some(vec![0.0; nv]),
            1 => {
                let cap = rng.uniform(0.2, 2.0);
                f.lower_bounds = Some(vec![-cap; nv]);
                f.upper_bounds = Some(vec![cap; nv]);
            }
            2 => {
                f.upper_bounds = Some((0..nv).map(|_| rng.uniform(0.0, 2.5)).collect());
            }
            _ => {
                f.lower_bounds = Some(
                    (0..nv)
                        .map(|_| {
                            if rng.next_f64() < 0.3 {
                                0.0
                            } else {
                                -rng.uniform(0.0, 2.0)
                            }
                        })
                        .collect(),
                );
                f.upper_bounds = Some(
                    (0..nv)
                        .map(|_| {
                            if rng.next_f64() < 0.2 {
                                f64::INFINITY
                            } else {
                                rng.uniform(0.0, 3.0)
                            }
                        })
                        .collect(),
                );
            }
        }
    }
    if mask.prescribed && coin(&mut rng, 0.4) && !interior.is_empty() {
        let mut used = vec![false; nv];
        for _ in 0..1 + (rng.next_u64() % 2) {
            let v = interior[(rng.next_u64() as usize) % interior.len()];
            if !used[v] {
                used[v] = true;
                f.prescribed_cones.push(Cone {
                    vertex: v as u32,
                    angle: rng.uniform(-1.5, 1.5),
                });
            }
        }
    }
    f
}

/// The outcome of the generic convex solver.
#[derive(Debug, Clone)]
pub enum Reference {
    Solved {
        objective: f64,
        /// Per vertex, zero on the boundary.
        log_scale: Vec<f64>,
        /// Per vertex, zero on the boundary.
        measure: Vec<f64>,
    },
    Infeasible,
    /// The generic solver did not reach its tolerances (its status is given).
    Unreliable(String),
}

/// Solves the canonical problem directly as a quadratic program with a generic interior-point
/// solver, in the variables `(u, µ, s)` on the unknown vertices:
///
/// ```text
/// minimize    ½ Σ wE_i m_i u_i² + λ Σ wR_i s_i
/// subject to  L_II u − µ = −Ω,   −s ≤ µ ≤ s,   lo ≤ µ ≤ hi,
///             µ_i = 0 at forbidden vertices,   µ_i = a_i at prescribed cones.
/// ```
///
/// Nothing of the solver under test is used: the operators come from `ddg-core`, the problem
/// data straight from the public problem description.
pub fn reference_solve(mesh: &TestMesh, f: &Features) -> Reference {
    reference_solve_with_boundary(mesh, f, BoundaryTreatment::Fixed)
}

/// [`reference_solve`] with the boundary values as unknowns. At a boundary vertex `b` the
/// state equation is `(L u)_b − ν_b = −k_b` (`k_b = π −` angle sum), and `ν_b` is
/// `Free`: unconstrained at no cost; `Penalized`: charged `λ wR_b |ν_b|` (`+∞` forces
/// `ν_b = 0`); `Convex`: `ν_b ≥ 0` at no cost. Everything is written out here from the problem
/// statement, not taken from the solver.
pub fn reference_solve_with_boundary(
    mesh: &TestMesh,
    f: &Features,
    boundary: BoundaryTreatment,
) -> Reference {
    let nv = mesh.mesh.n_vertices();
    let unknowns: Vec<usize> = if boundary == BoundaryTreatment::Fixed {
        mesh.interior()
    } else {
        (0..nv).collect()
    };
    let on_boundary = |v: usize| mesh.mesh.is_boundary_vertex(v);
    let n = unknowns.len();
    if n == 0 {
        return Reference::Solved {
            objective: 0.0,
            log_scale: vec![0.0; nv],
            measure: vec![0.0; nv],
        };
    }
    let mut index = vec![usize::MAX; nv];
    for (i, &v) in unknowns.iter().enumerate() {
        index[v] = i;
    }
    let laplacian = cotan_laplacian(&mesh.mesh, &mesh.geometry).unwrap();
    let mut omega = angle_defects(&mesh.mesh, &mesh.geometry).unwrap();
    if boundary != BoundaryTreatment::Fixed {
        let k = boundary_geodesic_curvatures(&mesh.mesh, &mesh.geometry).unwrap();
        for (o, k) in omega.iter_mut().zip(k) {
            *o += k;
        }
    }
    let mass = lumped_mass_normalized(&mesh.mesh, &mesh.geometry).unwrap();
    let w_e = |v: usize| f.distortion_weights.as_ref().map_or(1.0, |w| w[v]);
    let w_r = |v: usize| {
        if on_boundary(v) {
            match boundary {
                BoundaryTreatment::Penalized => f.cone_costs.as_ref().map_or(1.0, |w| w[v]),
                _ => 0.0,
            }
        } else {
            f.cone_costs.as_ref().map_or(1.0, |w| w[v])
        }
    };
    let lo = |v: usize| {
        if on_boundary(v) {
            if boundary == BoundaryTreatment::Convex {
                0.0
            } else {
                f64::NEG_INFINITY
            }
        } else {
            f.lower_bounds.as_ref().map_or(f64::NEG_INFINITY, |b| b[v])
        }
    };
    let hi = |v: usize| {
        if on_boundary(v) {
            f64::INFINITY
        } else {
            f.upper_bounds.as_ref().map_or(f64::INFINITY, |b| b[v])
        }
    };
    let mut fixed: Vec<Option<f64>> = unknowns
        .iter()
        .map(|&v| (w_r(v) == f64::INFINITY).then_some(0.0))
        .collect();
    for cone in &f.prescribed_cones {
        fixed[index[cone.vertex as usize]] = Some(cone.angle);
    }

    // Variables: u (0..n), µ (n..2n), s (2n..3n).
    let nvar = 3 * n;
    let p = CscMatrix::new(
        nvar,
        nvar,
        (0..=nvar).map(|j| j.min(n)).collect(),
        (0..n).collect(),
        unknowns.iter().map(|&v| w_e(v) * mass[v]).collect(),
    );
    let mut q = vec![0.0; nvar];
    for (i, &v) in unknowns.iter().enumerate() {
        if fixed[i].is_none() {
            q[2 * n + i] = f.lambda * w_r(v);
        }
    }
    // Constraint rows as triplets; equalities first.
    let mut rows: Vec<usize> = Vec::new();
    let mut cols: Vec<usize> = Vec::new();
    let mut vals: Vec<f64> = Vec::new();
    let mut b: Vec<f64> = Vec::new();
    let mut push = |row: usize, col: usize, val: f64| {
        rows.push(row);
        cols.push(col);
        vals.push(val);
    };
    // State equation.
    for (r, c, value) in laplacian.triplets() {
        let (r, c) = (index[r as usize], index[c as usize]);
        if r != usize::MAX && c != usize::MAX {
            push(r, c, value);
        }
    }
    for (i, &v) in unknowns.iter().enumerate() {
        push(i, n + i, -1.0);
        b.push(-omega[v]);
    }
    // Fixed cone angles.
    for (i, angle) in fixed.iter().enumerate() {
        if let Some(angle) = angle {
            push(b.len(), n + i, 1.0);
            b.push(*angle);
        }
    }
    let n_equalities = b.len();
    // µ − s ≤ 0 and −µ − s ≤ 0.
    for i in 0..n {
        push(b.len(), n + i, 1.0);
        push(b.len(), 2 * n + i, -1.0);
        b.push(0.0);
        push(b.len(), n + i, -1.0);
        push(b.len(), 2 * n + i, -1.0);
        b.push(0.0);
    }
    // Bounds of the penalized vertices.
    for (i, &v) in unknowns.iter().enumerate() {
        if fixed[i].is_some() {
            continue;
        }
        if hi(v).is_finite() {
            push(b.len(), n + i, 1.0);
            b.push(hi(v));
        }
        if lo(v).is_finite() {
            push(b.len(), n + i, -1.0);
            b.push(-lo(v));
        }
    }
    let a = CscMatrix::new_from_triplets(b.len(), nvar, rows, cols, vals);
    let cones = [
        SupportedConeT::ZeroConeT(n_equalities),
        SupportedConeT::NonnegativeConeT(b.len() - n_equalities),
    ];
    let settings = DefaultSettings {
        verbose: false,
        max_iter: 500,
        tol_gap_abs: 1e-11,
        tol_gap_rel: 1e-11,
        tol_feas: 1e-11,
        ..DefaultSettings::default()
    };
    let mut solver = match DefaultSolver::new(&p, &q, &a, &b, &cones, settings) {
        Ok(solver) => solver,
        Err(e) => return Reference::Unreliable(format!("set-up failed: {e:?}")),
    };
    solver.solve();
    match solver.solution.status {
        SolverStatus::Solved => {
            let x = &solver.solution.x;
            let mut log_scale = vec![0.0; nv];
            let mut measure = vec![0.0; nv];
            for (i, &v) in unknowns.iter().enumerate() {
                log_scale[v] = x[i];
                measure[v] = x[n + i];
            }
            // Evaluate the objective from (u, µ) itself rather than from the slack `s`.
            let objective: f64 = unknowns
                .iter()
                .enumerate()
                .map(|(i, &v)| {
                    let cost = if fixed[i].is_some() {
                        0.0
                    } else {
                        f.lambda * w_r(v) * x[n + i].abs()
                    };
                    0.5 * w_e(v) * mass[v] * x[i] * x[i] + cost
                })
                .sum();
            Reference::Solved {
                objective,
                log_scale,
                measure,
            }
        }
        SolverStatus::PrimalInfeasible | SolverStatus::AlmostPrimalInfeasible => {
            Reference::Infeasible
        }
        other => Reference::Unreliable(format!("{other:?}")),
    }
}

pub fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
}

/// A compact one-paragraph description of what a solve did, for failure messages.
pub fn summary(solution: &ddg_optimal_cones::Solution) -> String {
    use ddg_optimal_cones::StageKind;
    let r = &solution.report;
    let stages: Vec<String> = r
        .stages
        .iter()
        .map(|s| {
            let kind = match s.kind {
                StageKind::Regularized { gamma } => format!("R{gamma:.0e}"),
                StageKind::Proximal { gamma } => format!("P{gamma:.0e}"),
                StageKind::Exact => "E".to_string(),
            };
            format!(
                "{kind}:{}/{}/{}:{:?}",
                s.newton_iterations, s.damped_steps, s.constant_shifts, s.stop
            )
        })
        .collect();
    let c = &r.certificate;
    format!(
        "status {:?}, objective {:.12e}, {} cones, safeguard rounds {}, LDLt {} LU {} \
         (fallbacks {}), residuals: state {:.2e} adjoint {:.2e} cone {:.2e} bounds {:.2e} dual \
         {:.2e} GB {:.2e}, gap {:.2e} / {:.2e}, message {:?}\n  stages \
         (kind:newton/damped/shifts:stop): {}",
        r.status,
        solution.objective,
        solution.cones.len(),
        r.safeguard_rounds,
        r.linear_solver.ldlt_factorizations,
        r.linear_solver.lu_factorizations,
        r.linear_solver.ldlt_fallbacks,
        c.state_residual,
        c.adjoint_residual,
        c.cone_residual,
        c.bound_violation,
        c.dual_infeasibility,
        c.gauss_bonnet_residual,
        c.duality_gap,
        c.complementarity_gap,
        r.message,
        stages.join(" ")
    )
}

/// The mesh glued to a mirror copy of itself along its boundary: a closed mesh. Vertex `v` of
/// the original keeps its index; `copy[v]` is its mirror image (`v` itself on the boundary).
///
/// Any quantity per vertex of the original is, on the double, the same at both copies; at a
/// seam vertex the Laplacian row, the angle defect (`2k`) and the cone angle (`2ν`) are twice
/// those of the boundary vertex of the original, and the mass is the same (the area is
/// doubled, and so is the star of a seam vertex). That is why the double is a route to the
/// optimized boundary values that does not use the boundary code: see
/// [`doubled_features`].
///
/// `None` if the double is not a valid mesh: a triangle with three boundary vertices would be
/// glued to a second copy of itself.
pub fn double_mesh(mesh: &TriMesh, geometry: &Geometry) -> Option<(TriMesh, Geometry, Vec<u32>)> {
    let nv = mesh.n_vertices();
    let mut copy: Vec<u32> = (0..nv as u32).collect();
    let mut next = nv as u32;
    for (v, slot) in copy.iter_mut().enumerate() {
        if !mesh.is_boundary_vertex(v) {
            *slot = next;
            next += 1;
        }
    }
    let mut faces = mesh.faces().to_vec();
    faces.extend(mesh.faces().iter().map(|f| {
        [
            copy[f[0] as usize],
            copy[f[2] as usize],
            copy[f[1] as usize],
        ]
    }));
    let doubled = TriMesh::new(next as usize, faces).ok()?;
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
            geometry.edge_lengths()[edge]
        })
        .collect();
    let geometry = Geometry::from_edge_lengths(&doubled, lengths).expect("the double's geometry");
    Some((doubled, geometry, copy))
}

/// The features of the problem on the double that corresponds to the problem `f` on the
/// original with the boundary treatment `boundary` (`Free`, `Penalized` or `Convex`): interior
/// features at both copies, `λ/2`, and at the seam vertices the site of the treatment (cost 0
/// for `Free` and `Convex`, `wR` for `Penalized`; angle bound `≥ 0` for `Convex`). Every
/// objective and every solution value of the original are those of the double (a seam cone
/// angle is `2ν`).
pub fn doubled_features(
    mesh: &TriMesh,
    copy: &[u32],
    f: &Features,
    boundary: BoundaryTreatment,
) -> Features {
    let nv = mesh.n_vertices();
    let n2 = copy.iter().copied().max().unwrap() as usize + 1;
    let spread = |values: &Option<Vec<f64>>, default: f64, seam: &dyn Fn(usize) -> f64| {
        values.as_ref().map(|_| ()).map(|_| {
            let mut out = vec![default; n2];
            for v in 0..nv {
                let x = if mesh.is_boundary_vertex(v) {
                    seam(v)
                } else {
                    values.as_ref().unwrap()[v]
                };
                out[v] = x;
                out[copy[v] as usize] = x;
            }
            out
        })
    };
    let distortion_weights = spread(&f.distortion_weights, 1.0, &|v| {
        f.distortion_weights.as_ref().unwrap()[v]
    });
    let seam_cost = |v: usize| match boundary {
        BoundaryTreatment::Penalized => f.cone_costs.as_ref().map_or(1.0, |c| c[v]),
        _ => 0.0,
    };
    // The seam cost differs from the default even if the original has no costs.
    let costs = {
        let mut out = vec![1.0; n2];
        for v in 0..nv {
            let x = if mesh.is_boundary_vertex(v) {
                seam_cost(v)
            } else {
                f.cone_costs.as_ref().map_or(1.0, |c| c[v])
            };
            out[v] = x;
            out[copy[v] as usize] = x;
        }
        Some(out)
    };
    let lower = {
        let mut out = vec![f64::NEG_INFINITY; n2];
        for v in 0..nv {
            let x = if mesh.is_boundary_vertex(v) {
                if boundary == BoundaryTreatment::Convex {
                    0.0
                } else {
                    f64::NEG_INFINITY
                }
            } else {
                f.lower_bounds.as_ref().map_or(f64::NEG_INFINITY, |b| b[v])
            };
            out[v] = x;
            out[copy[v] as usize] = x;
        }
        Some(out)
    };
    let upper = {
        let mut out = vec![f64::INFINITY; n2];
        for v in 0..nv {
            let x = if mesh.is_boundary_vertex(v) {
                f64::INFINITY
            } else {
                f.upper_bounds.as_ref().map_or(f64::INFINITY, |b| b[v])
            };
            out[v] = x;
            out[copy[v] as usize] = x;
        }
        Some(out)
    };
    let mut prescribed = Vec::new();
    for c in &f.prescribed_cones {
        prescribed.push(*c);
        prescribed.push(Cone {
            vertex: copy[c.vertex as usize],
            angle: c.angle,
        });
    }
    Features {
        lambda: 0.5 * f.lambda,
        distortion_weights,
        cone_costs: costs,
        lower_bounds: lower,
        upper_bounds: upper,
        prescribed_cones: prescribed,
    }
}
