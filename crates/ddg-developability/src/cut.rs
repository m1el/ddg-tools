//! Cut and flatten (SPEC §6; paper Sec. 4.4; D22, D23).
//!
//! Terminals are the vertices whose developability energy exceeds `epsilon` or whose angle
//! defect exceeds `defect_tolerance` (the latter is an addition of ours, SPEC D30). A cut graph
//! through them is planned with [`ddg_bff::plan_cut`] (Steiner tree approximation plus handle
//! loops), the surface is cut, and every chart (a disk) is flattened with Boundary First
//! Flattening with zero boundary scale factors, i.e. isometrically along the boundary, so that
//! cutting along the singular set of a developable surface yields (nearly) isometric charts.
//!
//! Out of scope (D23): gluing small cracks and re-flattening, automatic avoidance of global
//! overlap (the paper does both by hand); overlap is reported per chart.

use ddg_bff::quality::{log_scale_factor, signed_double_area};
use ddg_bff::{
    BffError, ChartReport, CutPlan, CutPlanInput, DiskBff, DiskOptions, Extension, component_topology, cut_mesh,
    plan_cut,
};
use ddg_core::{Geometry, TriMesh};
use rayon::prelude::*;

use crate::eig::sym_eig3;
use crate::energy::{Energy, EnergyOptions, evaluate_energy};
use crate::error::DevError;
use crate::geometry::{corner_angles, face_frame, v3};

/// Default of [`CutOptions::defect_tolerance`] in radians (SPEC D30).
pub const DEFAULT_DEFECT_TOLERANCE: f64 = 0.1;

/// Options of [`cut_and_flatten`] (SPEC §6).
#[derive(Debug, Clone, PartialEq)]
pub struct CutOptions {
    /// The energy whose per-vertex values select the terminals (default `Covariance`). It is
    /// evaluated with `omit_valence_three = false` and boundary vertices included.
    pub energy: Energy,
    /// A vertex is a terminal when its energy exceeds this (default `1e-3`; paper Fig. 20 uses
    /// 0.01 … 0.001). An energy that is NaN never makes a terminal.
    pub epsilon: f64,
    /// Additional terminal criterion: an interior vertex with `|Ω_v| > defect_tolerance`
    /// (angle defect `Ω_v = 2π − Σθ`, radians) is a terminal as well. Default
    /// `Some(`[`DEFAULT_DEFECT_TOLERANCE`]`)` = 0.1 rad; `None` is the paper's rule (energy
    /// threshold only), SPEC D30. A cut through a vertex lets the flattening absorb its defect
    /// isometrically; the energies cannot see a defect on a star whose normals are coplanar but
    /// which is not embedded (D28). [formula: angle-defect]
    pub defect_tolerance: Option<f64>,
    /// Additional terminals (vertex ids), e.g. user-chosen cone sites.
    pub extra_terminals: Vec<u32>,
    /// `c_len` of the edge weight `w_e = c_len ℓ_e/ℓ̄ + c_align (1 − ½(|ê·x_a| + |ê·x_b|))`
    /// (default 1). [formula: dev-cut-edge-weight]
    pub length_weight: f64,
    /// `c_align` (default 0: length only, SPEC §11 answer 2). `x_v` is the eigenvector of the
    /// smallest eigenvalue of the covariance matrix `A_v` (§3.4) at the ends of the edge.
    pub alignment_weight: f64,
    /// Per edge of the mesh: whether it may be cut (`None`: all). Polygon edges for polygon input.
    pub may_cut: Option<Vec<bool>>,
    /// Edges (ids of `TriMesh::edges()`) that must be cut (user seams).
    pub must_cut: Vec<u32>,
    /// How BFF extends the boundary curve (default `Harmonic`, which reproduces the boundary
    /// lengths exactly; `Conjugate` is BFF's automatic mode).
    pub extension: Extension,
    /// Sample density of BFF's overlap estimate (see `ddg_bff::DiskOptions`); `0` switches it
    /// off.
    pub overlap_samples_per_triangle: f64,
    /// Threads of the energy evaluation (see `EnergyOptions::threads`).
    pub threads: Option<usize>,
}

impl Default for CutOptions {
    fn default() -> Self {
        CutOptions {
            energy: Energy::Covariance,
            epsilon: 1e-3,
            defect_tolerance: Some(DEFAULT_DEFECT_TOLERANCE),
            extra_terminals: Vec::new(),
            length_weight: 1.0,
            alignment_weight: 0.0,
            may_cut: None,
            must_cut: Vec::new(),
            extension: Extension::Harmonic,
            overlap_samples_per_triangle: ddg_bff::DEFAULT_OVERLAP_SAMPLES,
            threads: None,
        }
    }
}

impl CutOptions {
    /// Checks the options (the lengths of `may_cut` and the ranges of ids are checked against
    /// the mesh by [`cut_and_flatten`]).
    ///
    /// # Errors
    /// [`DevError::InvalidOption`] naming the first violation.
    pub fn validate(&self) -> Result<(), DevError> {
        let bad = |name, value, requirement| Err(DevError::InvalidOption { name, value, requirement });
        if self.epsilon.is_nan() || self.epsilon < 0.0 {
            return bad("epsilon", self.epsilon, ">= 0");
        }
        if let Some(t) = self.defect_tolerance
            && (t.is_nan() || t < 0.0)
        {
            return bad("defect_tolerance", t, ">= 0 (and not NaN)");
        }
        if !(self.length_weight.is_finite() && self.length_weight >= 0.0) {
            return bad("length_weight", self.length_weight, "finite and >= 0");
        }
        if !(self.alignment_weight.is_finite() && self.alignment_weight >= 0.0) {
            return bad("alignment_weight", self.alignment_weight, "finite and >= 0");
        }
        if self.length_weight == 0.0 && self.alignment_weight == 0.0 {
            return bad("length_weight", 0.0, "> 0 unless alignment_weight > 0 (all edge weights would be zero)");
        }
        if !(self.overlap_samples_per_triangle.is_finite() && self.overlap_samples_per_triangle >= 0.0) {
            return bad("overlap_samples_per_triangle", self.overlap_samples_per_triangle, "finite and >= 0");
        }
        if self.threads == Some(0) {
            return bad("threads", 0.0, ">= 1 (or None)");
        }
        Ok(())
    }
}

/// Result of [`cut_and_flatten`] (SPEC §6).
#[derive(Debug, Clone, PartialEq)]
pub struct DevelopableFlattening {
    /// Edge ids (of the input mesh) that were cut, ascending (interior edges only).
    pub cut_edges: Vec<u32>,
    /// The vertex pair `[a, b]` of every cut edge, same order.
    pub cut_vertex_pairs: Vec<[u32; 2]>,
    /// UV of every triangle corner: entry `3 f + k` is the corner of face `f` at
    /// `faces[f][k]`. Charts are not packed: all start with their first boundary vertex at the
    /// origin and lengths are those of the surface (zero boundary scale).
    pub corner_uvs: Vec<[f64; 2]>,
    /// The chart of every face (components of the cut surface, numbered by lowest face).
    pub face_chart: Vec<u32>,
    /// `u_f = ½ ln(A_uv / A_f)` per face with `A_uv` the absolute plane area; `NaN` for a face
    /// with zero plane area. Flipped faces are counted in `charts`.
    pub face_log_scale: Vec<f64>,
    /// Area-weighted mean of `face_log_scale` over the finite faces at each vertex (`NaN` if none).
    pub vertex_log_scale: Vec<f64>,
    /// Area-weighted mean of `|face_log_scale|` over the finite faces: the distortion number
    /// reported by the studies (0 for an isometric development).
    pub mean_abs_log_scale: f64,
    /// Terminals from the energy threshold and `extra_terminals` (sorted, unique).
    pub terminals: Vec<u32>,
    /// Terminals added on closed genus-0 components that had fewer than two (sorted).
    pub added_terminals: Vec<u32>,
    /// The planned cut.
    pub plan: CutPlan,
    /// One BFF report per chart (`first_face` is a face of the input mesh).
    pub charts: Vec<ChartReport>,
    /// The per-vertex energies of `CutOptions::energy` (no valence-3 omission).
    pub vertex_energy: Vec<f64>,
    /// The angle defect `Ω_v` of every interior vertex of the input (0 on the boundary).
    pub vertex_defect: Vec<f64>,
    /// [`star_embeddedness`] of every vertex of the input.
    pub star_embeddedness: Vec<f64>,
}

impl DevelopableFlattening {
    /// Number of charts.
    pub fn n_charts(&self) -> usize {
        self.charts.len()
    }

    /// Total number of faces whose image in the plane is flipped, over all charts.
    pub fn flipped_triangles(&self) -> usize {
        self.charts.iter().map(|c| c.distortion.flipped_triangles).sum()
    }

    /// Whether no chart has flipped triangles, wedges over `2π` or boundary crossings.
    pub fn is_injective(&self) -> bool {
        self.charts.iter().all(|c| c.is_injective())
    }
}

/// A per-vertex embeddedness diagnostic of the stars (SPEC D28, D30): the smallest cosine
/// `min_σ n̂_σ · N_v` between the unit normal of a face of the star and the vertex normal `N_v`
/// (§1). A positive value means that the star projects injectively onto the tangent plane (it is
/// a graph over it, hence embedded); a value `<= 0` means that some face is vertical or
/// folded over (a spike, a needle or a double cover in the sense of the paper's Fig. 9, or a
/// sharp but embedded hinge: the test is conservative). `NaN` for a vertex with `S_v = 0`.
/// Boundary vertices are included.
pub fn star_embeddedness(mesh: &TriMesh, positions: &[[f64; 3]]) -> Result<Vec<f64>, DevError> {
    let nv = mesh.n_vertices();
    let mut s = vec![[0.0f64; 3]; nv];
    let mut frames = Vec::with_capacity(mesh.n_faces());
    for (f, tri) in mesh.faces().iter().enumerate() {
        let fr = face_frame(tri.map(|v| positions[v as usize])).ok_or(DevError::DegenerateFace { face: f as u32 })?;
        for &v in tri {
            s[v as usize] = v3::add(s[v as usize], fr.n);
        }
        frames.push(fr);
    }
    let mut out = vec![f64::INFINITY; nv];
    for (tri, fr) in mesh.faces().iter().zip(&frames) {
        for &v in tri {
            let sv = s[v as usize];
            let l = v3::dot(sv, sv).sqrt();
            let c = if l > 0.0 { v3::dot(fr.normal, sv) / l } else { f64::NAN };
            let o = &mut out[v as usize];
            if c.is_nan() || o.is_nan() {
                *o = f64::NAN;
            } else {
                *o = o.min(c);
            }
        }
    }
    for o in out.iter_mut() {
        if *o == f64::INFINITY {
            *o = f64::NAN; // isolated vertex
        }
    }
    Ok(out)
}

/// Edge weights `w_e = c_len ℓ_e/ℓ̄ + c_align (1 − ½(|ê·x_a| + |ê·x_b|))` (SPEC §6 step 2).
///
/// [formula: dev-cut-edge-weight]
fn edge_weights(
    mesh: &TriMesh,
    positions: &[[f64; 3]],
    lengths: &[f64],
    opts: &CutOptions,
) -> Result<Vec<f64>, DevError> {
    let lbar = lengths.iter().sum::<f64>() / lengths.len().max(1) as f64;
    let mut w: Vec<f64> = lengths.iter().map(|&l| opts.length_weight * l / lbar).collect();
    if opts.alignment_weight == 0.0 {
        return Ok(w);
    }
    // covariance matrices A_v = sum theta N N^T
    let nv = mesh.n_vertices();
    let mut a = vec![[[0.0f64; 3]; 3]; nv];
    for (f, tri) in mesh.faces().iter().enumerate() {
        let p = tri.map(|v| positions[v as usize]);
        let fr = face_frame(p).ok_or(DevError::DegenerateFace { face: f as u32 })?;
        let th = corner_angles(p);
        for k in 0..3 {
            let m = &mut a[tri[k] as usize];
            for (row, &ni) in m.iter_mut().zip(&fr.normal) {
                for (cell, &nj) in row.iter_mut().zip(&fr.normal) {
                    *cell += th[k] * ni * nj;
                }
            }
        }
    }
    let mut x = vec![[0.0; 3]; nv];
    for v in 0..nv {
        x[v] = sym_eig3(a[v])?.vectors[0];
    }
    for (e, &[p, q]) in mesh.edges().iter().enumerate() {
        let d = v3::sub(positions[q as usize], positions[p as usize]);
        let l = lengths[e];
        let align = if l > 0.0 {
            let (da, db) = (v3::dot(d, x[p as usize]).abs() / l, v3::dot(d, x[q as usize]).abs() / l);
            1.0 - 0.5 * (da + db)
        } else {
            0.0
        };
        w[e] += opts.alignment_weight * align.max(0.0);
    }
    Ok(w)
}

/// Adds the highest-energy vertex of `component` that is not yet a terminal (ties: lower id).
fn add_terminal(
    component: u32,
    vertex_component: &[u32],
    vertex_energy: &[f64],
    is_terminal: &mut [bool],
    added: &mut Vec<u32>,
) -> Result<(), DevError> {
    let mut best: Option<u32> = None;
    for v in 0..vertex_energy.len() {
        if vertex_component[v] != component || is_terminal[v] {
            continue;
        }
        let key = |i: usize| if vertex_energy[i].is_nan() { f64::NEG_INFINITY } else { vertex_energy[i] };
        if best.is_none_or(|b| key(v) > key(b as usize)) {
            best = Some(v as u32);
        }
    }
    let b = best.ok_or_else(|| DevError::Internal(format!("component {component} has no vertex to add as a terminal")))?;
    is_terminal[b as usize] = true;
    added.push(b);
    Ok(())
}

/// Cuts the surface through the high-energy vertices and flattens the charts (SPEC §6).
///
/// 1. Terminals: vertices with `e_v > epsilon` or `|Ω_v| > defect_tolerance` (D30), plus
///    `extra_terminals`.
/// 2. Edge weights ([`CutOptions::length_weight`], [`CutOptions::alignment_weight`]).
/// 3. [`ddg_bff::plan_cut`]; a closed genus-0 component with fewer than two terminals gets its
///    highest-energy vertices added (reported in `added_terminals`).
/// 4. Cut, split into charts (all disks).
/// 5. Per chart, BFF with zero boundary scale factors.
/// 6. Per-face and per-vertex log scale factor.
///
/// # Errors
/// Invalid options, wrong lengths, non-finite positions, degenerate faces, undefined energies,
/// ids out of range, and everything the cut planner and BFF report ([`DevError::Bff`]).
pub fn cut_and_flatten(
    mesh: &TriMesh,
    positions: &[[f64; 3]],
    opts: &CutOptions,
) -> Result<DevelopableFlattening, DevError> {
    opts.validate()?;
    let nv = mesh.n_vertices();
    if positions.len() != nv {
        return Err(DevError::LengthMismatch { what: "positions", expected: nv, got: positions.len() });
    }
    if let Some(v) = positions.iter().position(|p| !v3::is_finite(*p)) {
        return Err(DevError::NonFinitePosition { vertex: v as u32 });
    }
    if let Some(&t) = opts.extra_terminals.iter().find(|&&t| t as usize >= nv) {
        return Err(DevError::VertexOutOfRange { what: "extra_terminals", vertex: t, n_vertices: nv });
    }
    if let Some(m) = &opts.may_cut
        && m.len() != mesh.n_edges()
    {
        return Err(DevError::LengthMismatch { what: "may_cut", expected: mesh.n_edges(), got: m.len() });
    }

    // 1. Terminals.
    let eopts = EnergyOptions {
        energy: opts.energy,
        include_boundary: true,
        omit_valence_three: false,
        threads: opts.threads,
    };
    let (_, vertex_energy) = evaluate_energy(mesh, positions, &eopts)?;
    let geometry = Geometry::from_positions(mesh, positions)?;
    let vertex_defect = ddg_core::operators::angle_defects(mesh, &geometry)?;
    let mut is_terminal: Vec<bool> = vertex_energy.iter().map(|&e| e > opts.epsilon).collect();
    if let Some(tol) = opts.defect_tolerance {
        for (t, &d) in is_terminal.iter_mut().zip(&vertex_defect) {
            *t |= d.abs() > tol;
        }
    }
    for &t in &opts.extra_terminals {
        is_terminal[t as usize] = true;
    }
    let terminals: Vec<u32> = (0..nv as u32).filter(|&v| is_terminal[v as usize]).collect();

    // 2. Edge weights.
    let weights = edge_weights(mesh, positions, geometry.edge_lengths(), opts)?;

    // 3. Closed genus-0 components need two terminals.
    let vertex_component = mesh.vertex_components();
    let mut added: Vec<u32> = Vec::new();
    for (c, topo) in component_topology(mesh).iter().enumerate() {
        if topo.is_closed() && topo.genus == 0 {
            let have = (0..nv).filter(|&v| is_terminal[v] && vertex_component[v] == c as u32).count();
            for _ in have..2 {
                add_terminal(c as u32, &vertex_component, &vertex_energy, &mut is_terminal, &mut added)?;
            }
        }
    }
    let mut retries = 0;
    let plan = loop {
        let cones: Vec<u32> = (0..nv as u32).filter(|&v| is_terminal[v as usize]).collect();
        match plan_cut(&CutPlanInput {
            mesh,
            edge_lengths: &weights,
            cone_vertices: &cones,
            may_cut: opts.may_cut.as_deref(),
            must_cut: &opts.must_cut,
        }) {
            Ok(p) => break p,
            // a closed piece after the required cuts that the component test did not foresee
            Err(BffError::SphereNeedsTwoCones { first_face, .. }) if retries < 8 => {
                retries += 1;
                let c = mesh.face_components()[first_face as usize];
                add_terminal(c, &vertex_component, &vertex_energy, &mut is_terminal, &mut added)?;
            }
            Err(e) => return Err(e.into()),
        }
    };
    added.sort_unstable();

    // 4. Cut and split.
    let cut = cut_mesh(mesh, &plan.cut_edges)?;
    let charts = cut.mesh.split_components();
    for (i, sm) in charts.iter().enumerate() {
        if sm.mesh.euler_characteristic() != 1 || sm.mesh.boundary_loops().len() != 1 {
            return Err(DevError::Internal(format!(
                "chart {i} (face {}) is not a disk after cutting: Euler characteristic {}, {} boundary loops",
                sm.face_to_parent[0],
                sm.mesh.euler_characteristic(),
                sm.mesh.boundary_loops().len()
            )));
        }
    }

    // 5. Flatten every chart.
    struct Flat {
        uv: Vec<[f64; 2]>,
        report: ChartReport,
        face_area: Vec<f64>,
    }
    let flats: Vec<Flat> = charts
        .par_iter()
        .map(|sm| -> Result<Flat, DevError> {
            let pos: Vec<[f64; 3]> =
                sm.vertex_to_parent.iter().map(|&v| positions[cut.vertex_to_original[v as usize] as usize]).collect();
            let g = Geometry::from_positions(&sm.mesh, &pos)?;
            let d = DiskOptions { overlap_samples_per_triangle: opts.overlap_samples_per_triangle, ..Default::default() };
            let mut bff = DiskBff::new(&sm.mesh, &g, d)?;
            let zeros = vec![0.0; bff.boundary().len()];
            let flat = bff.flatten_with_scale_factors(&zeros, opts.extension)?;
            let mut report = flat.report;
            report.first_face = sm.face_to_parent[0];
            Ok(Flat { uv: flat.uv, report, face_area: g.face_areas().to_vec() })
        })
        .collect::<Result<_, _>>()?;

    // 6. Assemble.
    let nf = mesh.n_faces();
    let mut corner_uvs = vec![[0.0; 2]; 3 * nf];
    let mut face_chart = vec![0u32; nf];
    let mut face_log_scale = vec![f64::NAN; nf];
    let mut face_area = vec![0.0; nf];
    for (ci, (sm, fl)) in charts.iter().zip(&flats).enumerate() {
        for (fl_i, tri) in sm.mesh.faces().iter().enumerate() {
            let f = sm.face_to_parent[fl_i] as usize;
            face_chart[f] = ci as u32;
            let uv = tri.map(|v| fl.uv[v as usize]);
            for k in 0..3 {
                corner_uvs[3 * f + k] = uv[k];
            }
            let plane = 0.5 * signed_double_area(uv[0], uv[1], uv[2]).abs();
            face_area[f] = fl.face_area[fl_i];
            if plane > 0.0 {
                face_log_scale[f] = log_scale_factor(fl.face_area[fl_i], plane);
            }
        }
    }
    let mut num = vec![0.0; nv];
    let mut den = vec![0.0; nv];
    let (mut abs_sum, mut area_sum) = (0.0, 0.0);
    for (f, tri) in mesh.faces().iter().enumerate() {
        let u = face_log_scale[f];
        if u.is_finite() {
            for &v in tri {
                num[v as usize] += face_area[f] * u;
                den[v as usize] += face_area[f];
            }
            abs_sum += face_area[f] * u.abs();
            area_sum += face_area[f];
        }
    }
    let vertex_log_scale: Vec<f64> =
        (0..nv).map(|v| if den[v] > 0.0 { num[v] / den[v] } else { f64::NAN }).collect();
    let mean_abs_log_scale = if area_sum > 0.0 { abs_sum / area_sum } else { f64::NAN };

    Ok(DevelopableFlattening {
        cut_vertex_pairs: plan.cut_edges.iter().map(|&e| mesh.edges()[e as usize]).collect(),
        cut_edges: plan.cut_edges.clone(),
        corner_uvs,
        face_chart,
        face_log_scale,
        vertex_log_scale,
        mean_abs_log_scale,
        terminals,
        added_terminals: added,
        plan,
        charts: flats.iter().map(|f| f.report).collect(),
        vertex_energy,
        vertex_defect,
        star_embeddedness: star_embeddedness(mesh, positions)?,
    })
}
