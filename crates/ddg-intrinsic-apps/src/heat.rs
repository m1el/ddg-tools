//! Geodesic distance by the heat method on an intrinsic triangulation (paper Sec. 5.3,
//! Fig. 18).
//!
//! The heat method of Crane, Weischedel and Wardetzky (2013; 2017): (I) heat flow from the
//! sources for a short time `t` by one backward Euler step, `(M + t L) u = δ`; (II) the unit
//! field `X = −∇u / |∇u|` in every triangle; (III) the potential of `X`, `L φ = −div X` (with
//! the workspace's positive semidefinite `L`; the paper writes `Δφ = ∇·X`), shifted to vanish at
//! the sources. `M` is the lumped mass, `L` the cotan Laplacian, `δ` the indicator of the
//! sources (spread over the corners of the triangle a source point lies in by its barycentric
//! coordinates). The time step is `t = m h²` with `h` the mean edge length of the triangulation
//! and `m = 1` by default (the paper's estimate, Sec. 3.2 of the 2017 version).
//!
//! On an intrinsic Delaunay triangulation the cotan weights are not negative and the method
//! behaves (paper Sec. 5.3: "by just using the iDT the error is dramatically reduced; further
//! refinement with iDR reduces the error to less than 1%"; measured in this crate's README).
//!
//! [`HeatSolver::new`] factors both systems once; [`HeatSolver::distance`] then costs three
//! back substitutions and two sparse passes per call, for any number of sources. Values come
//! out at every vertex of the triangulation (the input vertices first, with their ids); the
//! place of an inserted vertex on the input mesh is
//! [`IntrinsicTriangulation::vertex_location`], and [`HeatSolver::value_at`] interpolates the
//! result at any point of the input mesh through a point query (paper Alg. 11).
//!
//! # Boundary
//!
//! Neumann conditions (the natural ones, default), zero Dirichlet conditions for the heat
//! step, or the mean of the two heat solutions, which the 2017 paper advocates for smoothed
//! distances (its Sec. 3.4); the Poisson step always has natural conditions. A source with
//! weight at a boundary vertex is refused under the Dirichlet conditions, which would hold its
//! heat at zero.
//!
//! # References
//!
//! K. Crane, C. Weischedel, M. Wardetzky. *Geodesics in heat: A new approach to computing
//! distance based on heat flow.* ACM Transactions on Graphics 32(5), Article 152, 2013.
//! <https://doi.org/10.1145/2516971.2516977>
//!
//! K. Crane, C. Weischedel, M. Wardetzky. *The heat method for distance computation.*
//! Communications of the ACM 60(11), 90–99, 2017. <https://doi.org/10.1145/3131280>
//!
//! The source term and the shift that makes the result vanish at the sources follow
//! geometry-central's `heat_method_distance.cpp` (© 2017–2019 Nicholas Sharp and the
//! geometry-central contributors, MIT; read, not translated; notice in `LICENSES/`).

use std::fmt;
use std::time::Instant;

use crate::representation::InputPoints;
use ddg_core::SparseError;
use ddg_core::geometry::corner_cotan;
use ddg_core::sparse::{Factorization, Method, Parallelism, Symbolic};
use ddg_intrinsic::formulas::{corner_distances, triangle_layout};
use ddg_intrinsic::{IntrinsicError, IntrinsicTriangulation, SurfacePoint};

use crate::fem::{FemError, PinnedLaplacian};

/// Boundary conditions of the heat step. See the [module documentation](self).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum HeatBoundary {
    /// Zero flux (natural conditions).
    #[default]
    Neumann,
    /// Zero heat at the boundary vertices.
    Dirichlet,
    /// The mean of the Neumann and the Dirichlet solutions.
    Average,
}

/// Options of [`HeatSolver::new`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HeatOptions {
    /// `m` in the time step `t = m h²`, `h` the mean edge length of the triangulation.
    /// Default 1 (the paper's choice).
    pub time_factor: f64,
    /// The time step itself; overrides `time_factor` when given.
    pub time: Option<f64>,
    /// Boundary conditions of the heat step. Default Neumann.
    pub boundary: HeatBoundary,
}

impl Default for HeatOptions {
    fn default() -> Self {
        HeatOptions {
            time_factor: 1.0,
            time: None,
            boundary: HeatBoundary::Neumann,
        }
    }
}

/// A source of [`HeatSolver::distance`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum HeatSource {
    /// A vertex of the triangulation (an input vertex has its input id).
    Vertex(usize),
    /// A point of the input mesh, found on the triangulation by a point query.
    Point(SurfacePoint),
}

/// Why the heat method failed.
#[derive(Debug, Clone, PartialEq)]
pub enum HeatError {
    /// An option is out of range.
    InvalidOptions(String),
    /// No source was given.
    NoSource,
    /// A source vertex is not a vertex of the triangulation.
    SourceOutOfRange {
        /// The vertex.
        vertex: usize,
        /// The number of vertices.
        n_vertices: usize,
    },
    /// A source has weight at a boundary vertex, where the Dirichlet conditions hold the heat
    /// at zero.
    SourceOnDirichletBoundary {
        /// The boundary vertex.
        vertex: usize,
    },
    /// A point query or a point failed.
    Intrinsic(IntrinsicError),
    /// A factorization or a solve failed.
    Sparse(SparseError),
    /// The Poisson solve failed.
    Fem(FemError),
    /// The values given to [`HeatSolver::value_at`] are not one per vertex.
    ValueCount {
        /// Values given.
        got: usize,
        /// Vertices.
        expected: usize,
    },
}

impl fmt::Display for HeatError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            HeatError::InvalidOptions(m) => write!(f, "invalid heat options: {m}"),
            HeatError::NoSource => write!(f, "no source was given"),
            HeatError::SourceOutOfRange { vertex, n_vertices } => write!(
                f,
                "source vertex {vertex} does not exist (the triangulation has {n_vertices} vertices)"
            ),
            HeatError::SourceOnDirichletBoundary { vertex } => write!(
                f,
                "a source has weight at boundary vertex {vertex}, where the Dirichlet conditions \
                 hold the heat at zero; use Neumann conditions"
            ),
            HeatError::Intrinsic(e) => write!(f, "{e}"),
            HeatError::Sparse(e) => write!(f, "{e}"),
            HeatError::Fem(e) => write!(f, "{e}"),
            HeatError::ValueCount { got, expected } => {
                write!(f, "{got} values given for {expected} vertices")
            }
        }
    }
}

impl std::error::Error for HeatError {}

impl From<IntrinsicError> for HeatError {
    fn from(e: IntrinsicError) -> Self {
        HeatError::Intrinsic(e)
    }
}
impl From<SparseError> for HeatError {
    fn from(e: SparseError) -> Self {
        HeatError::Sparse(e)
    }
}
impl From<FemError> for HeatError {
    fn from(e: FemError) -> Self {
        HeatError::Fem(e)
    }
}

/// What [`HeatSolver::new`] set up.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HeatSetup {
    /// The time step `t`.
    pub time: f64,
    /// The mean edge length `h` of the triangulation.
    pub mean_edge_length: f64,
    /// Vertices of the triangulation.
    pub n_vertices: usize,
    /// Faces of the triangulation.
    pub n_faces: usize,
    /// Wall-clock seconds to assemble and factor.
    pub seconds: f64,
}

/// What one call of [`HeatSolver::distance`] did.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HeatReport {
    /// Number of sources.
    pub n_sources: usize,
    /// Faces in which the heat had no gradient (it underflowed to a constant), where `X` was
    /// set to zero.
    pub n_flat_faces: usize,
    /// The largest sum of the divergence over a component; zero up to rounding.
    pub inconsistency: f64,
    /// The constant subtracted so that the distance vanishes at the sources.
    pub shift: f64,
    /// Wall-clock seconds of the call.
    pub seconds: f64,
}

/// The distance and the report of one call.
#[derive(Debug, Clone, PartialEq)]
pub struct HeatDistance {
    /// The distance at every vertex of the triangulation (input vertices first).
    pub distances: Vec<f64>,
    /// What the call did.
    pub report: HeatReport,
}

/// The heat method on one intrinsic triangulation, with both systems factored.
pub struct HeatSolver<T: InputPoints + Clone = IntrinsicTriangulation> {
    tri: T,
    boundary: HeatBoundary,
    neumann: Option<Factorization>,
    /// The Dirichlet system on the interior vertices, and the row of every vertex
    /// (`u32::MAX` on the boundary).
    dirichlet: Option<(Factorization, Vec<u32>)>,
    poisson: PinnedLaplacian,
    setup: HeatSetup,
}

impl<T: InputPoints + Clone> HeatSolver<T> {
    /// Assembles and factors `M + t L` (for the boundary conditions asked for) and `L`.
    ///
    /// # Errors
    /// [`HeatError::InvalidOptions`] for a time or time factor that is not a positive number,
    /// [`HeatError::Sparse`], [`HeatError::Fem`] if a factorization fails.
    pub fn new(tri: &T, options: &HeatOptions) -> Result<HeatSolver<T>, HeatError> {
        let start = Instant::now();
        let mean = tri.edge_lengths().iter().sum::<f64>() / tri.complex().n_edges().max(1) as f64;
        let time = match options.time {
            Some(t) => t,
            None => {
                if !(options.time_factor > 0.0 && options.time_factor.is_finite()) {
                    return Err(HeatError::InvalidOptions(format!(
                        "time_factor must be a positive number, not {}",
                        options.time_factor
                    )));
                }
                options.time_factor * mean * mean
            }
        };
        if !(time > 0.0 && time.is_finite()) {
            return Err(HeatError::InvalidOptions(format!(
                "the time step must be a positive number, not {time}"
            )));
        }
        let n = tri.complex().n_vertices();
        let l = ddg_intrinsic::generic_ops::cotan_laplacian(tri);
        let mass = ddg_intrinsic::generic_ops::lumped_mass(tri);
        let mut heat = l.clone();
        heat.values_mut().iter_mut().for_each(|x| *x *= time);
        heat.add_to_diagonal(&mass)?;
        let neumann = if options.boundary != HeatBoundary::Dirichlet {
            Some(Symbolic::analyze(&heat, Method::Llt)?.factor(&heat, Parallelism::Sequential)?)
        } else {
            None
        };
        let dirichlet = if options.boundary != HeatBoundary::Neumann {
            let mut row = vec![u32::MAX; n];
            let mut keep = Vec::with_capacity(n);
            for (v, r) in row.iter_mut().enumerate() {
                if !tri.complex().is_boundary_vertex(v) {
                    *r = keep.len() as u32;
                    keep.push(v as u32);
                }
            }
            let a = heat.submatrix(&keep, &keep)?;
            let factor = Symbolic::analyze(&a, Method::Llt)?.factor(&a, Parallelism::Sequential)?;
            Some((factor, row))
        } else {
            None
        };
        let poisson = PinnedLaplacian::new(&l)?;
        Ok(HeatSolver {
            tri: tri.clone(),
            boundary: options.boundary,
            neumann,
            dirichlet,
            poisson,
            setup: HeatSetup {
                time,
                mean_edge_length: mean,
                n_vertices: n,
                n_faces: tri.complex().n_faces(),
                seconds: start.elapsed().as_secs_f64(),
            },
        })
    }

    /// The triangulation the solver works on.
    pub fn triangulation(&self) -> &T {
        &self.tri
    }

    /// The time step, the mean edge length and the time it took to factor.
    pub fn setup(&self) -> HeatSetup {
        self.setup
    }

    /// A source as a face of the triangulation and barycentric coordinates in it.
    fn locate(&self, source: &HeatSource) -> Result<(usize, [f64; 3]), HeatError> {
        let point = match *source {
            HeatSource::Vertex(v) => {
                if v >= self.tri.complex().n_vertices() {
                    return Err(HeatError::SourceOutOfRange {
                        vertex: v,
                        n_vertices: self.tri.complex().n_vertices(),
                    });
                }
                SurfacePoint::Vertex(v as u32)
            }
            HeatSource::Point(p) => self.tri.locate_on_triangulation(p)?,
        };
        Ok(match point {
            SurfacePoint::Vertex(v) => {
                let h = self.tri.complex().outgoing_halfedge(v as usize);
                let mut bary = [0.0; 3];
                bary[h % 3] = 1.0;
                (h / 3, bary)
            }
            SurfacePoint::Edge { edge, t } => {
                let h = self.tri.complex().edge_halfedge(edge as usize);
                let mut bary = [0.0; 3];
                bary[h % 3] = 1.0 - t;
                bary[(h + 1) % 3] = t;
                (h / 3, bary)
            }
            SurfacePoint::Face { face, bary } => (face as usize, bary),
        })
    }

    /// The geodesic distance from the nearest of `sources` at every vertex of the
    /// triangulation. With several sources the shift makes the weighted mean of the values at
    /// the sources zero, as the reference does.
    ///
    /// # Errors
    /// [`HeatError::NoSource`], [`HeatError::SourceOutOfRange`],
    /// [`HeatError::SourceOnDirichletBoundary`], [`HeatError::Intrinsic`] when a point query
    /// fails, and the errors of the solves.
    pub fn distance(&self, sources: &[HeatSource]) -> Result<HeatDistance, HeatError> {
        let start = Instant::now();
        if sources.is_empty() {
            return Err(HeatError::NoSource);
        }
        let tri = &self.tri;
        let n = tri.complex().n_vertices();
        let located: Vec<(usize, [f64; 3])> = sources
            .iter()
            .map(|s| self.locate(s))
            .collect::<Result<_, _>>()?;
        let mut delta = vec![0.0; n];
        for &(f, bary) in &located {
            for k in 0..3 {
                delta[tri.complex().faces()[f][k] as usize] += bary[k];
            }
        }
        let u = self.heat(&delta)?;
        let (div, n_flat) = self.divergence(&u);
        let inconsistency = self.poisson.inconsistency(&div);
        let mut phi: Vec<f64> = div.iter().map(|x| -x).collect();
        self.poisson.solve_in_place(&mut phi)?;
        // Shift: the values at the corners of every source face, minus the distances of the
        // source from those corners, weighted by the barycentric coordinates.
        let mut sum = 0.0;
        let mut weight = 0.0;
        for &(f, bary) in &located {
            let d = corner_distances(tri.face_lengths(f), bary)
                .expect("faces of a triangulation are triangles");
            for k in 0..3 {
                sum += bary[k] * (phi[tri.complex().faces()[f][k] as usize] - d[k]);
                weight += bary[k];
            }
        }
        let shift = sum / weight;
        phi.iter_mut().for_each(|x| *x -= shift);
        Ok(HeatDistance {
            distances: phi,
            report: HeatReport {
                n_sources: sources.len(),
                n_flat_faces: n_flat,
                inconsistency,
                shift,
                seconds: start.elapsed().as_secs_f64(),
            },
        })
    }

    /// Step I: the heat after one backward Euler step from `delta`.
    fn heat(&self, delta: &[f64]) -> Result<Vec<f64>, HeatError> {
        let n = delta.len();
        let neumann = match &self.neumann {
            Some(factor) => {
                let mut u = delta.to_vec();
                factor.solve_in_place(&mut u, Parallelism::Sequential)?;
                Some(u)
            }
            None => None,
        };
        let dirichlet = match &self.dirichlet {
            Some((factor, row)) => {
                if let Some(v) = (0..n).find(|&v| row[v] == u32::MAX && delta[v] != 0.0) {
                    return Err(HeatError::SourceOnDirichletBoundary { vertex: v });
                }
                let mut rhs: Vec<f64> = (0..n)
                    .filter(|&v| row[v] != u32::MAX)
                    .map(|v| delta[v])
                    .collect();
                factor.solve_in_place(&mut rhs, Parallelism::Sequential)?;
                Some(
                    (0..n)
                        .map(|v| {
                            if row[v] == u32::MAX {
                                0.0
                            } else {
                                rhs[row[v] as usize]
                            }
                        })
                        .collect::<Vec<f64>>(),
                )
            }
            None => None,
        };
        Ok(match (self.boundary, neumann, dirichlet) {
            (HeatBoundary::Average, Some(a), Some(b)) => {
                a.iter().zip(&b).map(|(x, y)| 0.5 * (x + y)).collect()
            }
            (_, Some(a), None) => a,
            (_, None, Some(b)) => b,
            _ => unreachable!("the factorizations match the boundary conditions"),
        })
    }

    /// Steps II and III's right-hand side: the integrated divergence of `X = −∇u/|∇u|` at every
    /// vertex, and the number of faces where `∇u` vanished.
    ///
    /// [formula: app-heat-gradient]
    /// [formula: app-heat-divergence]
    fn divergence(&self, u: &[f64]) -> (Vec<f64>, usize) {
        let tri = &self.tri;
        let mut div = vec![0.0; tri.complex().n_vertices()];
        let mut n_flat = 0;
        for (f, face) in tri.complex().faces().iter().enumerate() {
            let l = tri.face_lengths(f);
            let p = triangle_layout(l).expect("faces of a triangulation are triangles");
            let x = unit_descent(&p, [0, 1, 2].map(|k| u[face[k] as usize]));
            let Some(x) = x else {
                n_flat += 1;
                continue;
            };
            for k in 0..3 {
                let (k1, k2) = ((k + 1) % 3, (k + 2) % 3);
                // Edges from corner k; the angle opposite the edge k -> k1 is at k2.
                let e1 = [p[k1][0] - p[k][0], p[k1][1] - p[k][1]];
                let e2 = [p[k2][0] - p[k][0], p[k2][1] - p[k][1]];
                // l[k] is the side k -> k+1, l[k1] the side k1 -> k2, l[k2] the side k2 -> k.
                let cot_at_k2 = corner_cotan(l[k], l[k1], l[k2]).expect("a triangle");
                let cot_at_k1 = corner_cotan(l[k2], l[k], l[k1]).expect("a triangle");
                div[face[k] as usize] += 0.5
                    * (cot_at_k2 * (e1[0] * x[0] + e1[1] * x[1])
                        + cot_at_k1 * (e2[0] * x[0] + e2[1] * x[1]));
            }
        }
        (div, n_flat)
    }

    /// The value of the vertex function `values` (one per vertex of the triangulation, e.g. a
    /// result of [`distance`](Self::distance)) at a point of the input mesh: the point is found
    /// on the triangulation (paper Alg. 11) and the values of the corners of its face are
    /// interpolated linearly.
    ///
    /// # Errors
    /// [`HeatError::ValueCount`], and [`HeatError::Intrinsic`] when the point query fails.
    pub fn value_at(&self, values: &[f64], point: SurfacePoint) -> Result<f64, HeatError> {
        if values.len() != self.tri.complex().n_vertices() {
            return Err(HeatError::ValueCount {
                got: values.len(),
                expected: self.tri.complex().n_vertices(),
            });
        }
        let (f, bary) = self.locate(&HeatSource::Point(point))?;
        Ok((0..3)
            .map(|k| bary[k] * values[self.tri.complex().faces()[f][k] as usize])
            .sum())
    }
}

/// `−∇u / |∇u|` in the triangle with corners `p` (counter-clockwise) and values `u` at them, or
/// `None` if the gradient vanishes or is not finite. `∇u = Σ_k u_k J e_k / (2A)` with `e_k` the
/// side opposite corner `k` (counter-clockwise) and `J` the rotation by a quarter turn; the
/// factor `1/(2A)` does not change the direction and is left out.
///
/// [formula: app-heat-gradient]
fn unit_descent(p: &[[f64; 2]; 3], u: [f64; 3]) -> Option<[f64; 2]> {
    let mut g = [0.0; 2];
    for k in 0..3 {
        let (a, b) = (p[(k + 1) % 3], p[(k + 2) % 3]);
        let e = [b[0] - a[0], b[1] - a[1]];
        g[0] += u[k] * -e[1];
        g[1] += u[k] * e[0];
    }
    let norm = g[0].hypot(g[1]);
    if norm > 0.0 && norm.is_finite() {
        Some([-g[0] / norm, -g[1] / norm])
    } else {
        None
    }
}

/// The geodesic distance from `sources` by the heat method on `tri`: [`HeatSolver::new`]
/// followed by one [`HeatSolver::distance`]. Build the solver once to reuse the factorizations
/// for further sources.
///
/// # Errors
/// Those of [`HeatSolver::new`] and [`HeatSolver::distance`].
pub fn geodesic_distance<T: InputPoints + Clone>(
    tri: &T,
    sources: &[HeatSource],
    options: &HeatOptions,
) -> Result<HeatDistance, HeatError> {
    HeatSolver::new(tri, options)?.distance(sources)
}
