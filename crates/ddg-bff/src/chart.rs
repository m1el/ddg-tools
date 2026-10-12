//! Machinery shared by the disk flattenings and the cone flattening: Poisson solves on a mesh
//! whose corner data is given per face, and the quality report of a layout.

use std::f64::consts::PI;

use ddg_core::geometry::Geometry;
use ddg_core::sparse::Factorization;
use ddg_core::topology::{NONE, TriMesh};

use crate::curve::ClosureReport;
use crate::error::BffError;
use crate::laplace::{factor_spd, laplacian_apply, laplacian_block, relative_norm, solve_spd};
use crate::quality::{
    Coverage, DistortionAccumulator, DistortionReport, coverage, polygon_self_crossings,
    signed_corner_angle,
};

/// A mesh together with intrinsic data per face and corner.
///
/// For a cut mesh the per-face arrays are those of the mesh it was cut from: faces and corners
/// keep their ids under cutting.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Surface<'a> {
    pub mesh: &'a TriMesh,
    /// Corner cotangents, `cotans[f][k]` at corner `3 f + k`.
    pub cotans: &'a [[f64; 3]],
    /// Corner angles.
    pub angles: &'a [[f64; 3]],
    /// Face areas.
    pub areas: &'a [f64],
    /// Length of every halfedge `3 f + k`.
    pub halfedge_lengths: &'a [f64],
}

/// The length of every halfedge of `mesh` (which must be the mesh `geometry` was built for).
pub(crate) fn halfedge_lengths(mesh: &TriMesh, geometry: &Geometry) -> Vec<f64> {
    let lengths = geometry.edge_lengths();
    (0..mesh.n_halfedges())
        .map(|h| lengths[mesh.edge(h)])
        .collect()
}

impl Surface<'_> {
    /// `out = L x` with the cotan Laplacian of this surface.
    pub(crate) fn apply_laplacian(&self, x: &[f64], out: &mut [f64]) {
        laplacian_apply(self.mesh.faces(), self.cotans, x, out);
    }

    /// Per vertex, the sum of the corner angles at it.
    pub(crate) fn angle_sums(&self) -> Vec<f64> {
        let mut sums = vec![0.0; self.mesh.n_vertices()];
        for (face, angles) in self.mesh.faces().iter().zip(self.angles) {
            for k in 0..3 {
                sums[face[k] as usize] += angles[k];
            }
        }
        sums
    }

    /// The lumped vertex masses `m_i = ⅓ Σ A_f`.
    ///
    /// [formula: lumped-mass]
    pub(crate) fn masses(&self) -> Vec<f64> {
        let mut mass = vec![0.0; self.mesh.n_vertices()];
        for (face, &area) in self.mesh.faces().iter().zip(self.areas) {
            for &v in face {
                mass[v as usize] += area / 3.0;
            }
        }
        mass
    }
}

/// The Poisson problem `L_FF a_F = φ_F − L_FX a_X` on a set `F` of free vertices, with the
/// values `a_X` on the other vertices given: factorized once, solved for many right-hand sides.
#[derive(Debug)]
pub(crate) struct PoissonSystem {
    /// Per vertex: its row, or [`NONE`] if its value is prescribed.
    index: Vec<u32>,
    n_free: usize,
    factorization: Option<Factorization>,
    name: &'static str,
}

impl PoissonSystem {
    /// Factorizes `L_FF + shift · I` where `F` is the set of vertices with `fixed[v] == false`.
    pub(crate) fn new(
        name: &'static str,
        surface: &Surface<'_>,
        fixed: &[bool],
        shift: f64,
    ) -> Result<Self, BffError> {
        let mut index = vec![NONE; fixed.len()];
        let mut n_free = 0usize;
        for (v, _) in fixed.iter().enumerate().filter(|&(_, &f)| !f) {
            index[v] = n_free as u32;
            n_free += 1;
        }
        let factorization = if n_free == 0 {
            None
        } else {
            let block =
                laplacian_block(surface.mesh.faces(), surface.cotans, &index, n_free, shift);
            Some(factor_spd(name, &block)?)
        };
        Ok(PoissonSystem {
            index,
            n_free,
            factorization,
            name,
        })
    }

    /// Number of unknowns.
    pub(crate) fn n_free(&self) -> usize {
        self.n_free
    }

    /// Solves for `nrhs` columns stored one after another in `values`, each with one entry
    /// per vertex. On entry a column holds the prescribed value at fixed vertices and the
    /// source `φ` at free vertices; on exit it holds the solution at every vertex.
    ///
    /// Returns the largest relative residual `‖L_FF a_F − rhs‖ / ‖rhs‖` over the columns,
    /// evaluated with the unshifted Laplacian.
    ///
    /// [formula: bff-harmonic-extension]
    pub(crate) fn solve(
        &self,
        surface: &Surface<'_>,
        values: &mut [f64],
        nrhs: usize,
    ) -> Result<f64, BffError> {
        let n = self.index.len();
        debug_assert_eq!(values.len(), n * nrhs);
        let Some(factorization) = &self.factorization else {
            return Ok(0.0);
        };
        let mut rhs = vec![0.0; self.n_free * nrhs];
        let mut lifted = vec![0.0; n];
        let mut product = vec![0.0; n];
        for c in 0..nrhs {
            let column = &values[c * n..(c + 1) * n];
            for v in 0..n {
                lifted[v] = if self.index[v] == NONE {
                    column[v]
                } else {
                    0.0
                };
            }
            surface.apply_laplacian(&lifted, &mut product);
            for v in 0..n {
                if self.index[v] != NONE {
                    rhs[c * self.n_free + self.index[v] as usize] = column[v] - product[v];
                }
            }
        }
        let b = rhs.clone();
        solve_spd(self.name, factorization, &mut rhs, nrhs)?;
        let mut worst = 0.0f64;
        for c in 0..nrhs {
            let column = &mut values[c * n..(c + 1) * n];
            for v in 0..n {
                lifted[v] = if self.index[v] == NONE {
                    0.0
                } else {
                    rhs[c * self.n_free + self.index[v] as usize]
                };
            }
            surface.apply_laplacian(&lifted, &mut product);
            let b_column = &b[c * self.n_free..(c + 1) * self.n_free];
            let residual = relative_norm(
                (0..n)
                    .filter(|&v| self.index[v] != NONE)
                    .map(|v| product[v] - b_column[self.index[v] as usize]),
                b_column.iter().copied(),
            );
            worst = worst.max(residual);
            for v in 0..n {
                if self.index[v] != NONE {
                    column[v] = lifted[v];
                }
            }
        }
        Ok(worst)
    }
}

/// A boundary vertex counts as having more than a full turn around it if the angles of the
/// layout at it sum to more than `2π` plus this (radians). A flat vertex at the end of a cut
/// has exactly `2π`, up to the rounding of the solve.
pub(crate) const FULL_TURN_TOLERANCE: f64 = 1e-9;

/// Overlap of one chart's layout with itself: two exact signs of it, and an estimate of how
/// much area is concerned.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OverlapReport {
    /// Number of boundary vertices of the chart around which the layout covers more than a
    /// full turn: the signed corner angles of the layout at the vertex sum to more than `2π`
    /// (by more than `1e-9`).
    /// The triangles around such a vertex overlap each other. This is what happens at a cone
    /// of negative angle that is the end of a cut, or through which the cut passes unevenly.
    pub wedges_over_2pi: usize,
    /// Number of pairs of non-adjacent boundary edges of the layout that cross each other:
    /// two parts of the chart lie on top of each other.
    pub boundary_crossings: usize,
    /// The area of the plane that the chart covers, and the part of it that is covered more
    /// than once, estimated by point sampling ([`crate::quality::coverage`]). `None` if the
    /// estimate was switched off (`overlap_samples_per_triangle == 0`) or the chart has no
    /// area in the plane.
    pub coverage: Option<Coverage>,
}

impl OverlapReport {
    /// The estimated fraction of the chart's footprint that is covered more than once, or
    /// `None` if it was not estimated.
    pub fn overlapped_fraction(&self) -> Option<f64> {
        self.coverage.map(|c| c.overlapped_fraction())
    }
}

/// What is reported about the layout of one chart.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChartReport {
    /// Lowest-numbered face of the chart.
    pub first_face: u32,
    /// Number of vertices of the chart (after cutting).
    pub n_vertices: usize,
    /// Number of vertices on its boundary.
    pub n_boundary_vertices: usize,
    /// How its boundary curve was closed.
    pub closure: ClosureReport,
    /// Angle and area distortion, flipped and degenerate triangles.
    pub distortion: DistortionReport,
    /// Overlap of the layout with itself.
    pub overlap: OverlapReport,
}

impl ChartReport {
    /// Whether the layout of the chart is free of overlap: no flipped or degenerate triangle,
    /// no boundary vertex with more than a full turn around it, and a boundary polygon that
    /// does not cross itself.
    ///
    /// Why this is enough: without flipped or degenerate triangles the map preserves
    /// orientation, so the number of times it covers a point is the winding number of the
    /// boundary polygon around that point; a polygon without crossings and without a vertex
    /// that turns by more than a full turn winds at most once. The boundary may still *touch*
    /// itself — the two sides of a cut that ends at a flat vertex coincide — and coincident
    /// boundary pieces that are traversed in the same direction are not detected; the sampled
    /// [`OverlapReport::coverage`] is an independent check for that.
    ///
    /// [formula: disk-injectivity-criterion]
    pub fn is_injective(&self) -> bool {
        self.distortion.flipped_triangles == 0
            && self.distortion.degenerate_triangles == 0
            && self.overlap.wedges_over_2pi == 0
            && self.overlap.boundary_crossings == 0
    }
}

/// What [`layout_quality`] measures.
pub(crate) struct LayoutQuality {
    /// Distortion and overlap of every component, in component order.
    pub charts: Vec<(DistortionReport, OverlapReport)>,
    /// Per vertex, the sum of the signed corner angles of the layout at that vertex.
    pub angle_sums: Vec<f64>,
}

/// The quality of a layout `uv` (one point per vertex) of a surface whose components are
/// disks: one report per component, in component order, without the closure part. Also
/// returns, per vertex, the sum of the signed corner angles of the layout at that vertex.
pub(crate) fn layout_quality(
    surface: &Surface<'_>,
    uv: &[[f64; 2]],
    overlap_samples_per_triangle: f64,
) -> Result<LayoutQuality, BffError> {
    let mesh = surface.mesh;
    let component = mesh.face_components();
    let mut distortion: Vec<DistortionAccumulator> = (0..mesh.n_components())
        .map(|_| DistortionAccumulator::new())
        .collect();
    let mut angle_sums = vec![0.0; mesh.n_vertices()];
    for (f, face) in mesh.faces().iter().enumerate() {
        let q = [
            uv[face[0] as usize],
            uv[face[1] as usize],
            uv[face[2] as usize],
        ];
        let lengths = [
            surface.halfedge_lengths[3 * f],
            surface.halfedge_lengths[3 * f + 1],
            surface.halfedge_lengths[3 * f + 2],
        ];
        distortion[component[f] as usize].add(lengths, surface.areas[f], q);
        for k in 0..3 {
            angle_sums[face[k] as usize] +=
                signed_corner_angle(q[k], q[(k + 1) % 3], q[(k + 2) % 3]);
        }
    }
    let mut overlap = vec![
        OverlapReport {
            wedges_over_2pi: 0,
            boundary_crossings: 0,
            coverage: None,
        };
        mesh.n_components()
    ];
    let mut polygon: Vec<[f64; 2]> = Vec::new();
    for boundary in mesh.boundary_loops() {
        let c = component[mesh.outgoing_halfedge(boundary[0] as usize) / 3] as usize;
        polygon.clear();
        polygon.extend(boundary.iter().map(|&v| uv[v as usize]));
        overlap[c].boundary_crossings += polygon_self_crossings(&polygon);
        overlap[c].wedges_over_2pi += boundary
            .iter()
            .filter(|&&v| angle_sums[v as usize] > 2.0 * PI + FULL_TURN_TOLERANCE)
            .count();
    }
    if overlap_samples_per_triangle > 0.0 {
        // The triangles of every chart, chart after chart.
        let mut start = vec![0usize; mesh.n_components() + 1];
        for &c in component {
            start[c as usize + 1] += 1;
        }
        for c in 0..mesh.n_components() {
            start[c + 1] += start[c];
        }
        let mut next = start.clone();
        let mut triangles = vec![[[0.0; 2]; 3]; mesh.n_faces()];
        for (f, face) in mesh.faces().iter().enumerate() {
            let slot = &mut next[component[f] as usize];
            triangles[*slot] = face.map(|v| uv[v as usize]);
            *slot += 1;
        }
        for (c, report) in overlap.iter_mut().enumerate() {
            report.coverage = coverage(
                &triangles[start[c]..start[c + 1]],
                overlap_samples_per_triangle,
            );
        }
    }
    let mut reports = Vec::with_capacity(mesh.n_components());
    for (chart, (d, o)) in distortion.into_iter().zip(overlap).enumerate() {
        let d = d.finish();
        if d.degenerate_triangles == d.n_faces {
            return Err(BffError::DegenerateLayout {
                chart: chart as u32,
                n_faces: d.n_faces,
            });
        }
        if !d.plane_area.is_finite() {
            // Coordinates beyond about 1e154: finite, but their products are not.
            return Err(BffError::NonFiniteResult {
                what: "the area of the layout",
                chart: Some(chart as u32),
            });
        }
        reports.push((d, o));
    }
    Ok(LayoutQuality {
        charts: reports,
        angle_sums,
    })
}

/// Checks that every point of the boundary curve of a chart is finite.
pub(crate) fn check_finite_curve(chart: u32, positions: &[[f64; 2]]) -> Result<(), BffError> {
    if positions
        .iter()
        .all(|p| p[0].is_finite() && p[1].is_finite())
    {
        Ok(())
    } else {
        Err(BffError::NonFiniteResult {
            what: "the boundary curve",
            chart: Some(chart),
        })
    }
}

/// Checks that every coordinate of a layout is finite.
pub(crate) fn check_finite_layout(mesh: &TriMesh, uv: &[[f64; 2]]) -> Result<(), BffError> {
    match uv
        .iter()
        .position(|p| !(p[0].is_finite() && p[1].is_finite()))
    {
        None => Ok(()),
        Some(v) => Err(BffError::NonFiniteResult {
            what: "the layout",
            chart: Some(mesh.face_components()[mesh.outgoing_halfedge(v) / 3]),
        }),
    }
}
