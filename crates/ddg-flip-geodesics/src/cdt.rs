//! The intrinsic constrained Delaunay triangulation of a network (paper Sec. 6.3, Fig. 25;
//! SPEC §6.3, F20).
//!
//! After shortening, every path is an edge path of `T`. Marking every edge that carries a
//! segment as constrained and flipping the other edges to Delaunay (IT N6, IC §4.3: the shared
//! [`Triangulation::flip_to_delaunay`], which never flips a constrained edge, R3) gives the
//! intrinsic Delaunay triangulation constrained to the network, with the same vertex set:
//! the reference's `makeDelaunay` (`gc: fg.cpp:1035-1045`). Edges that carry segments are not
//! flipped, so the segments, their layering and every path stay exactly as they were.
//!
//! The marks are left on the triangulation (as the reference leaves them): the result is
//! Delaunay with respect to them, and a later refinement can use them. Marks that were there
//! before are kept. A later [`GeodesicNetwork::shorten`] that moves a path off a marked edge
//! cannot flip that edge any more (the triangulation refuses, a refused flip, F6); clear the
//! marks first with [`GeodesicNetwork::clear_constraints`] if the network is to change again.

use ddg_intrinsic::{DelaunayCriterion, DelaunayOptions, DelaunayReport, Triangulation};

use crate::error::FlipError;
use crate::network::GeodesicNetwork;

/// What [`GeodesicNetwork::make_delaunay`] did.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CdtReport {
    /// The Delaunay flipping. Its statistics leave the constrained edges out; its
    /// `n_constrained` counts every constrained edge at the end (those of the network and
    /// any marked before).
    pub delaunay: DelaunayReport,
    /// Edges that carry a segment, all constrained now.
    pub n_network_edges: usize,
    /// Of those, the ones that were not constrained before.
    pub n_newly_constrained: usize,
}

impl<T: Triangulation> GeodesicNetwork<T> {
    /// The intrinsic constrained Delaunay triangulation of the network (SPEC §6.3; module
    /// [`crate::cdt`]): every edge that carries a segment is marked constrained, then the
    /// triangulation is flipped to Delaunay with `options`, constrained edges left alone. The
    /// network does not change; every path is an edge path of the result.
    ///
    /// # Errors
    /// [`FlipError::InvalidOptions`] for a criterion tolerance that is negative or not finite
    /// (nothing is changed then); [`FlipError::Internal`] if the flipping changed an edge
    /// carrying a segment (a bug).
    pub fn make_delaunay(&mut self, options: &DelaunayOptions) -> Result<CdtReport, FlipError> {
        self.constrained_delaunay(options, true)
    }

    /// [`make_delaunay`](Self::make_delaunay); with `keep_marks = false` the marks it set are
    /// cleared again afterwards (Bézier's retry of a refused split, [`crate::bezier`]).
    pub(crate) fn constrained_delaunay(
        &mut self,
        options: &DelaunayOptions,
        keep_marks: bool,
    ) -> Result<CdtReport, FlipError> {
        let tolerance = match options.criterion {
            DelaunayCriterion::AngleSum { tolerance }
            | DelaunayCriterion::CotanWeight { tolerance } => tolerance,
        };
        if !(tolerance.is_finite() && tolerance >= 0.0) {
            return Err(FlipError::InvalidOptions {
                what: format!(
                    "the tolerance of the Delaunay criterion must be finite and >= 0, got \
                     {tolerance}"
                ),
            });
        }
        let n_edges = self.tri.complex().n_edges();
        let mut n_network_edges = 0;
        let mut n_newly_constrained = 0;
        let mut ends = Vec::new();
        let mut marked = Vec::new();
        for e in 0..n_edges {
            if self.network.carries_segment(e) {
                n_network_edges += 1;
                ends.push((e, self.tri.complex().edge_vertices(e)));
                if !self.tri.is_constrained(e) {
                    n_newly_constrained += 1;
                    self.tri.set_constrained(e, true)?;
                    marked.push(e);
                }
            }
        }
        let delaunay = self.tri.flip_to_delaunay(options);
        if !keep_marks {
            for &e in &marked {
                self.tri.set_constrained(e, false)?;
            }
        }
        let delaunay = delaunay?;
        let complex = self.tri.complex();
        for (e, v) in ends {
            if complex.edge_vertices(e) != v {
                return Err(FlipError::Internal {
                    what: format!(
                        "Delaunay flipping changed edge {e} ({v:?} → {:?}), which carries a \
                         segment and is constrained",
                        complex.edge_vertices(e)
                    ),
                });
            }
        }
        Ok(CdtReport {
            delaunay,
            n_network_edges,
            n_newly_constrained,
        })
    }

    /// Clears the constrained mark of every edge of the triangulation (those set by
    /// [`make_delaunay`](Self::make_delaunay) and any other). Returns how many were marked.
    /// Nothing else changes.
    pub fn clear_constraints(&mut self) -> usize {
        let marked: Vec<usize> = (0..self.tri.complex().n_edges())
            .filter(|&e| self.tri.is_constrained(e))
            .collect();
        // In range by construction, so `set_constrained` cannot fail.
        marked
            .into_iter()
            .filter(|&e| self.tri.set_constrained(e, false).is_ok())
            .count()
    }
}
