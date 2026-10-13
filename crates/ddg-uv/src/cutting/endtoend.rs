//! End to end: insert the cut, then unwrap with the old seams and the cut as seams.

use super::describe::CutOnPolygons;
use super::error::{CuttingError, internal};
use super::insert::{Inserted, insert};
use super::session::Session;
use crate::input::{MeshInput, PolygonMesh};
use crate::options::UnwrapOptions;
use crate::pipeline::{UvMap, unwrap};
use ddg_core::io::Cone;

/// Whether `unwrap` may place cone singularities after the cut.
///
/// The existing options of [`UnwrapOptions`] have no switch for cones as such; cones are
/// steered per vertex: `forbidden_vertices` of [`MeshInput`] gives a vertex an infinite cone
/// cost, `pinned_vertices` and `prescribed_cones` ask for cones. [`ConePlacement::Off`] uses
/// exactly that mechanism, for every vertex of the new mesh, so that the charts are the
/// regions of the cut, flattened without cones (the pipeline still cuts a chart that is not a
/// disk, along polygon edges, and a closed chart needs the policy of
/// [`UnwrapOptions::closed_surfaces`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ConePlacement {
    /// `unwrap` as it is: cones where the cone solver puts them, steered by the pinned,
    /// prescribed and forbidden vertices of the input. Vertices keep their numbers across the
    /// insertion, so that steering of the old mesh stays valid.
    #[default]
    Optimal,
    /// No cone anywhere: every vertex is forbidden. The input must not pin a vertex or
    /// prescribe a cone (an error names the vertex).
    Off,
}

/// The options of [`unwrap_with_cut`].
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct CutUnwrapOptions {
    /// The options of `unwrap`, as they are.
    pub unwrap: UnwrapOptions,
    /// Whether cones may be placed.
    pub cones: ConePlacement,
}

/// The result of [`unwrap_with_cut`].
#[derive(Debug, Clone, PartialEq)]
pub struct CutUnwrap {
    /// The insertion: the new polygon mesh and the maps.
    pub inserted: Inserted,
    /// The seams that were given to `unwrap`: the old seams (an old seam edge that received
    /// points is replaced by the pieces of its chain) and the cut.
    pub seams: Vec<[u32; 2]>,
    /// The UV map of the new mesh.
    pub map: UvMap,
}

/// Inserts `cut` into the mesh of `input` and unwraps the result with the seams of `input`
/// plus the cut as seams.
///
/// # Errors
/// As [`insert`]; [`CuttingError::InvalidOption`] for [`ConePlacement::Off`] with pinned or
/// prescribed vertices; [`CuttingError::Mesh`] with the error of `unwrap` (it names polygons
/// and vertices of the new mesh).
pub fn unwrap_with_cut(
    input: &MeshInput<'_>,
    cut: &CutOnPolygons,
    options: &CutUnwrapOptions,
) -> Result<CutUnwrap, CuttingError> {
    let mesh = PolygonMesh {
        positions: input.positions.to_vec(),
        corner_vertices: input.corner_vertices.to_vec(),
        polygon_offsets: input.polygon_offsets.to_vec(),
        triangles: input.triangles.to_vec(),
    };
    unwrap_inserted(
        &mesh,
        input.seam_edges,
        input.pinned_vertices,
        input.prescribed_cones,
        input.forbidden_vertices,
        cut,
        options,
    )
}

pub(super) fn unwrap_inserted(
    mesh: &PolygonMesh,
    seams: &[[u32; 2]],
    pinned: &[u32],
    prescribed: &[Cone],
    forbidden: &[u32],
    cut: &CutOnPolygons,
    options: &CutUnwrapOptions,
) -> Result<CutUnwrap, CuttingError> {
    let inserted = insert(mesh, cut)?;
    let all_seams = inserted.seams_after(seams);
    let forbidden_all: Vec<u32>;
    let forbidden: &[u32] = match options.cones {
        ConePlacement::Optimal => forbidden,
        ConePlacement::Off => {
            if let Some(&v) = pinned.first() {
                return Err(CuttingError::InvalidOption {
                    name: "cones",
                    reason: format!("cone placement is off but vertex {v} is pinned"),
                });
            }
            if let Some(c) = prescribed.first() {
                return Err(CuttingError::InvalidOption {
                    name: "cones",
                    reason: format!(
                        "cone placement is off but a cone of {} rad is prescribed at vertex {}",
                        c.angle, c.vertex
                    ),
                });
            }
            forbidden_all = (0..inserted.mesh.positions.len() as u32).collect();
            &forbidden_all
        }
    };
    let new_input = MeshInput {
        seam_edges: &all_seams,
        pinned_vertices: pinned,
        prescribed_cones: prescribed,
        forbidden_vertices: forbidden,
        ..inserted.mesh.input()
    };
    let map = unwrap(&new_input, &options.unwrap).map_err(CuttingError::Mesh)?;
    if map.uvs.len() != inserted.mesh.corner_vertices.len() {
        return Err(internal("the UV map does not have one UV per corner of the new mesh"));
    }
    Ok(CutUnwrap {
        inserted,
        seams: all_seams,
        map,
    })
}

impl Session {
    /// Describes the cut of the current state ([`Session::describe_cut`]), inserts it into the
    /// session's polygon mesh and unwraps the result: [`unwrap_with_cut`] on the session's
    /// mesh, seams and steering. Returns the description too.
    ///
    /// # Errors
    /// As [`Session::describe_cut`] and [`unwrap_with_cut`].
    pub fn unwrap(
        &self,
        min_crossing_offset: f64,
        options: &CutUnwrapOptions,
    ) -> Result<(CutOnPolygons, CutUnwrap), CuttingError> {
        let cut = self.describe_cut(min_crossing_offset)?;
        let result = unwrap_inserted(
            self.mesh(),
            self.seams(),
            &self.pinned,
            &self.prescribed,
            &self.forbidden,
            &cut,
            options,
        )?;
        Ok((cut, result))
    }
}
