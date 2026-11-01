//! OBJ input and output of triangle meshes through `ddg_core::io` (`SPEC.md` §9.1).
//!
//! Only triangles are accepted; a polygon with another number of corners is
//! [`RsError::PolygonFace`] naming the face. Texture coordinates, normals, groups and
//! materials are ignored on input and not written. Positions are read as they are; they are
//! validated (finite, no degenerate face) by the functions that use them
//! ([`Surface::validate_positions`]).

use std::path::Path;

use ddg_core::io::ObjMesh;

use crate::error::RsError;
use crate::surface::Surface;

fn surface_of(mesh: ObjMesh) -> Result<(Surface, Vec<[f64; 3]>), RsError> {
    let mut faces = Vec::with_capacity(mesh.n_faces());
    for f in 0..mesh.n_faces() {
        let r = mesh.face_corners(f);
        if r.len() != 3 {
            return Err(RsError::PolygonFace {
                face: f as u32,
                corners: r.len(),
            });
        }
        let c = &mesh.corner_verts[r];
        faces.push([c[0], c[1], c[2]]);
    }
    let surface = Surface::new(mesh.positions.len(), faces)?;
    Ok((surface, mesh.positions))
}

/// Reads a triangle mesh from an OBJ file and validates its topology ([`Surface::new`]).
///
/// # Errors
/// [`RsError::Io`], [`RsError::Parse`] (with the line), [`RsError::PolygonFace`], and the
/// topology errors of [`Surface::new`].
pub fn read_obj(path: impl AsRef<Path>) -> Result<(Surface, Vec<[f64; 3]>), RsError> {
    surface_of(ddg_core::io::read_obj(path)?)
}

/// As [`read_obj`], from OBJ text.
///
/// # Errors
/// As [`read_obj`].
pub fn parse_obj(text: &str) -> Result<(Surface, Vec<[f64; 3]>), RsError> {
    surface_of(ddg_core::io::parse_obj(text)?)
}

/// The OBJ text of `surface` with `positions` (`v` and `f` lines, faces in the caller's
/// vertex order).
///
/// # Errors
/// [`RsError::PositionCountMismatch`]; [`RsError::NonFinitePosition`].
pub fn obj_to_string(surface: &Surface, positions: &[[f64; 3]]) -> Result<String, RsError> {
    check_positions(surface, positions)?;
    let (corners, starts) = flat_faces(surface);
    Ok(ddg_core::io::obj_to_string(
        positions, &corners, &starts, None,
    )?)
}

/// Writes `surface` with `positions` as an OBJ file.
///
/// # Errors
/// As [`obj_to_string`]; [`RsError::Io`] if the file cannot be written.
pub fn write_obj(
    path: impl AsRef<Path>,
    surface: &Surface,
    positions: &[[f64; 3]],
) -> Result<(), RsError> {
    check_positions(surface, positions)?;
    let (corners, starts) = flat_faces(surface);
    Ok(ddg_core::io::write_obj(
        path, positions, &corners, &starts, None,
    )?)
}

fn check_positions(surface: &Surface, positions: &[[f64; 3]]) -> Result<(), RsError> {
    if positions.len() != surface.n_vertices() {
        return Err(RsError::PositionCountMismatch {
            expected: surface.n_vertices(),
            found: positions.len(),
        });
    }
    if let Some(v) = positions
        .iter()
        .position(|p| !p.iter().all(|x| x.is_finite()))
    {
        return Err(RsError::NonFinitePosition { vertex: v as u32 });
    }
    Ok(())
}

fn flat_faces(surface: &Surface) -> (Vec<u32>, Vec<u32>) {
    let corners: Vec<u32> = surface.faces().iter().flatten().copied().collect();
    let starts: Vec<u32> = (0..surface.n_faces() as u32).map(|f| 3 * f).collect();
    (corners, starts)
}
