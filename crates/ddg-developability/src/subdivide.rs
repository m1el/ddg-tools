//! 4-1 midpoint subdivision (SPEC §5.2, D20).
//!
//! Every edge, boundary edges included, gets a new vertex at its midpoint; every face is replaced
//! by four. No smoothing: the surface and the face normals are unchanged, hence every developability
//! energy is preserved (old stars keep their normals, angles and vector-area direction; new stars
//! are hinges).

use ddg_core::topology::TriMesh;

use crate::error::DevError;

/// Result of [`subdivide`].
#[derive(Debug, Clone, PartialEq)]
pub struct Subdivided {
    /// The refined mesh: `|V| + |E|` vertices, `4|F|` faces.
    pub mesh: TriMesh,
    /// Positions of the refined mesh; old vertices keep their ids and positions.
    pub positions: Vec<[f64; 3]>,
    /// For each edge of the parent mesh (`TriMesh::edges()` order), the id of its midpoint vertex
    /// (`|V| + e`).
    pub edge_midpoint: Vec<u32>,
    /// For each new face, the parent face it lies in (`t` for new faces `4t … 4t+3`).
    pub face_parent: Vec<u32>,
}

impl Subdivided {
    /// Pinned vertex ids of the refined mesh: old pinned vertices stay pinned, a midpoint is pinned
    /// iff both ends of its edge are (SPEC §5.2). `old_mesh` must be the mesh that was subdivided;
    /// ids out of range are an error.
    ///
    /// # Errors
    /// [`DevError::VertexOutOfRange`] for a pinned id `>= |V|`; [`DevError::Internal`] if
    /// `old_mesh` does not match.
    pub fn pinned_after(&self, old_mesh: &TriMesh, pinned: &[u32]) -> Result<Vec<u32>, DevError> {
        let nv = old_mesh.n_vertices();
        if self.edge_midpoint.len() != old_mesh.n_edges() {
            return Err(DevError::Internal("pinned_after: mesh does not match the subdivision".into()));
        }
        let mut is = vec![false; nv];
        for &p in pinned {
            if p as usize >= nv {
                return Err(DevError::VertexOutOfRange { what: "pinned", vertex: p, n_vertices: nv });
            }
            is[p as usize] = true;
        }
        let mut out: Vec<u32> = (0..nv as u32).filter(|&v| is[v as usize]).collect();
        for (e, &[a, b]) in old_mesh.edges().iter().enumerate() {
            if is[a as usize] && is[b as usize] {
                out.push(self.edge_midpoint[e]);
            }
        }
        Ok(out)
    }
}

/// Midpoint 4-1 subdivision of a triangle mesh.
///
/// New vertex of edge `e`: id `|V| + e`, position `(f_a + f_b)/2`. Face `t = (a, b, c)` with
/// midpoints `m_ab, m_bc, m_ca` becomes faces `4t … 4t+3`: `(a, m_ab, m_ca)`, `(b, m_bc, m_ab)`,
/// `(c, m_ca, m_bc)`, `(m_ab, m_bc, m_ca)`; orientation is kept.
///
/// [formula: dev-midpoint-subdivision]
///
/// # Errors
/// [`DevError::LengthMismatch`] if `positions.len() != |V|`, [`DevError::NonFinitePosition`], and
/// [`DevError::Mesh`] if the refined face list is rejected (cannot happen for a valid mesh).
pub fn subdivide(mesh: &TriMesh, positions: &[[f64; 3]]) -> Result<Subdivided, DevError> {
    let nv = mesh.n_vertices();
    if positions.len() != nv {
        return Err(DevError::LengthMismatch { what: "positions", expected: nv, got: positions.len() });
    }
    if let Some(v) = positions.iter().position(|p| p.iter().any(|c| !c.is_finite())) {
        return Err(DevError::NonFinitePosition { vertex: v as u32 });
    }
    let mut new_pos = Vec::with_capacity(nv + mesh.n_edges());
    new_pos.extend_from_slice(positions);
    let mut edge_midpoint = Vec::with_capacity(mesh.n_edges());
    for (e, &[a, b]) in mesh.edges().iter().enumerate() {
        let (pa, pb) = (positions[a as usize], positions[b as usize]);
        new_pos.push([0, 1, 2].map(|k| 0.5 * (pa[k] + pb[k])));
        edge_midpoint.push((nv + e) as u32);
    }
    let mid = |u: u32, v: u32| -> Result<u32, DevError> {
        mesh.find_edge(u as usize, v as usize)
            .map(|e| (nv + e) as u32)
            .ok_or_else(|| DevError::Internal(format!("edge ({u}, {v}) not found while subdividing")))
    };
    let mut faces = Vec::with_capacity(4 * mesh.n_faces());
    let mut face_parent = Vec::with_capacity(4 * mesh.n_faces());
    for (t, &[a, b, c]) in mesh.faces().iter().enumerate() {
        let (mab, mbc, mca) = (mid(a, b)?, mid(b, c)?, mid(c, a)?);
        faces.extend_from_slice(&[[a, mab, mca], [b, mbc, mab], [c, mca, mbc], [mab, mbc, mca]]);
        face_parent.extend_from_slice(&[t as u32; 4]);
    }
    let refined = TriMesh::new(new_pos.len(), faces)?;
    Ok(Subdivided { mesh: refined, positions: new_pos, edge_midpoint, face_parent })
}
