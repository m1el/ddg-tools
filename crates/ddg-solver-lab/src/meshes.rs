//! Test meshes: dataset files, generated meshes, and a dataset mesh with a hole cut in it.

use ddg_core::generators::{bumpy_grid_disk, icosphere, perturb};
use ddg_core::io::read_obj;
use ddg_core::{Geometry, TriMesh};

pub struct Case {
    pub name: String,
    pub mesh: TriMesh,
    pub geometry: Geometry,
}

pub fn data_dir() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../data/optimal-cones/optimal-cones-data/data")
}

fn build(name: &str, positions: Vec<[f64; 3]>, faces: Vec<[u32; 3]>) -> Result<Case, String> {
    let mesh = TriMesh::new(positions.len(), faces).map_err(|e| format!("{name}: {e}"))?;
    let geometry = Geometry::from_positions(&mesh, &positions).map_err(|e| format!("{name}: {e}"))?;
    Ok(Case { name: name.to_string(), mesh, geometry })
}

/// `NAME` (dataset), `NAME+hole` (the dataset mesh with the faces within a few rings of vertex 0
/// removed, so that it has a boundary), `bump:N` (bumpy disk with N x N cells), `sphere:L`
/// (perturbed icosphere of level L).
pub fn load(spec: &str) -> Result<Case, String> {
    if let Some(n) = spec.strip_prefix("bump:") {
        let n: usize = n.parse().map_err(|_| "bump:N".to_string())?;
        let g = bumpy_grid_disk(n, n, 0.3).map_err(|e| e.to_string())?;
        let mut p = g.positions;
        perturb(&mut p, [0.2 / n as f64, 0.2 / n as f64, 0.02], 7);
        return build(spec, p, g.faces);
    }
    if let Some(l) = spec.strip_prefix("sphere:") {
        let l: u32 = l.parse().map_err(|_| "sphere:L".to_string())?;
        let g = icosphere(l).map_err(|e| e.to_string())?;
        let mut p = g.positions;
        perturb(&mut p, [0.02, 0.02, 0.02], 3);
        return build(spec, p, g.faces);
    }
    let (base, hole) = match spec.strip_suffix("+hole") {
        Some(b) => (b, true),
        None => (spec, false),
    };
    let obj = read_obj(data_dir().join(format!("{base}.obj"))).map_err(|e| format!("{base}: {e}"))?;
    let faces = obj.triangles().ok_or_else(|| format!("{base}: not triangles"))?;
    if !hole {
        return build(spec, obj.positions, faces);
    }
    let full = TriMesh::new(obj.positions.len(), faces).map_err(|e| e.to_string())?;
    // remove faces within 3 rings of vertex 0
    let mut dist = vec![u32::MAX; full.n_vertices()];
    let mut queue = std::collections::VecDeque::new();
    dist[0] = 0;
    queue.push_back(0usize);
    while let Some(v) = queue.pop_front() {
        if dist[v] >= 3 {
            continue;
        }
        for w in full.vertex_neighbors(v) {
            if dist[w] == u32::MAX {
                dist[w] = dist[v] + 1;
                queue.push_back(w);
            }
        }
    }
    let keep: Vec<u32> = (0..full.n_faces() as u32)
        .filter(|&f| full.faces()[f as usize].iter().all(|&v| dist[v as usize] > 2))
        .collect();
    let sub = full.submesh(&keep).map_err(|e| e.to_string())?;
    let positions: Vec<[f64; 3]> = sub.vertex_to_parent.iter().map(|&v| obj.positions[v as usize]).collect();
    let geometry = Geometry::from_positions(&sub.mesh, &positions).map_err(|e| e.to_string())?;
    Ok(Case { name: spec.to_string(), mesh: sub.mesh, geometry })
}
