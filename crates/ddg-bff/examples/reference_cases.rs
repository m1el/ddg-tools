//! Writes generated test cases (OBJ files and cone files) for the differential test against
//! the reference C++ implementation.
//!
//! ```text
//! cargo run --release -p ddg-bff --example reference_cases -- OUT_DIR
//! ```
//!
//! For every case `OUT_DIR/<name>.obj` and `OUT_DIR/<name>-cones.txt` are written. Feed them
//! to the driver around the reference library and compare its output with
//! `examples/reference_diff.rs`. All cases are spheres or disks with one boundary loop, the
//! topologies the reference handles without hole filling.

use std::f64::consts::PI;
use std::path::Path;

use ddg_core::generators::{self, GeneratedMesh, SplitMix64};
use ddg_core::io::write_obj;
use ddg_core::operators::angle_defects;
use ddg_core::{Geometry, TriMesh};

fn write(dir: &Path, name: &str, g: &GeneratedMesh, cones: &[(u32, f64)]) {
    let corner_verts: Vec<u32> = g.faces.iter().flatten().copied().collect();
    let face_starts: Vec<u32> = (0..g.faces.len() as u32).map(|f| 3 * f).collect();
    write_obj(
        dir.join(format!("{name}.obj")),
        &g.positions,
        &corner_verts,
        &face_starts,
        None,
    )
    .unwrap_or_else(|e| panic!("{e}"));
    let text: String = cones
        .iter()
        .map(|&(v, angle)| format!("{} {:?}\n", v + 1, angle))
        .collect();
    std::fs::write(dir.join(format!("{name}-cones.txt")), text).unwrap();
    println!(
        "{name}: {} vertices, {} faces, {} cones",
        g.positions.len(),
        g.faces.len(),
        cones.len()
    );
}

/// `k` distinct interior vertices.
fn pick(mesh: &TriMesh, k: usize, rng: &mut SplitMix64) -> Vec<u32> {
    let mut interior = mesh.interior_vertices();
    for i in (1..interior.len()).rev() {
        interior.swap(i, (rng.next_u64() % (i as u64 + 1)) as usize);
    }
    interior.truncate(k);
    interior
}

fn main() {
    let dir = std::env::args()
        .nth(1)
        .expect("usage: reference_cases OUT_DIR");
    let dir = Path::new(&dir);
    std::fs::create_dir_all(dir).unwrap();

    // Sphere with the twelve cones of the icosahedron.
    let g = generators::icosphere(4).unwrap();
    let cones: Vec<(u32, f64)> = (0..12).map(|v| (v, PI / 3.0)).collect();
    write(dir, "sphere12", &g, &cones);

    // Bumpy spheres with random positive cones that sum to 4π.
    for (name, level, k, seed) in [("sphere5", 4u32, 5usize, 11u64), ("sphere9", 5, 9, 12)] {
        let mut g = generators::icosphere(level).unwrap();
        let mut rng = SplitMix64::new(seed);
        for p in &mut g.positions {
            let s = 1.0 + 0.15 * (3.0 * p[0]).sin() * (2.0 * p[1]).cos() + 0.1 * p[2];
            *p = [s * p[0], s * p[1], s * p[2]];
        }
        let mesh = g.topology().unwrap();
        let weights: Vec<f64> = (0..k).map(|_| rng.uniform(0.6, 1.4)).collect();
        let total: f64 = weights.iter().sum();
        let cones: Vec<(u32, f64)> = pick(&mesh, k, &mut rng)
            .into_iter()
            .zip(&weights)
            .map(|(v, w)| (v, 4.0 * PI * w / total))
            .collect();
        write(dir, name, &g, &cones);
    }

    // A flat octahedron with its own angle defects as cones: the unfolding is an isometry.
    {
        let mut positions = vec![
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
        ];
        let mut faces: Vec<[u32; 3]> = vec![
            [0, 2, 4],
            [2, 1, 4],
            [1, 3, 4],
            [3, 0, 4],
            [2, 0, 5],
            [1, 2, 5],
            [3, 1, 5],
            [0, 3, 5],
        ];
        for _ in 0..4 {
            let mut midpoint = std::collections::HashMap::new();
            let mut next = Vec::with_capacity(4 * faces.len());
            for &[a, b, c] in &faces {
                let mut mid = |p: u32, q: u32| {
                    *midpoint.entry((p.min(q), p.max(q))).or_insert_with(|| {
                        let (x, y) = (positions[p as usize], positions[q as usize]);
                        positions.push([
                            0.5 * (x[0] + y[0]),
                            0.5 * (x[1] + y[1]),
                            0.5 * (x[2] + y[2]),
                        ]);
                        positions.len() as u32 - 1
                    })
                };
                let (ab, bc, ca) = (mid(a, b), mid(b, c), mid(c, a));
                next.extend([[a, ab, ca], [ab, b, bc], [ca, bc, c], [ab, bc, ca]]);
            }
            faces = next;
        }
        let g = GeneratedMesh { positions, faces };
        let mesh = g.topology().unwrap();
        let geometry = Geometry::from_positions(&mesh, &g.positions).unwrap();
        let omega = angle_defects(&mesh, &geometry).unwrap();
        let cones: Vec<(u32, f64)> = (0..mesh.n_vertices())
            .filter(|&v| omega[v].abs() > 1e-9)
            .map(|v| (v as u32, omega[v]))
            .collect();
        write(dir, "octahedron", &g, &cones);
    }

    // Bumpy disks: without cones, and with cones of both signs.
    let g = generators::bumpy_grid_disk(40, 30, 0.15).unwrap();
    write(dir, "disk0", &g, &[]);
    let mesh = g.topology().unwrap();
    let mut rng = SplitMix64::new(21);
    let cones: Vec<(u32, f64)> = pick(&mesh, 6, &mut rng)
        .into_iter()
        .map(|v| (v, rng.uniform(-1.2, 1.6)))
        .collect();
    write(dir, "disk6", &g, &cones);
    let g = generators::bumpy_grid_disk(90, 70, 0.1).unwrap();
    let mesh = g.topology().unwrap();
    let cones: Vec<(u32, f64)> = pick(&mesh, 25, &mut rng)
        .into_iter()
        .map(|v| (v, rng.uniform(-0.8, 1.0)))
        .collect();
    write(dir, "disk25", &g, &cones);
}
