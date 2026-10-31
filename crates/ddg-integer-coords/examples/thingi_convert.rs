//! WP10, SPEC §13 Q1: converts the raw Thingi10k meshes to manifold triangle meshes.
//!
//! ```text
//! cargo run --release -p ddg-integer-coords --example thingi_convert -- \
//!     [--raw DIR] [--out DIR] [--tsv FILE] [--threads N] [--only ID,ID,...]
//! ```
//!
//! Defaults: `data/thingi10k/raw/Thingi10K/raw_meshes/`, `data/thingi10k/converted/`,
//! `data/thingi10k/converted.tsv`, 56 threads.
//!
//! The rule (the paper converted the raw files to PLY with MeshLab; ours is explicit):
//!
//! 1. read the file (STL ASCII and binary, OBJ through `ddg_core::io::parse_obj`, OFF, PLY ASCII
//!    and binary; the readers for the formats `ddg_core` lacks are in this file, minimal, and
//!    refuse what they do not understand);
//! 2. merge vertices whose three coordinates are exactly equal (`-0.0 == 0.0`), nothing else:
//!    no tolerance, no removal of degenerate or duplicate faces;
//! 3. triangulate polygons by a fan from their first corner;
//! 4. drop vertices that no face uses (STL has none; an OBJ may);
//! 5. keep the mesh if `ddg_core::TriMesh::new` accepts the faces (every edge on one or two
//!    faces, consistent orientation, one fan per vertex, no repeated vertex in a face, no
//!    duplicate face). Nothing is repaired; a mesh that `orient_faces` could repair is excluded
//!    and marked in the column `orient_repairable`.
//!
//! Kept meshes are written as OBJ (`v` with the shortest round-trip `f64` text, `f` 1-based).
//! The TSV has one row per file: id, format, status (`kept` / `excluded`), reason, sizes, the
//! metadata flags of Thingi10k (`geometry_data.csv`) when the CSV is found next to the data.

use std::collections::HashMap;
use std::fmt::Write as _;
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Condvar, Mutex};
use std::time::Instant;

use ddg_core::TriMesh;
use ddg_core::io::parse_obj;
use ddg_core::topology::{MeshError, compact_vertices, orient_faces};

/// Raw polygons: positions and faces of any size.
struct Raw {
    positions: Vec<[f64; 3]>,
    polygons: Vec<Vec<u32>>,
}

fn read_stl(bytes: &[u8]) -> Result<Raw, String> {
    // Binary: 80-byte header, u32 count, 50 bytes per triangle.
    if bytes.len() >= 84 {
        let n = u32::from_le_bytes(bytes[80..84].try_into().unwrap()) as usize;
        if bytes.len() == 84 + 50 * n {
            let mut positions = Vec::with_capacity(3 * n);
            let mut polygons = Vec::with_capacity(n);
            for t in 0..n {
                let base = 84 + 50 * t + 12;
                let mut poly = Vec::with_capacity(3);
                for k in 0..3 {
                    let o = base + 12 * k;
                    let c = |i: usize| {
                        f32::from_le_bytes(bytes[o + 4 * i..o + 4 * i + 4].try_into().unwrap())
                            as f64
                    };
                    poly.push(positions.len() as u32);
                    positions.push([c(0), c(1), c(2)]);
                }
                polygons.push(poly);
            }
            return Ok(Raw {
                positions,
                polygons,
            });
        }
    }
    // ASCII.
    let text = String::from_utf8_lossy(bytes);
    let mut positions = Vec::new();
    let mut polygons: Vec<Vec<u32>> = Vec::new();
    let mut current: Vec<u32> = Vec::new();
    let mut seen_solid = false;
    for line in text.lines() {
        let mut it = line.split_whitespace();
        match it.next() {
            Some(w) if w.eq_ignore_ascii_case("solid") => seen_solid = true,
            Some(w) if w.eq_ignore_ascii_case("vertex") => {
                let mut c = [0.0; 3];
                for slot in &mut c {
                    *slot = it
                        .next()
                        .and_then(|s| s.parse::<f64>().ok())
                        .ok_or_else(|| format!("bad vertex line {line:?}"))?;
                }
                current.push(positions.len() as u32);
                positions.push(c);
            }
            Some(w) if w.eq_ignore_ascii_case("endfacet") => {
                if current.len() != 3 {
                    return Err(format!("facet with {} vertices", current.len()));
                }
                polygons.push(std::mem::take(&mut current));
            }
            _ => {}
        }
    }
    if !seen_solid && polygons.is_empty() {
        return Err("neither a binary STL of the right size nor ASCII STL".into());
    }
    if !current.is_empty() {
        return Err("unterminated facet".into());
    }
    Ok(Raw {
        positions,
        polygons,
    })
}

fn read_obj_raw(bytes: &[u8]) -> Result<Raw, String> {
    let text = String::from_utf8_lossy(bytes);
    let obj = parse_obj(&text).map_err(|e| e.to_string())?;
    let polygons = (0..obj.n_faces())
        .map(|f| obj.corner_verts[obj.face_corners(f)].to_vec())
        .collect();
    Ok(Raw {
        positions: obj.positions,
        polygons,
    })
}

fn read_off(bytes: &[u8]) -> Result<Raw, String> {
    let text = String::from_utf8_lossy(bytes);
    let mut tokens = text
        .lines()
        .map(|l| l.split('#').next().unwrap_or(""))
        .flat_map(|l| l.split_whitespace());
    let magic = tokens.next().ok_or("empty file")?;
    if magic != "OFF" {
        return Err(format!("OFF magic expected, found {magic:?}"));
    }
    let mut num = |what: &str| -> Result<f64, String> {
        tokens
            .next()
            .and_then(|s| s.parse::<f64>().ok())
            .ok_or_else(|| format!("missing or bad number ({what})"))
    };
    let nv = num("vertex count")? as usize;
    let nf = num("face count")? as usize;
    let _ne = num("edge count")?;
    let mut positions = Vec::with_capacity(nv);
    for _ in 0..nv {
        positions.push([num("x")?, num("y")?, num("z")?]);
    }
    let mut polygons = Vec::with_capacity(nf);
    for _ in 0..nf {
        let k = num("face size")? as usize;
        let mut poly = Vec::with_capacity(k);
        for _ in 0..k {
            poly.push(num("index")? as u32);
        }
        polygons.push(poly);
    }
    Ok(Raw {
        positions,
        polygons,
    })
}

fn read_ply(bytes: &[u8]) -> Result<Raw, String> {
    // Header.
    let end = bytes
        .windows(10)
        .position(|w| w == b"end_header")
        .ok_or("no end_header")?;
    let mut body = end + 10;
    while body < bytes.len() && bytes[body] != b'\n' {
        body += 1;
    }
    body += 1;
    let header = String::from_utf8_lossy(&bytes[..end]).to_string();
    let mut format = String::new();
    struct Element {
        name: String,
        count: usize,
        props: Vec<(String, Option<(String, String)>, String)>, // (name, list(count,item), type)
    }
    let mut elements: Vec<Element> = Vec::new();
    for line in header.lines() {
        let w: Vec<&str> = line.split_whitespace().collect();
        match w.as_slice() {
            ["format", f, _] => format = (*f).to_string(),
            ["element", name, n] => elements.push(Element {
                name: (*name).to_string(),
                count: n.parse().map_err(|_| "bad element count")?,
                props: Vec::new(),
            }),
            ["property", "list", c, t, name] => {
                if let Some(e) = elements.last_mut() {
                    e.props.push((
                        (*name).to_string(),
                        Some(((*c).to_string(), (*t).to_string())),
                        String::new(),
                    ));
                }
            }
            ["property", t, name] => {
                if let Some(e) = elements.last_mut() {
                    e.props.push(((*name).to_string(), None, (*t).to_string()));
                }
            }
            _ => {}
        }
    }
    let size = |t: &str| -> Result<usize, String> {
        Ok(match t {
            "char" | "uchar" | "int8" | "uint8" => 1,
            "short" | "ushort" | "int16" | "uint16" => 2,
            "int" | "uint" | "float" | "int32" | "uint32" | "float32" => 4,
            "double" | "float64" => 8,
            _ => return Err(format!("unknown PLY type {t}")),
        })
    };
    let ascii = format == "ascii";
    let big = format == "binary_big_endian";
    if !ascii && format != "binary_little_endian" && !big {
        return Err(format!("unknown PLY format {format:?}"));
    }
    // A scalar reader over the body.
    let text = if ascii {
        Some(String::from_utf8_lossy(&bytes[body..]).to_string())
    } else {
        None
    };
    let mut toks = text.as_ref().map(|t| t.split_whitespace());
    let mut pos = body;
    let mut read = |t: &str| -> Result<f64, String> {
        if let Some(tk) = toks.as_mut() {
            return tk
                .next()
                .and_then(|s| s.parse::<f64>().ok())
                .ok_or_else(|| "bad ascii number".to_string());
        }
        let n = size(t)?;
        let slice = bytes.get(pos..pos + n).ok_or("truncated PLY")?;
        pos += n;
        let mut b = [0u8; 8];
        b[..n].copy_from_slice(slice);
        if big {
            b[..n].reverse();
        }
        Ok(match t {
            "char" | "int8" => b[0] as i8 as f64,
            "uchar" | "uint8" => b[0] as f64,
            "short" | "int16" => i16::from_le_bytes([b[0], b[1]]) as f64,
            "ushort" | "uint16" => u16::from_le_bytes([b[0], b[1]]) as f64,
            "int" | "int32" => i32::from_le_bytes([b[0], b[1], b[2], b[3]]) as f64,
            "uint" | "uint32" => u32::from_le_bytes([b[0], b[1], b[2], b[3]]) as f64,
            "float" | "float32" => f32::from_le_bytes([b[0], b[1], b[2], b[3]]) as f64,
            _ => f64::from_le_bytes(b),
        })
    };
    let mut positions = Vec::new();
    let mut polygons = Vec::new();
    for e in &elements {
        for _ in 0..e.count {
            let mut xyz = [0.0; 3];
            let mut poly: Option<Vec<u32>> = None;
            for (name, list, ty) in &e.props {
                if let Some((ct, it)) = list {
                    let n = read(ct)? as usize;
                    let mut items = Vec::with_capacity(n);
                    for _ in 0..n {
                        items.push(read(it)? as u32);
                    }
                    if e.name == "face" && (name == "vertex_indices" || name == "vertex_index") {
                        poly = Some(items);
                    }
                } else {
                    let v = read(ty)?;
                    match (e.name.as_str(), name.as_str()) {
                        ("vertex", "x") => xyz[0] = v,
                        ("vertex", "y") => xyz[1] = v,
                        ("vertex", "z") => xyz[2] = v,
                        _ => {}
                    }
                }
            }
            match e.name.as_str() {
                "vertex" => positions.push(xyz),
                "face" => polygons.push(poly.ok_or("face without vertex_indices")?),
                _ => {}
            }
        }
    }
    Ok(Raw {
        positions,
        polygons,
    })
}

/// The outcome for one file.
#[derive(Default)]
struct Row {
    id: String,
    format: String,
    status: &'static str,
    reason: String,
    detail: String,
    n_raw_vertices: usize,
    n_raw_polygons: usize,
    n_vertices: usize,
    n_faces: usize,
    n_boundary_edges: usize,
    n_components: usize,
    chi: i64,
    orient_repairable: &'static str,
    seconds: f64,
}

fn mesh_error_reason(e: &MeshError) -> &'static str {
    match e {
        MeshError::NoFaces => "no_faces",
        MeshError::TooLarge { .. } => "too_large",
        MeshError::VertexIndexOutOfRange { .. } => "index_out_of_range",
        MeshError::RepeatedVertexInFace { .. } => "degenerate_face_after_merge",
        MeshError::DuplicateFace { .. } => "duplicate_face",
        MeshError::NonManifoldEdge { .. } => "non_manifold_edge",
        MeshError::InconsistentOrientation { .. } => "inconsistent_orientation",
        MeshError::NonManifoldVertex { .. } => "non_manifold_vertex",
        MeshError::UnreferencedVertex { .. } => "unreferenced_vertex",
        _ => "other",
    }
}

fn convert(path: &Path, out_dir: &Path) -> Row {
    let t0 = Instant::now();
    let mut row = Row {
        id: path
            .file_stem()
            .map(|s| s.to_string_lossy().to_string())
            .unwrap_or_default(),
        format: path
            .extension()
            .map(|s| s.to_string_lossy().to_lowercase())
            .unwrap_or_default(),
        status: "excluded",
        orient_repairable: "-",
        ..Row::default()
    };
    let fail = |mut row: Row, reason: &str, detail: String, t0: Instant| {
        row.reason = reason.to_string();
        row.detail = detail.replace(['\t', '\n'], " ");
        row.seconds = t0.elapsed().as_secs_f64();
        row
    };
    let bytes = match std::fs::read(path) {
        Ok(b) => b,
        Err(e) => return fail(row, "io_error", e.to_string(), t0),
    };
    let raw = match row.format.as_str() {
        "stl" => read_stl(&bytes),
        "obj" => read_obj_raw(&bytes),
        "off" => read_off(&bytes),
        "ply" => read_ply(&bytes),
        other => Err(format!("unsupported format {other}")),
    };
    drop(bytes);
    let raw = match raw {
        Ok(r) => r,
        Err(e) => return fail(row, "parse_error", e, t0),
    };
    row.n_raw_vertices = raw.positions.len();
    row.n_raw_polygons = raw.polygons.len();
    if let Some(i) = raw
        .positions
        .iter()
        .position(|p| p.iter().any(|x| !x.is_finite()))
    {
        return fail(row, "non_finite_coordinate", format!("vertex {i}"), t0);
    }
    // Merge exactly coincident vertices.
    let mut map: HashMap<[u64; 3], u32> = HashMap::with_capacity(raw.positions.len());
    let mut positions: Vec<[f64; 3]> = Vec::new();
    let mut remap = Vec::with_capacity(raw.positions.len());
    for p in &raw.positions {
        let key = p.map(|x| if x == 0.0 { 0u64 } else { x.to_bits() });
        let id = *map.entry(key).or_insert_with(|| {
            positions.push(*p);
            (positions.len() - 1) as u32
        });
        remap.push(id);
    }
    drop(map);
    // Triangulate.
    let mut faces: Vec<[u32; 3]> = Vec::with_capacity(raw.polygons.len());
    for poly in &raw.polygons {
        if poly.len() < 3 {
            return fail(row, "degenerate_polygon", "fewer than 3 corners".into(), t0);
        }
        if let Some(&bad) = poly.iter().find(|&&i| i as usize >= remap.len()) {
            return fail(row, "index_out_of_range", format!("index {bad}"), t0);
        }
        for k in 1..poly.len() - 1 {
            faces.push([
                remap[poly[0] as usize],
                remap[poly[k] as usize],
                remap[poly[k + 1] as usize],
            ]);
        }
    }
    drop(raw);
    // Drop unreferenced vertices.
    let compact = match compact_vertices(positions.len(), &faces) {
        Ok(c) => c,
        Err(e) => return fail(row, mesh_error_reason(&e), e.to_string(), t0),
    };
    let positions: Vec<[f64; 3]> = compact
        .new_to_old
        .iter()
        .map(|&o| positions[o as usize])
        .collect();
    let faces = compact.faces;
    row.n_vertices = positions.len();
    row.n_faces = faces.len();
    let mesh = match TriMesh::new(positions.len(), faces.clone()) {
        Ok(m) => m,
        Err(e) => {
            row.orient_repairable = if matches!(e, MeshError::InconsistentOrientation { .. }) {
                match orient_faces(positions.len(), &faces)
                    .map(|o| TriMesh::new(positions.len(), o.faces))
                {
                    Ok(Ok(_)) => "yes",
                    _ => "no",
                }
            } else {
                "-"
            };
            let reason = mesh_error_reason(&e);
            return fail(row, reason, e.to_string(), t0);
        }
    };
    row.n_boundary_edges = mesh.n_boundary_edges();
    row.n_components = mesh.n_components();
    row.chi = mesh.euler_characteristic();
    // Write the OBJ.
    let mut s = String::with_capacity(positions.len() * 60 + faces.len() * 24);
    for p in &positions {
        let _ = writeln!(s, "v {} {} {}", p[0], p[1], p[2]);
    }
    for f in mesh.faces() {
        let _ = writeln!(s, "f {} {} {}", f[0] + 1, f[1] + 1, f[2] + 1);
    }
    let target = out_dir.join(format!("{}.obj", row.id));
    let tmp = out_dir.join(format!("{}.obj.tmp", row.id));
    if let Err(e) = std::fs::File::create(&tmp)
        .and_then(|mut f| f.write_all(s.as_bytes()))
        .and_then(|()| std::fs::rename(&tmp, &target))
    {
        return fail(row, "io_error", e.to_string(), t0);
    }
    row.status = "kept";
    row.seconds = t0.elapsed().as_secs_f64();
    row
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/thingi10k");
    let mut raw_dir = root.join("raw/Thingi10K/raw_meshes");
    let mut out_dir = root.join("converted");
    let mut tsv = root.join("converted.tsv");
    let mut threads = 56usize;
    let mut only: Option<Vec<String>> = None;
    let mut i = 0;
    while i < args.len() {
        let v = args.get(i + 1).cloned().unwrap_or_default();
        match args[i].as_str() {
            "--raw" => raw_dir = v.into(),
            "--out" => out_dir = v.into(),
            "--tsv" => tsv = v.into(),
            "--threads" => threads = v.parse().expect("--threads N"),
            "--only" => only = Some(v.split(',').map(str::to_string).collect()),
            other => panic!("unknown argument {other}"),
        }
        i += 2;
    }
    std::fs::create_dir_all(&out_dir).expect("create output directory");
    let mut files: Vec<(PathBuf, u64)> = std::fs::read_dir(&raw_dir)
        .expect("read raw directory")
        .filter_map(|e| e.ok())
        .map(|e| (e.path(), e.metadata().map_or(0, |m| m.len())))
        .filter(|(p, _)| {
            only.as_ref().is_none_or(|ids| {
                ids.iter()
                    .any(|id| p.file_stem().is_some_and(|s| s.to_string_lossy() == *id))
            })
        })
        .collect();
    // Largest first, so that the long jobs start early; a byte budget bounds the memory.
    files.sort_by_key(|f| std::cmp::Reverse(f.1));
    let next = AtomicUsize::new(0);
    let rows: Mutex<Vec<Row>> = Mutex::new(Vec::new());
    let budget = (20u64) << 30; // bytes of raw files in flight (peak memory is ~10x the file)
    let in_flight = (Mutex::new(0u64), Condvar::new());
    let t0 = Instant::now();
    std::thread::scope(|s| {
        for _ in 0..threads {
            s.spawn(|| {
                loop {
                    let k = next.fetch_add(1, Ordering::SeqCst);
                    let Some((path, size)) = files.get(k) else {
                        break;
                    };
                    let cost = (*size).max(1 << 20);
                    {
                        let mut g = in_flight.0.lock().unwrap();
                        while *g > 0 && *g + cost > budget {
                            g = in_flight.1.wait(g).unwrap();
                        }
                        *g += cost;
                    }
                    let row = convert(path, &out_dir);
                    {
                        *in_flight.0.lock().unwrap() -= cost;
                        in_flight.1.notify_all();
                    }
                    let mut r = rows.lock().unwrap();
                    r.push(row);
                    if r.len() % 500 == 0 {
                        eprintln!("{} files, {:.0} s", r.len(), t0.elapsed().as_secs_f64());
                    }
                }
            });
        }
    });
    let mut rows = rows.into_inner().unwrap();
    rows.sort_by_key(|r| r.id.parse::<u64>().unwrap_or(u64::MAX));
    let mut out = String::from(
        "id\tformat\tstatus\treason\tn_raw_vertices\tn_raw_polygons\tn_vertices\tn_faces\tn_boundary_edges\tn_components\tchi\torient_repairable\tseconds\tdetail\n",
    );
    for r in &rows {
        let _ = writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{}",
            r.id,
            r.format,
            r.status,
            if r.status == "kept" { "-" } else { &r.reason },
            r.n_raw_vertices,
            r.n_raw_polygons,
            r.n_vertices,
            r.n_faces,
            r.n_boundary_edges,
            r.n_components,
            r.chi,
            r.orient_repairable,
            r.seconds,
            r.detail
        );
    }
    std::fs::write(&tsv, out).expect("write tsv");
    let kept = rows.iter().filter(|r| r.status == "kept").count();
    eprintln!(
        "{} files, {kept} kept, {:.0} s; table in {}",
        rows.len(),
        t0.elapsed().as_secs_f64(),
        tsv.display()
    );
}
