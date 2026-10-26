//! Automatic splitting of charts that overlap themselves ([`SplitOptions`]).
//!
//! The pipeline flattens every chart with the cones it chose. A chart that is nearly
//! developable gets no cones, nothing keeps its flattening from lying on top of itself, and
//! the overlap is global: two far apart parts of the surface land on the same place of the
//! plane. This module finds, for such a chart, a cut that separates the two layers, adds it to
//! the seams and runs the pipeline again.
//!
//! **Finding the cut.** The overlap estimate says which triangles land on shared area and which
//! pairs of triangles share it. The pair `(a, b)` with the largest `overlapped samples ×
//! surface distance` (among the most overlapped triangles) names the two layers; neighbours
//! that fold over each other are not separable and have a small distance. With `d_a`, `d_b` the
//! surface distances from `a` and `b` (Dijkstra over the dual graph, a step across an edge costs
//! the mean height of its two triangles over it), the faces with `d_a - d_b < -band · D` and
//! with `d_a - d_b > band · D`, `D = d_a(b)`, are the two sides, and the cut is the minimum cut
//! (Dinic) between them in the dual graph, where cutting a polygon edge costs `ℓ / (1 + crease
//! bias · dihedral angle)`, cutting a diagonal of the triangulation is impossible, and the
//! boundary of the chart is free. The cut is therefore a curve from boundary to boundary in the
//! band between the layers, and prefers short stretches and creases.
//!
//! **Keeping it.** After the pipeline ran again with the cuts as seams, the overlapped area of
//! the parts of a chart is compared with that of the chart; the cut stays only if it is lower by
//! `min_gain` and no triangle is mirrored that was not before. Otherwise the cut is undone and
//! the chart is not tried again. The parts of a kept cut are split again if they still overlap,
//! up to `max_splits` kept cuts. Everything is deterministic.

use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::time::Instant;

use ddg_bff::quality::signed_double_area;

use crate::cones::solve_cones;
use crate::feedback::cone_feedback;
use crate::error::UvError;
use crate::flatten::flatten;
use crate::input::MeshInput;
use crate::options::{SplitOptions, UnwrapOptions};
use crate::pipeline::{UvMap, assemble};
use crate::prepare::{Piece, prepare};
use crate::report::{SplitInfo, SplitOutcome, SplitReport};

const INFINITE_CAPACITY: f64 = 1e30;
const MAX_CELLS: f64 = 16_777_216.0;

/// Which triangles of a layout share area with another one.
struct Overlap {
    /// Per triangle: number of its sample points that another triangle also covers.
    per_face: Vec<u32>,
    /// Pairs `(a, b)`, `a < b`, with the number of sample points both cover, ascending.
    pairs: Vec<(u32, u32, u32)>,
}

/// The sampling of `ddg_bff::quality::coverage`, keeping who covers what.
fn sample_overlap(triangles: &[[[f64; 2]; 3]], samples_per_triangle: f64) -> Option<Overlap> {
    if triangles.is_empty() || samples_per_triangle.is_nan() || samples_per_triangle <= 0.0 {
        return None;
    }
    let (mut lo, mut hi) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
    let mut total = 0.0;
    for t in triangles {
        for p in t {
            if !(p[0].is_finite() && p[1].is_finite()) {
                return None;
            }
            for k in 0..2 {
                lo[k] = lo[k].min(p[k]);
                hi[k] = hi[k].max(p[k]);
            }
        }
        total += 0.5 * signed_double_area(t[0], t[1], t[2]).abs();
    }
    let (width, height) = (hi[0] - lo[0], hi[1] - lo[1]);
    if !(total > 0.0 && total.is_finite() && width > 0.0 && height > 0.0) {
        return None;
    }
    let mut h = (total / (samples_per_triangle * triangles.len() as f64)).sqrt();
    h = h.max((width * height / MAX_CELLS).sqrt());
    if !(h > 0.0 && h.is_finite()) {
        return None;
    }
    let nx = (width / h).ceil().max(1.0) as usize;
    let ny = (height / h).ceil().max(1.0) as usize;
    let mut owner = vec![u32::MAX; nx * ny];
    let mut per_face = vec![0u32; triangles.len()];
    let mut keys: Vec<u64> = Vec::new();
    let owns = |dx: f64, dy: f64| dy < 0.0 || (dy == 0.0 && dx < 0.0);
    for (index, t) in triangles.iter().enumerate() {
        let double_area = signed_double_area(t[0], t[1], t[2]);
        if double_area == 0.0 {
            continue;
        }
        let [a, b, c] = if double_area > 0.0 {
            *t
        } else {
            [t[0], t[2], t[1]]
        };
        let first = |min: f64, origin: f64, n: usize| {
            (((min - origin) / h - 0.5).ceil().max(0.0) as usize).min(n)
        };
        let end = |max: f64, origin: f64, n: usize| {
            ((((max - origin) / h - 0.5).floor() + 1.0).max(0.0) as usize).min(n)
        };
        let (x0, x1) = (
            first(a[0].min(b[0]).min(c[0]), lo[0], nx),
            end(a[0].max(b[0]).max(c[0]), lo[0], nx),
        );
        let (y0, y1) = (
            first(a[1].min(b[1]).min(c[1]), lo[1], ny),
            end(a[1].max(b[1]).max(c[1]), lo[1], ny),
        );
        let edges = [(a, b), (b, c), (c, a)];
        for j in y0..y1 {
            let y = lo[1] + (j as f64 + 0.5) * h;
            for i in x0..x1 {
                let p = [lo[0] + (i as f64 + 0.5) * h, y];
                let inside = edges.iter().all(|&(u, v)| {
                    let side = signed_double_area(u, v, p);
                    side > 0.0 || (side == 0.0 && owns(v[0] - u[0], v[1] - u[1]))
                });
                if inside {
                    let slot = &mut owner[j * nx + i];
                    if *slot == u32::MAX {
                        *slot = index as u32;
                    } else {
                        let first = *slot;
                        per_face[first as usize] += 1;
                        per_face[index] += 1;
                        let (lo_id, hi_id) = (first.min(index as u32), first.max(index as u32));
                        keys.push(((lo_id as u64) << 32) | hi_id as u64);
                    }
                }
            }
        }
    }
    keys.sort_unstable();
    let mut pairs: Vec<(u32, u32, u32)> = Vec::new();
    for key in keys {
        let (a, b) = ((key >> 32) as u32, key as u32);
        match pairs.last_mut() {
            Some(last) if last.0 == a && last.1 == b => last.2 += 1,
            _ => pairs.push((a, b, 1)),
        }
    }
    Some(Overlap { per_face, pairs })
}

#[derive(PartialEq)]
struct Entry(f64, usize);
impl Eq for Entry {}
impl Ord for Entry {
    fn cmp(&self, other: &Self) -> Ordering {
        // A min-heap on (distance, index).
        other
            .0
            .total_cmp(&self.0)
            .then_with(|| other.1.cmp(&self.1))
    }
}
impl PartialOrd for Entry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// The dual graph of a piece: `nbr[f][k]` is the face across the side `k` of face `f` (the
/// halfedge `3f + k`), `usize::MAX` on the boundary.
struct Dual {
    nbr: Vec<[usize; 3]>,
    /// Surface distance between the centres of the two faces.
    step: Vec<[f64; 3]>,
}

impl Dual {
    fn new(piece: &Piece) -> Dual {
        let n = piece.mesh.n_faces();
        let (areas, lengths) = (piece.geometry.face_areas(), piece.geometry.edge_lengths());
        let mut nbr = vec![[usize::MAX; 3]; n];
        let mut step = vec![[0.0; 3]; n];
        for f in 0..n {
            for k in 0..3 {
                let h = 3 * f + k;
                if let Some(t) = piece.mesh.twin(h) {
                    let g = piece.mesh.face(t);
                    let l = lengths[piece.mesh.edge(h)];
                    nbr[f][k] = g;
                    step[f][k] = if l > 0.0 {
                        (areas[f] + areas[g]) / l
                    } else {
                        0.0
                    };
                }
            }
        }
        Dual { nbr, step }
    }

    fn distances(&self, source: usize) -> Vec<f64> {
        let mut dist = vec![f64::INFINITY; self.nbr.len()];
        let mut heap = BinaryHeap::new();
        dist[source] = 0.0;
        heap.push(Entry(0.0, source));
        while let Some(Entry(d, f)) = heap.pop() {
            if d > dist[f] {
                continue;
            }
            for k in 0..3 {
                let g = self.nbr[f][k];
                if g == usize::MAX {
                    continue;
                }
                let nd = d + self.step[f][k];
                if nd < dist[g] {
                    dist[g] = nd;
                    heap.push(Entry(nd, g));
                }
            }
        }
        dist
    }
}

/// Dinic's maximum flow on an undirected capacitated graph.
struct Flow {
    head: Vec<usize>,
    next: Vec<usize>,
    to: Vec<usize>,
    cap: Vec<f64>,
    level: Vec<i32>,
    cursor: Vec<usize>,
}

impl Flow {
    const NONE: usize = usize::MAX;
    const EPS: f64 = 1e-12;

    fn new(n: usize) -> Flow {
        Flow {
            head: vec![Self::NONE; n],
            next: Vec::new(),
            to: Vec::new(),
            cap: Vec::new(),
            level: vec![-1; n],
            cursor: vec![Self::NONE; n],
        }
    }

    /// An arc `u -> v` of capacity `forward` and its reverse of capacity `backward`.
    fn add(&mut self, u: usize, v: usize, forward: f64, backward: f64) {
        for (a, b, c) in [(u, v, forward), (v, u, backward)] {
            self.next.push(self.head[a]);
            self.to.push(b);
            self.cap.push(c);
            self.head[a] = self.to.len() - 1;
        }
    }

    fn bfs(&mut self, s: usize, t: usize) -> bool {
        self.level.fill(-1);
        let mut queue = std::collections::VecDeque::new();
        self.level[s] = 0;
        queue.push_back(s);
        while let Some(u) = queue.pop_front() {
            let mut a = self.head[u];
            while a != Self::NONE {
                let v = self.to[a];
                if self.cap[a] > Self::EPS && self.level[v] < 0 {
                    self.level[v] = self.level[u] + 1;
                    queue.push_back(v);
                }
                a = self.next[a];
            }
        }
        self.level[t] >= 0
    }

    fn blocking(&mut self, s: usize, t: usize) -> f64 {
        self.cursor.copy_from_slice(&self.head);
        let mut total = 0.0;
        let mut path: Vec<usize> = Vec::new();
        let mut u = s;
        loop {
            if u == t {
                let bottleneck = path.iter().map(|&a| self.cap[a]).fold(f64::INFINITY, f64::min);
                for &a in &path {
                    self.cap[a] -= bottleneck;
                    self.cap[a ^ 1] += bottleneck;
                }
                total += bottleneck;
                path.clear();
                u = s;
                continue;
            }
            let mut found = Self::NONE;
            while self.cursor[u] != Self::NONE {
                let a = self.cursor[u];
                let v = self.to[a];
                if self.cap[a] > Self::EPS && self.level[v] == self.level[u] + 1 {
                    found = a;
                    break;
                }
                self.cursor[u] = self.next[a];
            }
            if found != Self::NONE {
                path.push(found);
                u = self.to[found];
            } else if let Some(a) = path.pop() {
                // Dead end: this node is useless in this phase.
                self.level[u] = -1;
                u = self.to[a ^ 1];
                self.cursor[u] = self.next[self.cursor[u]];
            } else {
                return total;
            }
        }
    }

    /// Runs the flow; returns its value and the nodes on the source side of a minimum cut.
    fn run(&mut self, s: usize, t: usize) -> (f64, Vec<bool>) {
        let mut total = 0.0;
        while self.bfs(s, t) {
            total += self.blocking(s, t);
            if total >= INFINITE_CAPACITY {
                break;
            }
        }
        // After the last BFS the reachable nodes are the source side.
        self.bfs(s, t);
        (total, self.level.iter().map(|&l| l >= 0).collect())
    }
}

/// A cut found for a chart.
struct Cut {
    /// Edges of the piece.
    edges: Vec<u32>,
    length: f64,
}

fn normal(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> [f64; 3] {
    let (u, v) = (
        [b[0] - a[0], b[1] - a[1], b[2] - a[2]],
        [c[0] - a[0], c[1] - a[1], c[2] - a[2]],
    );
    let n = [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ];
    let l = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
    if l > 0.0 {
        [n[0] / l, n[1] / l, n[2] / l]
    } else {
        [0.0; 3]
    }
}

/// How many of the most overlapped triangles are tried as the first layer.
const SEED_CANDIDATES: usize = 8;

/// Finds the cut that separates the two layers of the overlap of the piece `on`, whose corner
/// UVs are `uvs`; see the module documentation. `None` if the layers are neighbours or no
/// finite cut exists.
fn find_cut(
    positions: &[[f64; 3]],
    on: &Piece,
    uvs: &[[f64; 2]],
    samples: f64,
    options: &SplitOptions,
) -> Option<Cut> {
    let n = on.mesh.n_faces();
    let triangles: Vec<[[f64; 2]; 3]> = (0..n)
        .map(|f| [uvs[3 * f], uvs[3 * f + 1], uvs[3 * f + 2]])
        .collect();
    let overlap = sample_overlap(&triangles, samples)?;
    let dual = Dual::new(on);
    let total_area: f64 = on.geometry.face_areas().iter().sum();
    let min_distance = 3.0 * (total_area / n as f64).sqrt();

    let mut partners: Vec<Vec<(u32, u32)>> = vec![Vec::new(); n];
    for &(a, b, c) in &overlap.pairs {
        partners[a as usize].push((b, c));
        partners[b as usize].push((a, c));
    }
    let mut order: Vec<usize> = (0..n).filter(|&f| overlap.per_face[f] > 0).collect();
    order.sort_by(|&x, &y| overlap.per_face[y].cmp(&overlap.per_face[x]).then(x.cmp(&y)));

    // The two layers.
    let mut chosen = None;
    for &a in order.iter().take(SEED_CANDIDATES) {
        let da = dual.distances(a);
        let best = partners[a]
            .iter()
            .filter(|&&(b, _)| da[b as usize] >= min_distance)
            .max_by(|x, y| {
                let (wx, wy) = (
                    x.1 as f64 * da[x.0 as usize],
                    y.1 as f64 * da[y.0 as usize],
                );
                wx.total_cmp(&wy).then(y.0.cmp(&x.0))
            });
        if let Some(&(b, _)) = best {
            chosen = Some((a, b as usize, da));
            break;
        }
    }
    let (a, b, da) = chosen?;
    let db = dual.distances(b);
    let band = options.band * da[b];

    // Normals for the dihedral angles.
    let p = |v: u32| positions[on.vertex_input[v as usize] as usize];
    let normals: Vec<[f64; 3]> = on
        .mesh
        .faces()
        .iter()
        .map(|f| normal(p(f[0]), p(f[1]), p(f[2])))
        .collect();
    let lengths = on.geometry.edge_lengths();

    let (source, sink) = (n, n + 1);
    let mut flow = Flow::new(n + 2);
    let mut interior: Vec<(usize, usize, u32, bool)> = Vec::new();
    for h in 0..on.mesh.n_halfedges() {
        let Some(t) = on.mesh.twin(h) else { continue };
        if h > t {
            continue;
        }
        let (f, g) = (on.mesh.face(h), on.mesh.face(t));
        let e = on.mesh.edge(h);
        let cuttable = on.is_polygon_edge[e];
        let cap = if cuttable {
            let cos = (0..3)
                .map(|k| normals[f][k] * normals[g][k])
                .sum::<f64>()
                .clamp(-1.0, 1.0);
            lengths[e] / (1.0 + options.crease_bias * cos.acos())
        } else {
            INFINITE_CAPACITY
        };
        flow.add(f, g, cap, cap);
        interior.push((f, g, e as u32, cuttable));
    }
    let (mut n_a, mut n_b) = (0, 0);
    for f in 0..n {
        let d = da[f] - db[f];
        if d < -band {
            flow.add(source, f, INFINITE_CAPACITY, 0.0);
            n_a += 1;
        } else if d > band {
            flow.add(f, sink, INFINITE_CAPACITY, 0.0);
            n_b += 1;
        }
    }
    if n_a == 0 || n_b == 0 || a == b {
        return None;
    }
    let (value, side) = flow.run(source, sink);
    if value.is_nan() || value >= 0.5 * INFINITE_CAPACITY {
        return None;
    }
    let mut edges = Vec::new();
    let mut length = 0.0;
    for &(f, g, e, cuttable) in &interior {
        if side[f] != side[g] {
            if !cuttable {
                return None;
            }
            edges.push(e);
            length += lengths[e as usize];
        }
    }
    (!edges.is_empty()).then_some(Cut { edges, length })
}

/// A cut that was made and awaits the verdict of the next run.
struct Pending {
    /// The polygons of the chart before the cut, ascending.
    polygons: Vec<u32>,
    /// The edges added to the seams.
    edges: Vec<[u32; 2]>,
    before_area: f64,
    before_flipped: usize,
    /// Index into the report's `splits`.
    info: usize,
}

fn sorted_pair(a: u32, b: u32) -> [u32; 2] {
    [a.min(b), a.max(b)]
}

/// [`crate::unwrap`] with the splitting of overlapping charts switched on.
///
/// # Errors
/// As [`crate::unwrap`].
pub(crate) fn unwrap_splitting(
    input: &MeshInput<'_>,
    options: &UnwrapOptions,
    split: &SplitOptions,
) -> Result<UvMap, UvError> {
    let start = Instant::now();
    options.validate()?;
    let mut extra: Vec<[u32; 2]> = Vec::new();
    let mut pending: Vec<Pending> = Vec::new();
    let mut frozen: Vec<Vec<u32>> = Vec::new();
    let mut report = SplitReport::default();
    let max_runs = 2 * split.max_splits + 4;
    let mut kept = 0usize;
    let mut seams: Vec<[u32; 2]>;
    loop {
        seams = input.seam_edges.to_vec();
        seams.extend_from_slice(&extra);
        let run_input = MeshInput {
            seam_edges: &seams,
            ..*input
        };
        let prepared = prepare(&run_input, options)?;
        let mut cones = solve_cones(&prepared, options)?;
        let mut flat = flatten(&prepared, &cones, options)?;
        cone_feedback(&prepared, &mut cones, &mut flat, options)?;
        report.runs += 1;
        let n_pieces = prepared.pieces.len();
        let on = |i: usize| cones.pieces[i].slit_piece.as_ref().unwrap_or(&prepared.pieces[i]);
        let piece_polygons: Vec<Vec<u32>> = (0..n_pieces)
            .map(|i| {
                let mut p: Vec<u32> = on(i)
                    .face_triangle
                    .iter()
                    .map(|&t| prepared.triangle_polygon[t as usize])
                    .collect();
                p.sort_unstable();
                p.dedup();
                p
            })
            .collect();
        let area = |i: usize| flat.pieces[i].coverage.map_or((0.0, 0.0), |c| c);
        let fraction = |i: usize| {
            let (o, c) = area(i);
            if c > 0.0 { o / c } else { 0.0 }
        };
        if report.runs == 1 {
            let (o, c) = (0..n_pieces)
                .map(area)
                .fold((0.0, 0.0), |s, x| (s.0 + x.0, s.1 + x.1));
            report.overlap_before = (c > 0.0).then(|| o / c);
            report.max_chart_overlap_before = (0..n_pieces)
                .filter(|&i| flat.pieces[i].coverage.is_some())
                .map(fraction)
                .fold(None, |b: Option<f64>, x| Some(b.map_or(x, |b| b.max(x))));
        }
        let mut piece_of_polygon = vec![u32::MAX; prepared.n_polygons];
        for (i, polygons) in piece_polygons.iter().enumerate() {
            for &p in polygons {
                piece_of_polygon[p as usize] = i as u32;
            }
        }

        // The verdict on the cuts of the last run.
        let mut changed = false;
        for cut in pending.drain(..) {
            let mut parts: Vec<u32> = cut
                .polygons
                .iter()
                .map(|&p| piece_of_polygon[p as usize])
                .collect();
            parts.sort_unstable();
            parts.dedup();
            let (mut over, mut covered, mut flipped) = (0.0, 0.0, 0);
            for &i in &parts {
                let (o, c) = area(i as usize);
                over += o;
                covered += c;
                flipped += flat.pieces[i as usize].report.flipped_triangles;
            }
            let info = &mut report.splits[cut.info];
            info.n_parts = parts.len();
            info.overlap_after = if covered > 0.0 { over / covered } else { 0.0 };
            if over <= (1.0 - split.min_gain) * cut.before_area && flipped <= cut.before_flipped {
                info.outcome = SplitOutcome::Kept;
                report.n_kept += 1;
                report.cut_length += info.cut_length;
                kept += 1;
            } else {
                info.outcome = SplitOutcome::Rejected;
                report.n_rejected += 1;
                extra.retain(|e| !cut.edges.contains(e));
                frozen.push(cut.polygons);
                changed = true;
            }
        }

        // New cuts.
        let mut exhausted = false;
        if report.runs < max_runs {
            #[allow(clippy::needless_range_loop)]
            for i in 0..n_pieces {
                let fraction_i = fraction(i);
                if flat.pieces[i].coverage.is_none() || fraction_i <= split.min_overlap_fraction {
                    continue;
                }
                if frozen.contains(&piece_polygons[i]) {
                    continue;
                }
                if kept + pending.len() >= split.max_splits {
                    exhausted = true;
                    break;
                }
                let piece = on(i);
                let found = find_cut(
                    input.positions,
                    piece,
                    &flat.pieces[i].corner_uvs,
                    options.overlap_samples_per_triangle,
                    split,
                );
                let first_polygon = piece_polygons[i][0];
                match found {
                    None => {
                        report.n_no_cut += 1;
                        report.splits.push(SplitInfo {
                            first_polygon,
                            round: report.runs - 1,
                            n_cut_edges: 0,
                            cut_length: 0.0,
                            overlap_before: fraction_i,
                            overlap_after: 0.0,
                            n_parts: 1,
                            outcome: SplitOutcome::NoCut,
                        });
                        frozen.push(piece_polygons[i].clone());
                    }
                    Some(cut) => {
                        let edges: Vec<[u32; 2]> = cut
                            .edges
                            .iter()
                            .map(|&e| {
                                let [a, b] = piece.mesh.edge_vertices(e as usize);
                                sorted_pair(piece.vertex_input[a], piece.vertex_input[b])
                            })
                            .collect();
                        report.splits.push(SplitInfo {
                            first_polygon,
                            round: report.runs - 1,
                            n_cut_edges: edges.len(),
                            cut_length: cut.length,
                            overlap_before: fraction_i,
                            overlap_after: 0.0,
                            n_parts: 1,
                            outcome: SplitOutcome::Rejected,
                        });
                        extra.extend_from_slice(&edges);
                        pending.push(Pending {
                            polygons: piece_polygons[i].clone(),
                            edges,
                            before_area: area(i).0,
                            before_flipped: flat.pieces[i].report.flipped_triangles,
                            info: report.splits.len() - 1,
                        });
                        changed = true;
                    }
                }
            }
        }
        if !changed {
            // Overlapping charts that were left alone for want of budget.
            report.exhausted = exhausted
                || (0..n_pieces).any(|i| {
                    flat.pieces[i].coverage.is_some()
                        && fraction(i) > split.min_overlap_fraction
                        && !frozen.contains(&piece_polygons[i])
                });
            extra.sort_unstable();
            extra.dedup();
            let mut map = assemble(&run_input, &prepared, &cones, &flat, options)?;
            if report.exhausted {
                map.report.warnings.push(format!(
                    "overlapping charts remain after {} kept cuts (split_overlaps.max_splits \
                     {})",
                    report.n_kept, split.max_splits
                ));
            }
            map.report.splitting = Some(report);
            map.report.timings.total = start.elapsed().as_secs_f64();
            return Ok(map);
        }
        extra.sort_unstable();
        extra.dedup();
    }
}
