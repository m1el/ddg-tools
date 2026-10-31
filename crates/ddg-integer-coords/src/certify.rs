//! The certificate (SPEC §2.4, §8 layer 1): checks that the integers are those of a
//! realization of the input edges over the triangulation.
//!
//! [`IntegerTriangulation::certify`] checks, in this order, and stops at the first violation:
//!
//! 0. **Ranges** of the stored data (§2.2): array sizes, `n ≥ −1`, `n = −1` on boundary edges,
//!    a roundabout in `[0, deg₀)` exactly on halfedges whose tail is an input vertex, input
//!    vertices at themselves and inserted vertices at a point of an input edge or face.
//! 1. **Item 1**: in every face the numerators of Eq. 1 are even (`ic-face-validity`).
//!    Item 2 (at most one corner of a face emanates) follows from Eq. 2 (Lean
//!    `emanating_exclusive`) and is not checked separately.
//! 2. **Item 3**: every edge with `n = −1` joins two distinct vertices, and runs along one input
//!    edge: the roundabouts of its halfedges that leave input vertices name it, its input-vertex
//!    ends are end points of it, its inserted ends lie on it.
//! 3. **Item 4**: around every input vertex the roundabouts agree with the emanating counts
//!    ([`IntegerTriangulation::input_directions`]); every input edge, followed from its tail
//!    (Alg. 4 with Alg. 1 and the continuation through inserted vertices), arrives at its head
//!    by the end that the roundabouts there give the opposite direction, passing only through
//!    inserted vertices that lie on it (in order along it, with exactly two curve ends each);
//!    every crossing and every edge with `n = −1` is used by exactly one input edge.
//!
//! Cost: linear in `|E¹| + Σ n⁺ + Σ_{V⁰} deg₀` (every crossing is visited once; a bitmap of
//! `Σ n⁺` bits marks them). Measured in the crate's README (SPEC I4).

use ddg_intrinsic::SurfacePoint;

use crate::coords::{IntegerTriangulation, NO_ROUNDABOUT};
use crate::error::{IntegerError, Place, ValidityItem};
use crate::formulas::{face_valid, n_plus};
use crate::trace::{Crossing, CurveEnd};

/// What [`IntegerTriangulation::certify`] checked.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CertifyReport {
    /// Number of faces checked (item 1).
    pub n_faces: usize,
    /// Number of edges with `n = −1` (item 3).
    pub n_along: usize,
    /// `Σ n⁺`: the number of crossings, each visited once (item 4).
    pub n_crossings: u64,
    /// Number of input edges extracted (item 4).
    pub n_input_edges: usize,
    /// Number of inserted vertices that lie on input edges (passed by extraction).
    pub n_inserted_on_edges: usize,
    /// The largest number of crossings on one input edge.
    pub max_crossings_per_input_edge: u64,
}

/// One bit per crossing and per edge with `n = −1`.
struct Marks {
    /// Per edge: the first crossing's bit.
    offsets: Vec<u64>,
    bits: Vec<u64>,
    along: Vec<bool>,
}

impl Marks {
    fn new(t: &IntegerTriangulation) -> Marks {
        let mut offsets = Vec::with_capacity(t.normal.len() + 1);
        let mut total = 0u64;
        for &n in &t.normal {
            offsets.push(total);
            total += n_plus(n) as u64;
        }
        offsets.push(total);
        Marks {
            offsets,
            bits: vec![0; total.div_ceil(64) as usize],
            along: vec![false; t.normal.len()],
        }
    }

    fn total(&self) -> u64 {
        *self
            .offsets
            .last()
            .expect("one offset per edge and one more")
    }

    /// Marks crossing `z`; an error if it was marked already.
    fn mark(&mut self, t: &IntegerTriangulation, z: Crossing) -> Result<(), IntegerError> {
        let h = z.halfedge as usize;
        let e = t.complex.edge(h);
        let n = n_plus(t.normal[e]);
        let p = if t.complex.is_forward(h) {
            z.index
        } else {
            n - 1 - z.index
        };
        let bit = self.offsets[e] + p as u64;
        let (w, b) = ((bit / 64) as usize, bit % 64);
        if self.bits[w] & (1 << b) != 0 {
            return Err(IntegerError::invalid(
                ValidityItem::Extraction,
                z.place(),
                "the crossing belongs to two input-edge curves (or one curve passes it twice)"
                    .to_string(),
            ));
        }
        self.bits[w] |= 1 << b;
        Ok(())
    }

    fn mark_along(&mut self, e: u32) -> Result<(), IntegerError> {
        if std::mem::replace(&mut self.along[e as usize], true) {
            return Err(IntegerError::invalid(
                ValidityItem::Extraction,
                Place::Edge(e),
                "the edge with n = −1 lies along two input-edge curves".to_string(),
            ));
        }
        Ok(())
    }

    /// The first crossing not marked, as `(edge, position from its halfedge's tail)`.
    fn first_unmarked(&self) -> Option<(usize, u64)> {
        let total = self.total();
        for (w, &word) in self.bits.iter().enumerate() {
            if word != u64::MAX {
                let bit = 64 * w as u64 + (!word).trailing_zeros() as u64;
                if bit >= total {
                    return None;
                }
                let e = self.offsets.partition_point(|&o| o <= bit) - 1;
                return Some((e, bit - self.offsets[e]));
            }
        }
        None
    }
}

impl IntegerTriangulation {
    /// Checks that the integer data is valid (SPEC §2.4, all four items; see the
    /// [module documentation](self) for what exactly), and says what failed and where if not.
    ///
    /// # Errors
    /// [`IntegerError::InvalidCoordinates`] with the item of SPEC §2.4 that fails and the place
    /// (face, edge, halfedge, vertex, crossing or input edge).
    pub fn certify(&self) -> Result<CertifyReport, IntegerError> {
        self.check_ranges()?;
        self.check_faces()?;
        let n_along = self.check_along()?;
        let mut report = self.check_extraction()?;
        report.n_faces = self.complex.n_faces();
        report.n_along = n_along;
        Ok(report)
    }

    fn check_ranges(&self) -> Result<(), IntegerError> {
        let ranges = |place: Place, detail: String| {
            IntegerError::invalid(ValidityItem::Ranges, place, detail)
        };
        let c = &self.complex;
        let n_input = self.n_input_vertices();
        let mesh = self.input.mesh();
        if self.normal.len() != c.n_edges()
            || self.lengths.len() != c.n_edges()
            || self.roundabouts.len() != c.n_halfedges()
            || self.locations.len() != c.n_vertices()
            || c.n_vertices() < n_input
        {
            return Err(ranges(
                Place::Vertex(0),
                "array lengths do not match the element counts".to_string(),
            ));
        }
        for (e, &n) in self.normal.iter().enumerate() {
            if n < -1 {
                return Err(ranges(
                    Place::Edge(e as u32),
                    format!("n = {n} is below −1"),
                ));
            }
            if c.is_boundary_edge(e) && n != -1 {
                return Err(ranges(
                    Place::Edge(e as u32),
                    format!("a boundary edge has n = {n}; boundary edges run along the input's"),
                ));
            }
        }
        for h in 0..c.n_halfedges() {
            let a = c.from_vertex(h);
            let r = self.roundabouts[h];
            if a < n_input {
                let deg = self.local.degree(a);
                if r == NO_ROUNDABOUT || r >= deg {
                    return Err(ranges(
                        Place::Halfedge(h as u32),
                        format!(
                            "roundabout {} at input vertex {a} is not in [0, {deg})",
                            if r == NO_ROUNDABOUT {
                                "none".to_string()
                            } else {
                                r.to_string()
                            }
                        ),
                    ));
                }
            } else if r != NO_ROUNDABOUT {
                return Err(ranges(
                    Place::Halfedge(h as u32),
                    format!("inserted vertex {a} has no input halfedges, but roundabout {r}"),
                ));
            }
        }
        for (v, &q) in self.locations.iter().enumerate() {
            let ok = match q {
                SurfacePoint::Vertex(u) => v < n_input && u as usize == v,
                SurfacePoint::Edge { edge, t } => {
                    v >= n_input && (edge as usize) < mesh.n_edges() && t > 0.0 && t < 1.0
                }
                SurfacePoint::Face { face, bary } => {
                    v >= n_input
                        && (face as usize) < mesh.n_faces()
                        && bary.iter().all(|&b| b > 0.0 && b < 1.0)
                }
            };
            if !ok {
                return Err(ranges(
                    Place::Vertex(v as u32),
                    format!(
                        "{} vertex placed at {q}",
                        if v < n_input { "input" } else { "inserted" }
                    ),
                ));
            }
        }
        for a in 0..n_input {
            if self.local.is_boundary(a) != c.is_boundary_vertex(a) {
                return Err(ranges(
                    Place::Vertex(a as u32),
                    "the vertex is on the boundary of one triangulation and not of the other"
                        .to_string(),
                ));
            }
        }
        Ok(())
    }

    /// Item 1.
    ///
    /// [formula: ic-face-validity]
    fn check_faces(&self) -> Result<(), IntegerError> {
        for f in 0..self.complex.n_faces() {
            if !face_valid(self.face_normal(f)) {
                return Err(self.face_coordinates(f).unwrap_err());
            }
        }
        Ok(())
    }

    /// Item 3; returns the number of edges with `n = −1`.
    fn check_along(&self) -> Result<usize, IntegerError> {
        let c = &self.complex;
        let mesh = self.input.mesh();
        let mut count = 0;
        for e in 0..c.n_edges() {
            if self.normal[e] != -1 {
                continue;
            }
            count += 1;
            let bad = |detail: String| {
                IntegerError::invalid(ValidityItem::AlongInputEdge, Place::Edge(e as u32), detail)
            };
            let h = c.edge_halfedge(e);
            let [x, y] = [c.from_vertex(h), c.to_vertex(h)];
            if x == y {
                return Err(bad(format!(
                    "a self-edge at vertex {x} cannot run along an input edge"
                )));
            }
            // The input edge named at each end: by the roundabout at an input vertex (the
            // last direction for the head of a boundary edge), by the place of an inserted one.
            let mut named: Option<usize> = None;
            for (v, out) in [(x, Some(h)), (y, c.twin(h))] {
                let input_edge = if self.is_input_vertex(v) {
                    let k = match out {
                        Some(g) => self.roundabouts[g],
                        None => self.local.degree(v) - 1,
                    };
                    let ie = self.local.edge(v, k);
                    let [p, q] = mesh.edge_vertices(ie);
                    if p != v && q != v {
                        return Err(bad(format!("input edge {ie} does not end at vertex {v}")));
                    }
                    ie
                } else {
                    match self.locations[v] {
                        SurfacePoint::Edge { edge, .. } => edge as usize,
                        other => {
                            return Err(bad(format!(
                                "inserted vertex {v} at its end lies at {other}, not on an input edge"
                            )));
                        }
                    }
                };
                match named {
                    None => named = Some(input_edge),
                    Some(ie) if ie != input_edge => {
                        return Err(bad(format!(
                            "its ends name different input edges: {ie} at vertex {x}, \
                             {input_edge} at vertex {v}"
                        )));
                    }
                    _ => {}
                }
            }
            // Both ends input vertices: the input edge joins them.
            if self.is_input_vertex(x) && self.is_input_vertex(y) {
                let ie = named.expect("both ends were looked at");
                let [p, q] = mesh.edge_vertices(ie);
                if !((p == x && q == y) || (p == y && q == x)) {
                    return Err(bad(format!(
                        "input edge {ie} joins {p} and {q}, not the ends {x} and {y}"
                    )));
                }
            }
        }
        Ok(count)
    }

    /// Item 4.
    fn check_extraction(&self) -> Result<CertifyReport, IntegerError> {
        let mesh = self.input.mesh();
        let n_input = self.n_input_vertices();
        let directions: Vec<Vec<CurveEnd>> = (0..n_input)
            .map(|a| self.input_directions(a))
            .collect::<Result<_, _>>()?;
        let mut marks = Marks::new(self);
        let mut n_inserted_on_edges = 0usize;
        let mut max_per_edge = 0u64;
        let mut passed = vec![false; self.complex.n_vertices() - n_input];
        for ie in 0..mesh.n_edges() {
            let h0 = mesh.edge_halfedge(ie);
            let (a, b) = (mesh.from_vertex(h0), mesh.to_vertex(h0));
            let k_b = self
                .local
                .of_edge_at(mesh, ie, b)
                .expect("an input edge ends at the head of its halfedge");
            let expected = directions[b][k_b as usize];
            let (mut v, mut end) = (a, directions[a][self.local.of_halfedge(h0) as usize]);
            let mut last_t = None;
            let mut on_edge = 0u64;
            loop {
                if let CurveEnd::Along { edge } = end {
                    marks.mark_along(edge)?;
                }
                let mut on_crossing = |z: Crossing| {
                    on_edge += 1;
                    marks.mark(self, z)
                };
                let (x, arrived) = self.follow(v, end, None, &mut on_crossing, None)?;
                if self.is_input_vertex(x) {
                    if x != b || arrived != expected {
                        return Err(IntegerError::invalid(
                            ValidityItem::Extraction,
                            Place::InputEdge(ie as u32),
                            format!(
                                "the curve from vertex {a} arrives at vertex {x} by {arrived:?}, \
                                 but the roundabouts at its head {b} give {expected:?}"
                            ),
                        ));
                    }
                    break;
                }
                let m = x - n_input;
                if std::mem::replace(&mut passed[m], true) {
                    return Err(IntegerError::invalid(
                        ValidityItem::Extraction,
                        Place::Vertex(x as u32),
                        format!("the curve of input edge {ie} passes this vertex a second time"),
                    ));
                }
                n_inserted_on_edges += 1;
                end = self.pass_inserted(x, ie, arrived)?;
                // In order along the input edge (its parameter runs from its lower vertex).
                if let SurfacePoint::Edge { t, .. } = self.locations[x] {
                    let t = if a == mesh.edge_vertices(ie)[0] {
                        t
                    } else {
                        1.0 - t
                    };
                    if last_t.is_some_and(|s| t <= s) {
                        return Err(IntegerError::invalid(
                            ValidityItem::Extraction,
                            Place::Vertex(x as u32),
                            format!(
                                "inserted vertices on input edge {ie} are out of order along it"
                            ),
                        ));
                    }
                    last_t = Some(t);
                }
                v = x;
            }
            max_per_edge = max_per_edge.max(on_edge);
        }
        if let Some((e, p)) = marks.first_unmarked() {
            return Err(IntegerError::invalid(
                ValidityItem::Extraction,
                Place::Crossing {
                    halfedge: self.complex.edge_halfedge(e) as u32,
                    index: p as i64,
                    face: (self.complex.edge_halfedge(e) / 3) as u32,
                },
                "the crossing belongs to no input edge (a closed curve or a stray arc)".to_string(),
            ));
        }
        for (e, &n) in self.normal.iter().enumerate() {
            if n == -1 && !marks.along[e] {
                return Err(IntegerError::invalid(
                    ValidityItem::Extraction,
                    Place::Edge(e as u32),
                    "the edge with n = −1 lies along no input edge".to_string(),
                ));
            }
        }
        // Inserted vertices: on an input edge exactly when a curve passes them.
        for (m, &was_passed) in passed.iter().enumerate() {
            let v = m + n_input;
            let on_edge = matches!(self.locations[v], SurfacePoint::Edge { .. });
            if on_edge != was_passed {
                return Err(IntegerError::invalid(
                    ValidityItem::Extraction,
                    Place::Vertex(v as u32),
                    format!(
                        "the vertex lies at {} but {} input edge passes through it",
                        self.locations[v],
                        if was_passed { "an" } else { "no" }
                    ),
                ));
            }
        }
        Ok(CertifyReport {
            n_faces: 0,
            n_along: 0,
            n_crossings: marks.total(),
            n_input_edges: mesh.n_edges(),
            n_inserted_on_edges,
            max_crossings_per_input_edge: max_per_edge,
        })
    }
}
