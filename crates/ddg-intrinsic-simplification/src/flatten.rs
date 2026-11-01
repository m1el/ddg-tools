//! Flattening a vertex (paper Sec. 4.1, Eq. 3, 4; SPEC §3.4, S9): the conformal factor `u`
//! of one vertex that makes its angle sum `2π` (inside) or `π` (boundary), found by Newton's
//! method on the extended angles of Springborn, Schröder and Pinkall (2008), safeguarded by a
//! bracket.

use crate::formulas::{face_angle_sum_and_derivative, projective_interpolation_signed};
use crate::mesh::{CoarseMesh, Location, face_lengths};

/// The largest `|u|` the bracket may expand to: lengths change by at most `e^{|u|/2}`
/// (`e^{30} ≈ 1e13`), beyond which nothing is described by the lengths any more.
pub const FLATTEN_MAX_FACTOR: f64 = 60.0;

/// A face of the star of a vertex, once per face: its id and which slots are the vertex.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct StarFace {
    pub(crate) face: usize,
    pub(crate) at: [bool; 3],
}

/// What a flattening solve found.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Solution {
    pub(crate) u: f64,
    pub(crate) iterations: usize,
    pub(crate) bisections: usize,
    pub(crate) residual: f64,
    /// Whether the star satisfies the strict triangle inequality at `u`.
    pub(crate) valid: bool,
}

/// What the flattening did (one per removal).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct FlattenReport {
    /// The conformal factor of the vertex (`ℓ̃ = e^{u/2} ℓ` on its edges, Eq. 3).
    pub u: f64,
    /// Steps (Newton or bisection) over all attempts.
    pub iterations: usize,
    /// Of these, bisection steps (Newton left the bracket or the slope was not negative).
    pub bisections: usize,
    /// Solves made: one plus the number of fallback flips.
    pub attempts: usize,
    /// `|Θ̃ − Θ|` at the end, in radians.
    pub residual: f64,
}

/// Why a solve failed: the best residual and the steps taken.
pub(crate) struct SolveFailure {
    pub(crate) residual: f64,
    pub(crate) iterations: usize,
}

impl CoarseMesh {
    /// The faces of the star of complex vertex `v`, each once (at the first of its slots that
    /// is `v`), in the order of the fan from the reference halfedge.
    pub(crate) fn star_faces(&self, v: usize) -> Vec<StarFace> {
        let c = &self.complex;
        let mut out = Vec::new();
        for h in self.fan(v) {
            let f = h / 3;
            let face = c.faces()[f];
            let at = face.map(|x| x as usize == v);
            let first = (0..3).find(|&k| at[k]).unwrap_or(0);
            if first == h % 3 {
                out.push(StarFace { face: f, at });
            }
        }
        out
    }

    /// The lengths of a star face scaled by the conformal factor `u` of the star's vertex
    /// (Eq. 3 with `u_j = 0`): `e^{u/2}` per end point of a side that is the vertex.
    #[inline]
    pub(crate) fn scaled_face_lengths(&self, sf: &StarFace, half: f64, full: f64) -> [f64; 3] {
        let l = face_lengths(&self.complex, &self.length, sf.face);
        [0, 1, 2].map(|k| match (sf.at[k], sf.at[(k + 1) % 3]) {
            (true, true) => l[k] * full,
            (true, false) | (false, true) => l[k] * half,
            (false, false) => l[k],
        })
    }

    /// `Θ(u)`, `dΘ/du` and whether every face is valid, for the star `faces`.
    fn star_angle_sum(&self, faces: &[StarFace], u: f64) -> Option<(f64, f64, bool)> {
        let (half, full) = ((0.5 * u).exp(), u.exp());
        let mut theta = 0.0;
        let mut d = 0.0;
        let mut valid = true;
        for sf in faces {
            let l = self.scaled_face_lengths(sf, half, full);
            let (t, dt, ok) = face_angle_sum_and_derivative(l, sf.at)?;
            theta += t;
            d += dt;
            valid &= ok;
        }
        Some((theta, d, valid))
    }

    /// Solves `Θ(u) = target` for the star `faces` (SPEC §3.4, S9): Newton's step (Eq. 4)
    /// from `u = 0`, kept inside a bracket `[lo, hi]` with `Θ(lo) ≥ target ≥ Θ(hi)` that grows
    /// geometrically while one side is open and is bisected whenever Newton leaves it or the
    /// slope is not negative. Stops at `|Θ − target| ≤ tol · target`; a bracket that shrinks
    /// to adjacent doubles is accepted if the residual is within `1000 tol · target`.
    pub(crate) fn solve_flattening(
        &self,
        faces: &[StarFace],
        target: f64,
    ) -> Result<Solution, SolveFailure> {
        let tol = self.options.flatten_tolerance * target;
        let max_iterations = self.options.flatten_max_iterations;
        let (mut lo, mut hi) = (f64::NEG_INFINITY, f64::INFINITY);
        let mut u = 0.0f64;
        let mut step = 1.0f64;
        let mut bisections = 0;
        let mut best = (f64::INFINITY, 0.0f64, false);
        let fail = |residual: f64, iterations: usize| SolveFailure {
            residual,
            iterations,
        };
        for iteration in 0..=max_iterations {
            let Some((theta, d, valid)) = self.star_angle_sum(faces, u) else {
                return Err(fail(best.0, iteration));
            };
            let g = theta - target;
            if g.abs() < best.0 {
                best = (g.abs(), u, valid);
            }
            if g.abs() <= tol {
                return Ok(Solution {
                    u,
                    iterations: iteration,
                    bisections,
                    residual: g.abs(),
                    valid,
                });
            }
            if iteration == max_iterations {
                break;
            }
            // Θ is non-increasing in u: a positive g needs a larger u.
            if g > 0.0 {
                lo = u;
            } else {
                hi = u;
            }
            if lo.is_finite()
                && hi.is_finite()
                && hi - lo <= 2.0 * f64::EPSILON * lo.abs().max(hi.abs()).max(1.0)
            {
                // The bracket is as small as doubles allow.
                if best.0 <= 1000.0 * tol {
                    let (_, _, valid) = self
                        .star_angle_sum(faces, best.1)
                        .ok_or_else(|| fail(best.0, iteration))?;
                    return Ok(Solution {
                        u: best.1,
                        iterations: iteration,
                        bisections,
                        residual: best.0,
                        valid,
                    });
                }
                return Err(fail(best.0, iteration));
            }
            let newton = if d < 0.0 { u - g / d } else { f64::NAN };
            let next = if newton > lo && newton < hi {
                newton
            } else if lo.is_finite() && hi.is_finite() {
                bisections += 1;
                0.5 * (lo + hi)
            } else if hi.is_infinite() {
                let x = lo + step;
                step *= 2.0;
                x
            } else {
                let x = hi - step;
                step *= 2.0;
                x
            };
            if next.is_nan() || next.abs() > FLATTEN_MAX_FACTOR {
                return Err(fail(best.0, iteration + 1));
            }
            u = next;
        }
        Err(fail(best.0, max_iterations))
    }

    /// The edges Springborn's fallback may flip for the star `faces` at the factor `u` (SPEC
    /// §3.4): in every scaled face that violates the triangle inequality, the longest side
    /// (the one opposite the angle `π`); spokes first, then outer edges, each group by edge id.
    pub(crate) fn fallback_candidates(&self, faces: &[StarFace], u: f64) -> Vec<usize> {
        let (half, full) = ((0.5 * u).exp(), u.exp());
        let mut spokes = Vec::new();
        let mut outer = Vec::new();
        for sf in faces {
            let l = self.scaled_face_lengths(sf, half, full);
            if ddg_core::geometry::triangle_area(l[0], l[1], l[2]).is_some() {
                continue;
            }
            let k = (0..3).max_by(|&a, &b| l[a].total_cmp(&l[b])).unwrap_or(0);
            let e = self.complex.edge(3 * sf.face + k);
            if sf.at[k] || sf.at[(k + 1) % 3] {
                spokes.push(e);
            } else if !self.complex.is_boundary_edge(e) {
                outer.push(e);
            }
        }
        for list in [&mut spokes, &mut outer] {
            list.sort_unstable();
            list.dedup();
        }
        spokes.retain(|e| !outer.contains(e));
        spokes.extend(outer);
        spokes
    }

    /// Applies the factor `u` to the star `faces` of complex vertex `v` (SPEC §3.4, "Commit"):
    /// scales its edges (Eq. 3), recomputes the angle sums and angular coordinates of `v` and
    /// of every vertex of the star from their references (S4), and, if `track`, moves the
    /// points of the star by Eq. 10. Writes are recorded.
    pub(crate) fn commit_flattening(&mut self, v: usize, faces: &[StarFace], u: f64, track: bool) {
        let (half, full) = ((0.5 * u).exp(), u.exp());
        let epoch = self.next_epoch();
        let mut vertices = Vec::new();
        for sf in faces {
            for k in 0..3 {
                let e = self.complex.edge(3 * sf.face + k);
                let x = self.complex.faces()[sf.face][k] as usize;
                if self.vertex_mark[x] != epoch {
                    self.vertex_mark[x] = epoch;
                    vertices.push(x);
                }
                if self.edge_mark[e] == epoch {
                    continue;
                }
                let factor = match (sf.at[k], sf.at[(k + 1) % 3]) {
                    (true, true) => full,
                    (true, false) | (false, true) => half,
                    (false, false) => continue,
                };
                self.edge_mark[e] = epoch;
                let l = self.length[e] * factor;
                self.set_length(e, l);
            }
        }
        debug_assert!(vertices.contains(&v));
        for x in vertices {
            self.refresh_vertex(x);
        }
        if track {
            for sf in faces {
                let n = self.face_points[sf.face].len();
                for slot in 0..n {
                    let p = self.face_points[sf.face][slot];
                    if let Location::Face { bary, .. } = self.points[p as usize] {
                        let b = projective_interpolation_signed(
                            bary,
                            sf.at,
                            u,
                            self.options.projective_sign,
                        );
                        self.set_bary(p, b);
                    }
                }
            }
        }
    }
}
