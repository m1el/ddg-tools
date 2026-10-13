//! The linear system of the cut mesh on a fixed set of vertices of the input mesh.
//!
//! The Laplacian of the cut mesh, restricted to vertices of the input mesh, has an entry for
//! `(i, j)` only if `ij` is an edge of the input mesh that is not crossed (SPEC §3.4, remark).
//! Its pattern is therefore a subset of the pattern of the input mesh's Laplacian, whatever the
//! cut: one symbolic analysis serves every cut, and only the values change. Entries of crossed
//! edges are stored as explicit zeros.

use ddg_core::sparse::{Factorization, Method, Parallelism, Symbolic};
use ddg_core::topology::NONE;
use ddg_core::{SparseMatrix, TriMesh};

use crate::cut_mesh::CutTriangles;
use crate::error::CutError;

/// The cotan Laplacian of a cut mesh on a subset `U` of the vertices of the input mesh (the
/// *unknowns*), with every other vertex of the cut mesh held at a given value: the matrix
/// `L_UU`, its symbolic analysis, and its current numeric factorization.
#[derive(Debug)]
pub(crate) struct LaplaceSystem {
    /// For every vertex of the input mesh, its index among the unknowns, or `NONE`.
    unknown_of_vertex: Vec<u32>,
    /// For every unknown, its vertex.
    vertex_of_unknown: Vec<u32>,
    /// For every unknown, the position of its diagonal entry in the matrix values.
    diagonal_slot: Vec<u32>,
    /// For every edge of the input mesh, the position of its entry in the matrix values, or
    /// `NONE` unless both ends are unknowns.
    edge_slot: Vec<u32>,
    /// The lower triangle of `L_UU` (the symmetric factorizations read nothing else).
    matrix: SparseMatrix,
    /// `None` if there are no unknowns.
    symbolic: Option<Symbolic>,
    /// The factorization of the values assembled last, once [`LaplaceSystem::factor`] ran.
    factorization: Option<Factorization>,
    /// Whether `factorization` belongs to the current values of `matrix`.
    factorized: bool,
}

impl LaplaceSystem {
    /// Analyses the pattern for the unknowns selected by `is_unknown`.
    pub(crate) fn new(
        mesh: &TriMesh,
        is_unknown: impl Fn(usize) -> bool,
    ) -> Result<LaplaceSystem, CutError> {
        let nv = mesh.n_vertices();
        let mut unknown_of_vertex = vec![NONE; nv];
        let mut vertex_of_unknown = Vec::new();
        for (v, slot) in unknown_of_vertex.iter_mut().enumerate() {
            if is_unknown(v) {
                *slot = vertex_of_unknown.len() as u32;
                vertex_of_unknown.push(v as u32);
            }
        }
        let n = vertex_of_unknown.len();
        let mut triplets: Vec<(u32, u32, f64)> = Vec::with_capacity(n + mesh.n_edges());
        for i in 0..n as u32 {
            triplets.push((i, i, 0.0));
        }
        for &[a, b] in mesh.edges() {
            let (i, j) = (unknown_of_vertex[a as usize], unknown_of_vertex[b as usize]);
            if i != NONE && j != NONE {
                triplets.push((i.max(j), i.min(j), 0.0));
            }
        }
        let matrix = SparseMatrix::from_triplets(n, n, &triplets)?;
        let slot = |row: u32, col: u32| -> u32 {
            matrix
                .find(row as usize, col as usize)
                .expect("the entry was inserted above") as u32
        };
        let diagonal_slot = (0..n as u32).map(|i| slot(i, i)).collect();
        let edge_slot = mesh
            .edges()
            .iter()
            .map(|&[a, b]| {
                let (i, j) = (unknown_of_vertex[a as usize], unknown_of_vertex[b as usize]);
                if i != NONE && j != NONE {
                    slot(i.max(j), i.min(j))
                } else {
                    NONE
                }
            })
            .collect();
        let symbolic = if n > 0 {
            Some(Symbolic::analyze(&matrix, Method::Llt)?)
        } else {
            None
        };
        Ok(LaplaceSystem {
            unknown_of_vertex,
            vertex_of_unknown,
            diagonal_slot,
            edge_slot,
            matrix,
            symbolic,
            factorization: None,
            factorized: false,
        })
    }

    /// Number of unknowns.
    #[inline]
    pub(crate) fn n_unknowns(&self) -> usize {
        self.vertex_of_unknown.len()
    }

    /// For every vertex of the input mesh, its index among the unknowns, or `NONE`.
    #[inline]
    pub(crate) fn unknown_of_vertex(&self) -> &[u32] {
        &self.unknown_of_vertex
    }

    /// For every unknown, its vertex of the input mesh.
    #[inline]
    pub(crate) fn vertex_of_unknown(&self) -> &[u32] {
        &self.vertex_of_unknown
    }

    /// Whether vertex `v` of the cut mesh is an unknown.
    #[inline]
    pub(crate) fn is_unknown(&self, v: usize) -> bool {
        v < self.unknown_of_vertex.len() && self.unknown_of_vertex[v] != NONE
    }

    /// The symbolic analysis (`None` if there are no unknowns).
    #[inline]
    pub(crate) fn symbolic(&self) -> Option<&Symbolic> {
        self.symbolic.as_ref()
    }

    /// The factorization of the matrix assembled last, if [`LaplaceSystem::factor`] succeeded
    /// since.
    #[inline]
    pub(crate) fn factorization(&self) -> Option<&Factorization> {
        if self.factorized {
            self.factorization.as_ref()
        } else {
            None
        }
    }

    /// Assembles `L_UU` of the cut mesh `tri`: the cotan Laplacian restricted to the unknowns
    /// (SPEC §3.4). `pinned` flags vertices of the input mesh whose unknown is taken out of the
    /// system: its row and column become those of the identity. A pinned vertex must not share
    /// a triangle with an unknown that is not pinned (it is used for whole components).
    ///
    /// [formula: cut-yamabe-system]
    pub(crate) fn assemble(&mut self, mesh: &TriMesh, tri: &CutTriangles, pinned: &[bool]) {
        self.factorized = false;
        let nv = tri.n_original;
        let unknown = |v: u32| -> u32 {
            if (v as usize) < nv && !pinned[v as usize] {
                self.unknown_of_vertex[v as usize]
            } else {
                NONE
            }
        };
        let mut values = vec![0.0f64; self.matrix.nnz()];
        for (t, (face, cot)) in tri.faces.iter().zip(&tri.cot).enumerate() {
            // The edge opposite corner k runs from corner k + 1 to corner k + 2.
            for k in 0..3 {
                let (i, j) = (face[(k + 1) % 3], face[(k + 2) % 3]);
                let (ui, uj) = (unknown(i), unknown(j));
                if ui == NONE && uj == NONE {
                    continue;
                }
                let w = 0.5 * cot[k];
                if ui != NONE {
                    values[self.diagonal_slot[ui as usize] as usize] += w;
                }
                if uj != NONE {
                    values[self.diagonal_slot[uj as usize] as usize] += w;
                }
                if ui != NONE && uj != NONE {
                    let h = tri
                        .original_halfedge(mesh, t, i, j)
                        .expect("two input vertices of a cut-mesh triangle span an input edge");
                    values[self.edge_slot[mesh.edge(h)] as usize] -= w;
                }
            }
        }
        for (v, &is_pinned) in pinned.iter().enumerate() {
            if is_pinned && self.unknown_of_vertex[v] != NONE {
                values[self.diagonal_slot[self.unknown_of_vertex[v] as usize] as usize] = 1.0;
            }
        }
        self.matrix.values_mut().copy_from_slice(&values);
    }

    /// Factorizes the assembled matrix, reusing the symbolic analysis and the memory of the
    /// previous factorization.
    pub(crate) fn factor(&mut self, par: Parallelism<'_>) -> Result<(), CutError> {
        let Some(symbolic) = &self.symbolic else {
            self.factorized = true;
            return Ok(());
        };
        match &mut self.factorization {
            Some(factorization) => factorization.refactor(&self.matrix, par)?,
            None => self.factorization = Some(symbolic.factor(&self.matrix, par)?),
        }
        self.factorized = true;
        Ok(())
    }

    /// Solves `L_UU x = b` with the current factorization. `rhs` gives `b` per vertex of the cut
    /// mesh (only unknowns are read); the solution is written to the unknown entries of `x`,
    /// which is indexed the same way. Other entries of `x` are left alone.
    pub(crate) fn solve(
        &self,
        rhs: impl Fn(usize) -> f64,
        x: &mut [f64],
        par: Parallelism<'_>,
    ) -> Result<(), CutError> {
        let Some(factorization) = self.factorization() else {
            return if self.n_unknowns() == 0 {
                Ok(())
            } else {
                Err(CutError::Sparse(ddg_core::SparseError::NotFactorized))
            };
        };
        let mut b: Vec<f64> = self
            .vertex_of_unknown
            .iter()
            .map(|&v| rhs(v as usize))
            .collect();
        factorization.solve_in_place(&mut b, par)?;
        for (&v, value) in self.vertex_of_unknown.iter().zip(b) {
            x[v as usize] = value;
        }
        Ok(())
    }
}
