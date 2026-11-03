//! A simplicial sparse `LDLᵀ` of the preconditioner's `L′` with triangular solves in a fixed,
//! documented operation order (GPU-2, `NOTES-GPU2.md`), so that the CPU and the GPU
//! (`src/gpu/kernels.cu`, `ldl_*` kernels) compute the same solves bit for bit.
//!
//! **Factorization.** faer's sparse `LDLᵀ` without pivoting (as `ddg_core::sparse` uses it) of
//! `P A Pᵀ`, `P` the caller's elimination order followed by a postorder of its elimination tree
//! (same fill, same factor up to round-off; with the postorder a node with a single child has that
//! child just before it, so the chains of the tree are index ranges). The factor is read out as the
//! strictly lower `L` by columns (CSC, rows increasing) and by rows (CSR, columns increasing), and
//! `D`.
//!
//! **Solves** (`X` is `n × m`, row-major in the permuted numbering; column `c` of `X` is one
//! right-hand side and depends only on itself):
//!
//! ```text
//! forward    for j = 0 … n−1:   x_j = b_j − S_j,   S_j = Σ_{k ∈ row j of L} L_jk x_k
//! diagonal   x_j = x_j / d_j
//! backward   for j = n−1 … 0:   x_j = x_j − T_j,   T_j = Σ_{i ∈ column j of L} L_ij x_i
//! ```
//!
//! The terms of a sum are taken in a fixed sequence (forward: `k` increasing; backward: `i`
//! decreasing) and cut into consecutive chunks of [`CHUNK`]; each chunk is summed from `0.0` in
//! sequence, `p = 0; p += l·x`, and the chunk sums are subtracted in sequence, `s = b; s −= p`.
//! Every value is a fixed function of the factor and the right-hand side whatever the schedule:
//! the GPU computes the chunks of rows whose terms are known in parallel and reproduces it
//! (`src/gpu/kernels.cu`, `ldl_*`).

use ddg_core::SparseMatrix;

use crate::error::RsError;

/// Terms per chunk of a triangular-solve sum (module docs).
pub const CHUNK: usize = 32;

const NONE: u32 = u32::MAX;

/// A simplicial `LDLᵀ` factor (module docs).
#[derive(Debug, Clone)]
pub struct SimplicialLdlt {
    n: usize,
    /// `perm[new] = old`.
    perm: Vec<u32>,
    /// Elimination tree in the new numbering (`u32::MAX`: a root).
    parent: Vec<u32>,
    /// `L` by columns.
    lp: Vec<usize>,
    li: Vec<u32>,
    lx: Vec<f64>,
    /// `L` by rows.
    rp: Vec<usize>,
    rj: Vec<u32>,
    rx: Vec<f64>,
    d: Vec<f64>,
}

/// Facts about the factor's structure (for the notes and the GPU schedule).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LdltStats {
    /// Dimension.
    pub n: usize,
    /// Entries of the strictly lower `L`.
    pub nnz_l: usize,
    /// Height of the elimination tree (nodes on the longest leaf-to-root path).
    pub height: usize,
    /// Longest row of `L`.
    pub max_row: usize,
    /// Longest column of `L`.
    pub max_col: usize,
}

fn zero_pivot(step: usize) -> RsError {
    RsError::SingularSaddle { pivot: step }
}

type Factor = (Vec<usize>, Vec<u32>, Vec<f64>, Vec<f64>);

/// faer's sparse `LDLᵀ` of `P A Pᵀ` (`perm[new] = old`, the order is used as given), read out
/// as the strictly lower `L` by columns (rows increasing, every stored entry kept, explicit zeros
/// of relaxed supernodes too) and `D`.
fn faer_ldlt(a: &SparseMatrix, perm: &[u32], pinv: &[u32]) -> Result<Factor, RsError> {
    use faer::dyn_stack::{MemBuffer, MemStack};
    use faer::linalg::cholesky::ldlt::factor::LdltRegularization;
    use faer::sparse::linalg::cholesky as fchol;
    use faer::sparse::{SparseColMatRef, SymbolicSparseColMatRef};
    let n = a.nrows();
    let oom = |_| RsError::InvalidOption {
        name: "preconditioner factorization",
        value: format!("dimension {n}"),
        requirement: "enough memory for the sparse factor",
    };
    let col_ptr = a.col_ptr().to_vec();
    let row_idx: Vec<usize> = a.row_idx().iter().map(|&r| r as usize).collect();
    // SAFETY: `SparseMatrix` keeps monotone column pointers and sorted, in-range, unique rows.
    let pattern = unsafe { SymbolicSparseColMatRef::<usize>::new_unchecked(n, n, &col_ptr, None, &row_idx) };
    let fwd: Vec<usize> = perm.iter().map(|&i| i as usize).collect();
    let inv: Vec<usize> = pinv.iter().map(|&i| i as usize).collect();
    // SAFETY: `perm` and `pinv` are inverse permutations of 0..n (checked by the caller).
    let ord = fchol::SymmetricOrdering::Custom(unsafe { faer::perm::PermRef::new_unchecked(&fwd, &inv, n) });
    let sym = fchol::factorize_symbolic_cholesky(pattern, faer::Side::Lower, ord, Default::default())
        .map_err(|e| RsError::InvalidOption {
            name: "preconditioner analysis",
            value: format!("{e:?}"),
            requirement: "a factorizable preconditioner Laplacian",
        })?;
    let mut values = vec![0.0f64; sym.len_val()];
    let mut mem = MemBuffer::try_new(sym.factorize_numeric_ldlt_scratch::<f64>(faer::Par::Seq, Default::default()))
        .map_err(oom)?;
    sym.factorize_numeric_ldlt::<f64>(
        &mut values,
        SparseColMatRef::<usize, f64>::new(pattern, a.values()),
        faer::Side::Lower,
        LdltRegularization::default(),
        faer::Par::Seq,
        MemStack::new(&mut mem),
        Default::default(),
    )
    .map_err(|e| match e {
        faer::linalg::cholesky::ldlt::factor::LdltError::ZeroPivot { index } => zero_pivot(index),
    })?;
    let mut d = vec![0.0f64; n];
    let mut lp = vec![0usize; n + 1];
    let mut li: Vec<u32>;
    let mut lx: Vec<f64>;
    match sym.raw() {
        fchol::SymbolicCholeskyRaw::Simplicial(s) => {
            let (cp, ri) = (s.col_ptr(), s.row_idx());
            for j in 0..n {
                lp[j + 1] = lp[j] + ri[cp[j]..cp[j + 1]].iter().filter(|&&i| i > j).count();
            }
            li = Vec::with_capacity(lp[n]);
            lx = Vec::with_capacity(lp[n]);
            for j in 0..n {
                for p in cp[j]..cp[j + 1] {
                    let i = ri[p];
                    if i == j {
                        d[j] = values[p];
                    } else if i > j {
                        li.push(i as u32);
                        lx.push(values[p]);
                    }
                }
            }
        }
        fchol::SymbolicCholeskyRaw::Supernodal(s) => {
            let (begin, end) = (s.supernode_begin(), s.supernode_end());
            let (cpr, cpv, ri) = (s.col_ptr_for_row_idx(), s.col_ptr_for_val(), s.row_idx());
            for k in 0..s.n_supernodes() {
                let (s0, s1) = (begin[k], end[k]);
                let np = cpr[k + 1] - cpr[k];
                for j in s0..s1 {
                    lp[j + 1] = lp[j] + (s1 - j - 1) + np;
                }
            }
            li = Vec::with_capacity(lp[n]);
            lx = Vec::with_capacity(lp[n]);
            for k in 0..s.n_supernodes() {
                let (s0, s1) = (begin[k], end[k]);
                let pat = &ri[cpr[k]..cpr[k + 1]];
                let w = s1 - s0;
                let rows = w + pat.len();
                let block = &values[cpv[k]..cpv[k + 1]];
                for c in 0..w {
                    let col = &block[c * rows..(c + 1) * rows];
                    d[s0 + c] = col[c];
                    for r in c + 1..w {
                        li.push((s0 + r) as u32);
                        lx.push(col[r]);
                    }
                    for (q, &i) in pat.iter().enumerate() {
                        li.push(i as u32);
                        lx.push(col[w + q]);
                    }
                }
            }
        }
    }
    // Rows increasing within every column (sort a column only if needed).
    for j in 0..n {
        let r = lp[j]..lp[j + 1];
        if !li[r.clone()].is_sorted() {
            let mut e: Vec<(u32, f64)> = li[r.clone()].iter().copied().zip(lx[r.clone()].iter().copied()).collect();
            e.sort_unstable_by_key(|t| t.0);
            for (q, (i, x)) in r.zip(e) {
                li[q] = i;
                lx[q] = x;
            }
        }
    }
    if let Some(j) = d.iter().position(|v| *v == 0.0 || !v.is_finite()) {
        return Err(zero_pivot(j));
    }
    Ok((lp, li, lx, d))
}

impl SimplicialLdlt {
    /// Factors the symmetric `a` (both triangles stored; only entries with permuted
    /// `row ≤ col` are read) in the elimination order `order` (`order[step] = index`), refined by
    /// a postorder of the elimination tree.
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] for a non-square `a` or an order that is not a permutation;
    /// [`RsError::SingularSaddle`] for a zero or non-finite pivot.
    pub fn factor(a: &SparseMatrix, order: &[u32]) -> Result<SimplicialLdlt, RsError> {
        let n = a.nrows();
        let bad = |what: String| RsError::InvalidOption {
            name: "LDLT elimination order",
            value: what,
            requirement: "a square matrix and a permutation of its indices",
        };
        if a.ncols() != n || order.len() != n {
            return Err(bad(format!("{} x {} matrix, {} order entries", n, a.ncols(), order.len())));
        }
        let mut pinv = vec![NONE; n];
        for (s, &i) in order.iter().enumerate() {
            if i as usize >= n || pinv[i as usize] != NONE {
                return Err(bad(format!("entry {s} = {i}")));
            }
            pinv[i as usize] = s as u32;
        }
        // Elimination tree of the given order, then its postorder (children in increasing order).
        let parent0 = etree(a, order, &pinv);
        let post = postorder(&parent0);
        let perm: Vec<u32> = post.iter().map(|&s| order[s as usize]).collect();
        for (s, &i) in perm.iter().enumerate() {
            pinv[i as usize] = s as u32;
        }
        let (lp, li, lx, d) = faer_ldlt(a, &perm, &pinv)?;
        let nnz = li.len();
        // The elimination tree of the stored pattern: the first sub-diagonal row of a column.
        let parent: Vec<u32> = (0..n).map(|j| if lp[j] < lp[j + 1] { li[lp[j]] } else { NONE }).collect();
        // Rows of L (CSR, columns increasing).
        let mut rp = vec![0usize; n + 1];
        for &i in &li {
            rp[i as usize + 1] += 1;
        }
        for k in 0..n {
            rp[k + 1] += rp[k];
        }
        let mut next = rp.clone();
        let mut rj = vec![0u32; nnz];
        let mut rx = vec![0.0f64; nnz];
        for j in 0..n {
            for p in lp[j]..lp[j + 1] {
                let i = li[p] as usize;
                rj[next[i]] = j as u32;
                rx[next[i]] = lx[p];
                next[i] += 1;
            }
        }
        Ok(SimplicialLdlt {
            n,
            perm,
            parent,
            lp,
            li,
            lx,
            rp,
            rj,
            rx,
            d,
        })
    }

    /// Dimension.
    pub fn dim(&self) -> usize {
        self.n
    }

    /// `perm[new] = old`.
    pub fn perm(&self) -> &[u32] {
        &self.perm
    }

    /// The elimination tree (new numbering; `u32::MAX` for roots).
    pub fn parent(&self) -> &[u32] {
        &self.parent
    }

    /// `L` by columns: pointers, rows, values.
    pub fn columns(&self) -> (&[usize], &[u32], &[f64]) {
        (&self.lp, &self.li, &self.lx)
    }

    /// `L` by rows: pointers, columns, values.
    pub fn rows(&self) -> (&[usize], &[u32], &[f64]) {
        (&self.rp, &self.rj, &self.rx)
    }

    /// `D`.
    pub fn diagonal(&self) -> &[f64] {
        &self.d
    }

    /// Structure facts.
    pub fn stats(&self) -> LdltStats {
        let n = self.n;
        let mut h = vec![1usize; n];
        let mut height = if n > 0 { 1 } else { 0 };
        for j in 0..n {
            let p = self.parent[j];
            if p != NONE {
                h[p as usize] = h[p as usize].max(h[j] + 1);
                height = height.max(h[p as usize]);
            }
        }
        let row = |i: usize| self.rp[i + 1] - self.rp[i];
        let col = |j: usize| self.lp[j + 1] - self.lp[j];
        LdltStats {
            n,
            nnz_l: self.li.len(),
            height,
            max_row: (0..n).map(row).max().unwrap_or(0),
            max_col: (0..n).map(col).max().unwrap_or(0),
        }
    }

    /// The solve of the columns `c0 .. c0 + W` of the `n × m` row-major `x` (module docs).
    fn solve_group<const W: usize>(&self, x: &mut [f64], m: usize, c0: usize) {
        let n = self.n;
        // Forward, by rows (terms k increasing).
        for j in 0..n {
            let (p0, p1) = (self.rp[j], self.rp[j + 1]);
            let mut s: [f64; W] = std::array::from_fn(|c| x[j * m + c0 + c]);
            for (ks, ls) in self.rj[p0..p1].chunks(CHUNK).zip(self.rx[p0..p1].chunks(CHUNK)) {
                let mut cur = [0.0f64; W];
                for (&k, &l) in ks.iter().zip(ls) {
                    let base = k as usize * m + c0;
                    for c in 0..W {
                        cur[c] += l * x[base + c];
                    }
                }
                for c in 0..W {
                    s[c] -= cur[c];
                }
            }
            x[j * m + c0..j * m + c0 + W].copy_from_slice(&s);
        }
        for j in 0..n {
            let dj = self.d[j];
            for c in 0..W {
                x[j * m + c0 + c] /= dj;
            }
        }
        // Backward, by columns (terms i decreasing: chunks from the end of the column).
        for j in (0..n).rev() {
            let (p0, p1) = (self.lp[j], self.lp[j + 1]);
            let mut s: [f64; W] = std::array::from_fn(|c| x[j * m + c0 + c]);
            for (is, ls) in self.li[p0..p1].rchunks(CHUNK).zip(self.lx[p0..p1].rchunks(CHUNK)) {
                let mut cur = [0.0f64; W];
                for (&i, &l) in is.iter().zip(ls).rev() {
                    let base = i as usize * m + c0;
                    for c in 0..W {
                        cur[c] += l * x[base + c];
                    }
                }
                for c in 0..W {
                    s[c] -= cur[c];
                }
            }
            x[j * m + c0..j * m + c0 + W].copy_from_slice(&s);
        }
    }

    /// Solves `P A Pᵀ X = B` in place for `m` right-hand sides (`x`: `n × m` row-major, permuted
    /// numbering), in the operation order of the module docs.
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] for a wrong length; [`RsError::SingularSaddle`] for a
    /// non-finite solution.
    pub fn solve_rows(&self, x: &mut [f64], m: usize) -> Result<(), RsError> {
        let n = self.n;
        if x.len() != n * m {
            return Err(RsError::InvalidOption {
                name: "LDLT right-hand side",
                value: format!("{} entries for {n} x {m}", x.len()),
                requirement: "n x m entries",
            });
        }
        let mut c0 = 0;
        while c0 < m {
            match m - c0 {
                1 => self.solve_group::<1>(x, m, c0),
                2 => self.solve_group::<2>(x, m, c0),
                3 => self.solve_group::<3>(x, m, c0),
                _ => self.solve_group::<4>(x, m, c0),
            }
            c0 += 4;
        }
        if let Some(i) = x.iter().position(|v| !v.is_finite()) {
            return Err(zero_pivot(i / m.max(1)));
        }
        Ok(())
    }
}

/// Elimination tree of `P A Pᵀ` (`pinv[old] = new`, `order[new] = old`).
fn etree(a: &SparseMatrix, order: &[u32], pinv: &[u32]) -> Vec<u32> {
    let n = a.nrows();
    let (cp, ri) = (a.col_ptr(), a.row_idx());
    let mut parent = vec![NONE; n];
    let mut ancestor = vec![NONE; n];
    for k in 0..n {
        let col = order[k] as usize;
        for p in cp[col]..cp[col + 1] {
            let mut i = pinv[ri[p] as usize] as usize;
            if i >= k {
                continue;
            }
            // Path compression (Liu).
            while ancestor[i] != NONE && ancestor[i] != k as u32 {
                let next = ancestor[i] as usize;
                ancestor[i] = k as u32;
                i = next;
            }
            if ancestor[i] == NONE {
                ancestor[i] = k as u32;
                parent[i] = k as u32;
            }
        }
    }
    parent
}

/// A postorder of the forest `parent` (`post[new] = old`), children in increasing order, roots in
/// increasing order.
fn postorder(parent: &[u32]) -> Vec<u32> {
    let n = parent.len();
    let mut head = vec![NONE; n];
    let mut next = vec![NONE; n];
    // Children lists in increasing order: insert in decreasing order at the head.
    for j in (0..n).rev() {
        let p = parent[j];
        if p != NONE {
            next[j] = head[p as usize];
            head[p as usize] = j as u32;
        }
    }
    let mut post = Vec::with_capacity(n);
    let mut stack: Vec<u32> = Vec::new();
    for root in 0..n {
        if parent[root] != NONE {
            continue;
        }
        stack.push(root as u32);
        while let Some(&top) = stack.last() {
            let child = head[top as usize];
            if child == NONE {
                stack.pop();
                post.push(top);
            } else {
                head[top as usize] = next[child as usize];
                stack.push(child);
            }
        }
    }
    post
}

/// Longest chain of the GPU schedule (longer chains are cut into pieces; `kernels.cu` `LDL_MAXC`).
pub const MAX_CHAIN: usize = 1024;
/// Rows per thread of a chain block (`kernels.cu` `LDL_ITEMS`).
pub const CHAIN_ITEMS: usize = 4;

/// One sweep (forward or backward) of [`LdltSchedule`].
#[derive(Debug, Clone, Default)]
pub struct Sweep {
    /// Chains grouped by level (`chains[level_ptr[l]..level_ptr[l + 1]]`), in processing order.
    pub level_ptr: Vec<usize>,
    /// Chain ids.
    pub chains: Vec<u32>,
    /// Threads per block of every level (a multiple of 32, at most 256).
    pub level_block: Vec<u32>,
    /// The chunk slots of every level (`slot_*[slot_ptr[l]..slot_ptr[l + 1]]`).
    pub slot_ptr: Vec<usize>,
    /// Row (forward) or column (backward) of a slot.
    pub slot_row: Vec<u32>,
    /// Chunk index of a slot within its row.
    pub slot_chunk: Vec<u32>,
    /// Per row: the number of leading chunks whose terms all lie outside the row's chain.
    pub nfull: Vec<u32>,
    /// Per row: the slot of its first chunk.
    pub slot_base: Vec<u32>,
}

/// The GPU schedule of the solves of a [`SimplicialLdlt`] (`NOTES-GPU2.md`): the elimination tree
/// cut into chains (index ranges `j, j + 1, …` with `parent(j) = j + 1` and `j + 1` having no other
/// child; at most [`MAX_CHAIN`] long), grouped into levels (forward: by height, every chain after
/// the chains it reads; backward: by depth). Per level the chunks of terms from earlier levels are
/// summed in parallel (`ldl_partials`), then one block per chain runs the chain in sequence
/// (`ldl_chains`), adding the in-chain terms as they become known. The result is that of
/// [`SimplicialLdlt::solve_rows`] bit for bit ([`SimplicialLdlt::solve_scheduled`] runs the same
/// schedule on the CPU).
#[derive(Debug, Clone)]
pub struct LdltSchedule {
    /// First index of every chain.
    pub chain_start: Vec<u32>,
    /// One past the last index of every chain.
    pub chain_end: Vec<u32>,
    /// Forward sweep (rows of `L`).
    pub fwd: Sweep,
    /// Backward sweep (columns of `L`, rows decreasing: [`LdltSchedule::ci`], [`LdltSchedule::cx`]).
    pub bwd: Sweep,
    /// Rows of every column of `L`, decreasing (pointers: [`SimplicialLdlt::columns`]).
    pub ci: Vec<u32>,
    /// The values of [`LdltSchedule::ci`].
    pub cx: Vec<f64>,
    /// Number of chunk slots of the forward and backward sweeps.
    pub n_slots: (usize, usize),
}

fn group(levels: &[usize], nlev: usize, chain_len: &[usize]) -> (Vec<usize>, Vec<u32>, Vec<u32>) {
    let mut ptr = vec![0usize; nlev + 1];
    for &l in levels {
        ptr[l + 1] += 1;
    }
    for l in 0..nlev {
        ptr[l + 1] += ptr[l];
    }
    let mut next = ptr.clone();
    let mut chains = vec![0u32; levels.len()];
    for (c, &l) in levels.iter().enumerate() {
        chains[next[l]] = c as u32;
        next[l] += 1;
    }
    let block = (0..nlev)
        .map(|l| {
            let longest = chains[ptr[l]..ptr[l + 1]].iter().map(|&c| chain_len[c as usize]).max().unwrap_or(1);
            (longest.div_ceil(CHAIN_ITEMS).div_ceil(32) * 32).clamp(32, 256) as u32
        })
        .collect();
    (ptr, chains, block)
}

impl SimplicialLdlt {
    /// The GPU schedule of the solves ([`LdltSchedule`]).
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] if an entry of `L` does not respect the tree (a factor whose
    /// pattern is not closed; not expected from a Cholesky factorization).
    pub fn schedule(&self) -> Result<LdltSchedule, RsError> {
        let n = self.n;
        let mut nchild = vec![0u32; n];
        for &p in &self.parent {
            if p != NONE {
                nchild[p as usize] += 1;
            }
        }
        let mut chain_start: Vec<u32> = Vec::new();
        let mut chain_of = vec![0u32; n];
        for j in 0..n {
            let joined = j > 0
                && self.parent[j - 1] == j as u32
                && nchild[j] == 1
                && j - (*chain_start.last().expect("a chain") as usize) < MAX_CHAIN;
            if !joined {
                chain_start.push(j as u32);
            }
            chain_of[j] = (chain_start.len() - 1) as u32;
        }
        let nch = chain_start.len();
        let chain_end: Vec<u32> = (0..nch)
            .map(|c| chain_start.get(c + 1).copied().unwrap_or(n as u32))
            .collect();
        let chain_len: Vec<usize> = (0..nch).map(|c| (chain_end[c] - chain_start[c]) as usize).collect();
        let parent_chain = |c: usize| -> Option<usize> {
            let p = self.parent[chain_end[c] as usize - 1];
            (p != NONE).then(|| chain_of[p as usize] as usize)
        };
        let mut height = vec![0usize; nch];
        for c in 0..nch {
            if let Some(pc) = parent_chain(c) {
                height[pc] = height[pc].max(height[c] + 1);
            }
        }
        let mut depth = vec![0usize; nch];
        for c in (0..nch).rev() {
            if let Some(pc) = parent_chain(c) {
                depth[c] = depth[pc] + 1;
            }
        }
        for j in 0..n {
            for p in self.lp[j]..self.lp[j + 1] {
                let i = self.li[p] as usize;
                let (ci, cj) = (chain_of[i] as usize, chain_of[j] as usize);
                if ci != cj && !(height[cj] < height[ci] && depth[ci] < depth[cj]) {
                    return Err(RsError::InvalidOption {
                        name: "LDLT schedule",
                        value: format!("entry ({i}, {j}) of L"),
                        requirement: "a factor pattern closed under its elimination tree",
                    });
                }
            }
        }
        let nlev_h = height.iter().max().map_or(0, |h| h + 1);
        let nlev_d = depth.iter().max().map_or(0, |d| d + 1);
        let mut sweeps = Vec::new();
        for (levels, nlev, fwd) in [(&height, nlev_h, true), (&depth, nlev_d, false)] {
            let (level_ptr, chains, level_block) = group(levels, nlev, &chain_len);
            let mut nfull = vec![0u32; n];
            let mut slot_base = vec![0u32; n];
            let mut slot_ptr = vec![0usize];
            let mut slot_row = Vec::new();
            let mut slot_chunk = Vec::new();
            for l in 0..nlev {
                for &c in &chains[level_ptr[l]..level_ptr[l + 1]] {
                    let (c0, c1) = (chain_start[c as usize], chain_end[c as usize]);
                    for i in c0..c1 {
                        let iu = i as usize;
                        let n_ext = if fwd {
                            self.rj[self.rp[iu]..self.rp[iu + 1]].iter().filter(|&&k| k < c0).count()
                        } else {
                            self.li[self.lp[iu]..self.lp[iu + 1]].iter().filter(|&&k| k >= c1).count()
                        };
                        let nf = n_ext / CHUNK;
                        nfull[iu] = nf as u32;
                        slot_base[iu] = slot_row.len() as u32;
                        for q in 0..nf {
                            slot_row.push(i);
                            slot_chunk.push(q as u32);
                        }
                    }
                }
                slot_ptr.push(slot_row.len());
            }
            sweeps.push(Sweep {
                level_ptr,
                chains,
                level_block,
                slot_ptr,
                slot_row,
                slot_chunk,
                nfull,
                slot_base,
            });
        }
        let bwd = sweeps.pop().expect("two sweeps");
        let fwd = sweeps.pop().expect("two sweeps");
        let mut ci = Vec::with_capacity(self.li.len());
        let mut cx = Vec::with_capacity(self.lx.len());
        for j in 0..n {
            for p in (self.lp[j]..self.lp[j + 1]).rev() {
                ci.push(self.li[p]);
                cx.push(self.lx[p]);
            }
        }
        if u32::try_from(self.li.len()).is_err() || u32::try_from(fwd.slot_row.len().max(bwd.slot_row.len())).is_err() {
            return Err(RsError::InvalidOption {
                name: "LDLT schedule",
                value: format!("{} entries of L", self.li.len()),
                requirement: "fewer than 2^32 entries",
            });
        }
        let n_slots = (fwd.slot_row.len(), bwd.slot_row.len());
        Ok(LdltSchedule {
            chain_start,
            chain_end,
            fwd,
            bwd,
            ci,
            cx,
            n_slots,
        })
    }

    /// [`SimplicialLdlt::solve_rows`] by the GPU schedule `s`, on the CPU (the reference of the
    /// kernels `ldl_partials` and `ldl_chains`; bitwise the same result as `solve_rows`).
    ///
    /// # Errors
    /// As [`SimplicialLdlt::solve_rows`].
    pub fn solve_scheduled(&self, s: &LdltSchedule, x: &mut [f64], m: usize) -> Result<(), RsError> {
        let n = self.n;
        if x.len() != n * m {
            return Err(RsError::InvalidOption {
                name: "LDLT right-hand side",
                value: format!("{} entries for {n} x {m}", x.len()),
                requirement: "n x m entries",
            });
        }
        let mut w = x.to_vec();
        let mut z = vec![0.0; n * m];
        for fwd in [true, false] {
            let (sw, ptr, idx, val) = if fwd {
                (&s.fwd, &self.rp, &self.rj, &self.rx)
            } else {
                (&s.bwd, &self.lp, &s.ci, &s.cx)
            };
            let mut part = vec![0.0; sw.slot_row.len() * m];
            for l in 0..sw.level_ptr.len() - 1 {
                let data = &mut w;
                // ldl_partials
                for slot in sw.slot_ptr[l]..sw.slot_ptr[l + 1] {
                    let row = sw.slot_row[slot] as usize;
                    let p0 = ptr[row] + sw.slot_chunk[slot] as usize * CHUNK;
                    for c in 0..m {
                        let mut cur = 0.0;
                        for p in p0..p0 + CHUNK {
                            cur += val[p] * data[idx[p] as usize * m + c];
                        }
                        part[slot * m + c] = cur;
                    }
                }
                // ldl_chains
                for &ch in &sw.chains[sw.level_ptr[l]..sw.level_ptr[l + 1]] {
                    let (c0, c1) = (s.chain_start[ch as usize] as usize, s.chain_end[ch as usize] as usize);
                    for c in 0..m {
                        let mut st: Vec<(f64, f64, usize, usize, usize)> = (c0..c1)
                            .map(|i| {
                                let mut sv = data[i * m + c];
                                let (nf, base) = (sw.nfull[i] as usize, sw.slot_base[i] as usize);
                                for q in 0..nf {
                                    sv -= part[(base + q) * m + c];
                                }
                                let (mut p, e) = (ptr[i] + nf * CHUNK, ptr[i + 1]);
                                let (mut cu, mut cn) = (0.0, 0);
                                while p < e && {
                                    let k = idx[p] as usize;
                                    if fwd { k < c0 } else { k >= c1 }
                                } {
                                    cu += val[p] * data[idx[p] as usize * m + c];
                                    cn += 1;
                                    p += 1;
                                }
                                (sv, cu, cn, p, e)
                            })
                            .collect();
                        let order: Vec<usize> = if fwd { (c0..c1).collect() } else { (c0..c1).rev().collect() };
                        for &j in &order {
                            let r = j - c0;
                            let (mut sv, cu, cn, _, _) = st[r];
                            if cn > 0 {
                                sv -= cu;
                            }
                            data[j * m + c] = sv;
                            if fwd {
                                z[j * m + c] = sv / self.d[j];
                            }
                            for e in st.iter_mut() {
                                let (sv2, cu2, cn2, p2, e2) = e;
                                if *p2 < *e2 && idx[*p2] as usize == j {
                                    *cu2 += val[*p2] * sv;
                                    *cn2 += 1;
                                    *p2 += 1;
                                    if *cn2 == CHUNK {
                                        *sv2 -= *cu2;
                                        *cu2 = 0.0;
                                        *cn2 = 0;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if fwd {
                std::mem::swap(&mut w, &mut z);
            }
        }
        x.copy_from_slice(&w);
        if let Some(i) = x.iter().position(|v| !v.is_finite()) {
            return Err(zero_pivot(i / m.max(1)));
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn laplacian_1d(n: usize) -> SparseMatrix {
        let mut t = Vec::new();
        for i in 0..n {
            t.push((i as u32, i as u32, 2.5));
            if i + 1 < n {
                t.push((i as u32, i as u32 + 1, -1.0));
                t.push((i as u32 + 1, i as u32, -1.0));
            }
        }
        SparseMatrix::from_triplets(n, n, &t).unwrap()
    }

    #[test]
    fn solves_saddle() {
        // A 1-D Laplacian block with a dense constraint row (a long row) and a short one.
        let nv = 1500;
        let k = laplacian_1d(nv);
        let mut t: Vec<(u32, u32, f64)> = k.triplets().collect();
        for v in 0..nv as u32 {
            t.push((nv as u32, v, 1.0));
            t.push((v, nv as u32, 1.0));
        }
        t.push((nv as u32 + 1, 7, 1.0));
        t.push((7, nv as u32 + 1, 1.0));
        let n = nv + 2;
        let a = SparseMatrix::from_triplets(n, n, &t).unwrap();
        let order: Vec<u32> = (0..n as u32).rev().filter(|&i| i < nv as u32).chain([nv as u32, nv as u32 + 1]).collect();
        let f = SimplicialLdlt::factor(&a, &order).unwrap();
        assert_eq!(f.stats().max_row, nv);
        let m = 2;
        let xs: Vec<f64> = (0..n * m).map(|i| ((i * 31) % 17) as f64 - 8.0).collect();
        // b = A x in the original numbering, then permuted rows.
        let mut b = vec![0.0; n * m];
        for c in 0..m {
            let xc: Vec<f64> = (0..n).map(|i| xs[i * m + c]).collect();
            let mut y = vec![0.0; n];
            a.mul_vec(&xc, &mut y).unwrap();
            for (new, &old) in f.perm().iter().enumerate() {
                b[new * m + c] = y[old as usize];
            }
        }
        let mut bs = b.clone();
        f.solve_rows(&mut b, m).unwrap();
        let sched = f.schedule().unwrap();
        assert!(sched.chain_start.len() > 1);
        f.solve_scheduled(&sched, &mut bs, m).unwrap();
        assert!(b.iter().zip(&bs).all(|(u, v)| u.to_bits() == v.to_bits()));
        for (new, &old) in f.perm().iter().enumerate() {
            for c in 0..m {
                let e = (b[new * m + c] - xs[old as usize * m + c]).abs();
                assert!(e < 1e-8, "{new} {c}: {e}");
            }
        }
    }

    #[test]
    fn schedule_matches_sequential_on_a_grid() {
        // A 2-D grid Laplacian (branching elimination tree) with two dense rows.
        let (w, h) = (60usize, 45usize);
        let nv = w * h;
        let mut t = Vec::new();
        for y in 0..h {
            for x in 0..w {
                let v = (y * w + x) as u32;
                t.push((v, v, 4.0 + 1e-3 * (v % 7) as f64));
                if x + 1 < w {
                    t.push((v, v + 1, -1.0));
                    t.push((v + 1, v, -1.0));
                }
                if y + 1 < h {
                    t.push((v, v + w as u32, -1.0 - 1e-2 * (v % 3) as f64));
                    t.push((v + w as u32, v, -1.0 - 1e-2 * (v % 3) as f64));
                }
                t.push((nv as u32, v, 1.0 / nv as f64));
                t.push((v, nv as u32, 1.0 / nv as f64));
            }
        }
        t.push((nv as u32 + 1, 100, 1.0));
        t.push((100, nv as u32 + 1, 1.0));
        let n = nv + 2;
        let a = SparseMatrix::from_triplets(n, n, &t).unwrap();
        let kb: Vec<u32> = (0..nv as u32).collect();
        let amd = ddg_core::sparse::fill_reducing_ordering(&a.submatrix(&kb, &kb).unwrap()).unwrap();
        let order: Vec<u32> = amd.into_iter().chain([nv as u32, nv as u32 + 1]).collect();
        let f = SimplicialLdlt::factor(&a, &order).unwrap();
        let sched = f.schedule().unwrap();
        assert!(sched.fwd.level_ptr.len() > 5);
        for m in [1usize, 3, 6] {
            let b: Vec<f64> = (0..n * m).map(|i| ((i * 7919) % 1013) as f64 / 500.0 - 1.0).collect();
            let (mut x1, mut x2) = (b.clone(), b.clone());
            f.solve_rows(&mut x1, m).unwrap();
            f.solve_scheduled(&sched, &mut x2, m).unwrap();
            assert!(x1.iter().zip(&x2).all(|(u, v)| u.to_bits() == v.to_bits()), "m = {m}");
        }
    }
}
