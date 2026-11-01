//! The undo log of a removal (SPEC §3.6, S5).
//!
//! While a removal is pending, every write to the mesh's lengths, angular coordinates,
//! reference halfedges, angle sums, connectivity (flips) and tracked points records the value
//! it overwrites. Reverting replays the records backwards, which restores every array bitwise
//! (including the order of the per-face point lists). Outside a removal nothing is recorded.

use ddg_intrinsic::complex::ComplexFlipUndo;

/// One overwritten value.
#[derive(Debug, Clone, Copy)]
pub(crate) enum Undo {
    /// The length of an edge.
    Length(u32, f64),
    /// The angular coordinate of a halfedge.
    Phi(u32, f64),
    /// The reference halfedge of a (complex) vertex.
    Reference(u32, u32),
    /// The cached angle sum of a (complex) vertex.
    AngleSum(u32, f64),
    /// A flip of the complex: index into [`UndoLog::flips`].
    Flip(u32),
    /// The barycentric coordinates of a point that stayed in its face.
    Bary(u32, [f64; 3]),
    /// A point that moved from `face` (where it was at `slot` of the list) to the end of
    /// another face's list, with its old coordinates.
    Moved {
        point: u32,
        face: u32,
        slot: u32,
        bary: [f64; 3],
    },
}

/// The records of the pending removal.
#[derive(Debug, Clone, Default)]
pub(crate) struct UndoLog {
    /// Whether writes are recorded.
    pub(crate) active: bool,
    pub(crate) entries: Vec<Undo>,
    pub(crate) flips: Vec<ComplexFlipUndo>,
}

impl UndoLog {
    /// Starts recording (forgetting earlier records).
    pub(crate) fn begin(&mut self) {
        self.active = true;
        self.entries.clear();
        self.flips.clear();
    }

    /// Stops recording and forgets the records.
    pub(crate) fn end(&mut self) {
        self.active = false;
        self.entries.clear();
        self.flips.clear();
    }

    #[inline]
    pub(crate) fn push(&mut self, entry: Undo) {
        if self.active {
            self.entries.push(entry);
        }
    }

    #[inline]
    pub(crate) fn push_flip(&mut self, undo: ComplexFlipUndo) {
        if self.active {
            self.entries.push(Undo::Flip(self.flips.len() as u32));
            self.flips.push(undo);
        }
    }
}
