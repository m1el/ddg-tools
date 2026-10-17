//! Applications of intrinsic triangulations, from Sec. 5 of
//!
//! > Nicholas Sharp, Yousuf Soliman, Keenan Crane. *Navigating Intrinsic Triangulations.*
//! > ACM Transactions on Graphics 38(4), 2019. <https://doi.org/10.1145/3306346.3322979>
//!
//! built on the signpost data structure of `ddg-intrinsic`. The decisions the paper leaves
//! open are in `SPEC.md` next to this crate's `Cargo.toml`; the measurements in `README.md`.
//!
//! * [`amr`] (Sec. 5.4): adaptive mesh refinement driven by local a posteriori error
//!   estimates (Morin, Nochetto, Siebert 2002; Mekchay, Nochetto 2005), for the harmonic
//!   Green's function, the short-time heat kernel and a harmonic map of a disk.
//! * [`vector`] (Sec. 5.5): connection Laplacians on an intrinsic triangulation (finite
//!   differences, Sharp, Soliman, Crane 2019 "The Vector Heat Method", Sec. 5.3; finite
//!   elements, Knöppel, Crane, Pinkall, Schröder 2013, Sec. 6.1.1), harmonic vector fields
//!   with boundary data and the vector maximum principle, smoothest `n`-direction fields.

#![warn(missing_docs)]

pub mod amr;
pub mod vector;
