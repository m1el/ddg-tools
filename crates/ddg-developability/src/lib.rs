//! Developability of triangle meshes.
//!
//! An implementation, from the text of the paper, of
//!
//! > Oded Stein, Eitan Grinspun, Keenan Crane. *Developability of Triangle Meshes.*
//! > ACM Transactions on Graphics 37(4), 2018. <https://doi.org/10.1145/3197517.3201303>
//!
//! Contract: `SPEC.md` next to this crate's `Cargo.toml`. This crate is MIT licensed; it contains
//! no code from the authors' GPL implementation or from the LGPL C code of Kopp's eigensolvers.
//!
//! The eigenvalue energy `E_λ` of the paper (Sec. 4.3) needs an accurate 3x3 symmetric
//! eigensolver (the paper cites J. Kopp, *Efficient numerical diagonalization of hermitian 3x3
//! matrices*, Int. J. Mod. Phys. C 19 (2008) 523-548, arXiv:physics/0610206). [`eig`] provides it,
//! implemented here with cyclic Jacobi rotations (a textbook method).

pub mod eig;
pub mod energy;
pub mod error;
pub mod geometry;
pub mod optimize;
pub mod remesh;
pub mod subdivide;

pub use energy::{Energy, EnergyOptions, Evaluation, evaluate, evaluate_energy};
pub use error::DevError;
