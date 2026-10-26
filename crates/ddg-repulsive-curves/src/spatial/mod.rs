//! Spatial data structures: bounding boxes, BVH, block cluster tree (SPEC §8).
//!
//! [`aabb`] and [`bvh`] (work package WP4) are generic in the dimension and the element type and
//! contain no curve types; [`bct`] is filled by WP5.

pub mod aabb;
pub mod bct;
pub mod bvh;

pub use aabb::Aabb;
pub use bvh::{Bvh, Element, Node};
