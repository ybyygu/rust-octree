#[cfg(test)]
#[macro_use]
extern crate approx;

mod octant;
mod query;
mod octree;

// set up aliases for convenience
pub use self::octree::Octree;
