//! SimInf-core: Pure Rust implementation of the Gillespie SSA solver
//! for stochastic disease spread simulations.

mod model;
mod solver;
mod events;
mod trajectory;

pub use model::{Model, ModelBuilder, Compartment, CompartmentId, Transition, PropensityFn};
pub use trajectory::TrajectoryResult;
pub use solver::Solver;