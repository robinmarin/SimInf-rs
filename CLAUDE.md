# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Git

Never add `AGENTS.md` or `CLAUDE.md` to git commits.

## Commands

```bash
# Build
cargo build
cargo build --release

# Run all tests (including validation against R fixtures)
cargo test --workspace

# Run a single test by name
cargo test --package siminf-core test_sir_conservation_of_individuals

# Run the SIR example (writes trajectory.csv to repo root)
cargo run --example sir --package siminf-core

# Build WebAssembly bindings (requires wasm-pack)
cd siminf-wasm && wasm-pack build --target web
```

## Architecture

This is a Cargo workspace (`resolver = "2"`) with three crates:

- **`siminf-core`** — pure Rust solver, no FFI. The only crate with real implementation.
- **`siminf-wasm`** — `wasm-bindgen` wrapper over `siminf-core` (scaffold, not yet functional).
- **`siminf-py`** — PyO3 wrapper over `siminf-core` (scaffold, not yet functional).

### Core solver data flow

1. **`ModelBuilder`** (`model.rs`) — fluent builder that assembles compartments, transitions, global parameters (`gdata`), tspan, optional scheduled events, and a seed into a `Model`.
2. **`Solver::new(model)`** (`solver.rs`) — precomputes two CSC sparse matrices:
   - **S matrix** — stoichiometry: which compartments change and by how much per transition.
   - **G matrix** — dependency graph: which transitions are affected when a given transition fires (used to selectively recompute rates rather than recomputing all rates after every event).
3. **`Solver::run()`** — executes the Gillespie SSA per node (no parallelism yet despite `rayon` being a dependency). Stores compartment counts at every timepoint in a flat `Vec<i32>` indexed as `[timepoint][node][compartment]`.
4. **`TrajectoryResult`** (`trajectory.rs`) — output type. Key methods: `mean_compartments()` (returns mean across nodes per timepoint) and `to_csv()`.

### Propensity functions

`PropensityFn` is currently an enum with two hardcoded variants (`SToI`, `IToR`) plus a `Custom` closure variant. The hardcoded variants hard-code SIR-specific indexing (`u[0]` = S, `u[1]` = I, `u[2]` = R) and parameter positions in `gdata` (`gdata[0]` = beta, `gdata[1]` = gamma). This is a known limitation — the `mparse`-style string DSL is on the roadmap to replace this.

### Scheduled events

`ScheduledEvent` supports four event types: Exit, Enter, InternalTransfer, ExternalTransfer. The `EventProcessor` struct in `events.rs` is currently unused; event processing is inlined directly in `Solver::run()` and only handles ExternalTransfer (node-to-node transfers).

### Validation

Reference output is in `tests/validation/`. The key fixture is `trajectory_means_r.csv` — mean S/I/R per timepoint generated from the original R SimInf package with seed 42. The test `test_sir_validation_against_r_means` in `siminf-core/tests/sir_validation.rs` checks relative error < 1e-6 at timepoints 1, 50, 100, 150. The validation test silently skips if the fixture file is missing.

The reference implementation for all solver behaviour is the original SimInf C source (`stewid/SimInf/src/`).
