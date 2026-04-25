# SimInf-rs

A Rust port of [SimInf](https://github.com/stewid/SimInf) — a framework for data-driven stochastic
disease spread simulations. This project reimplements the SimInf simulation engine in Rust, with
the goal of broader platform reach: native Rust, Python, and browser (WebAssembly) targets from a
single core codebase.

> **Status:** Core solver is working and validated against R reference output. WASM and Python bindings are scaffolded. Do not use in production.

---

## Motivation

The original SimInf R package is an excellent, well-validated framework for epidemiological
simulation. It uses the Gillespie stochastic simulation algorithm (SSA) with OpenMP parallelism
over subpopulations (nodes), and has been applied to national-scale livestock disease modelling.

However, its R/C architecture limits who can use it:

- Requires R and a C toolchain to install
- Cannot run in a browser
- Cannot be called from Python without subprocess hacks
- Not available as a standalone library for systems software

This port aims to change that by implementing the simulation engine in Rust, which allows the same
validated core to be compiled to native binaries, Python wheels (via PyO3), and WebAssembly (via
wasm-bindgen) — making SimInf-style simulations accessible from any environment.

---

## Architecture

This repository is a Cargo workspace with three crates:

```
siminf-rs/
├── crates/
│   ├── siminf-core/    # Pure Rust solver — no FFI, no runtime dependencies
│   ├── siminf-wasm/    # WebAssembly bindings via wasm-bindgen
│   └── siminf-py/      # Python bindings via PyO3 (scaffolded)
├── examples/
│   └── sir.rs
└── tests/
    └── fixtures/       # CSV outputs from the R package used for validation
```

### `siminf-core`

The heart of the project. Implements:

- Gillespie SSA inner loop
- Node-level parallelism via `rayon` (replaces OpenMP)
- Scheduled events: births, deaths, transfers between nodes
- A typed `Model` builder API
- `Trajectory` output type (compartment counts per node per timepoint)

`siminf-core` has no FFI dependencies and is intended to be usable in `no_std` environments
where possible.

### `siminf-wasm`

**Scaffold.** Thin `wasm-bindgen` layer over `siminf-core`. Enables SimInf simulations to run in a browser tab
with no server required — suitable for outbreak dashboards, policy tools, and interactive
educational interfaces built in React, Vue, or plain TypeScript.

> Note: the WASM build is single-threaded (Rayon parallelism is disabled) due to browser
> threading constraints. This is acceptable for exploratory, smaller-scale models typical of
> browser use.

### `siminf-py`

**Scaffold.** PyO3 bindings exposing `siminf-core` as a Python package. Intended to integrate naturally with
NumPy, pandas, and the broader scientific Python ecosystem. Currently scaffolded; not yet
functional.

---

## Relationship to the original SimInf

This project is an independent port, not a fork. It aims to produce numerically equivalent output
to the original R package and uses the original C source as the reference specification for solver
behaviour.

The original package should be cited if you use this software in research:

> Widgren S, Bauer P, Eriksson R, Engblom S (2019). SimInf: An R Package for Data-Driven
> Stochastic Disease Spread Simulations. *Journal of Statistical Software*, 91(12), 1–42.
> https://doi.org/10.18637/jss.v091.i12

> Bauer P, Engblom S, Widgren S (2016). Fast event-based epidemiological simulations on national
> scales. *International Journal of High Performance Computing Applications*, 30(4), 438–453.
> https://doi.org/10.1177/1094342016635723

---

## Quickstart (native Rust)

> Prerequisites: Rust 1.75+ via [rustup](https://rustup.rs)

```bash
git clone https://github.com/your-username/SimInf-rs.git
cd SimInf-rs
cargo run --example sir
```

The SIR example replicates the canonical SimInf README scenario: 1000 nodes, S=99, I=5, R=0,
β=0.16, γ=0.077, over 150 time steps.

```rust
use siminf_core::{Model, Solver};

let model = Model::builder()
    .compartments(["S", "I", "R"])
    .transition("S -> beta*S*I/N -> I")
    .transition("I -> gamma*I -> R")
    .gdata([("beta", 0.16), ("gamma", 0.077)])
    .u0(/* 1000 nodes, S=99 I=5 R=0 */)
    .tspan(1..=150)
    .build()?;

let result = Solver::new().seed(42).run(&model)?;
let trajectory = result.trajectory(); // returns a DataFrame-like structure
```

> API shown above is illustrative and subject to change during early development.

---

## WebAssembly / TypeScript

```bash
cd crates/siminf-wasm
wasm-pack build --target web
```

This produces a `pkg/` directory with a `.wasm` binary, JavaScript glue, and TypeScript type
definitions. Import it in any web project:

```typescript
import init, { SirModel, run } from './pkg/siminf_wasm';

await init();

const model = new SirModel({ beta: 0.16, gamma: 0.077, nodes: 1000 });
const result = run(model, { tspan: [1, 150] });
```

---

## Validation

Correctness is validated against output from the original R package. Reference fixtures are
generated by running SimInf with fixed seeds and exporting results to CSV:

```bash
# Requires R and the SimInf package
Rscript tests/fixtures/generate.R
```

The Rust test suite then asserts that aggregate statistics (mean compartment counts per timepoint
across nodes) match the R output within a tolerance of 1e-6:

```bash
cargo test --workspace
```

---

## Roadmap

- [x] `siminf-core`: Gillespie SSA solver
- [x] `siminf-core`: validation against R package fixtures
- [ ] `siminf-core`: scheduled events (births, deaths, node transfers) — **scaffold**
- [ ] `siminf-wasm`: wasm-bindgen bindings — **scaffold**
- [ ] `siminf-wasm`: TypeScript example (React dashboard)
- [ ] `siminf-py`: PyO3 bindings — **scaffold**
- [ ] `siminf-py`: pandas-native trajectory output
- [ ] `mparse`-style string DSL for model definition
- [ ] Continuous compartments (beyond integer counts)
- [ ] Publish `siminf-core` to crates.io

---

## Contributing

The project is in early development and the architecture is still being established. If you are
interested in contributing, opening an issue to discuss before submitting a PR is strongly
recommended.

The original SimInf R package source is the reference for all solver behaviour. When in doubt
about a numerical detail, the C source in `stewid/SimInf/src/` is the spec.

---

## License

Licensed under [GPL-3.0](LICENSE), consistent with the original SimInf package.
