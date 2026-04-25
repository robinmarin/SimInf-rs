# SimInf-rs Code Review: Improvements

Grouped into sessions:
- **Group 1**: Sampling fixes (`events.rs`)
- **Group 2**: Event processing (`solver.rs`)
- **Group 3**: Data structures / indexing (`model.rs`, `solver.rs`)
- **Group 4**: Model field management (`solver.rs`)
- **Group 5**: PropensityFn redesign (big refactor)
- **Group 6**: Continuous state variables (separate feature)
- **Group 7**: Small isolated fix

---

Compared against the reference C implementation in `~/Gits/SimInf/src/`.

---

## 1. Critical: Time-stepping logic diverges from the C solver

**Status: DONE**

**File:** `solver.rs:123-231`

The C solver advances time in **integer unit steps** (`next_unit_of_time += 1.0`), processes events at unit boundaries, and stores output whenever `tt > tspan[U_it]`. The Rust solver instead jumps from one `tspan` output point to the next (`t = next_output_time` on line 217), running only a single Gillespie burst per output interval.

This means:
- If `tspan = [1, 50, 100, 150]`, the Rust solver runs a single Gillespie inner loop from t=1 to t=50 with `next_unit_of_time = t.floor() + 1.0` (line 130), but then jumps `t` straight to 50. The inner loop only simulates up to `next_unit_of_time` (line 143: `while node_t < next_unit_of_time`), which is at most `t.floor() + 1.0 = 2.0`. So **the simulation only advances 1 time unit per output interval**, skipping all intermediate time between tspan points.
- The C solver loops through every integer time unit regardless of output spacing, running SSA within each unit, processing events at each unit boundary, and only storing output at tspan points.

This is a fundamental correctness bug for non-unit-spaced tspan values. It currently passes validation only because the test uses `tspan = 1..=150` (unit spacing), so the outer loop happens to advance one unit at a time.

**Fix:** Add an inner loop over integer time units between output points, matching the C solver's `while(!done)` structure. Process events at each unit boundary, not just at output times.

---

## 2. Critical: `gdata` ordering is non-deterministic

**Status: DONE**

**File:** `solver.rs:41`

```rust
let gdata: Vec<f64> = model.gdata.values.values().cloned().collect();
```

`GlobalData` uses a `HashMap<String, f64>`. HashMap iteration order is non-deterministic. The hardcoded `PropensityFn::SToI` reads `gdata[0]` as beta and `gdata[1]` as gamma. If the HashMap happens to iterate in a different order, the parameters are swapped silently.

This currently works by luck. The `Custom` propensity variant receives the same shuffled `gdata` slice, so users of the Custom API also cannot rely on parameter ordering.

**Fix:** Either use a `BTreeMap` (deterministic order) or, better, switch `gdata` to a flat `Vec<f64>` with named indices (matching the C solver's flat `gdata` array with known positions). The HashMap-based `GlobalData` is fundamentally at odds with positional indexing.

---

## 3. Critical: `PropensityFn::Clone` silently produces a broken closure

**Status: DONE**

**File:** `model.rs:62-73`

```rust
impl Clone for PropensityFn {
    fn clone(&self) -> Self {
        match self {
            PropensityFn::Custom { name, .. } => PropensityFn::Custom {
                name: name.clone(),
                eval: Box::new(|_, _, _, _| 0.0),  // <-- always returns 0
            },
            ...
        }
    }
}
```

Cloning a `Custom` propensity replaces the actual rate function with a no-op that returns 0.0. This is a silent correctness trap. Any code path that clones a `Model` (which contains `Transition` which contains `PropensityFn`) will produce a model that silently computes zero rates for all custom transitions.

`Solver::new()` currently takes ownership so it doesn't clone, but the `Solver` struct itself derives `Clone`, meaning `solver.clone().run()` would silently produce garbage.

**Fix:** Use `Arc<dyn Fn(...)>` instead of `Box<dyn Fn(...)>` so the closure can be cheaply cloned by sharing. Alternatively, remove `Clone` from `PropensityFn` and `Transition` entirely and let the compiler enforce that they aren't cloned.

---

## 4. Critical: Dependency graph (`G` matrix) is incomplete

**Status: DONE**

**File:** `solver.rs:302-335`

The `build_g_matrix` function only checks if a transition's `to` compartments overlap with another transition's `from` compartments:

```rust
for from in &transitions[tr_idx].to {
    for other_from in &transitions[other_idx].from {
        if from.0 == other_from.0 { return true; }
    }
}
```

This misses the case where a transition's `from` compartments overlap with another transition's `from` compartments. In SIR: when S->I fires, both S and I change. The S->I transition itself depends on both S and I, so it should be in its own dependency set. Currently it only appears if the `to` (I) matches another transition's `from` -- which it does for I->R, but not for S->I itself (S is consumed, not produced, by S->I).

In the C solver, the G matrix is built externally (by the R layer) with full knowledge of which compartments are read by each propensity function, not just which compartments appear in `from`/`to`. The dependency is: "if transition `j` fires and changes compartment `c`, and transition `k`'s rate depends on compartment `c`, then `k` must be in column `j` of G."

For SIR this happens to work because:
- S->I fires: S decreases, I increases. G column for S->I should list both S->I (depends on S,I) and I->R (depends on I). The current code finds I->R (I is in `to` of S->I, and `from` of I->R) but also needs S->I itself (I is in `to` of S->I, and `from` of S->I -- wait, `from` of S->I is S, and `to` is I, so actually the code finds S->I in the affected set if I is in another transition's from... let me re-check).

Actually, re-checking: when S->I fires (tr_idx=0), the code checks `transitions[0].to = [I]` against `transitions[other].from`. For other=0: `transitions[0].from = [S]`, so I != S, no match. For other=1: `transitions[1].from = [I]`, so I == I, match. So column 0 of G only contains transition 1, not transition 0. This means after S->I fires, the S->I rate itself is NOT recalculated via the dependency graph.

However, looking at solver.rs:169, the code recalculates `rates[tr]` unconditionally after firing:
```rust
rates[tr] = self.compute_rate(node_state, tr, node_t);
```

This patches over the incomplete G matrix, but it's a divergence from the C solver's approach and won't generalize correctly to more complex models where the `from` compartments of the firing transition also affect other transitions not captured by the `to`-only dependency check.

**Fix:** Build the G matrix from the full set of compartments that each transition reads (from + to), not just `to`. Or better, derive dependencies from the S matrix columns: if any compartment changed by transition `j` (i.e., any row with non-zero entry in S column `j`) is also a compartment that transition `k` reads, then `k` depends on `j`.

---

## 5. Major: Sampling uses wrong distribution

**Status: DONE**

**File:** `events.rs:124-147` (`determine_n`)

The C solver uses `gsl_ran_binomial(rng, proportion, Nindividuals)` to determine how many individuals to sample when `n == 0`. The Rust port uses:

```rust
let p = (rng.gen::<f64>() * (total as f64) * proportion).round() as usize;
```

This is `round(U * total * proportion)` where U ~ Uniform(0,1). This is NOT a binomial distribution. It produces a different distribution with different mean and variance. The binomial has mean `n*p` and variance `n*p*(1-p)`; this custom formula has mean `total * proportion / 2` (half the expected value) and a uniform-ish distribution shape.

**Fix:** Use a proper binomial distribution. The `rand_distr` crate provides `Binomial`.

---

## 6. Major: Sampling uses wrong algorithm (not hypergeometric)

**Status: DONE**

**File:** `events.rs:224-278`

The C solver uses `gsl_ran_hypergeometric` for the equal-weight case (multivariate hypergeometric distribution, Gentle 2003 p.206). The Rust port uses an ad-hoc probability-based sampling:

```rust
let prob = remaining as f64 / available as f64;
if rng.gen::<f64>() < prob { remaining } else { 0 }
```

This doesn't implement the hypergeometric distribution. The hypergeometric models drawing without replacement from a finite population with known composition. The Rust code's approach doesn't correctly account for the changing population composition as individuals are drawn.

**Fix:** Implement multivariate hypergeometric sampling matching the C solver's algorithm, or use a crate that provides it.

---

## 7. Major: Enter events use wrong sampling function

**Status: DONE**

**File:** `events.rs:341-363`

The C solver has two distinct sampling functions: `SimInf_sample_select` (for Exit/InternalTransfer/ExternalTransfer) and `SimInf_sample_select_enter` (for Enter). The Enter variant samples individuals proportionally to weights in the select matrix, distributing them into compartments based on those weights. It does NOT use the current node state for sampling -- it distributes new arrivals according to weights.

The Rust port uses the same `sample_from_select` for all event types, including Enter. This means Enter events incorrectly use the current node state to distribute incoming individuals rather than the select matrix weights.

**Fix:** Implement a separate `sample_from_select_enter` matching the C solver's `SimInf_sample_select_enter`.

---

## 8. Major: E2 events incorrectly require shift >= 0

**Status: DONE**

**File:** `events.rs:94-96`

```rust
if event.shift < 0 {
    return Err(EventError::InvalidShift);
}
```

The C solver allows `shift < 0` for ExternalTransfer events. When shift is negative, individuals are transferred to the same compartment indices in the destination node (no shift applied). The Rust code errors on this valid case.

**Fix:** Match the C behavior: when `shift < 0`, transfer individuals to the same compartment indices without shifting.

---

## 9. Major: No rate recalculation after events

**Status: DONE (architectural note)**

**File:** `solver.rs:178-215`

After processing E1 and E2 events, the C solver marks affected nodes in `update_node[]` and recalculates their transition rates before advancing to the next time unit. The Rust architecture does not store persistent per-node rate arrays — rates are recomputed from scratch at the start of each unit-of-time iteration via `compute_all_rates()`. This is functionally equivalent to the C solver's selective recalculation (since `compute_rate` is a pure function of `node_state` and `gdata`), but the optimization is not replicated.

Note: A naive implementation that computes rates after events and discards them (no-op) was attempted but removed as it served no purpose.

---

## 10. Major: No negative state checking after transitions

**Status: DONE**

**File:** `solver.rs:167`

The C solver checks for negative compartment counts after every transition:
```c
if (m.u[node * m.Nc + m.irS[j]] < 0) {
    m.error = SIMINF_ERR_NEGATIVE_STATE;
}
```

The Rust solver applies transitions without any validation:
```rust
apply_transition(node_state, &self.S, tr);
```

A negative compartment count indicates a bug in the propensity function (it should have returned 0 to prevent the transition). Silent negative counts corrupt the simulation.

**Fix:** Check for negative state after `apply_transition` and either return an error or log a diagnostic.

---

## 11. Major: No rate validity checking

**Status: DONE**

**File:** `solver.rs:247-273`

The C solver checks every computed rate for `!R_FINITE(rate) || rate < 0.0` and sets an error flag. The Rust solver does no validation on computed rates. A NaN or negative rate from a buggy propensity function will corrupt the simulation silently (NaN propagates through the cumulative sum, tau becomes NaN, etc.).

**Fix:** Validate rates after computation. At minimum check for NaN/infinity and negative values.

---

## 12. Major: Floating-point edge case in transition selection

**Status: DONE**

**File:** `solver.rs:156-165`

The C solver has an elaborate floating-point fix after selecting a transition:
```c
if (tr >= m.Nt) tr = m.Nt - 1;
if (m.t_rate[node * m.Nt + tr] == 0.0) {
    // Go backwards to find first nonzero rate
    for (; tr > 0 && m.t_rate[node * m.Nt + tr] == 0.0; tr--);
    if (m.t_rate[node * m.Nt + tr] == 0.0) {
        m.sum_t_rate[node] = 0.0;
        break; // nil event
    }
}
```

The Rust solver has none of this. If floating-point drift causes `u_rand` to exceed `cum` for all transitions (because `sum_t_rate` was computed incrementally via deltas and has drifted), the loop exits with `tr = 0` (the initial value, since the break was never hit). This silently fires the wrong transition.

**Fix:** Add the clamp-and-backtrack logic from the C solver.

---

## 13. Moderate: `rng.gen::<f64>()` can return 0.0

**Status: DONE**

**File:** `solver.rs:149`

```rust
let tau = -rng.gen::<f64>().ln() / sum_t_rate;
```

`rng.gen::<f64>()` returns values in `[0, 1)`, so it can return exactly 0.0. `(0.0f64).ln()` is `-inf`, making `tau = inf`, which means the break on line 150 fires and the loop exits. This is benign in practice but differs from the C solver which uses `gsl_rng_uniform_pos` (returns `(0, 1)` exclusive of 0). Same issue on line 156.

**Fix:** Use `rng.gen_range(f64::EPSILON..1.0)` or a wrapper that excludes 0.0, matching `gsl_rng_uniform_pos` semantics.

---

## 14. Moderate: `SparseMatrix` uses `i32` indices, limiting scale

**Status: DONE**

**File:** `model.rs:76-119`

The `SparseMatrix` struct is now generic: `SparseMatrix<T>` with `ir: Vec<i32>`, `jc: Vec<i32>`, `pr: Vec<T>`. The `usize` index concern remains (still uses `i32` for compatibility with C sparse matrix format), but the dual `pr`/`pr_f64` union design is eliminated. `SparseMatrix<f64>` is used for select_matrix (weights) and `SparseMatrix<i32>` for S, G, and shift matrices (stoichiometry/dependencies). The old `get_f64`/`get_i32` ambiguity is replaced by explicit typed methods.

---

## 15. Moderate: `Solver` stores both `model` and copies of its fields

**Status: DONE**

**File:** `solver.rs:21-40`

The `model: Model` field was removed from `Solver`. The struct now stores only the extracted fields: `transitions`, `compartment_names`, `initial_counts`, `num_compartments`, `num_transitions`, `num_nodes`, `S`, `G`, `gdata`, `ldata`, `u0`, `tspan`, `events`, `select_matrix`, `shift_matrix`, `e2_rng`, `seed`. The model is consumed in `Solver::new()` and only the needed fields are retained, eliminating duplication.

---

## 16. Moderate: Per-node initial state is identical

**Status: DONE**

**File:** `model.rs:215`, `solver.rs:128-136`

`Model` now has an optional `u0: Option<Vec<i32>>` field — a flat Vec indexed as `[node * Nc + compartment]`. If `Some`, it provides per-node initial compartment counts. If `None`, all nodes use `Compartment.initial_count`. The `ModelBuilder::u0()` method accepts the flat Vec. Validation ensures `u0.len() == num_nodes * num_compartments` if provided.

---

## 17. Moderate: `ldata` (local data) is never passed to propensity functions

**Status: DONE**

**File:** `model.rs:199`, `solver.rs:30, 65-68, 332-360`

`ldata: Vec<LocalData>` (one entry per node, each with `values: Vec<f64>`) is now flattened into `Option<Arc<Vec<f64>>>` in `Solver` during construction, stored as `[node * Nld + ld_idx]`. `compute_rate` now takes a `node_idx` parameter and passes the per-node `ldata` slice to `Custom` propensity closures. `ModelBuilder::build()` validates all `LocalData` entries have the same number of values so the flat indexing is well-defined.

---

## 18. Moderate: No continuous state variables (`v`)

**Status: DONE (Group 6)**

**Files:** `model.rs:167-221`, `solver.rs:22-89, 30-38, 196-209, 361-380`

The Rust solver now supports continuous state variables (`v`, `v_new`) alongside discrete compartment counts (`u`).

**Changes made:**

1. **`PropensityFn::Custom` signature** updated to include `&v` (5th parameter):
   ```rust
   Fn(&[i32] /* u */, &[f64] /* v */, &[f64] /* ldata */, &[f64] /* gdata */, f64 /* t */) -> f64
   ```

2. **`Model` / `ModelBuilder`** now have:
   - `nd: usize` — number of continuous variables per node
   - `v0: Option<Vec<f64>>` — initial continuous state (flattened `[node * Nd + d]`)
   - `pts_fun: Option<Arc<dyn Fn(&mut [f64], &[i32], &[f64], &[f64], &[f64], usize, f64)>>` — post-timestep callback with signature `(v_new, u, v, ldata, gdata, node_idx, t)`
   - Builder methods: `.nd()`, `.v0()`, `.pts_fun()`
   - `build()` validates `v0.len() == num_nodes * nd` when `nd > 0`

3. **`Solver`** now stores:
   - `nd: usize`, `v: Vec<f64>`, `v_new: Vec<f64>`, `pts_fun: Option<Arc<...>>`
   - Initialized from `model.nd`, `model.v0` in `Solver::new()`

4. **SSA loop** (`solver.rs:196-209`): Each node's `v` slice passed to `compute_all_rates`/`compute_rate`, which passes it to `Custom` propensity closures.

5. **Post-timestep** (`solver.rs:361-380`): After E1/E2 events and before advancing `t`, for each node `pts_fun(v_new_slice, u_slice, v_slice, ldata, gdata, node_idx, next_unit_of_time)` is called, then `v`/`v_new` are swapped. Skipped when `nd == 0`.

---

## 19. Moderate: `select_matrix` default is wrong shape

**Status: DONE (Group 7)**

**File:** `solver.rs:70-82`, `model.rs:328-330`

Two changes made:

1. **Buggy `jc` construction fixed** (`solver.rs:73-78`): The old loop `jc.push(jc.len() as i32 + 1)` produced `[0, 2, 3, 4, ...]` because `jc.len()` grew with each push. Replaced with `jc.push((i + 1) as i32)` which correctly produces `[0, 1, 2, 3, ...]`.

2. **Validation added** (`model.rs:328-330`): `ModelBuilder::build()` now returns an error if events are present but `select_matrix` is not provided. The old implicit diagonal default was semantically wrong (it created columns indexed by compartment number, but `event.select` is an indirect index into selection groups). Users with events must now explicitly provide a correctly-structured `select_matrix`.

---

## 20. Minor: `apply_transition` has a redundant cast

**Status: DONE**

**File:** `solver.rs:344`

```rust
u[row] = (u[row] as i32 + val) as i32;
```

`u[row]` is already `i32` and `val` is already `i32`. The casts are no-ops.

---

## 21. Minor: `to_csv` truncates time to integer

**Status: DONE**

**File:** `trajectory.rs:48`

```rust
s.push_str(&format!("{},{}", node + 1, t as i32));
```

The `t as i32` cast silently truncates fractional tspan values. If someone uses non-integer timepoints, the CSV output will be wrong.

**Fix:** Format `t` as a float.

---

## 22. Minor: Event errors are silently ignored

**Status: DONE**

**File:** `solver.rs:195, 212`

```rust
let _ = event_processor.process_e1_events(...);
// ...
if let Err(_) = e2_processor.apply_e2_event(...) {}
```

Event processing errors are explicitly discarded. This hides bugs in event configuration. The C solver sets an error flag and aborts the simulation on event errors.

**Fix:** Propagate event errors. Change `run()` to return `Result<TrajectoryResult, SimInfError>`.

---

## 23. Minor: `TrajectoryView::time()` takes `self` by value

**Status: DONE**

**File:** `trajectory.rs:68`

```rust
pub fn time(mut self, idx: usize) -> Self {
```

This consumes the view. The `time_index` field set here is also never used by `all_nodes_at_time` (which takes its own `t_idx` parameter), making the `time()` method only useful for `node_state()`. The API is inconsistent.

**Fix:** Either make `time()` take `&mut self` or redesign the view API to be consistent.

---

## 24. Minor: No `#[must_use]` on `ModelBuilder::build()`

**Status: DONE**

**File:** `model.rs:309`

The `build()` method returns `Result<Model, String>`. Without `#[must_use]`, a user could call `.build()` and ignore the result. Also, using `String` for errors instead of a proper error type prevents structured error handling.

**Fix:** Add `#[must_use]` and consider a proper error enum.

---

## 25. Design: `PropensityFn` enum with hardcoded variants doesn't scale

**Status: DONE (Group 5)**

**File:** `model.rs:52-60`

The `SToI` and `IToR` variants were removed. `PropensityFn` is now a simple `#[derive(Clone)]` enum with only one variant:

```rust
#[derive(Clone)]
pub enum PropensityFn {
    Custom {
        name: String,
        eval: Arc<dyn Fn(&[i32], &[f64], &[f64], f64) -> f64 + Send + Sync>,
    },
}
```

The `Arc` wrapping the closure makes cloning cheap (fixes issue #3). The SIR-specific propensities (beta, gamma lookup, SIR rate formulas) were moved to `examples/sir.rs` and `tests/sir_validation.rs` as explicit closures. `compute_rate` in `solver.rs` now handles only the `Custom` variant.

---

## 26. Design: Events are filtered/cloned repeatedly per output interval

**Status: DONE**

**File:** `solver.rs:51, 140-148`

Events are now sorted by time during `Solver::new()`. The main loop in `run()` uses a cursor (`events_index`) that advances through the sorted events, replacing the O(events * nodes * timepoints) filter+clone pattern with O(1) amortized cursor increment.

---

## 27. Design: Node state as `Vec<Vec<i32>>` instead of flat array

**Status: DONE**

**File:** `solver.rs:111-127, 151, 256-269`

`u_current` is now a flat `Vec<i32>` indexed as `u_current[node * Nc + compartment]`. The parallel SSA loop uses `par_chunks_mut(num_compartments)` to get per-node mutable slices, replacing `par_iter_mut()` over `Vec<Vec<i32>>`. `store_solution` works directly with flat slices. E1 events are processed per-node by passing `&mut u_current[node_offset..node_offset+Nc]` slices. E2 (ExternalTransfer) events are handled inline in the solver loop with direct flat indexing, avoiding the two-mutable-borrow issue that would arise from passing the flat `u_current` to a function that needs two node slices simultaneously.