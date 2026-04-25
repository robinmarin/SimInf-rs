use rand::Rng;
use rand::rngs::StdRng;
use rand::SeedableRng;
use rayon::prelude::*;
use std::sync::Arc;

use crate::events::EventError;
use crate::model::{
    Transition, PropensityFn, Model, SparseMatrix,
};
use crate::events::EventProcessor;
use crate::model::EventType;

#[derive(Debug, Clone)]
pub enum SimInfError {
    NegativeState(usize),
    InvalidRate(usize),
    EventError(EventError),
}

#[derive(Clone)]
#[allow(non_snake_case)]
pub struct Solver {
    transitions: Arc<Vec<Transition>>,
    compartment_names: Vec<String>,
    initial_counts: Vec<i32>,
    num_compartments: usize,
    num_transitions: usize,
    num_nodes: usize,
    nd: usize,
    S: Arc<SparseMatrix<i32>>,
    G: Arc<SparseMatrix<i32>>,
    gdata: Arc<Vec<f64>>,
    ldata: Option<Arc<Vec<f64>>>,
    u0: Option<Vec<i32>>,
    v: Vec<f64>,
    v_new: Vec<f64>,
    pts_fun: Option<Arc<dyn Fn(&mut [f64], &[i32], &[f64], &[f64], &[f64], usize, f64) + Send + Sync>>,
    tspan: Arc<Vec<f64>>,
    events: Arc<Vec<crate::model::ScheduledEvent>>,
    select_matrix: Arc<SparseMatrix<f64>>,
    shift_matrix: Arc<SparseMatrix<i32>>,
    e2_rng: StdRng,
    seed: u64,
}

impl std::fmt::Debug for Solver {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Solver")
            .field("transitions", &self.transitions)
            .field("compartment_names", &self.compartment_names)
            .field("initial_counts", &self.initial_counts)
            .field("num_compartments", &self.num_compartments)
            .field("num_transitions", &self.num_transitions)
            .field("num_nodes", &self.num_nodes)
            .field("nd", &self.nd)
            .field("S", &self.S)
            .field("G", &self.G)
            .field("gdata", &self.gdata)
            .field("ldata", &self.ldata)
            .field("u0", &self.u0)
            .field("v", &self.v)
            .field("v_new", &self.v_new)
            .field("pts_fun", &"Option<Arc<dyn Fn(...)>>")
            .field("tspan", &self.tspan)
            .field("events", &self.events)
            .field("select_matrix", &self.select_matrix)
            .field("shift_matrix", &self.shift_matrix)
            .field("seed", &self.seed)
            .finish()
    }
}

#[allow(non_snake_case)]
impl Solver {
    pub fn new(model: Model) -> Self {
        let num_compartments = model.compartments.len();
        let num_transitions = model.transitions.len();
        let num_nodes = model.num_nodes;

        let s = build_s_matrix(&model.transitions, num_compartments);
        let g = build_g_matrix(&model.transitions, num_transitions);

        let gdata: Vec<f64> = model.gdata.to_vec();
        let seed = model.seed.unwrap_or(42);
        let mut events = model.events.clone();
        let tspan = model.tspan.clone();
        let transitions = model.transitions.clone();
        let compartment_names: Vec<String> = model.compartments.iter().map(|c| c.name.clone()).collect();
        let initial_counts: Vec<i32> = model.compartments.iter().map(|c| c.initial_count).collect();

        let ldata: Option<Arc<Vec<f64>>> = if model.ldata.is_empty() {
            None
        } else {
            let flattened: Vec<f64> = model.ldata.iter().map(|ld| ld.values.clone()).collect::<Vec<_>>().into_iter().flatten().collect();
            Some(Arc::new(flattened))
        };

        events.sort_by(|a, b| a.time.partial_cmp(&b.time).unwrap_or(std::cmp::Ordering::Equal));

        let select_matrix = model.select_matrix.clone().unwrap_or_else(|| {
            SparseMatrix::new(
                (0..num_compartments).map(|i| i as i32).collect(),
                {
                    let mut jc = vec![0];
                    for i in 0..num_compartments {
                        jc.push((i + 1) as i32);
                    }
                    jc
                },
                vec![1.0; num_compartments],
            )
        });
        let shift_matrix = model.shift_matrix.clone().unwrap_or_else(|| {
            SparseMatrix::new(
                vec![],
                vec![0; num_compartments + 1],
                vec![],
            )
        });

        let e2_rng = StdRng::seed_from_u64(seed.wrapping_add(0xE2E2E2E2E2E2E2E2u64));

        let u0 = model.u0;

        let nd = model.nd;
        let pts_fun = model.pts_fun;

        let v = if nd > 0 {
            model.v0.clone().unwrap_or_else(|| vec![0.0; num_nodes * nd])
        } else {
            vec![]
        };
        let v_new = if nd > 0 {
            v.clone()
        } else {
            vec![]
        };

        Self {
            transitions: Arc::new(transitions),
            compartment_names,
            initial_counts,
            num_compartments,
            num_transitions,
            num_nodes,
            nd,
            S: Arc::new(s),
            G: Arc::new(g),
            gdata: Arc::new(gdata),
            ldata,
            u0,
            v,
            v_new,
            pts_fun,
            tspan: Arc::new(tspan),
            events: Arc::new(events),
            select_matrix: Arc::new(select_matrix),
            shift_matrix: Arc::new(shift_matrix),
            e2_rng,
            seed,
        }
    }

    pub fn run(&mut self) -> Result<super::TrajectoryResult, SimInfError> {
        // Each node gets its own RNG seeded from the master seed so that node
        // execution order does not affect reproducibility when run in parallel.
        let mut node_rngs: Vec<StdRng> = (0..self.num_nodes)
            .map(|node| StdRng::seed_from_u64(self.seed.wrapping_add(node as u64)))
            .collect();

        let tspan = &self.tspan;
        let num_timepoints = tspan.len();
        let t_end = tspan[tspan.len() - 1];

        let mut u = vec![0i32; self.num_nodes * self.num_compartments * num_timepoints];

        for node in 0..self.num_nodes {
            let src = if let Some(ref u0) = self.u0 {
                let offset = node * self.num_compartments;
                &u0[offset..offset + self.num_compartments]
            } else {
                &self.initial_counts
            };
            let dst_offset = node * self.num_compartments;
            u[dst_offset..dst_offset + self.num_compartments].copy_from_slice(src);
        }

        let mut u_current: Vec<i32> = vec![0i32; self.num_nodes * self.num_compartments];
        for node in 0..self.num_nodes {
            let offset = node * self.num_compartments;
            u_current[offset..offset + self.num_compartments]
                .copy_from_slice(&u[offset..offset + self.num_compartments]);
        }

        let mut t = tspan[0];
        let mut output_idx = 0;
        let mut events_index = 0usize;

        store_solution(&mut u, &u_current, output_idx, self.num_nodes, self.num_compartments);

        let mut next_unit_of_time = t.floor() + 1.0;

        while t < t_end {
            let next_output_time = if output_idx + 1 < num_timepoints {
                tspan[output_idx + 1]
            } else {
                t_end + 1.0
            };

            let events_in_interval: Vec<_> = self.events[events_index..]
                .iter()
                .take_while(|e| e.time > t && e.time <= next_output_time)
                .cloned()
                .collect();
            events_index += events_in_interval.len();

            u_current.par_chunks_mut(self.num_compartments)
                .zip(node_rngs.par_iter_mut())
                .enumerate()
                .for_each(|(node_idx, (node_state, rng))| {
                    let node_v = if self.nd > 0 {
                        let offset = node_idx * self.nd;
                        &self.v[offset..offset + self.nd]
                    } else {
                        &[]
                    };
                    let mut node_t = t;
                    let mut rates = self.compute_all_rates(node_state, node_v, node_idx, node_t);

                    while node_t < next_unit_of_time {
                        let sum_t_rate: f64 = rates.iter().sum();
                        if sum_t_rate <= 0.0 {
                            break;
                        }

                        let tau = -rng.gen_range(f64::EPSILON..1.0).ln() / sum_t_rate;
                        if node_t + tau >= next_unit_of_time {
                            break;
                        }

                        node_t += tau;

                        let u_rand = rng.gen_range(0.0..1.0) * sum_t_rate;
                        let mut cum = 0.0;
                        let mut tr = 0;
                        for (i, rate) in rates.iter().enumerate() {
                            cum += *rate;
                            if u_rand <= cum {
                                tr = i;
                                break;
                            }
                        }

                        if tr >= self.num_transitions {
                            tr = self.num_transitions.saturating_sub(1);
                        }
                        if rates[tr] == 0.0 {
                            while tr > 0 && rates[tr] == 0.0 {
                                tr -= 1;
                            }
                            if rates[tr] == 0.0 {
                                break;
                            }
                        }

                        apply_transition(node_state, &self.S, tr).expect("Negative compartment count detected");

                        rates[tr] = self.compute_rate(node_state, node_v, node_idx, tr, node_t);

                        for j in self.G.jc[tr] as usize..self.G.jc[tr + 1] as usize {
                            let affected_tr = self.G.ir[j] as usize;
                            rates[affected_tr] = self.compute_rate(node_state, node_v, node_idx, affected_tr, node_t);
                        }
                    }
                });

            let e1_events: Vec<_> = events_in_interval.iter()
                .filter(|e| e.event_type != EventType::ExternalTransfer)
                .cloned()
                .collect();

            for node_idx in 0..self.num_nodes {
                let node_events: Vec<_> = e1_events.iter()
                    .filter(|e| e.node == node_idx)
                    .cloned()
                    .collect();

                if !node_events.is_empty() {
                    let node_offset = node_idx * self.num_compartments;
                    let event_processor = EventProcessor::new(
                        self.select_matrix.as_ref(),
                        self.shift_matrix.as_ref(),
                        self.num_compartments,
                    );
                    let _ = event_processor.process_e1_events(
                        &node_events,
                        node_idx,
                        &mut u_current[node_offset..node_offset + self.num_compartments],
                        &mut node_rngs[node_idx],
                    );
                }
            }

            drop(e1_events);

            for event in &events_in_interval {
                if event.event_type == EventType::ExternalTransfer {
                    let node_idx = event.node;
                    let dest = match event.dest {
                        Some(d) => d,
                        None => continue,
                    };

                    let src_offset = node_idx * self.num_compartments;
                    let dst_offset = dest * self.num_compartments;

                    let n = {
                        let src_state = &u_current[src_offset..src_offset + self.num_compartments];
                        let event_processor = EventProcessor::new(
                            self.select_matrix.as_ref(),
                            self.shift_matrix.as_ref(),
                            self.num_compartments,
                        );
                        event_processor.determine_n(event, src_state, &mut self.e2_rng)
                            .map_err(SimInfError::EventError)?
                    };

                    let individuals = {
                        let src_state = &u_current[src_offset..src_offset + self.num_compartments];
                        let event_processor = EventProcessor::new(
                            self.select_matrix.as_ref(),
                            self.shift_matrix.as_ref(),
                            self.num_compartments,
                        );
                        event_processor.sample_from_select(event.select as usize, src_state, n, &mut self.e2_rng)
                            .map_err(SimInfError::EventError)?
                    };

                    let shift_idx = if event.shift < 0 { 0 } else { event.shift as usize };

                    let select_col = event.select as usize;
                    for i in 0..self.select_matrix.num_rows_in_col(select_col) {
                        let jj = self.select_matrix.ir[self.select_matrix.jc[select_col] as usize + i] as usize;
                        let ll = self.shift_matrix.get_i32(shift_idx, i) as isize;
                        let dst_idx = jj as isize + ll;

                        if dst_idx < 0 || dst_idx >= self.num_compartments as isize {
                            continue;
                        }

                        let src_idx = jj;
                        let count = individuals[src_idx];

                        u_current[src_offset + src_idx] -= count;
                        u_current[dst_offset + dst_idx as usize] += count;
                    }
                }
            }

            if self.nd > 0 {
                let node_ldata = self.ldata.as_ref();
                for node_idx in 0..self.num_nodes {
                    let v_new_slice = &mut self.v_new[node_idx * self.nd..node_idx * self.nd + self.nd];
                    let u_slice = &u_current[node_idx * self.num_compartments..node_idx * self.num_compartments + self.num_compartments];
                    let v_slice = &self.v[node_idx * self.nd..node_idx * self.nd + self.nd];
                    let ldata_slice = if let Some(ref ldata) = node_ldata {
                        let nld = ldata.len() / self.num_nodes;
                        let offset = node_idx * nld;
                        &ldata[offset..offset + nld]
                    } else {
                        &[]
                    };
                    if let Some(ref pts_fun) = self.pts_fun {
                        pts_fun(v_new_slice, u_slice, v_slice, ldata_slice, &self.gdata, node_idx, next_unit_of_time);
                    }
                }
                std::mem::swap(&mut self.v, &mut self.v_new);
            }

            t = next_unit_of_time;

            while output_idx < num_timepoints && t >= tspan[output_idx] {
                store_solution(&mut u, &u_current, output_idx, self.num_nodes, self.num_compartments);
                output_idx += 1;
            }

            if t >= t_end {
                break;
            }

            next_unit_of_time = (t.floor() + 1.0).max(t + 1.0);
        }

        let final_offset = self.num_nodes * self.num_compartments;
        u[..final_offset].copy_from_slice(&u_current);

        Ok(super::TrajectoryResult {
            u: u,
            num_nodes: self.num_nodes,
            num_compartments: self.num_compartments,
            tspan: (*self.tspan).clone(),
            compartment_names: self.compartment_names.clone(),
        })
    }

    fn compute_rate(&self, u: &[i32], node_v: &[f64], node_idx: usize, transition_idx: usize, t: f64) -> f64 {
        let node_ldata = if let Some(ref ldata) = self.ldata {
            let nld = ldata.len() / self.num_nodes;
            let offset = node_idx * nld;
            &ldata[offset..offset + nld]
        } else {
            &[]
        };
        let rate = match &self.transitions[transition_idx].propensity_fn {
            PropensityFn::Custom { name: _, eval } => {
                eval(u, node_v, node_ldata, &self.gdata, t)
            }
        };
        if !rate.is_finite() || rate < 0.0 {
            0.0
        } else {
            rate
        }
    }

    fn compute_all_rates(&self, u: &[i32], node_v: &[f64], node_idx: usize, t: f64) -> Vec<f64> {
        (0..self.num_transitions)
            .map(|i| self.compute_rate(u, node_v, node_idx, i, t))
            .collect()
    }
}

fn build_s_matrix(transitions: &[Transition], _num_compartments: usize) -> SparseMatrix<i32> {
    let mut ir = Vec::new();
    let mut jc = Vec::new();
    let mut pr = Vec::new();

    jc.push(0);

    for tr in transitions {
        for (from, to) in tr.from.iter().zip(tr.to.iter()) {
            ir.push(from.0 as i32);
            pr.push(-1);
            ir.push(to.0 as i32);
            pr.push(1);
        }
        jc.push(ir.len() as i32);
    }

    SparseMatrix::new(ir, jc, pr)
}

fn build_g_matrix(transitions: &[Transition], num_transitions: usize) -> SparseMatrix<i32> {
    let mut ir = Vec::new();
    let mut jc = Vec::new();
    let mut pr = Vec::new();

    jc.push(0);

    for tr_idx in 0..num_transitions {
        let affected: Vec<usize> = (0..num_transitions)
            .filter(|&other_idx| {
                let tr_compartments = transitions[tr_idx].from.iter().chain(transitions[tr_idx].to.iter());
                for changed_compartment in tr_compartments {
                    for other_from in &transitions[other_idx].from {
                        if changed_compartment.0 == other_from.0 {
                            return true;
                        }
                    }
                }
                false
            })
            .collect();

        for affected_idx in &affected {
            ir.push(*affected_idx as i32);
            pr.push(1);
        }
        jc.push(ir.len() as i32);
    }

    while jc.len() < num_transitions + 1 {
        jc.push(ir.len() as i32);
    }

    SparseMatrix::new(ir, jc, pr)
}

fn apply_transition(u: &mut [i32], s: &SparseMatrix<i32>, tr: usize) -> Result<(), &'static str> {
    let start = s.jc[tr] as usize;
    let end = s.jc[tr + 1] as usize;

    for j in start..end {
        let row = s.ir[j] as usize;
        let val = s.pr[j];
        u[row] = u[row] + val;
        if u[row] < 0 {
            return Err("Negative compartment count after transition");
        }
    }
    Ok(())
}

fn store_solution(u: &mut [i32], node_states: &[i32], time_idx: usize, num_nodes: usize, nc: usize) {
    let src_len = num_nodes * nc;
    let dst_offset = time_idx * src_len;
    u[dst_offset..dst_offset + src_len].copy_from_slice(&node_states[..src_len]);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::CompartmentId;

    #[test]
    fn test_s_matrix_construction() {
        let sto_i = PropensityFn::Custom {
            name: "S_to_I".to_string(),
            eval: Arc::new(|_: &[i32], _: &[f64], _: &[f64], _: &[f64], _: f64| 0.0),
        };
        let ito_r = PropensityFn::Custom {
            name: "I_to_R".to_string(),
            eval: Arc::new(|_: &[i32], _: &[f64], _: &[f64], _: &[f64], _: f64| 0.0),
        };
        let transitions = vec![
            Transition {
                name: "S_to_I".to_string(),
                from: vec![CompartmentId(0)],
                to: vec![CompartmentId(1)],
                propensity_fn: sto_i,
            },
            Transition {
                name: "I_to_R".to_string(),
                from: vec![CompartmentId(1)],
                to: vec![CompartmentId(2)],
                propensity_fn: ito_r,
            },
        ];

        let S = build_s_matrix(&transitions, 3);

        assert_eq!(S.jc.len(), 3);
        assert_eq!(S.ir.len(), 4);
        assert_eq!(S.pr.len(), 4);
    }

    #[test]
    fn test_apply_transition() {
        let S = SparseMatrix::new(
            vec![0, 1, 1, 2],
            vec![0, 2, 4, 4],
            vec![-1, 1, -1, 1],
        );

        let mut u = vec![10, 5, 0];

        apply_transition(&mut u, &S, 0);
        assert_eq!(u[0], 9);
        assert_eq!(u[1], 6);

        apply_transition(&mut u, &S, 1);
        assert_eq!(u[1], 5);
        assert_eq!(u[2], 1);
    }
}