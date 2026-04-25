use rand::Rng;
use rand::rngs::StdRng;
use rand::SeedableRng;
use rayon::prelude::*;
use std::sync::Arc;

use crate::model::{
    Transition, PropensityFn, Model, SparseMatrix,
};
use crate::events::EventProcessor;

#[derive(Debug, Clone)]
#[allow(non_snake_case)]
pub struct Solver {
    model: Model,
    num_compartments: usize,
    num_transitions: usize,
    num_nodes: usize,
    S: Arc<SparseMatrix>,
    G: Arc<SparseMatrix>,
    gdata: Arc<Vec<f64>>,
    tspan: Arc<Vec<f64>>,
    events: Arc<Vec<crate::model::ScheduledEvent>>,
    select_matrix: Arc<SparseMatrix>,
    shift_matrix: Arc<SparseMatrix>,
    e2_rng: StdRng,
    seed: u64,
}

#[allow(non_snake_case)]
impl Solver {
    pub fn new(model: Model) -> Self {
        let model = model;
        let num_compartments = model.compartments.len();
        let num_transitions = model.transitions.len();
        let num_nodes = model.num_nodes;

        let s = build_s_matrix(&model.transitions, num_compartments);
        let g = build_g_matrix(&model.transitions, num_transitions);

        let gdata: Vec<f64> = model.gdata.values.values().cloned().collect();
        let seed = model.seed.unwrap_or(42);
        let events = model.events.clone();
        let tspan = model.tspan.clone();

        let select_matrix = model.select_matrix.clone().unwrap_or_else(|| {
            SparseMatrix::new(
                (0..num_compartments).map(|i| i as i32).collect(),
                {
                    let mut jc = vec![0];
                    for _ in 0..num_compartments {
                        jc.push(jc.len() as i32 + 1);
                    }
                    jc
                },
                vec![1; num_compartments],
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

        Self {
            model,
            num_compartments,
            num_transitions,
            num_nodes,
            S: Arc::new(s),
            G: Arc::new(g),
            gdata: Arc::new(gdata),
            tspan: Arc::new(tspan),
            events: Arc::new(events),
            select_matrix: Arc::new(select_matrix),
            shift_matrix: Arc::new(shift_matrix),
            e2_rng,
            seed,
        }
    }

    pub fn run(&mut self) -> super::TrajectoryResult {
        // Each node gets its own RNG seeded from the master seed so that node
        // execution order does not affect reproducibility when run in parallel.
        let mut node_rngs: Vec<StdRng> = (0..self.num_nodes)
            .map(|node| StdRng::seed_from_u64(self.seed.wrapping_add(node as u64)))
            .collect();

        let tspan = &self.tspan;
        let num_timepoints = tspan.len();
        let t_end = tspan[tspan.len() - 1];

        let mut u = vec![0i32; self.num_nodes * self.num_compartments * num_timepoints];

        let initial_state: Vec<i32> = self.model.compartments
            .iter()
            .map(|c| c.initial_count)
            .collect();

        for node in 0..self.num_nodes {
            for comp in 0..self.num_compartments {
                u[node * self.num_compartments + comp] = initial_state[comp];
            }
        }

        let mut u_current: Vec<Vec<i32>> = (0..self.num_nodes)
            .map(|node| {
                let offset = node * self.num_compartments;
                u[offset..offset + self.num_compartments].to_vec()
            })
            .collect();

        let mut t = tspan[0];
        let mut output_idx = 0;

        store_solution(&mut u, &u_current, output_idx, self.num_nodes, self.num_compartments);

        while t < t_end {
            let next_output_time = if output_idx + 1 < num_timepoints {
                tspan[output_idx + 1]
            } else {
                t_end + 1.0
            };

            let next_unit_of_time = (t.floor() + 1.0).max(t);

            let events_in_interval: Vec<_> = self.events.iter()
                .filter(|e| e.time > t && e.time <= next_output_time)
                .cloned()
                .collect();

            u_current.par_iter_mut()
                .zip(node_rngs.par_iter_mut())
                .for_each(|(node_state, rng)| {
                    let mut node_t = t;
                    let mut rates = self.compute_all_rates(node_state, node_t);

                    while node_t < next_unit_of_time {
                        let sum_t_rate: f64 = rates.iter().sum();
                        if sum_t_rate <= 0.0 {
                            break;
                        }

                        let tau = -rng.gen::<f64>().ln() / sum_t_rate;
                        if node_t + tau >= next_unit_of_time {
                            break;
                        }

                        node_t += tau;

                        let u_rand = rng.gen::<f64>() * sum_t_rate;
                        let mut cum = 0.0;
                        let mut tr = 0;
                        for (i, rate) in rates.iter().enumerate() {
                            cum += *rate;
                            if u_rand <= cum {
                                tr = i;
                                break;
                            }
                        }

                        apply_transition(node_state, &self.S, tr);

                        rates[tr] = self.compute_rate(node_state, tr, node_t);

                        for j in self.G.jc[tr] as usize..self.G.jc[tr + 1] as usize {
                            let affected_tr = self.G.ir[j] as usize;
                            rates[affected_tr] = self.compute_rate(node_state, affected_tr, node_t);
                        }
                    }
                });

            let e1_events: Vec<_> = events_in_interval.iter()
                .filter(|e| e.event_type != crate::model::EventType::ExternalTransfer)
                .cloned()
                .collect();

            for node_idx in 0..self.num_nodes {
                let node_events: Vec<_> = e1_events.iter()
                    .filter(|e| e.node == node_idx)
                    .cloned()
                    .collect();

                if !node_events.is_empty() {
                    let event_processor = EventProcessor::new(
                        self.select_matrix.as_ref(),
                        self.shift_matrix.as_ref(),
                        self.num_compartments,
                    );
                    let _ = event_processor.process_e1_events(
                        &node_events,
                        &mut u_current,
                        &mut node_rngs[node_idx],
                    );
                }
            }

            drop(e1_events);

            let e2_processor = EventProcessor::new(
                        self.select_matrix.as_ref(),
                        self.shift_matrix.as_ref(),
                        self.num_compartments,
                    );
                    for event in &events_in_interval {
                        if event.event_type == crate::model::EventType::ExternalTransfer {
                            if let Err(_) = e2_processor.apply_e2_event(event, &mut u_current, &mut self.e2_rng) {
                            }
                        }
                    }

            t = next_output_time;
            output_idx += 1;
            if output_idx < num_timepoints {
                store_solution(&mut u, &u_current, output_idx, self.num_nodes, self.num_compartments);
            }

            if t >= t_end {
                break;
            }

            while output_idx < num_timepoints && t > tspan[output_idx] {
                store_solution(&mut u, &u_current, output_idx, self.num_nodes, self.num_compartments);
                output_idx += 1;
            }
        }

        for node in 0..self.num_nodes {
            let offset = node * self.num_compartments;
            u[offset..offset + self.num_compartments].copy_from_slice(&u_current[node]);
        }

        super::TrajectoryResult {
            u: u,
            num_nodes: self.num_nodes,
            num_compartments: self.num_compartments,
            tspan: (*self.tspan).clone(),
            compartment_names: self.model.compartments.iter().map(|c| c.name.clone()).collect(),
        }
    }

    fn compute_rate(&self, u: &[i32], transition_idx: usize, _t: f64) -> f64 {
        match &self.model.transitions[transition_idx].propensity_fn {
            PropensityFn::SToI => {
                let beta = self.gdata[0];
                let s = u[0] as f64;
                let i = u[1] as f64;
                let n = (u[0] + u[1] + u[2]) as f64;
                if n > 0.0 && s > 0.0 && i > 0.0 {
                    (beta * s * i) / n
                } else {
                    0.0
                }
            }
            PropensityFn::IToR => {
                let gamma = self.gdata[1];
                let i = u[1] as f64;
                if i > 0.0 {
                    gamma * i
                } else {
                    0.0
                }
            }
            PropensityFn::Custom { name: _, eval } => {
                eval(u, &[], &self.gdata, _t)
            }
        }
    }

    fn compute_all_rates(&self, u: &[i32], t: f64) -> Vec<f64> {
        (0..self.num_transitions)
            .map(|i| self.compute_rate(u, i, t))
            .collect()
    }
}

fn build_s_matrix(transitions: &[Transition], _num_compartments: usize) -> SparseMatrix {
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

fn build_g_matrix(transitions: &[Transition], num_transitions: usize) -> SparseMatrix {
    let mut ir = Vec::new();
    let mut jc = Vec::new();
    let mut pr = Vec::new();

    jc.push(0);

    for tr_idx in 0..num_transitions {
        let affected: Vec<usize> = (0..num_transitions)
            .filter(|&other_idx| {
                for from in &transitions[tr_idx].to {
                    for other_from in &transitions[other_idx].from {
                        if from.0 == other_from.0 {
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

fn apply_transition(u: &mut [i32], s: &SparseMatrix, tr: usize) {
    let start = s.jc[tr] as usize;
    let end = s.jc[tr + 1] as usize;

    for j in start..end {
        let row = s.ir[j] as usize;
        let val = s.pr[j];
        u[row] = (u[row] as i32 + val) as i32;
    }
}

fn store_solution(u: &mut [i32], node_states: &[Vec<i32>], time_idx: usize, num_nodes: usize, nc: usize) {
    for (node, state) in node_states.iter().enumerate() {
        let offset = time_idx * num_nodes * nc + node * nc;
        u[offset..offset + nc].copy_from_slice(state);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::CompartmentId;

    #[test]
    fn test_s_matrix_construction() {
        let transitions = vec![
            Transition {
                name: "S_to_I".to_string(),
                from: vec![CompartmentId(0)],
                to: vec![CompartmentId(1)],
                propensity_fn: PropensityFn::SToI,
            },
            Transition {
                name: "I_to_R".to_string(),
                from: vec![CompartmentId(1)],
                to: vec![CompartmentId(2)],
                propensity_fn: PropensityFn::IToR,
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