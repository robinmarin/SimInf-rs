use serde::{Deserialize, Serialize};
use std::fmt;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrajectoryResult {
    pub u: Vec<i32>,
    pub num_nodes: usize,
    pub num_compartments: usize,
    pub tspan: Vec<f64>,
    pub compartment_names: Vec<String>,
}

impl TrajectoryResult {
    pub fn trajectory(&self) -> TrajectoryView<'_> {
        TrajectoryView {
            result: self,
            time_index: 0,
        }
    }

    pub fn mean_compartments(&self) -> Vec<Vec<f64>> {
        let num_timepoints = self.tspan.len();
        let mut means = vec![vec![0.0; self.num_compartments]; num_timepoints];

        for t in 0..num_timepoints {
            for node in 0..self.num_nodes {
                for c in 0..self.num_compartments {
                    let idx = t * self.num_nodes * self.num_compartments + node * self.num_compartments + c;
                    means[t][c] += self.u[idx] as f64;
                }
            }
            for c in 0..self.num_compartments {
                means[t][c] /= self.num_nodes as f64;
            }
        }

        means
    }

    pub fn to_csv(&self) -> String {
        let mut s = String::new();
        s.push_str("node,time,");
        s.push_str(&self.compartment_names.join(","));
        s.push('\n');

        for (t_idx, &t) in self.tspan.iter().enumerate() {
            for node in 0..self.num_nodes {
                s.push_str(&format!("{},{}", node + 1, t as i32));
                for c in 0..self.num_compartments {
                    let idx = t_idx * self.num_nodes * self.num_compartments + node * self.num_compartments + c;
                    s.push_str(&format!(",{}", self.u[idx]));
                }
                s.push('\n');
            }
        }

        s
    }
}

#[derive(Debug)]
pub struct TrajectoryView<'a> {
    result: &'a TrajectoryResult,
    time_index: usize,
}

impl<'a> TrajectoryView<'a> {
    pub fn time(mut self, idx: usize) -> Self {
        self.time_index = idx;
        self
    }

    pub fn node_state(&self, node: usize) -> Vec<i32> {
        let offset = self.time_index * self.result.num_nodes * self.result.num_compartments
            + node * self.result.num_compartments;
        self.result.u[offset..offset + self.result.num_compartments].to_vec()
    }

    pub fn all_nodes_at_time(&self, t_idx: usize) -> Vec<Vec<i32>> {
        let offset = t_idx * self.result.num_nodes * self.result.num_compartments;
        (0..self.result.num_nodes)
            .map(|node| {
                let off = offset + node * self.result.num_compartments;
                self.result.u[off..off + self.result.num_compartments].to_vec()
            })
            .collect()
    }
}

impl fmt::Display for TrajectoryResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "TrajectoryResult {{")?;
        writeln!(f, "  num_nodes: {}", self.num_nodes)?;
        writeln!(f, "  num_compartments: {}", self.num_compartments)?;
        writeln!(f, "  num_timepoints: {}", self.tspan.len())?;
        writeln!(f, "  tspan: {:?}", self.tspan)?;
        writeln!(f, "  compartment_names: {:?}", self.compartment_names)?;
        writeln!(f, "}}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_trajectory_result_creation() {
        let result = TrajectoryResult {
            u: vec![99, 5, 0, 98, 6, 0, 97, 7, 0, 96, 8, 0],
            num_nodes: 2,
            num_compartments: 3,
            tspan: vec![1.0, 2.0],
            compartment_names: vec!["S".to_string(), "I".to_string(), "R".to_string()],
        };

        let means = result.mean_compartments();
        assert_eq!(means.len(), 2);
        assert!((means[0][0] - 98.5).abs() < 1e-6);
    }
}