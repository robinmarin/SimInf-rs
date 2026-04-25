use siminf_core::{Model, CompartmentId, Transition, PropensityFn, Solver};
use std::path::Path;

fn build_sir_model() -> Model {
    Model::builder()
        .compartments(&["S", "I", "R"], &[99, 5, 0])
        .transition(Transition {
            name: "S_to_I".to_string(),
            from: vec![CompartmentId(0)],
            to: vec![CompartmentId(1)],
            propensity_fn: PropensityFn::SToI,
        })
        .transition(Transition {
            name: "I_to_R".to_string(),
            from: vec![CompartmentId(1)],
            to: vec![CompartmentId(2)],
            propensity_fn: PropensityFn::IToR,
        })
        .global_data("beta", 0.16)
        .global_data("gamma", 0.077)
        .num_nodes(1000)
        .tspan(1.0, 150.0)
        .seed(42)
        .build()
        .expect("Failed to build model")
}

#[test]
fn test_sir_model_builds_correctly() {
    let model = build_sir_model();
    assert_eq!(model.num_nodes, 1000);
    assert_eq!(model.compartments.len(), 3);
    assert_eq!(model.transitions.len(), 2);
    assert_eq!(model.gdata.get("beta"), 0.16);
    assert_eq!(model.gdata.get("gamma"), 0.077);
}

#[test]
fn test_sir_simulation_runs() {
    let model = build_sir_model();
    let mut solver = Solver::new(model);
    let result = solver.run();

    assert_eq!(result.num_nodes, 1000);
    assert_eq!(result.num_compartments, 3);
    assert_eq!(result.tspan.len(), 150);
    assert!(!result.u.is_empty());
}

#[test]
fn test_sir_conservation_of_individuals() {
    let model = build_sir_model();
    let mut solver = Solver::new(model);
    let result = solver.run();

    for t_idx in 0..result.tspan.len() {
        for node in 0..result.num_nodes {
            let offset = t_idx * result.num_nodes * result.num_compartments
                + node * result.num_compartments;
            let total: i32 = result.u[offset..offset + 3].iter().sum();
            assert_eq!(total, 104, "Total individuals should be conserved at time {}", t_idx);
        }
    }
}

#[test]
fn test_sir_initial_condition() {
    let model = build_sir_model();
    let mut solver = Solver::new(model);
    let result = solver.run();

    let offset = 0;
    assert_eq!(result.u[offset] + result.u[offset + 1] + result.u[offset + 2], 104);
}

fn load_reference_csv(path: &Path) -> Vec<(usize, i32, i32, i32, i32)> {
    let content = std::fs::read_to_string(path).expect("Failed to read reference file");
    let mut data = Vec::new();

    for line in content.lines().skip(1) {
        let line = line.trim_matches('"');
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() >= 5 {
            let node: usize = parts[0].parse().unwrap();
            let time: i32 = parts[1].parse().unwrap();
            let S: i32 = parts[2].parse().unwrap();
            let I: i32 = parts[3].parse().unwrap();
            let R: i32 = parts[4].parse().unwrap();
            data.push((node, time, S, I, R));
        }
    }

    data
}

fn compute_mean_at_time(result: &siminf_core::TrajectoryResult, t_idx: usize) -> [f64; 3] {
    let mut sum = [0.0, 0.0, 0.0];
    for node in 0..result.num_nodes {
        let offset = t_idx * result.num_nodes * result.num_compartments
            + node * result.num_compartments;
        for c in 0..3 {
            sum[c] += result.u[offset + c] as f64;
        }
    }
    for c in 0..3 {
        sum[c] /= result.num_nodes as f64;
    }
    sum
}

fn load_means_csv(path: &Path) -> Vec<(i32, f64, f64, f64)> {
    let content = std::fs::read_to_string(path).expect("Failed to read reference file");
    let mut data = Vec::new();

    for line in content.lines().skip(1) {
        let line = line.trim_matches('"');
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() >= 4 {
            let time: i32 = parts[0].parse().unwrap();
            let S: f64 = parts[1].parse().unwrap();
            let I: f64 = parts[2].parse().unwrap();
            let R: f64 = parts[3].parse().unwrap();
            data.push((time, S, I, R));
        }
    }

    data
}

#[test]
fn test_sir_validation_against_r_means() {
    let reference_path = Path::new("tests/validation/trajectory_means_r.csv");

    if !reference_path.exists() {
        eprintln!("Warning: Reference file not found at {:?}, skipping validation", reference_path);
        return;
    }

    let model = build_sir_model();
    let mut solver = Solver::new(model);
    let result = solver.run();

    let reference = load_means_csv(reference_path);

    let tolerance = 1e-6;
    let check_times = [1, 50, 100, 150];

    for t in check_times {
        let t_idx = (t - 1) as usize;
        let our_means = compute_mean_at_time(&result, t_idx);

        let ref_row = reference.iter().find(|(time, _, _, _)| *time == t as i32);

        if let Some((_, ref_S, ref_I, ref_R)) = ref_row {
            let err_S = ((our_means[0] - ref_S) / (ref_S + 1e-10)).abs();
            let err_I = ((our_means[1] - ref_I) / (ref_I + 1e-10)).abs();
            let err_R = ((our_means[2] - ref_R) / (ref_R + 1e-10)).abs();

            assert!(err_S < tolerance, "Mean S at time {}: {} vs {}, error {}", t, our_means[0], ref_S, err_S);
            assert!(err_I < tolerance, "Mean I at time {}: {} vs {}, error {}", t, our_means[1], ref_I, err_I);
            assert!(err_R < tolerance, "Mean R at time {}: {} vs {}, error {}", t, our_means[2], ref_R, err_R);
        }
    }
}