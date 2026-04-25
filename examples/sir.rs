use rand::SeedableRng;
use rand::rngs::StdRng;
use siminf_core::{Model, CompartmentId, Transition, PropensityFn, Solver};

fn main() {
    let _rng = StdRng::seed_from_u64(42);

    let model = Model::builder()
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
        .build()
        .expect("Failed to build model");

    println!("Model: SimInf_model");
    println!("Number of nodes: {}", model.num_nodes);
    println!("Number of transitions: {}", model.transitions.len());
    println!("Number of scheduled events: {}", model.events.len());
    println!();
    println!("Global data");
    println!("-----------");
    println!("  Parameter Value");
    for (name, value) in &model.gdata.values {
        println!("  {:>8} {:.3}", name, value);
    }
    println!();
    println!("Compartments");
    println!("------------");

    let solver = Solver::new(model);
    println!("Running simulation...");
    let result = solver.run();

    println!("Simulation complete.");
    println!();
    println!("{}", result);

    let csv = result.to_csv();
    println!("Writing trajectory to trajectory.csv...");
    std::fs::write("trajectory.csv", csv).expect("Failed to write CSV");

    let means = result.mean_compartments();
    println!("\nMean compartment counts across nodes:");
    println!("{:>6} {:>8} {:>8} {:>8}", "time", "S", "I", "R");
    for (t_idx, t) in result.tspan.iter().enumerate() {
        println!("{:>6.0} {:>8.2} {:>8.2} {:>8.2}", 
            *t as i32, means[t_idx][0], means[t_idx][1], means[t_idx][2]);
    }
}