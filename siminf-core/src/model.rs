use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct CompartmentId(pub usize);

#[derive(Debug, Clone)]
pub struct Compartment {
    pub name: String,
    pub initial_count: i32,
}

pub struct Transition {
    pub name: String,
    pub from: Vec<CompartmentId>,
    pub to: Vec<CompartmentId>,
    pub propensity_fn: PropensityFn,
}

impl Clone for Transition {
    fn clone(&self) -> Self {
        Transition {
            name: self.name.clone(),
            from: self.from.clone(),
            to: self.to.clone(),
            propensity_fn: self.propensity_fn.clone(),
        }
    }
}

impl std::fmt::Debug for Transition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Transition")
            .field("name", &self.name)
            .field("from", &self.from)
            .field("to", &self.to)
            .field("propensity_fn", &self.propensity_fn)
            .finish()
    }
}

impl std::fmt::Debug for PropensityFn {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PropensityFn::SToI => write!(f, "PropensityFn::SToI"),
            PropensityFn::IToR => write!(f, "PropensityFn::IToR"),
            PropensityFn::Custom { name, .. } => write!(f, "PropensityFn::Custom({})", name),
        }
    }
}

pub enum PropensityFn {
    SToI,
    IToR,
    Custom {
        name: String,
        #[allow(dead_code)]
        eval: Box<dyn Fn(&[i32], &[f64], &[f64], f64) -> f64 + Send + Sync>,
    },
}

impl Clone for PropensityFn {
    fn clone(&self) -> Self {
        match self {
            PropensityFn::SToI => PropensityFn::SToI,
            PropensityFn::IToR => PropensityFn::IToR,
            PropensityFn::Custom { name, .. } => PropensityFn::Custom {
                name: name.clone(),
                eval: Box::new(|_, _, _,_| 0.0),
            },
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SparseMatrix {
    pub ir: Vec<i32>,
    pub jc: Vec<i32>,
    pub pr: Vec<i32>,
}

impl SparseMatrix {
    pub fn new(ir: Vec<i32>, jc: Vec<i32>, pr: Vec<i32>) -> Self {
        Self { ir, jc, pr }
    }
}

#[derive(Debug, Clone, Default)]
pub struct GlobalData {
    pub values: HashMap<String, f64>,
}

impl GlobalData {
    pub fn get(&self, name: &str) -> f64 {
        self.values.get(name).copied().unwrap_or(0.0)
    }
}

#[derive(Debug, Clone)]
pub struct LocalData {
    pub values: Vec<f64>,
}

#[derive(Debug, Clone)]
pub struct ScheduledEvent {
    pub event_type: EventType,
    pub time: f64,
    pub node: usize,
    pub dest: Option<usize>,
    pub n: usize,
    pub proportion: f64,
    pub select: i32,
    pub shift: i32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EventType {
    Exit,
    Enter,
    InternalTransfer,
    ExternalTransfer,
}

#[derive(Debug, Clone)]
pub struct Model {
    pub compartments: Vec<Compartment>,
    pub transitions: Vec<Transition>,
    pub gdata: GlobalData,
    pub ldata: Vec<LocalData>,
    pub num_nodes: usize,
    pub tspan: Vec<f64>,
    pub events: Vec<ScheduledEvent>,
    pub seed: Option<u64>,
}

impl Model {
    pub fn builder() -> ModelBuilder {
        ModelBuilder::new()
    }
}

#[derive(Debug, Default)]
pub struct ModelBuilder {
    compartments: Vec<Compartment>,
    transitions: Vec<Transition>,
    gdata: GlobalData,
    ldata: Vec<LocalData>,
    num_nodes: Option<usize>,
    tspan: Option<Vec<f64>>,
    events: Vec<ScheduledEvent>,
    seed: Option<u64>,
}

impl ModelBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn compartment(mut self, name: &str, initial_count: i32) -> Self {
        self.compartments.push(Compartment {
            name: name.to_string(),
            initial_count,
        });
        self
    }

    pub fn compartments(mut self, names: &[&str], initial_counts: &[i32]) -> Self {
        for (name, count) in names.iter().zip(initial_counts.iter()) {
            self.compartments.push(Compartment {
                name: name.to_string(),
                initial_count: *count,
            });
        }
        self
    }

    pub fn transition(mut self, transition: Transition) -> Self {
        self.transitions.push(transition);
        self
    }

    pub fn global_data(mut self, key: &str, value: f64) -> Self {
        self.gdata.values.insert(key.to_string(), value);
        self
    }

    pub fn local_data(mut self, data: Vec<f64>) -> Self {
        self.ldata.push(LocalData { values: data });
        self
    }

    pub fn num_nodes(mut self, n: usize) -> Self {
        self.num_nodes = Some(n);
        self
    }

    pub fn tspan(mut self, start: f64, end: f64) -> Self {
        let start = start as i32;
        let end = end as i32;
        self.tspan = Some((start..=end).map(|x| x as f64).collect());
        self
    }

    pub fn tspan_points(mut self, points: Vec<f64>) -> Self {
        self.tspan = Some(points);
        self
    }

    pub fn event(mut self, event: ScheduledEvent) -> Self {
        self.events.push(event);
        self
    }

    pub fn seed(mut self, seed: u64) -> Self {
        self.seed = Some(seed);
        self
    }

    pub fn build(self) -> Result<Model, String> {
        let num_nodes = self.num_nodes.ok_or("num_nodes is required")?;
        let tspan = self.tspan.ok_or("tspan is required")?;

        if self.compartments.is_empty() {
            return Err("At least one compartment is required".to_string());
        }

        if self.transitions.is_empty() {
            return Err("At least one transition is required".to_string());
        }

        Ok(Model {
            compartments: self.compartments,
            transitions: self.transitions,
            gdata: self.gdata,
            ldata: self.ldata,
            num_nodes,
            tspan,
            events: self.events,
            seed: self.seed,
        })
    }
}