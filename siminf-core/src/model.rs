use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::sync::Arc;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct CompartmentId(pub usize);

#[derive(Debug, Clone)]
pub struct Compartment {
    pub name: String,
    pub initial_count: i32,
}

#[derive(Clone)]
pub struct Transition {
    pub name: String,
    pub from: Vec<CompartmentId>,
    pub to: Vec<CompartmentId>,
    pub propensity_fn: PropensityFn,
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
            PropensityFn::Custom { name, .. } => write!(f, "PropensityFn::Custom({})", name),
        }
    }
}

#[derive(Clone)]
pub enum PropensityFn {
    Custom {
        name: String,
        eval: Arc<dyn Fn(&[i32], &[f64], &[f64], &[f64], f64) -> f64 + Send + Sync>,
    },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SparseMatrix<T> {
    pub ir: Vec<i32>,
    pub jc: Vec<i32>,
    pub pr: Vec<T>,
}

impl<T> SparseMatrix<T> {
    pub fn new(ir: Vec<i32>, jc: Vec<i32>, pr: Vec<T>) -> Self {
        Self { ir, jc, pr }
    }

    pub fn num_rows_in_col(&self, col: usize) -> usize {
        if col + 1 >= self.jc.len() {
            return 0;
        }
        (self.jc[col + 1] - self.jc[col]) as usize
    }
}

impl SparseMatrix<f64> {
    pub fn get_f64(&self, col: usize, row_in_col: usize) -> f64 {
        let start = self.jc[col] as usize;
        self.pr[start + row_in_col]
    }
}

impl SparseMatrix<i32> {
    pub fn get_i32(&self, col: usize, row_in_col: usize) -> i32 {
        let start = self.jc[col] as usize;
        self.pr[start + row_in_col]
    }
}

#[derive(Debug, Clone, Default)]
pub struct GlobalData {
    pub values: BTreeMap<String, f64>,
}

impl GlobalData {
    pub fn get(&self, name: &str) -> f64 {
        self.values.get(name).copied().unwrap_or(0.0)
    }

    pub fn to_vec(&self) -> Vec<f64> {
        self.values.values().cloned().collect()
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

impl ScheduledEvent {
    pub fn new(
        event_type: EventType,
        time: f64,
        node: usize,
    ) -> Self {
        Self {
            event_type,
            time,
            node,
            dest: None,
            n: 0,
            proportion: 0.0,
            select: 0,
            shift: -1,
        }
    }

    pub fn n(mut self, n: usize) -> Self {
        self.n = n;
        self
    }

    pub fn proportion(mut self, proportion: f64) -> Self {
        self.proportion = proportion;
        self
    }

    pub fn select(mut self, select: i32) -> Self {
        self.select = select;
        self
    }

    pub fn shift(mut self, shift: i32) -> Self {
        self.shift = shift;
        self
    }

    pub fn dest(mut self, dest: usize) -> Self {
        self.dest = Some(dest);
        self
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EventType {
    Exit,
    Enter,
    InternalTransfer,
    ExternalTransfer,
}

#[derive(Clone)]
pub struct Model {
    pub compartments: Vec<Compartment>,
    pub transitions: Vec<Transition>,
    pub gdata: GlobalData,
    pub ldata: Vec<LocalData>,
    pub num_nodes: usize,
    pub tspan: Vec<f64>,
    pub events: Vec<ScheduledEvent>,
    pub seed: Option<u64>,
    pub select_matrix: Option<SparseMatrix<f64>>,
    pub shift_matrix: Option<SparseMatrix<i32>>,
    pub u0: Option<Vec<i32>>,
    pub nd: usize,
    pub v0: Option<Vec<f64>>,
    pub pts_fun: Option<Arc<dyn Fn(&mut [f64], &[i32], &[f64], &[f64], &[f64], usize, f64) + Send + Sync>>,
}

impl std::fmt::Debug for Model {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Model")
            .field("compartments", &self.compartments)
            .field("transitions", &self.transitions)
            .field("gdata", &self.gdata)
            .field("ldata", &self.ldata)
            .field("num_nodes", &self.num_nodes)
            .field("tspan", &self.tspan)
            .field("events", &self.events)
            .field("seed", &self.seed)
            .field("select_matrix", &self.select_matrix)
            .field("shift_matrix", &self.shift_matrix)
            .field("u0", &self.u0)
            .field("nd", &self.nd)
            .field("v0", &self.v0)
            .field("pts_fun", &"Option<Arc<dyn Fn(...)>>")
            .finish()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ModelError {
    MissingField(&'static str),
    Validation(&'static str),
}

impl std::fmt::Display for ModelError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ModelError::MissingField(field) => write!(f, "Missing required field: {}", field),
            ModelError::Validation(msg) => write!(f, "Validation error: {}", msg),
        }
    }
}

impl std::error::Error for ModelError {}

impl Model {
    pub fn builder() -> ModelBuilder {
        ModelBuilder::new()
    }
}

#[derive(Default)]
pub struct ModelBuilder {
    compartments: Vec<Compartment>,
    transitions: Vec<Transition>,
    gdata: GlobalData,
    ldata: Vec<LocalData>,
    num_nodes: Option<usize>,
    tspan: Option<Vec<f64>>,
    events: Vec<ScheduledEvent>,
    seed: Option<u64>,
    select_matrix: Option<SparseMatrix<f64>>,
    shift_matrix: Option<SparseMatrix<i32>>,
    u0: Option<Vec<i32>>,
    nd: Option<usize>,
    v0: Option<Vec<f64>>,
    pts_fun: Option<Arc<dyn Fn(&mut [f64], &[i32], &[f64], &[f64], &[f64], usize, f64) + Send + Sync>>,
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

    pub fn gdata_struct(mut self, gd: GlobalData) -> Self {
        self.gdata = gd;
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

    pub fn select_matrix(mut self, matrix: SparseMatrix<f64>) -> Self {
        self.select_matrix = Some(matrix);
        self
    }

    pub fn shift_matrix(mut self, matrix: SparseMatrix<i32>) -> Self {
        self.shift_matrix = Some(matrix);
        self
    }

    pub fn u0(mut self, u0: Vec<i32>) -> Self {
        self.u0 = Some(u0);
        self
    }

    pub fn nd(mut self, nd: usize) -> Self {
        self.nd = Some(nd);
        self
    }

    pub fn v0(mut self, v0: Vec<f64>) -> Self {
        self.v0 = Some(v0);
        self
    }

    pub fn pts_fun<F>(mut self, pts_fun: F) -> Self
    where
        F: Fn(&mut [f64], &[i32], &[f64], &[f64], &[f64], usize, f64) + Send + Sync + 'static,
    {
        self.pts_fun = Some(Arc::new(pts_fun));
        self
    }

    #[must_use]
    pub fn build(self) -> Result<Model, ModelError> {
        let num_nodes = self.num_nodes.ok_or(ModelError::MissingField("num_nodes"))?;
        let tspan = self.tspan.ok_or(ModelError::MissingField("tspan"))?;

        if self.compartments.is_empty() {
            return Err(ModelError::Validation("At least one compartment is required"));
        }

        if self.transitions.is_empty() {
            return Err(ModelError::Validation("At least one transition is required"));
        }

        let num_compartments = self.compartments.len();
        if !self.ldata.is_empty() {
            let first_len = self.ldata[0].values.len();
            for ld in &self.ldata {
                if ld.values.len() != first_len {
                    return Err(ModelError::Validation("All LocalData entries must have the same number of values"));
                }
            }
        }

        if let Some(ref u0) = self.u0 {
            if u0.len() != num_nodes * num_compartments {
                return Err(ModelError::Validation("u0 length must be num_nodes * num_compartments"));
            }
        }

        if !self.events.is_empty() && self.select_matrix.is_none() {
            return Err(ModelError::Validation("select_matrix is required when events are present"));
        }

        let nd = self.nd.unwrap_or(0);

        if let Some(ref v0) = self.v0 {
            if v0.len() != num_nodes * nd {
                return Err(ModelError::Validation("v0 length must be num_nodes * nd"));
            }
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
            select_matrix: self.select_matrix,
            shift_matrix: self.shift_matrix,
            u0: self.u0,
            nd,
            v0: self.v0,
            pts_fun: self.pts_fun,
        })
    }
}