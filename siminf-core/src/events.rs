use rand::Rng;
use crate::model::{EventType, ScheduledEvent, SparseMatrix};

#[derive(Debug, Clone, PartialEq)]
pub enum EventError {
    NodeOutOfBounds,
    InvalidShift,
    ShiftOutOfBounds,
    NegativeState,
    SampleError,
}

pub struct EventProcessor<'a> {
    select_matrix: &'a SparseMatrix,
    shift_matrix: &'a SparseMatrix,
    num_compartments: usize,
}

impl<'a> EventProcessor<'a> {
    pub fn new(select_matrix: &'a SparseMatrix, shift_matrix: &'a SparseMatrix, num_compartments: usize) -> Self {
        Self {
            select_matrix,
            shift_matrix,
            num_compartments,
        }
    }

    pub fn process_e1_events(
        &self,
        events: &[ScheduledEvent],
        u_current: &mut [Vec<i32>],
        rng: &mut impl Rng,
    ) -> Result<(), EventError> {
        for event in events {
            match event.event_type {
                EventType::Exit | EventType::Enter | EventType::InternalTransfer => {
                    self.apply_e1_event(event, u_current, rng)?;
                }
                EventType::ExternalTransfer => {}
            }
        }
        Ok(())
    }

    fn apply_e1_event(
        &self,
        event: &ScheduledEvent,
        u_current: &mut [Vec<i32>],
        rng: &mut impl Rng,
    ) -> Result<(), EventError> {
        let node = event.node;
        if node >= u_current.len() {
            return Ok(());
        }
        let node_state = &mut u_current[node];

        let n = self.determine_n(event, node_state, rng)?;

        match event.event_type {
            EventType::Exit => {
                self.apply_exit(event, node_state, n, rng)?;
            }
            EventType::Enter => {
                self.apply_enter(event, node_state, n, rng)?;
            }
            EventType::InternalTransfer => {
                self.apply_internal_transfer(event, node_state, n, rng)?;
            }
            EventType::ExternalTransfer => {}
        }
        Ok(())
    }

    pub fn apply_e2_event(
        &self,
        event: &ScheduledEvent,
        u_current: &mut [Vec<i32>],
        rng: &mut impl Rng,
    ) -> Result<(), EventError> {
        if event.event_type != EventType::ExternalTransfer {
            return Ok(());
        }

        let node = event.node;
        let dest = match event.dest {
            Some(d) => d,
            None => return Ok(()),
        };

        if node >= u_current.len() || dest >= u_current.len() {
            return Err(EventError::NodeOutOfBounds);
        }

        if event.shift < 0 {
            return Err(EventError::InvalidShift);
        }

        let n = self.determine_n(event, &u_current[node], rng)?;
        let individuals = self.sample_from_select(event.select as usize, &u_current[node], n, rng)?;

        for i in 0..self.select_matrix.num_rows_in_col(event.select as usize) {
            let jj = self.select_matrix.ir[self.select_matrix.jc[event.select as usize] as usize + i] as usize;
            let ll = self.shift_matrix.get_i32(event.shift as usize, i) as isize;
            let dst_idx = jj as isize + ll;

            if dst_idx < 0 || dst_idx >= self.num_compartments as isize {
                return Err(EventError::ShiftOutOfBounds);
            }

            let src_idx = jj;
            let count = individuals[src_idx];

            if u_current[node][src_idx] < count {
                return Err(EventError::NegativeState);
            }

            u_current[node][src_idx] -= count;
            u_current[dest][dst_idx as usize] += count;
        }

        Ok(())
    }

    fn determine_n(&self, event: &ScheduledEvent, node_state: &[i32], rng: &mut impl Rng) -> Result<usize, EventError> {
        if event.n > 0 {
            return Ok(event.n);
        }

        let select_col = event.select as usize;
        let mut total = 0usize;
        for i in 0..self.select_matrix.num_rows_in_col(select_col) {
            let row = self.select_matrix.ir[self.select_matrix.jc[select_col] as usize + i] as usize;
            total += node_state[row].max(0) as usize;
        }

        if total == 0 {
            return Ok(0);
        }

        let proportion = event.proportion;
        if proportion > 1.0 {
            return Err(EventError::SampleError);
        }

        let p = (rng.gen::<f64>() * (total as f64) * proportion).round() as usize;
        Ok(p.min(total))
    }

    fn sample_from_select(
        &self,
        select: usize,
        node_state: &[i32],
        n: usize,
        rng: &mut impl Rng,
    ) -> Result<Vec<i32>, EventError> {
        let mut individuals = vec![0i32; self.num_compartments];
        let n = n as i32;

        if n == 0 {
            return Ok(individuals);
        }

        let start = self.select_matrix.jc[select] as usize;
        let end = self.select_matrix.jc[select + 1] as usize;
        let num_rows = end - start;

        let mut n_individuals = 0usize;
        let mut n_kinds = 0usize;

        for i in start..end {
            let row = self.select_matrix.ir[i] as usize;
            let count = node_state[row].max(0) as usize;
            if count > 0 {
                n_kinds += 1;
            }
            n_individuals += count;
        }

        if n_individuals == n as usize {
            for i in start..end {
                let row = self.select_matrix.ir[i] as usize;
                individuals[row] = node_state[row].max(0);
            }
            return Ok(individuals);
        }

        if num_rows == 1 {
            individuals[self.select_matrix.ir[start] as usize] = n;
            return Ok(individuals);
        }

        if n_kinds == 1 {
            for i in start..end {
                let row = self.select_matrix.ir[i] as usize;
                if node_state[row] > 0 {
                    individuals[row] = n;
                    break;
                }
            }
            return Ok(individuals);
        }

        let pr_f64 = self.select_matrix.pr_f64.as_deref();
        let all_weights_equal = if let Some(pr) = pr_f64 {
            if pr.is_empty() {
                true
            } else {
                let first_weight = pr[start];
                let mut equal = true;
                for i in (start + 1)..end {
                    if (pr[i] - first_weight).abs() > 1e-12 {
                        equal = false;
                        break;
                    }
                }
                equal
            }
        } else {
            true
        };

        let no_weights = pr_f64.map(|pr| pr.is_empty()).unwrap_or(true);

        if no_weights || all_weights_equal {
            let mut remaining = n;
            let total_available = n_individuals as i32;

            for i in start..end {
                if remaining == 0 {
                    break;
                }
                let row = self.select_matrix.ir[i] as usize;
                let count = node_state[row].max(0) as i32;

                if count == 0 {
                    continue;
                }

                let available = (total_available - individuals.iter().sum::<i32>()).max(0);
                if available == 0 {
                    break;
                }

                let sampled = if remaining >= count && count <= available {
                    count
                } else if remaining >= available {
                    available
                } else {
                    let prob = remaining as f64 / available as f64;
                    if rng.gen::<f64>() < prob {
                        remaining
                    } else {
                        0
                    }
                };

                individuals[row] = sampled;
                remaining -= sampled;
            }

            if remaining > 0 {
                for i in start..end {
                    let row = self.select_matrix.ir[i] as usize;
                    let space = node_state[row].max(0) - individuals[row];
                    if space > 0 {
                        let add = remaining.min(space);
                        individuals[row] += add;
                        remaining -= add;
                    }
                    if remaining == 0 {
                        break;
                    }
                }
            }
        } else {
            individuals = self.sample_biased_urn(select, node_state, n, rng)?;
        }

        Ok(individuals)
    }

    fn sample_biased_urn(&self, select: usize, node_state: &[i32], n: i32, rng: &mut impl Rng) -> Result<Vec<i32>, EventError> {
        let mut individuals = vec![0i32; self.num_compartments];
        let mut remaining = n;

        let start = self.select_matrix.jc[select] as usize;
        let end = self.select_matrix.jc[select + 1] as usize;
        let pr = self.select_matrix.pr_f64.as_deref().unwrap_or(&[]);

        while remaining > 0 {
            let mut cum = 0.0;
            for i in start..end {
                let row = self.select_matrix.ir[i] as usize;
                let weight = if pr.is_empty() { 1.0 } else { pr[i] };
                cum += weight * (node_state[row].max(0) as f64 - individuals[row] as f64);
            }

            if cum <= 0.0 {
                break;
            }

            let mut rand = rng.gen::<f64>() * cum;
            let mut i = start;
            while i < end {
                let row = self.select_matrix.ir[i] as usize;
                let weight = if pr.is_empty() { 1.0 } else { pr[i] };
                let w = weight * (node_state[row].max(0) as f64 - individuals[row] as f64);

                if rand <= w {
                    individuals[row] += 1;
                    remaining -= 1;
                    break;
                }
                rand -= w;
                i += 1;
            }

            if i >= end && start < end {
                let row = self.select_matrix.ir[end - 1] as usize;
                individuals[row] += 1;
                remaining -= 1;
            }
        }

        Ok(individuals)
    }

    fn apply_exit(&self, event: &ScheduledEvent, node_state: &mut [i32], n: usize, rng: &mut impl Rng) -> Result<(), EventError> {
        let individuals = self.sample_from_select(event.select as usize, node_state, n, rng)?;

        for i in 0..self.select_matrix.num_rows_in_col(event.select as usize) {
            let row = self.select_matrix.ir[self.select_matrix.jc[event.select as usize] as usize + i] as usize;
            if node_state[row] < individuals[row] {
                return Err(EventError::NegativeState);
            }
            node_state[row] -= individuals[row];
        }
        Ok(())
    }

    fn apply_enter(&self, event: &ScheduledEvent, node_state: &mut [i32], n: usize, rng: &mut impl Rng) -> Result<(), EventError> {
        let individuals = self.sample_from_select(event.select as usize, node_state, n, rng)?;

        let start = self.select_matrix.jc[event.select as usize] as usize;
        let end = self.select_matrix.jc[event.select as usize + 1] as usize;

        if event.shift < 0 {
            for i in start..end {
                let row = self.select_matrix.ir[i] as usize;
                node_state[row] += individuals[row];
            }
        } else {
            for i in start..end {
                let row = self.select_matrix.ir[i] as usize;
                let ll = self.shift_matrix.get_i32(event.shift as usize, i) as isize;
                let dst = row as isize + ll;
                if dst < 0 || dst >= self.num_compartments as isize {
                    return Err(EventError::ShiftOutOfBounds);
                }
                node_state[dst as usize] += individuals[row];
            }
        }
        Ok(())
    }

    fn apply_internal_transfer(&self, event: &ScheduledEvent, node_state: &mut [i32], n: usize, rng: &mut impl Rng) -> Result<(), EventError> {
        if event.shift < 0 {
            return Err(EventError::InvalidShift);
        }

        let individuals = self.sample_from_select(event.select as usize, node_state, n, rng)?;

        let start = self.select_matrix.jc[event.select as usize] as usize;
        let end = self.select_matrix.jc[event.select as usize + 1] as usize;

        for i in start..end {
            let row = self.select_matrix.ir[i] as usize;
            let ll = self.shift_matrix.get_i32(event.shift as usize, i) as isize;
            let dst = row as isize + ll;

            if dst < 0 || dst >= self.num_compartments as isize {
                return Err(EventError::ShiftOutOfBounds);
            }

            let count = individuals[row];
            if node_state[row] < count {
                return Err(EventError::NegativeState);
            }
            node_state[row] -= count;
            node_state[dst as usize] += count;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    fn make_select_matrix(ir: Vec<i32>, jc: Vec<i32>) -> SparseMatrix {
        SparseMatrix::new(ir, jc, vec![])
    }

    fn make_shift_matrix(ir: Vec<i32>, jc: Vec<i32>, pr: Vec<i32>) -> SparseMatrix {
        SparseMatrix::new(ir, jc, pr)
    }

    #[test]
    fn test_determine_n_with_absolute_count() {
        let select = make_select_matrix(vec![0, 1, 2], vec![0, 3, 3, 3]);
        let shift = make_shift_matrix(vec![], vec![0, 0], vec![]);
        let ep = EventProcessor::new(&select, &shift, 3);

        let node_state = vec![10, 20, 30];
        let mut rng = StdRng::seed_from_u64(42);

        let event = ScheduledEvent {
            event_type: EventType::Exit,
            time: 1.0,
            node: 0,
            dest: None,
            n: 5,
            proportion: 0.0,
            select: 0,
            shift: -1,
        };

        let n = ep.determine_n(&event, &node_state, &mut rng).unwrap();
        assert_eq!(n, 5);
    }

    #[test]
    fn test_determine_n_with_proportion() {
        let select = make_select_matrix(vec![0, 1, 2], vec![0, 3, 3, 3]);
        let shift = make_shift_matrix(vec![], vec![0, 0], vec![]);
        let ep = EventProcessor::new(&select, &shift, 3);

        let node_state = vec![10, 20, 30];
        let mut rng = StdRng::seed_from_u64(42);

        let event = ScheduledEvent {
            event_type: EventType::Exit,
            time: 1.0,
            node: 0,
            dest: None,
            n: 0,
            proportion: 0.5,
            select: 0,
            shift: -1,
        };

        let n = ep.determine_n(&event, &node_state, &mut rng).unwrap();
        let total: i32 = node_state.iter().sum();
        assert!(n <= (total as f64 * 0.5).ceil() as usize + 2);
        assert!(n <= total as usize);
    }

    #[test]
    fn test_exit_event_removes_individuals() {
        let select = make_select_matrix(vec![0, 1, 2], vec![0, 3, 3, 3]);
        let shift = make_shift_matrix(vec![], vec![0, 0], vec![]);
        let ep = EventProcessor::new(&select, &shift, 3);

        let mut node_state = vec![10, 20, 30];
        let mut rng = StdRng::seed_from_u64(42);

        let event = ScheduledEvent {
            event_type: EventType::Exit,
            time: 1.0,
            node: 0,
            dest: None,
            n: 5,
            proportion: 0.0,
            select: 0,
            shift: -1,
        };

        let n = ep.determine_n(&event, &node_state, &mut rng).unwrap();
        ep.apply_exit(&event, &mut node_state, n, &mut rng).unwrap();

        let total: i32 = node_state.iter().sum();
        assert_eq!(total, 60 - n as i32);
    }

    #[test]
    fn test_enter_event_adds_individuals() {
        let select = make_select_matrix(vec![0, 1, 2], vec![0, 3, 3, 3]);
        let shift = make_shift_matrix(vec![], vec![0, 0], vec![]);
        let ep = EventProcessor::new(&select, &shift, 3);

        let mut node_state = vec![10, 20, 30];
        let mut rng = StdRng::seed_from_u64(42);

        let event = ScheduledEvent {
            event_type: EventType::Enter,
            time: 1.0,
            node: 0,
            dest: None,
            n: 5,
            proportion: 0.0,
            select: 0,
            shift: -1,
        };

        let n = ep.determine_n(&event, &node_state, &mut rng).unwrap();
        ep.apply_enter(&event, &mut node_state, n, &mut rng).unwrap();

        let total: i32 = node_state.iter().sum();
        assert_eq!(total, 60 + n as i32);
    }

    #[test]
    fn test_internal_transfer_shifts_compartments() {
        let select = make_select_matrix(vec![0, 1, 2], vec![0, 2, 2, 2]);
        let shift = make_shift_matrix(vec![0, 1], vec![0, 2, 2], vec![1, 1]);
        let ep = EventProcessor::new(&select, &shift, 3);

        let mut node_state = vec![10, 20, 0];
        let mut rng = StdRng::seed_from_u64(42);

        let event = ScheduledEvent {
            event_type: EventType::InternalTransfer,
            time: 1.0,
            node: 0,
            dest: None,
            n: 5,
            proportion: 0.0,
            select: 0,
            shift: 0,
        };

        let n = ep.determine_n(&event, &node_state, &mut rng).unwrap();
        let before_total: i32 = node_state.iter().sum();
        ep.apply_internal_transfer(&event, &mut node_state, n, &mut rng).unwrap();
        let after_total: i32 = node_state.iter().sum();

        assert_eq!(before_total, after_total);
        assert_eq!(node_state[0], 10 - n as i32);
        assert_eq!(node_state[1], 20 + n as i32);
    }

    #[test]
    fn test_sample_from_select_all_individuals() {
        let select = make_select_matrix(vec![0, 1, 2], vec![0, 3, 3, 3]);
        let shift = make_shift_matrix(vec![], vec![0, 0], vec![]);
        let ep = EventProcessor::new(&select, &shift, 3);

        let node_state = vec![10, 20, 30];
        let mut rng = StdRng::seed_from_u64(42);

        let individuals = ep.sample_from_select(0, &node_state, 60, &mut rng).unwrap();

        assert_eq!(individuals[0], 10);
        assert_eq!(individuals[1], 20);
        assert_eq!(individuals[2], 30);
    }

    #[test]
    fn test_sample_from_select_single_compartment() {
        let select = make_select_matrix(vec![0, 1, 2], vec![0, 1, 1, 1]);
        let shift = make_shift_matrix(vec![], vec![0, 0], vec![]);
        let ep = EventProcessor::new(&select, &shift, 3);

        let node_state = vec![10, 20, 30];
        let mut rng = StdRng::seed_from_u64(42);

        let individuals = ep.sample_from_select(0, &node_state, 5, &mut rng).unwrap();

        assert_eq!(individuals[0], 5);
        assert_eq!(individuals[1], 0);
        assert_eq!(individuals[2], 0);
    }
}