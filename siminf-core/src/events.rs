use crate::model::{EventType, ScheduledEvent};

#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct EventProcessor {
    events: Vec<ScheduledEvent>,
}

#[allow(dead_code)]
impl EventProcessor {
    pub fn new(events: Vec<ScheduledEvent>) -> Self {
        Self { events }
    }

    pub fn process_events_at_time(
        &self,
        time: f64,
        u: &mut [i32],
        nc: usize,
        rng: &mut impl rand::Rng,
    ) {
        for event in &self.events {
            if (event.time - time).abs() < 1e-12 {
                self.apply_event(event, u, nc, rng);
            }
        }
    }

    fn apply_event(
        &self,
        event: &ScheduledEvent,
        u: &mut [i32],
        nc: usize,
        rng: &mut impl rand::Rng,
    ) {
        match event.event_type {
            EventType::Exit => {
                let _n = self.determine_n(event, u, nc, rng);
                let offset = event.node * nc;
                for i in 0..nc {
                    u[offset + i] = u[offset + i].saturating_sub(rng.gen_range(0..=u[offset + i]));
                }
            }
            EventType::Enter => {
                let n = self.determine_n(event, u, nc, rng);
                let offset = event.node * nc;
                for i in 0..nc {
                    u[offset + i] += rng.gen_range(0..=n);
                }
            }
            EventType::InternalTransfer => {
                let n = self.determine_n(event, u, nc, rng);
                let offset = event.node * nc;
                let shift = event.shift as usize;
                if shift < nc && offset + shift < u.len() {
                    let transfer = rng.gen_range(0..=n.min(u[offset]) as i32);
                    u[offset] = u[offset].saturating_sub(transfer);
                    u[offset + shift] = u[offset + shift].saturating_add(transfer);
                }
            }
            EventType::ExternalTransfer => {
                if let Some(dest) = event.dest {
                    let n = self.determine_n(event, u, nc, rng);
                    let src_offset = event.node * nc;
                    let dst_offset = dest * nc;
                    let transfer = rng.gen_range(0..=n.min(u[src_offset]));
                    u[src_offset] = u[src_offset].saturating_sub(transfer);
                    u[dst_offset] = u[dst_offset].saturating_add(transfer);
                }
            }
        }
    }

    fn determine_n(
        &self,
        event: &ScheduledEvent,
        _u: &[i32],
        _nc: usize,
        rng: &mut impl rand::Rng,
    ) -> i32 {
        if event.n > 0 {
            event.n as i32
        } else {
            let proportion = event.proportion;
            let base: f64 = rng.gen();
            (base * proportion) as i32
        }
    }
}

impl Default for EventProcessor {
    fn default() -> Self {
        Self { events: Vec::new() }
    }
}