use std::collections::VecDeque;
use universal_kmc::{events::{Event, EventModel}, Engine};
#[derive(Default)]
struct Queue { jobs: VecDeque<String>, completed: u64, next_id: u64 }
fn main() -> Result<(), String> {
    let mut model = EventModel::new(Queue::default(), |_: &Queue| {
        let arrival = Event::new(|_: &Queue| 2.0, |q: &mut Queue| {
            q.jobs.push_back(format!("job-{}", q.next_id)); q.next_id += 1;
        });
        let mut service = Event::new(|_: &Queue| 3.0, |q: &mut Queue| { q.jobs.pop_front(); q.completed += 1; });
        service.enabled = Box::new(|q| !q.jobs.is_empty());
        vec![arrival, service]
    });
    let mut engine = Engine::new(42)?;
    engine.run(&mut model, 10000, 100.0, |_| {})?;
    println!("{} {} {}", model.data.jobs.len(), model.data.completed, engine.steps);
    Ok(())
}
