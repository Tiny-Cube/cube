use universal_kmc::{io::loads, Engine};
fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mut model = loads(&std::fs::read_to_string(&args[1]).unwrap()).unwrap();
    let mut engine = Engine::new(args[3].parse().unwrap()).unwrap();
    let steps = args[2].parse().unwrap();
    let mut checksum: i64 = 0;
    let started = std::time::Instant::now();
    engine.run(model.as_mut(), steps, f64::INFINITY,
               |row| checksum += row.event + 1).unwrap();
    let elapsed = started.elapsed().as_secs_f64();
    print!("{:.17e} {} {:.17e} {}", elapsed, engine.steps, engine.time, checksum);
    for x in model.position() { print!(" {:.17e}", x); }
    for x in model.state() { print!(" {}", x); }
    println!();
}
