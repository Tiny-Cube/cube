use std::io::{self, Write};
use universal_kmc::{ensure, io::{integer, loads, number}, Engine, Result};
fn run() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    ensure((2..=5).contains(&args.len()), "usage: ukmc-rust MODEL [STEPS [SEED [UNTIL]]]")?;
    let steps = if args.len() > 2 { integer(&args[2])? } else { 1000 };
    let seed = if args.len() > 3 { integer(&args[3])? } else { 1 };
    let until = if args.len() > 4 { let t = number(&args[4])?; ensure(t.is_finite() && t >= 0.0, "invalid time limit")?; t } else { f64::INFINITY };
    let mut model = loads(&std::fs::read_to_string(&args[1]).map_err(|e| e.to_string())?)?;
    let mut engine = Engine::new(seed)?;
    let mut out = io::BufWriter::new(io::stdout().lock());
    write!(out, "step,time,event,x,y,z").map_err(|e| e.to_string())?;
    for i in 0..model.state().len() { write!(out, ",state{i}").map_err(|e| e.to_string())?; }
    writeln!(out).map_err(|e| e.to_string())?;
    let mut output_error = None;
    engine.run(model.as_mut(), steps, until, |r| {
        let result = (|| -> io::Result<()> {
            write!(out, "{},{:.17e},{}", r.step, r.time, r.event)?;
            for x in r.position { write!(out, ",{x:.17e}")?; }
            for n in &r.state { write!(out, ",{n}")?; }
            writeln!(out)
        })();
        if let Err(e) = result { output_error = Some(e); }
    })?;
    if let Some(e) = output_error { return Err(e.to_string()); }
    out.flush().map_err(|e| e.to_string())?;
    eprintln!("status={} events={} time={:.17e}", engine.status, engine.steps, engine.time);
    Ok(())
}
fn main() { if let Err(e) = run() { eprintln!("ukmc: {e}"); std::process::exit(2); } }
