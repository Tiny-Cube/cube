//! Exact direct SSA. See the repository's docs/algorithm.md and docs/format.md.
pub mod events;
pub mod io;
pub mod rates;
pub const MAX_COUNT: u64 = 9_007_199_254_740_991;
pub type Result<T> = std::result::Result<T, String>;
pub fn ensure(ok: bool, message: &str) -> Result<()> {
    if ok { Ok(()) } else { Err(message.into()) }
}
pub fn count(n: u64) -> Result<u64> { ensure(n <= MAX_COUNT, "count exceeds 2^53-1")?; Ok(n) }
pub fn nonnegative(a: f64) -> Result<f64> { ensure(a.is_finite() && a >= 0.0, "rate must be finite and nonnegative")?; Ok(a) }

pub trait RandomSource { fn uniform(&mut self) -> f64; }
pub struct PortableRng { pub state: u64 }
impl PortableRng {
    pub fn new(seed: u64) -> Result<Self> { Ok(Self { state: count(seed)? % 2_147_483_646 + 1 }) }
}
impl RandomSource for PortableRng {
    fn uniform(&mut self) -> f64 {
        self.state = (48271 * self.state) % 2_147_483_647;
        self.state as f64 / 2_147_483_647.0
    }
}
pub trait Model {
    fn propensities(&self) -> Result<Vec<f64>>;
    fn fire(&mut self, event: usize) -> Result<()>;
    fn state(&self) -> Vec<u64> { vec![] }
    fn position(&self) -> [f64; 3] { [0.0; 3] }
}
#[derive(Debug, Clone)]
pub struct Record { pub step: u64, pub time: f64, pub event: i64, pub position: [f64; 3], pub state: Vec<u64> }
pub struct Engine<R = PortableRng> {
    pub rng: R, pub time: f64, pub steps: u64, pub status: &'static str,
}
impl Engine<PortableRng> {
    pub fn new(seed: u64) -> Result<Self> { Ok(Self::with_rng(PortableRng::new(seed)?)) }
}
impl<R: RandomSource> Engine<R> {
    pub fn with_rng(rng: R) -> Self { Self { rng, time: 0.0, steps: 0, status: "ready" } }
    fn uniform(&mut self) -> Result<f64> {
        let u = self.rng.uniform();
        ensure(u.is_finite() && u > 0.0 && u < 1.0, "RNG must return a value in (0,1)")?;
        Ok(u)
    }
    pub fn snapshot(&self, model: &dyn Model, event: i64) -> Record {
        Record { step: self.steps, time: self.time, event, position: model.position(), state: model.state() }
    }
    pub fn step(&mut self, model: &mut dyn Model, until: f64) -> Result<Option<Record>> {
        ensure(!until.is_nan() && until >= self.time, "until must be >= current time")?;
        if self.time == until { self.status = "time_limit"; return Ok(None); }
        let rates = model.propensities()?;
        let mut total = 0.0;
        for a in &rates { total += nonnegative(*a)?; }
        ensure(total.is_finite(), "total propensity overflow")?;
        if total == 0.0 {
            self.status = "absorbed";
            if until.is_finite() { self.time = until; }
            return Ok(None);
        }
        let dt = -self.uniform()?.ln() / total;
        if dt > until - self.time { self.time = until; self.status = "time_limit"; return Ok(None); }
        ensure(dt.is_finite() && (self.time + dt).is_finite() && self.time + dt > self.time, "waiting time outside floating-point range")?;
        let threshold = self.uniform()? * total;
        let mut cumulative = 0.0;
        let mut chosen = 0;
        for (i, a) in rates.iter().enumerate() {
            if *a > 0.0 { chosen = i; }
            cumulative += a;
            if threshold < cumulative { break; }
        }
        model.fire(chosen)?;
        self.time += dt;
        self.steps += 1;
        self.status = "running";
        Ok(Some(self.snapshot(model, chosen as i64)))
    }
    pub fn run(&mut self, model: &mut dyn Model, max_steps: u64, until: f64, mut observe: impl FnMut(&Record)) -> Result<()> {
        count(max_steps)?;
        ensure(!until.is_nan() && until >= self.time, "until must be >= current time")?;
        observe(&self.snapshot(model, -1));
        for _ in 0..max_steps {
            let before = self.time;
            match self.step(model, until)? {
                Some(row) => observe(&row),
                None => {
                    if self.time != before { observe(&self.snapshot(model, -1)); }
                    return Ok(());
                }
            }
        }
        self.status = if self.time == until { "time_limit" } else { "step_limit" };
        Ok(())
    }
}
#[derive(Clone, Debug)]
pub struct Edge { pub source: u64, pub target: u64, pub rate: f64, pub displacement: [f64; 3] }
pub struct Graph { edges: Vec<Edge>, site: u64, xyz: [f64; 3] }
impl Graph {
    pub fn new(sites: u64, edges: Vec<Edge>, initial: u64) -> Result<Self> {
        ensure(count(sites)? > 0 && count(initial)? < sites, "invalid graph size or initial site")?;
        for e in &edges {
            ensure(count(e.source)? < sites && count(e.target)? < sites, "edge endpoint out of bounds")?;
            nonnegative(e.rate)?;
            ensure(e.displacement.iter().all(|x| x.is_finite()), "non-finite displacement")?;
        }
        Ok(Self { edges, site: initial, xyz: [0.0; 3] })
    }
}
impl Model for Graph {
    fn propensities(&self) -> Result<Vec<f64>> {
        Ok(self.edges.iter().map(|e| if e.source == self.site { e.rate } else { 0.0 }).collect())
    }
    fn fire(&mut self, event: usize) -> Result<()> {
        let e = self.edges.get(event).ok_or("event out of bounds")?;
        ensure(e.source == self.site && e.rate > 0.0, "inactive event")?;
        let mut next = self.xyz;
        for (x, dx) in next.iter_mut().zip(e.displacement) {
            *x += dx;
            ensure(x.is_finite(), "position overflow")?;
        }
        self.xyz = next; self.site = e.target;
        Ok(())
    }
    fn state(&self) -> Vec<u64> { vec![self.site] }
    fn position(&self) -> [f64; 3] { self.xyz }
}
#[derive(Clone, Debug)]
pub struct Reaction { pub rate: f64, pub reactants: Vec<u64>, pub products: Vec<u64> }
pub struct ReactionNetwork { populations: Vec<u64>, reactions: Vec<Reaction> }
impl ReactionNetwork {
    pub fn new(populations: Vec<u64>, reactions: Vec<Reaction>) -> Result<Self> {
        ensure(!populations.is_empty(), "at least one species required")?;
        for n in &populations { count(*n)?; }
        for r in &reactions {
            nonnegative(r.rate)?;
            ensure(r.reactants.len() == populations.len() && r.products.len() == populations.len(), "stoichiometry dimension mismatch")?;
            for n in r.reactants.iter().chain(&r.products) { count(*n)?; }
        }
        Ok(Self { populations, reactions })
    }
}
impl Model for ReactionNetwork {
    fn propensities(&self) -> Result<Vec<f64>> {
        let mut result = Vec::new();
        for r in &self.reactions {
            if r.rate == 0.0 || self.populations.iter().zip(&r.reactants).any(|(n, a)| n < a) {
                result.push(0.0); continue;
            }
            let mut a = r.rate;
            for (n, order) in self.populations.iter().zip(&r.reactants) {
                for k in 1..=(*order).min(n - order) {
                    a = nonnegative(a * ((n - k + 1) as f64 / k as f64))?;
                }
            }
            result.push(a);
        }
        Ok(result)
    }
    fn fire(&mut self, event: usize) -> Result<()> {
        let r = self.reactions.get(event).ok_or("event out of bounds")?;
        ensure(r.rate > 0.0, "inactive reaction")?;
        let mut next = Vec::new();
        for ((n, a), b) in self.populations.iter().zip(&r.reactants).zip(&r.products) {
            ensure(n >= a, "inactive reaction")?;
            next.push(count(n - a + b)?);
        }
        self.populations = next;
        Ok(())
    }
    fn state(&self) -> Vec<u64> { self.populations.clone() }
    fn position(&self) -> [f64; 3] { [0.0; 3] }
}
