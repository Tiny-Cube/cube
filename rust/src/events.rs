//! Domain-free events over arbitrary user-defined Rust state.
use std::cell::RefCell;
use crate::{ensure, nonnegative, Model, Result};
pub struct Event<S> {
    pub rate: Box<dyn Fn(&S) -> f64>,
    pub action: Box<dyn Fn(&mut S)>,
    pub enabled: Box<dyn Fn(&S) -> bool>,
    pub name: String,
}
impl<S> Event<S> {
    pub fn new(rate: impl Fn(&S) -> f64 + 'static, action: impl Fn(&mut S) + 'static) -> Self {
        Self { rate: Box::new(rate), action: Box::new(action), enabled: Box::new(|_| true), name: String::new() }
    }
}
pub type EventFactory<S> = Box<dyn Fn(&S) -> Vec<Event<S>>>;
pub type Encoder<S> = Box<dyn Fn(&S) -> Vec<u64>>;
pub struct EventModel<S> {
    pub data: S,
    factory: EventFactory<S>,
    cache: RefCell<(Vec<Event<S>>, Vec<f64>)>,
    pub encode: Encoder<S>,
}
impl<S> EventModel<S> {
    pub fn new(data: S, factory: impl Fn(&S) -> Vec<Event<S>> + 'static) -> Self {
        Self { data, factory: Box::new(factory), cache: RefCell::new((vec![], vec![])), encode: Box::new(|_| vec![]) }
    }
}
impl<S> Model for EventModel<S> {
    fn propensities(&self) -> Result<Vec<f64>> {
        let events = (self.factory)(&self.data);
        let mut rates = Vec::new();
        for e in &events { rates.push(if (e.enabled)(&self.data) { nonnegative((e.rate)(&self.data))? } else { 0.0 }); }
        *self.cache.borrow_mut() = (events, rates.clone());
        Ok(rates)
    }
    fn fire(&mut self, event: usize) -> Result<()> {
        let (events, rates) = self.cache.get_mut();
        ensure(event < rates.len() && rates[event] > 0.0, "inactive or unevaluated event")?;
        (events[event].action)(&mut self.data);
        events.clear(); rates.clear();
        Ok(())
    }
    fn state(&self) -> Vec<u64> { (self.encode)(&self.data) }
}
