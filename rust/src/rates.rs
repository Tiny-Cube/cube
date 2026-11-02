//! Energies: eV; temperature: K; returned rates: s^-1.
use crate::{ensure, nonnegative, Result};
use std::f64::consts::PI;
pub const KB: f64 = 8.617333262145e-5;
pub const HBAR: f64 = 6.582119569e-16;
fn thermal(t: f64) -> Result<f64> { ensure(t.is_finite() && t > 0.0, "positive finite temperature required")?; Ok(KB * t) }
pub fn marcus(coupling: f64, lambda: f64, dg: f64, t: f64) -> Result<f64> {
    let kt = thermal(t)?;
    ensure(coupling.is_finite() && lambda.is_finite() && lambda > 0.0 && dg.is_finite(), "invalid Marcus energies")?;
    if coupling == 0.0 { return Ok(0.0); }
    let x = (dg + lambda) / (lambda * kt).sqrt();
    nonnegative(((2.0*PI/HBAR).ln() + 2.0*coupling.abs().ln() - 0.5*(4.0*PI*lambda*kt).ln() - x*x/4.0).exp())
}
pub fn arrhenius(attempt: f64, barrier: f64, t: f64) -> Result<f64> {
    Ok(nonnegative(attempt)? * (-nonnegative(barrier)? / thermal(t)?).exp())
}
pub fn miller_abrahams(attempt: f64, distance: f64, localization: f64, dg: f64, t: f64) -> Result<f64> {
    let kt = thermal(t)?;
    nonnegative(attempt)?; nonnegative(distance)?;
    ensure(localization.is_finite() && localization > 0.0 && dg.is_finite(), "invalid hopping parameters")?;
    Ok(attempt * (-2.0*distance/localization - dg.max(0.0)/kt).exp())
}
