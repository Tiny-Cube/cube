"""Energies in eV, temperature in K, all returned rates in s^-1."""
import math
from .core import nonnegative

KB = 8.617333262145e-5
HBAR = 6.582119569e-16


def thermal(temperature: float) -> float:
    if not math.isfinite(temperature) or temperature <= 0:
        raise ValueError("temperature must be finite and positive")
    return KB * temperature


def marcus(coupling: float, reorganization: float, delta_g: float, temperature: float = 300.0) -> float:
    kt = thermal(temperature)
    if not all(math.isfinite(v) for v in (coupling, reorganization, delta_g)) or reorganization <= 0:
        raise ValueError("finite energies and positive reorganization required")
    if coupling == 0:
        return 0.0
    # Log form avoids squaring the electronic coupling before exponential damping.
    exponent = ((delta_g + reorganization) / math.sqrt(reorganization * kt)) ** 2 / 4
    log_rate = math.log(2 * math.pi / HBAR) + 2 * math.log(abs(coupling)) - 0.5 * math.log(4 * math.pi * reorganization * kt) - exponent
    return nonnegative(math.exp(log_rate))


def arrhenius(attempt: float, barrier: float, temperature: float = 300.0) -> float:
    return nonnegative(attempt) * math.exp(-nonnegative(barrier, "barrier") / thermal(temperature))


def miller_abrahams(attempt: float, distance: float, localization: float, delta_g: float, temperature: float = 300.0) -> float:
    kt = thermal(temperature)
    nonnegative(attempt)
    nonnegative(distance, "distance")
    if not math.isfinite(localization) or localization <= 0 or not math.isfinite(delta_g):
        raise ValueError("positive localization length and finite delta_g required")
    return attempt * math.exp(-2 * distance / localization - max(delta_g, 0.0) / kt)
