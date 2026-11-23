"""Exact direct SSA for state-dependent, time-homogeneous Markov jump processes."""
from dataclasses import dataclass
import math
from typing import Iterator, Protocol, Sequence

MAX_COUNT = 2**53 - 1  # common exact integer range for all four implementations


def count(value: int, name: str = "count") -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not 0 <= value <= MAX_COUNT:
        raise ValueError(f"{name} must be an integer in [0, {MAX_COUNT}]")
    return value


def nonnegative(value: float, name: str = "rate") -> float:
    if not math.isfinite(value) or value < 0:
        raise ValueError(f"{name} must be finite and nonnegative")
    return value


class RandomSource(Protocol):
    def uniform(self) -> float: ...


class PortableRNG:
    """MINSTD (48271); portable reference stream, not a parallel/HPC RNG."""
    def __init__(self, seed: int = 1):
        count(seed, "seed")
        self.state = seed % 2147483646 + 1

    def uniform(self) -> float:
        self.state = (48271 * self.state) % 2147483647
        return self.state / 2147483647.0


class Model(Protocol):
    """Event order must stay fixed between propensities() and fire()."""
    def propensities(self) -> Sequence[float]: ...
    def fire(self, event: int) -> None: ...


@dataclass(frozen=True)
class Record:
    step: int
    time: float
    event: int
    position: tuple[float, ...]
    state: tuple[int, ...]


class Engine:
    def __init__(self, model: Model, seed: int = 1, rng: RandomSource | None = None):
        self.model = model
        self.rng = rng if rng is not None else PortableRNG(seed)
        self.time = 0.0
        self.steps = 0
        self.status = "ready"

    def snapshot(self, event: int = -1) -> Record:
        position_fn = getattr(self.model, "position", None)
        state_fn = getattr(self.model, "state", None)
        position = position_fn() if callable(position_fn) else (0.0, 0.0, 0.0)
        state = state_fn() if callable(state_fn) else ()
        return Record(self.steps, self.time, event, tuple(position), tuple(state))

    def _uniform(self) -> float:
        u = self.rng.uniform()
        if not math.isfinite(u) or not 0 < u < 1:
            raise ValueError("RNG must return a finite value strictly inside (0, 1)")
        return u

    def step(self, until: float = math.inf) -> Record | None:
        if math.isnan(until) or until < self.time:
            raise ValueError("until must be >= current time")
        if self.time == until:
            self.status = "time_limit"
            return None
        rates = list(self.model.propensities())
        total = 0.0
        for rate in rates:
            total += nonnegative(rate)
        if not math.isfinite(total):
            raise ValueError("total propensity overflow")
        if total == 0:
            self.status = "absorbed"
            if math.isfinite(until):
                self.time = until
            return None
        dt = -math.log(self._uniform()) / total
        if dt > until - self.time:
            self.time = until
            self.status = "time_limit"
            return None
        if not math.isfinite(dt) or not math.isfinite(self.time + dt) or self.time + dt <= self.time:
            raise ValueError("waiting time is outside floating-point range")
        threshold = self._uniform() * total
        cumulative = 0.0
        chosen = -1
        for i, rate in enumerate(rates):
            if rate > 0:
                chosen = i  # safe fallback for rounding at the last interval
            cumulative += rate
            if threshold < cumulative:
                break
        self.model.fire(chosen)
        self.time += dt
        self.steps += 1
        self.status = "running"
        return self.snapshot(chosen)

    def run(self, max_steps: int = 1000, until: float = math.inf) -> Iterator[Record]:
        count(max_steps, "max_steps")
        if math.isnan(until) or until < self.time:
            raise ValueError("until must be >= current time")
        yield self.snapshot()
        for _ in range(max_steps):
            before = self.time
            record = self.step(until)
            if record is None:
                if self.time != before:
                    yield self.snapshot()  # terminal observation, not an event
                return
            yield record
        self.status = "time_limit" if self.time == until else "step_limit"
