"""Domain-free callback adapter. State can be any Python object."""
from dataclasses import dataclass
from typing import Callable, Generic, Sequence, TypeVar
from .core import nonnegative

S = TypeVar('S')


@dataclass(frozen=True)
class Event(Generic[S]):
    rate: Callable[[S], float]
    action: Callable[[S], None]
    enabled: Callable[[S], bool] = lambda _: True
    name: str = ''


class EventModel(Generic[S]):
    """Fixed or dynamically generated events acting on arbitrary mutable data.

    A factory is called once per propensity evaluation. Actions must mutate data
    in place. Do not mutate data externally between evaluating and firing events.
    The optional encode function is only for integer CSV diagnostics; observers
    may access arbitrary fields in data directly, including real values/objects.
    """
    def __init__(self, data: S,
                 events: Sequence[Event[S]] | Callable[[S], Sequence[Event[S]]],
                 encode: Callable[[S], Sequence[int]] = lambda _: ()):
        self.data = data
        self._factory = events if callable(events) else lambda _: events
        self._encode = encode
        self.events: tuple[Event[S], ...] = ()
        self._rates: tuple[float, ...] = ()

    def propensities(self) -> tuple[float, ...]:
        events = tuple(self._factory(self.data))
        rates = tuple(nonnegative(e.rate(self.data)) if e.enabled(self.data) else 0.0 for e in events)
        self.events, self._rates = events, rates
        return rates

    def fire(self, event: int) -> None:
        if not 0 <= event < len(self.events) or self._rates[event] <= 0:
            raise ValueError('inactive or unevaluated event')
        self.events[event].action(self.data)
        self._rates = ()
        self.events = ()

    def state(self) -> tuple[int, ...]:
        return tuple(self._encode(self.data))
