"""Queueing and Ising spin dynamics: neither needs the graph/reaction adapters."""
from dataclasses import dataclass, field
import math
from universal_kmc import Engine, Event, EventModel


@dataclass
class Queue:
    jobs: list[str] = field(default_factory=list)
    completed: int = 0
    next_id: int = 0


def queue_model() -> EventModel[Queue]:
    def arrive(q: Queue) -> None:
        q.jobs.append(f'job-{q.next_id}')
        q.next_id += 1

    def serve(q: Queue) -> None:
        q.jobs.pop(0)
        q.completed += 1

    return EventModel(Queue(), [
        Event(lambda q: 2.0, arrive, name='arrival'),
        Event(lambda q: 3.0, serve, enabled=lambda q: bool(q.jobs), name='service'),
    ], encode=lambda q: (len(q.jobs), q.completed))


@dataclass
class Spins:
    spins: list[int]
    beta: float = 0.6
    coupling: float = 1.0


def ising_model(size: int = 16) -> EventModel[Spins]:
    if size < 3:
        raise ValueError('ring must contain at least three spins')

    def events(s: Spins) -> list[Event[Spins]]:
        result = []
        for i in range(len(s.spins)):
            def rate(state: Spins, i: int = i) -> float:
                values = state.spins
                delta = 2 * state.coupling * values[i] * (values[i - 1] + values[(i + 1) % len(values)])
                return math.exp(-max(state.beta * delta, 0.0))

            def flip(state: Spins, i: int = i) -> None:
                state.spins[i] *= -1

            result.append(Event(rate, flip, name=f'flip-{i}'))
        return result

    return EventModel(Spins([1] * size), events, encode=lambda s: s.spins)


def main() -> None:
    for model in (queue_model(), ising_model()):
        engine = Engine(model, seed=42)
        for _ in engine.run(10000, until=100.0):
            pass
        print(type(model.data).__name__, engine.status, engine.time, model.data)


if __name__ == '__main__':
    main()
