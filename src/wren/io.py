"""Strict, whitespace-delimited KMC1 exchange format (see docs/format.md)."""
from pathlib import Path
import re
from .core import Model, count
from .models import Edge, Graph, Reaction, ReactionNetwork


def loads(text: str) -> Model:
    words = iter(' '.join(line.split('#', 1)[0] for line in text.splitlines()).split())

    def token() -> str:
        try:
            return next(words)
        except StopIteration as exc:
            raise ValueError("truncated KMC1 input") from exc

    def integer() -> int:
        word = token()
        if not re.fullmatch(r'[0-9]+', word):
            raise ValueError("expected unsigned decimal integer")
        return count(int(word))

    def number() -> float:
        word = token()
        if not re.fullmatch(r'[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?', word):
            raise ValueError("expected decimal floating-point number")
        return float(word)

    if token() != 'KMC1':
        raise ValueError("expected KMC1 header")
    kind = token()
    n, m = integer(), integer()
    if kind == 'graph':
        initial = integer()
        model = Graph(n, [Edge(integer(), integer(), number(), (number(), number(), number())) for _ in range(m)], initial)
    elif kind == 'reaction':
        populations = [integer() for _ in range(n)]
        model = ReactionNetwork(populations, [Reaction(number(), tuple(integer() for _ in range(n)), tuple(integer() for _ in range(n))) for _ in range(m)])
    else:
        raise ValueError("model kind must be graph or reaction")
    if next(words, None) is not None:
        raise ValueError("unexpected trailing tokens")
    return model


def load(path: str | Path) -> Model:
    return loads(Path(path).read_text())
