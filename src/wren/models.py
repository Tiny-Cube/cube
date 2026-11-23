"""Graph transport and stochastic mass-action models."""
from dataclasses import dataclass
import math
from .core import MAX_COUNT, count, nonnegative


@dataclass(frozen=True)
class Edge:
    source: int
    target: int
    rate: float
    displacement: tuple[float, float, float] = (0.0, 0.0, 0.0)


class Graph:
    """One carrier; explicit edge displacements support unwrapped periodic motion."""
    def __init__(self, sites: int, edges: list[Edge], initial: int = 0):
        count(sites, "sites")
        count(initial, "initial")
        if sites == 0 or initial >= sites:
            raise ValueError("invalid graph size or initial site")
        self.sites = sites
        self.edges = tuple(edges)
        for edge in self.edges:
            if count(edge.source) >= sites or count(edge.target) >= sites:
                raise ValueError("edge endpoint out of bounds")
            nonnegative(edge.rate)
            if len(edge.displacement) != 3 or not all(math.isfinite(x) for x in edge.displacement):
                raise ValueError("displacement must have three finite components")
        self.site = initial
        self.xyz = [0.0, 0.0, 0.0]

    def propensities(self) -> list[float]:
        return [e.rate if e.source == self.site else 0.0 for e in self.edges]

    def fire(self, event: int) -> None:
        if not 0 <= event < len(self.edges):
            raise ValueError("event index out of bounds")
        edge = self.edges[event]
        if edge.source != self.site or edge.rate <= 0:
            raise ValueError("inactive event")
        new = [x + dx for x, dx in zip(self.xyz, edge.displacement)]
        if not all(math.isfinite(x) for x in new):
            raise ValueError("position overflow")
        self.site, self.xyz = edge.target, new

    def state(self) -> tuple[int]:
        return (self.site,)

    def position(self) -> tuple[float, ...]:
        return tuple(self.xyz)


@dataclass(frozen=True)
class Reaction:
    rate: float
    reactants: tuple[int, ...]
    products: tuple[int, ...]


class ReactionNetwork:
    """a_j = c_j * product_i binomial(n_i, reactants_ji).

    c_j is a stochastic rate constant; volume/factorial conversions are caller-owned.
    """
    def __init__(self, populations: list[int], reactions: list[Reaction]):
        if not populations:
            raise ValueError("at least one species is required")
        self.populations = [count(n) for n in populations]
        self.reactions = tuple(reactions)
        for reaction in self.reactions:
            nonnegative(reaction.rate)
            if len(reaction.reactants) != len(populations) or len(reaction.products) != len(populations):
                raise ValueError("stoichiometry dimension mismatch")
            for n in (*reaction.reactants, *reaction.products):
                count(n, "stoichiometry")

    def propensities(self) -> list[float]:
        rates = []
        for reaction in self.reactions:
            if reaction.rate == 0 or any(n < r for n, r in zip(self.populations, reaction.reactants)):
                rates.append(0.0)
                continue
            a = reaction.rate
            for n, r in zip(self.populations, reaction.reactants):
                for k in range(1, min(r, n - r) + 1):
                    a *= (n - k + 1) / k
                    nonnegative(a, "propensity")
            rates.append(a)
        return rates

    def fire(self, event: int) -> None:
        if not 0 <= event < len(self.reactions):
            raise ValueError("event index out of bounds")
        r = self.reactions[event]
        if r.rate <= 0 or any(n < a for n, a in zip(self.populations, r.reactants)):
            raise ValueError("inactive reaction")
        new = [n - a + b for n, a, b in zip(self.populations, r.reactants, r.products)]
        if any(n > MAX_COUNT for n in new):
            raise ValueError("population overflow")
        self.populations = new

    def state(self) -> tuple[int, ...]:
        return tuple(self.populations)

    def position(self) -> tuple[float, float, float]:
        return (0.0, 0.0, 0.0)
