import argparse
import csv
import math
import sys
from .core import Engine
from .io import load


def main() -> int:
    parser = argparse.ArgumentParser(description="Continuous-time KMC; writes CSV to stdout")
    parser.add_argument('model')
    parser.add_argument('steps', nargs='?', type=int, default=1000)
    parser.add_argument('seed', nargs='?', type=int, default=1)
    parser.add_argument('until', nargs='?', type=float, default=math.inf)
    args = parser.parse_args()
    try:
        if len(sys.argv) > 4 and (not math.isfinite(args.until) or args.until < 0):
            raise ValueError("explicit time limit must be finite and nonnegative")
        model = load(args.model)
        engine = Engine(model, args.seed)
        writer = csv.writer(sys.stdout, lineterminator='\n')
        writer.writerow(['step', 'time', 'event', 'x', 'y', 'z', *[f'state{i}' for i in range(len(model.state()))]])
        for row in engine.run(args.steps, args.until):
            writer.writerow([row.step, format(row.time, '.17g'), row.event, *[format(x, '.17g') for x in row.position], *row.state])
        print(f"status={engine.status} events={engine.steps} time={engine.time:.17g}", file=sys.stderr)
    except (ValueError, OverflowError, OSError) as exc:
        print(f"ukmc: {exc}", file=sys.stderr)
        return 2
    return 0
