"""Time the public run API, including snapshots but excluding parsing and CSV."""
import sys
import time
from universal_kmc import Engine
from universal_kmc.io import load

model = load(sys.argv[1])
engine = Engine(model, int(sys.argv[3]))
checksum = 0
started = time.perf_counter()
for row in engine.run(int(sys.argv[2])):
    checksum += row.event + 1
elapsed = time.perf_counter() - started
print(elapsed, engine.steps, engine.time, checksum, *model.position(), *model.state())
