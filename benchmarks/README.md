# Four-language performance benchmark

Run from the repository root after `uv sync --locked`:

```bash
uv run python scripts/benchmark.py --output runs/performance_new
```

CMake 3.25+ configures a dedicated Release tree in `build/benchmark`. The
`wren_benchmarks` target builds the current C++, Fortran and Cargo kernel and
CLI programs. See its CMakeCache.txt for compiler paths and flags. The Rust
examples use the Cargo release profile. Fortran checks internal to the model
remain enabled in Release builds.
Python runs directly with the uv-created `.venv/bin/python` (no uv startup cost).
Use `--python /path/to/python` to compare another interpreter.

Models and work per sample:

- Two-state graph: 200,000 events.
- Birth/death reaction network: 200,000 events.
- Actual molecular hole graph: 20,000 events, 1,100 directed channels.
- Actual molecular electron graph: 20,000 events, 1,100 directed channels.

Every run uses seed 42 and the same existing direct-SSA implementation and model
input. The graph adapter scans all channels, including currently disabled ones;
this is not the sparse outgoing-edge adapter used for the molecular ensemble.

Two timings are reported separately:

1. `kernel`: monotonic wall time inside each benchmark driver, covering the public
   `Engine.run` API, including rate construction, validation, RNG, event execution,
   snapshots and a tiny event-index checksum observer. Model parsing and printing
   the final result are outside the timed region. Each language retains its
   existing allocation and dispatch behavior; these are implementation benchmarks,
   not universal comparisons of language performance.
2. `cli`: parent-process wall time around the existing CLI, including interpreter
   or runtime startup, parsing, computation and full per-event CSV formatting.
   Output is discarded via `/dev/null`; storage-device write speed is excluded.

On Linux, the runner pins itself and child processes to the first allowed logical
CPU. Every workload/mode/language receives one complete unrecorded warmup and
five recorded runs. Language order is randomly interleaved with a fixed ordering
seed. No benchmark subprocesses run concurrently. Compilation is outside the
measurements. CPU frequency, turbo boost and other system activity are not locked.

For every kernel run, the event count, final simulated time, event-index checksum,
position and model state must match the reference across all languages. CLI runs
must match the event count and simulated time. Existing trajectory conformance
tests provide more detailed event-by-event validation. Result JSON contains each
sample, medians, ranges, event throughput, versions, CPU, source hashes and final
signatures. An existing completed output directory is never overwritten.
