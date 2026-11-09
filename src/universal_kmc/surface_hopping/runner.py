"""Model ensembles and audited outputs; no atomistic/GPU claims."""
from dataclasses import asdict
import json
from pathlib import Path
import time

from universal_kmc.quantum.cli import _write_result
from .control import load_control, sha256


def run_control(filename, *, output=None, resume=False):
    config, directory, provenance, fingerprint = load_control(filename, output)
    result_path = directory/'result.json'
    trace_path = directory/'trajectories.jsonl'
    if result_path.exists():
        previous = json.loads(result_path.read_text())
        if not resume or previous.get('fingerprint') != fingerprint or previous.get('status') != 'completed':
            raise ValueError('existing surface-hopping result; use --resume only with identical completed inputs')
        if not trace_path.is_file() or sha256(trace_path) != previous['trajectory_sha256']:
            raise ValueError('cannot resume: trajectory output missing or changed')
        return {**previous, 'reused': True, 'output': str(result_path)}
    if trace_path.exists():
        raise ValueError('partial trajectory output exists; select a new output directory')
    import numpy as np
    from .engine import FSSH
    from .models import TullyModel, LinearVibronicModel
    if config.model == 'lvc':
        data = json.loads(Path(config.model_file).read_text())
        if set(data) != {'schema', 'units', 'h0', 'linear', 'harmonic'} or data['schema'] != 'cube-lvc-1' or data['units'] != 'atomic':
            raise ValueError('LVC requires cube-lvc-1, atomic units, h0, linear and harmonic only')
        provider = LinearVibronicModel(data['h0'], data['linear'], data['harmonic'])
    else:
        provider = TullyModel(int(config.model[-1]))
    engine = FSSH(provider, config.masses, dt=config.dt, electronic_substeps=config.substeps,
                  decoherence=config.decoherence, edc_constant=config.edc_constant)
    # Validate provider and initial state before creating scientific outputs.
    initial = engine.initialize(config.coordinates, config.momenta, config.active)
    n = len(initial.c)
    children = np.random.SeedSequence(config.seed).spawn(config.trajectories)
    directory.mkdir(parents=True, exist_ok=True)
    result = {'schema': 'cube-surface-hopping-1', 'status': 'running', 'scope': 'model_hamiltonian',
              'fingerprint': fingerprint, 'provenance': provenance, 'config': asdict(config),
              'numpy_version': np.__version__, 'rng': 'PCG64 with SeedSequence.spawn per trajectory',
              'limitations': ['real fixed orthonormal diabatic basis', 'no atomistic QC, MD replay or GPU backend',
                              'frustrated hops keep momentum; no velocity reversal',
                              'completed-result reuse only; no interrupted-trajectory restart',
                              'FSSH is approximate; convergence and observable validation required']}
    _write_result(result_path, result)
    counts = np.zeros(n, dtype=int)
    populations, diabatic = np.zeros(n), np.zeros(n)
    events = {'accepted': 0, 'frustrated': 0}
    exits = {'left': 0, 'right': 0, 'max_steps': 0}
    channels = {side: [0]*n for side in exits}
    max_drift = max_norm = max_jump = max_probability = 0.
    finals = []
    started = time.perf_counter()
    def snapshot(index, step, state):
        site = abs(state.frame.vectors @ state.c)**2
        return {'trajectory': index, 'step': step, 'time_au': state.time, 'q_bohr': state.q.tolist(),
                'p_au': state.p.tolist(), 'active_state': state.active, 'energy_hartree': engine.energy(state),
                'adiabatic_populations': (abs(state.c)**2).tolist(), 'diabatic_populations': site.tolist(),
                'diabatic_ipr': float(1/np.sum(site**2))}
    try:
        with trace_path.open('x') as trace:
            for index, child in enumerate(children):
                rng = np.random.Generator(np.random.PCG64(child))
                state = engine.initialize(config.coordinates, config.momenta, config.active)
                energy0 = engine.energy(state)
                trace.write(json.dumps(snapshot(index, 0, state), allow_nan=False)+'\n')
                reason = 'max_steps'
                for step in range(1, config.steps+1):
                    state, event = engine.step(state, rng)
                    if event['event'] in events:
                        events[event['event']] += 1
                    max_drift = max(max_drift, abs(engine.energy(state)-energy0))
                    max_norm = max(max_norm, abs(float(np.vdot(state.c, state.c).real)-1))
                    max_jump = max(max_jump, abs(event['jump_energy_error']))
                    max_probability = max(max_probability, event['probability_sum'])
                    stop = False
                    if config.boundaries:
                        if state.q[0] <= config.boundaries[0]:
                            reason, stop = 'left', True
                        elif state.q[0] >= config.boundaries[1]:
                            reason, stop = 'right', True
                    if step % config.stride == 0 or stop or step == config.steps:
                        trace.write(json.dumps(snapshot(index, step, state), allow_nan=False)+'\n')
                    if stop:
                        break
                counts[state.active] += 1
                populations += abs(state.c)**2
                diabatic += abs(state.frame.vectors @ state.c)**2
                exits[reason] += 1
                channels[reason][state.active] += 1
                finals.append({**snapshot(index, step, state), 'termination': reason})
                if (index+1) % max(1, config.trajectories//8) == 0:
                    print(f'cube: FSSH {index+1}/{config.trajectories} trajectories', flush=True)
        frequencies = counts/config.trajectories
        result.update(status='partial' if config.boundaries and exits['max_steps'] else 'completed', elapsed_seconds=time.perf_counter()-started,
                      active_counts=counts.tolist(), active_fractions=frequencies.tolist(),
                      active_fraction_standard_errors=np.sqrt(frequencies*(1-frequencies)/config.trajectories).tolist(),
                      mean_final_adiabatic_populations=(populations/config.trajectories).tolist(),
                      mean_final_diabatic_populations=(diabatic/config.trajectories).tolist(),
                      exits=exits, scattering_counts=channels, hops=events,
                      max_energy_drift_hartree=max_drift, max_norm_error=max_norm,
                      max_hop_energy_error_hartree=max_jump, max_probability_sum=max_probability,
                      final_states=finals, trajectories=str(trace_path), trajectory_sha256=sha256(trace_path))
        _write_result(result_path, result)
    except Exception as exc:
        result.update(status='failed', error={'type': type(exc).__name__, 'message': str(exc)})
        _write_result(result_path, result)
        raise
    return {**result, 'reused': False, 'output': str(result_path)}
