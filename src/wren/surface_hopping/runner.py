"""Model and molecular adiabatic ensembles, diagnostics and audited outputs."""
from dataclasses import asdict, replace
import json
from pathlib import Path
import time

from wren.quantum.cli import _write_result
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
    from .models import TullyModel, LinearVibronicModel, FragmentSitesModel
    quantum = config.model == 'pyscf-sa-casscf'
    if quantum:
        from .casscf_input import read_casscf_model
        from .pyscf_casscf import SACASSCF
        settings, _, _, sources = read_casscf_model(config.model_file)
        if any(provenance['sources'][p] != digest for p, digest in sources.items()):
            raise ValueError('SA-CASSCF model or XYZ changed while preparing dynamics')
        provider = SACASSCF(settings)
    elif config.model == 'lvc':
        data = json.loads(Path(config.model_file).read_text())
        if set(data) != {'schema', 'units', 'h0', 'linear', 'harmonic'} or data['schema'] != 'cube-lvc-1' or data['units'] != 'atomic':
            raise ValueError('LVC requires cube-lvc-1, atomic units, h0, linear and harmonic only')
        provider = LinearVibronicModel(data['h0'], data['linear'], data['harmonic'])
    elif config.model == 'fragment-sites':
        from .fragment_sites import read_fragment_model
        fragments, digest = read_fragment_model(config.model_file)
        if digest != provenance['sources'][config.model_file]:
            raise ValueError('fragment model changed while preparing dynamics')
        provider = FragmentSitesModel(fragments)
    else:
        provider = TullyModel(int(config.model[-1]))
    engine = FSSH(provider, config.masses, dt=config.dt, electronic_substeps=config.substeps,
                  decoherence=config.decoherence, edc_constant=config.edc_constant)
    result = {'schema': 'cube-surface-hopping-1', 'status': 'running', 'scope': 'model_hamiltonian',
              'fingerprint': fingerprint, 'provenance': provenance, 'config': asdict(config),
              'numpy_version': np.__version__, 'rng': 'PCG64 with SeedSequence.spawn per trajectory',
              'limitations': ['real fixed orthonormal diabatic basis', 'no atomistic QC, MD replay or GPU backend',
                              'frustrated hops keep momentum; no velocity reversal',
                              'completed-result reuse only; no interrupted-trajectory restart',
                              'FSSH is approximate; convergence and observable validation required']}
    initial_amplitudes = None
    if quantum:
        result.update(scope='on_the_fly_molecular_sa_casscf', basis='moving_adiabatic_many_electron_states',
                      atom_symbols=settings['symbols'], quantum_settings=settings,
                      limitations=['isolated molecule or finite cluster, equal-weight same-spin states, real wavefunctions, CPU',
                                   'active-space choice and initialization determine the accessible states',
                                   'energy-ordered roots with overlap checks; near-degeneracies or swaps stop the trajectory',
                                   'frustrated hops keep momentum; completed-result reuse only'])
        # The first geometry already requires orbital, CI and derivative solves.
        # Claim its output and preserve the cause if any of these fails.
        directory.mkdir(parents=True, exist_ok=True)
        with result_path.open('x') as initial_output:
            json.dump(result, initial_output, allow_nan=False)
    try:
        initial = engine.initialize(config.coordinates, config.momenta, config.active)
        if quantum and config.charge_localization:
            from .fragment_charges import localize_charge
            carrier, name = config.charge_localization
            initial_amplitudes, charge = localize_charge(initial.frame.diagnostics['fragment_charge_operators_e'], name, carrier)
            result.update(localization={'carrier': carrier, 'fragment': name, 'extremal_charge_e': charge,
                                        'amplitudes_adiabatic': initial_amplitudes.tolist()},
                          initial_active_probabilities=(abs(initial_amplitudes)**2).tolist())
    except Exception as exc:
        if quantum:
            result.update(status='failed', error={'type': type(exc).__name__, 'message': str(exc)})
            _write_result(result_path, result)
        raise
    n = len(initial.c)
    children = np.random.SeedSequence(config.seed).spawn(config.trajectories)
    directory.mkdir(parents=True, exist_ok=True)
    if config.model == 'fragment-sites':
        result['site_ids'] = list(fragments.names)
        result['basis'] = 'fixed_orthonormal_fragment_sites'
        if config.fragment_site is not None:
            initial_amplitudes = initial.frame.vectors[fragments.names.index(config.fragment_site), :].copy()
            result['initial_active_probabilities'] = (abs(initial_amplitudes)**2).tolist()
    _write_result(result_path, result)
    counts = np.zeros(n, dtype=int)
    populations, diabatic = np.zeros(n), np.zeros(n)
    events = {'accepted': 0, 'frustrated': 0}
    exits = {'left': 0, 'right': 0, 'max_steps': 0}
    channels = {side: [0]*n for side in exits}
    max_drift = max_norm = max_jump = max_probability = 0.
    if quantum:
        minimum_gap = initial.frame.diagnostics['minimum_gap_hartree']
        minimum_overlap = 1.
    finals = []
    started = time.perf_counter()
    def snapshot(index, step, state):
        data = {'trajectory': index, 'step': step, 'time_au': state.time, 'q_bohr': state.q.tolist(),
                'p_au': state.p.tolist(), 'active_state': state.active, 'energy_hartree': engine.energy(state),
                'adiabatic_populations': (abs(state.c)**2).tolist()}
        if quantum:
            data.update(energies_hartree=state.frame.energies.tolist(),
                        gradients_hartree_per_bohr=state.frame.gradients.tolist(),
                        derivative_couplings_per_bohr=state.frame.nac.tolist(), quantum=state.frame.diagnostics)
            if 'fragments' in settings:
                from .fragment_charges import fragment_observables
                obs = state.frame.diagnostics
                data.update(fragment_observables(state.c, obs['fragment_charge_operators_e'],
                                                 obs['fragment_nuclear_charges_e']))
                data.update(amplitudes_adiabatic_real=state.c.real.tolist(),
                            amplitudes_adiabatic_imag=state.c.imag.tolist())
        else:
            site = abs(state.frame.vectors @ state.c)**2
            data.update(diabatic_populations=site.tolist(), diabatic_ipr=float(1/np.sum(site**2)))
        return data
    try:
        with trace_path.open('x') as trace:
            for index, child in enumerate(children):
                rng = np.random.Generator(np.random.PCG64(child))
                active = config.active if initial_amplitudes is None else int(rng.choice(n, p=abs(initial_amplitudes)**2))
                # The initial molecular frame is identical for all trajectories.
                # step() produces fresh states and uses each trajectory's own preceding wavefunctions.
                if quantum:
                    state = initial if initial_amplitudes is None else replace(
                        initial, active=active, c=initial_amplitudes.astype(complex))
                else:
                    state = engine.initialize(config.coordinates, config.momenta, active, amplitudes=initial_amplitudes)
                energy0 = engine.energy(state)
                trace.write(json.dumps(snapshot(index, 0, state), allow_nan=False)+'\n')
                reason = 'max_steps'
                for step in range(1, config.steps+1):
                    state, event = engine.step(state, rng)
                    if quantum:
                        qc = state.frame.diagnostics
                        minimum_overlap = min(minimum_overlap, float(np.min(np.diag(qc['adjacent_wavefunction_overlaps']))))
                        minimum_gap = min(minimum_gap, qc['minimum_gap_hartree'])
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
                    if step % config.stride == 0 or stop or step == config.steps or (quantum and event['event'] != 'none'):
                        row = snapshot(index, step, state)
                        if quantum:
                            row['hop'] = event
                        trace.write(json.dumps(row, allow_nan=False)+'\n')
                    if stop:
                        break
                counts[state.active] += 1
                populations += abs(state.c)**2
                if not quantum:
                    diabatic += abs(state.frame.vectors @ state.c)**2
                exits[reason] += 1
                channels[reason][state.active] += 1
                finals.append({**snapshot(index, step, state), 'termination': reason})
                if (index+1) % max(1, config.trajectories//8) == 0:
                    print(f'wren: FSSH {index+1}/{config.trajectories} trajectories', flush=True)
        frequencies = counts/config.trajectories
        result.update(status='partial' if config.boundaries and exits['max_steps'] else 'completed', elapsed_seconds=time.perf_counter()-started,
                      active_counts=counts.tolist(), active_fractions=frequencies.tolist(),
                      active_fraction_standard_errors=np.sqrt(frequencies*(1-frequencies)/config.trajectories).tolist(),
                      mean_final_adiabatic_populations=(populations/config.trajectories).tolist(),
                      exits=exits, scattering_counts=channels, hops=events,
                      max_energy_drift_hartree=max_drift, max_norm_error=max_norm,
                      max_hop_energy_error_hartree=max_jump, max_probability_sum=max_probability,
                      final_states=finals, trajectories=str(trace_path), trajectory_sha256=sha256(trace_path))
        if quantum:
            result.update(minimum_adjacent_state_overlap=minimum_overlap, minimum_state_gap_hartree=minimum_gap)
            if 'fragments' in settings:
                result['mean_final_fragment_net_charges_e'] = {
                    name: float(np.mean([row['fragment_net_charges_e'][name] for row in finals]))
                    for name in settings['fragments']}
        else:
            result['mean_final_diabatic_populations'] = (diabatic/config.trajectories).tolist()
        _write_result(result_path, result)
    except Exception as exc:
        result.update(status='failed', error={'type': type(exc).__name__, 'message': str(exc)})
        _write_result(result_path, result)
        raise
    return {**result, 'reused': False, 'output': str(result_path)}
