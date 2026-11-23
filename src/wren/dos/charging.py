"""First vertical removal and attachment energies for isolated MD monomers."""
import csv

import numpy as np

from wren.workflow.scf import run_single_point
from .core import broaden, disorder_statistics


def calculate_frame(records, settings, directory):
    """Compute closed-shell singlet and two doublets without changing geometry."""
    from pyscf.data import elements
    from pyscf.data.nist import HARTREE2EV
    from pyscf.scf import chkfile
    from wren.quantum.geometry import normalize_atoms
    molecules = []
    for index, record in enumerate(records):
        atoms = normalize_atoms(zip(record['elements'], record['coordinates_angstrom'], strict=True))
        electrons = sum(elements.charge(symbol) for symbol, _ in atoms)
        if electrons < 2 or electrons % 2:
            raise ValueError(f"vertical-dscf requires an even-electron neutral singlet, molecule {record['id']}")
        energies, statuses = {}, {}
        for state, charge, spin in (('neutral', 0, 0), ('cation', 1, 1), ('anion', -1, 1)):
            prefix = directory/f'site{index:05d}'/state
            status = run_single_point(atoms, settings, prefix, charge=charge, spin=spin)
            occupation = np.asarray(chkfile.load(status['checkpoint'], 'scf/mo_occ'))
            if (not np.isclose(occupation.sum(), electrons - charge, rtol=0, atol=1e-8)
                    or (spin and (occupation.ndim != 2 or occupation.shape[0] != 2
                                  or abs(occupation[0].sum()-occupation[1].sum()-spin) > 1e-8))):
                raise ValueError(f'{state} SCF occupations disagree with requested electron count or spin; check basis')
            energies[state] = status['energy_hartree']
            statuses[state] = {
                'status_file': str(prefix.with_suffix('.status.json').relative_to(directory.parent)),
                **{key: status[key] for key in ('converged', 'backend', 'source_class', 'charge', 'spin_2s',
                                                'spin_square', 'checkpoint_sha256')}}
        ip = (energies['cation']-energies['neutral'])*HARTREE2EV
        ea = (energies['neutral']-energies['anion'])*HARTREE2EV
        molecules.append({'id': str(record['id']), 'energy_hartree': energies, 'scf': statuses,
                          'vertical_ip_ev': ip, 'vertical_ea_ev': ea,
                          'removal_ev': -ip, 'addition_ev': -ea})
    return molecules


def write_distributions(directory, frames, grid, settings):
    """Write each transition density on the same electron-energy grid."""
    channels = {}
    for channel in ('removal', 'addition'):
        levels = np.array([[m[channel+'_ev'] for m in frame['molecules']] for frame in frames])
        per_frame = [broaden(row, grid, settings['sigma'], weights=np.full(row.size, 1/row.size),
                             backend=settings['dos-backend']) for row in levels]
        density = np.mean(per_frame, axis=0)
        integral = float(np.trapezoid(density, grid))
        filename = channel+'.csv'
        with (directory/filename).open('w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['energy_ev', 'density_per_ev', *[f'frame_{i}' for i in range(len(frames))]])
            writer.writerows(zip(grid, density, *per_frame))
        channels[channel] = {'file': filename, 'electron_energy': '-IP' if channel == 'removal' else '-EA',
                             'expected_integral': 1., 'finite_grid_integral': integral,
                             'integral_absolute_error': abs(integral-1.),
                             'sample_statistics': disorder_statistics(levels)}
    with (directory/'vertical-energies.csv').open('w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['frame', 'time_ps', 'mol', 'neutral_hartree', 'cation_hartree', 'anion_hartree',
                         'vertical_ip_ev', 'vertical_ea_ev', 'removal_ev', 'addition_ev'])
        for index, frame in enumerate(frames):
            for mol in frame['molecules']:
                writer.writerow([frame['frame'] if frame['frame'] is not None else index,
                                 frame['time_ps'], mol['id'],
                                 *[mol['energy_hartree'][state] for state in ('neutral', 'cation', 'anion')],
                                 *[mol[key] for key in ('vertical_ip_ev', 'vertical_ea_ev', 'removal_ev', 'addition_ev')]])
    return {'scope': 'isolated_monomer_first_vertical_charge_transition_distributions',
            'normalization': 'one_transition_per_molecule_per_channel_equal_frame_weights',
            'energy_reference': 'electron_at_rest_in_vacuum', 'geometry': 'fixed_input_coordinates',
            'channels': channels, 'raw_energies': 'vertical-energies.csv',
            'nonpositive_vertical_ea_count': sum(m['vertical_ea_ev'] <= 0 for f in frames for m in f['molecules'])}
