"""Export completed four-point results without running quantum calculations."""
import csv
import math
import os
from pathlib import Path
import tempfile


ENERGY_COLUMNS = ('E0_R0', 'Eq_R0', 'Eq_Rq', 'E0_Rq')
RELAXATION_COLUMNS = ('charged_relaxation', 'neutral_relaxation', 'total')
COLUMNS = ('mol', *ENERGY_COLUMNS, *RELAXATION_COLUMNS)
UNITS = {**dict.fromkeys(ENERGY_COLUMNS, 'Hartree'),
         **dict.fromkeys(RELAXATION_COLUMNS, 'eV')}


def write_reorganization_csv(results, directory, carriers):
    """Write one table per requested carrier, preserving molecule order.

    ``results`` maps molecule IDs to completed reorganization.json records.
    Tables describe only this selection; missing/failed results are never
    represented as zero energies. Empty results produce header-only tables.
    """
    carriers = tuple(carriers)
    if not carriers or len(set(carriers)) != len(carriers) or set(carriers)-{'hole', 'electron'}:
        raise ValueError('CSV carriers must be unique hole/electron channels')
    tables = {carrier: [] for carrier in carriers}
    for mol, result in results.items():
        if result.get('schema') != 'cube-reorganization-1' or result.get('status') != 'completed':
            raise ValueError(f'molecule {mol}: CSV requires a completed four-point result')
        for carrier in carriers:
            try:
                block = result['results'][carrier]
                row = {'mol': str(mol),
                       **{key: block['energies_hartree'][key] for key in ENERGY_COLUMNS},
                       **{key: block['relaxation_ev'][key] for key in RELAXATION_COLUMNS}}
            except KeyError as exc:
                raise ValueError(f'molecule {mol}: missing {carrier} four-point CSV data') from exc
            if any(isinstance(row[key], bool) or not isinstance(row[key], (int, float))
                   or not math.isfinite(row[key]) for key in COLUMNS[1:]):
                raise ValueError(f'molecule {mol}: CSV energies must be finite numbers')
            tables[carrier].append(row)
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    paths = {}
    for carrier, rows in tables.items():
        path = directory/f'{carrier}.csv'
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', newline='',
                                             dir=directory, prefix=path.name+'.', delete=False) as file:
                temporary = Path(file.name)
                writer = csv.DictWriter(file, fieldnames=COLUMNS)
                writer.writeheader()
                writer.writerows(rows)
            os.replace(temporary, path)
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)
        paths[carrier] = str(path)
    return paths
