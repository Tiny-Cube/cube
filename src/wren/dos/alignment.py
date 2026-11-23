"""Define a common energy axis for periodic spectra from multiple configurations."""
import csv
import math


def mean_local_potential_ev(cell):
    """Cell-mean Hartree plus local ionic GTH potential in PySCF's neutral 3D convention."""
    from pyscf.data.nist import HARTREE2EV
    from pyscf.pbc.gto.pseudo.pp import get_alphas
    if (cell.dimension != 3 or cell.charge != 0 or
            any(cell.atom_symbol(i) not in cell._pseudo for i in range(cell.natm))):
        raise ValueError('cell-mean-local needs a neutral 3D cell with a GTH pseudopotential on every atom')
    value = float(get_alphas(cell).sum()/cell.vol)*HARTREE2EV
    if not math.isfinite(value):
        raise ValueError('nonfinite periodic mean local potential')
    return value


def reference_values(settings, datasets):
    mode = settings.get('periodic-reference')
    if len(datasets) > 1:
        if mode is None:
            raise ValueError('multi-frame periodic DOS needs explicit potential alignment or periodic-reference=frame-vbm')
        first = datasets[0]['monomers']
        topology = [(m['id'], m['elements']) for m in first]
        if any([(m['id'], m['elements']) for m in d['monomers']] != topology for d in datasets[1:]):
            raise ValueError('periodic frames require the same ordered molecules and atomic species')
    if mode != 'external':
        return {}
    expected = {d.get('frame'): d['_source_sha256'] for d in datasets}
    if any(type(frame) is not int or frame < 0 for frame in expected):
        raise ValueError('external alignment requires a nonnegative integer frame ID in each manifest')
    refs = {}
    with open(settings['reference-file'], newline='', encoding='utf-8-sig') as stream:
        reader = csv.DictReader(stream)
        names = reader.fieldnames or []
        fields = set(names)
        if (len(fields) != len(names) or not {'frame', 'reference_ev'} <= fields
                or fields - {'frame', 'reference_ev', 'manifest_sha256'}):
            raise ValueError('reference CSV requires frame,reference_ev and optionally manifest_sha256')
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError('reference CSV row has an incorrect number of fields')
            frame = int(row['frame'])
            value = float(row['reference_ev'])
            if frame in refs or frame not in expected or not math.isfinite(value):
                raise ValueError('reference CSV must have one finite energy per selected frame')
            if 'manifest_sha256' in row and row['manifest_sha256'].strip().lower() != expected[frame]:
                raise ValueError(f'reference CSV manifest hash mismatch for frame {frame}')
            refs[frame] = value
    if refs.keys() != expected.keys():
        raise ValueError('reference CSV does not cover every selected frame')
    return refs


def aligned_levels(settings, frame, levels, refs):
    mode = settings.get('periodic-reference')
    if mode is None:
        return levels
    edges = frame['band_edges']
    if mode == 'frame-vbm':
        if edges['vbm_ev'] is None or edges['partial_occupation'] or not edges['sampled_gap_ev']:
            raise ValueError('frame-vbm requires occupied and empty bands separated by a sampled gap')
        reference = edges['vbm_ev']
    elif mode == 'cell-mean-local':
        reference = frame['cell_mean_local_ev']
    else:
        reference = refs[frame['frame']]
    frame['energy_alignment'] = {'mode': mode, 'reference_ev': reference,
                                 'aligned_vbm_ev': None if edges['vbm_ev'] is None else edges['vbm_ev']-reference,
                                 'aligned_cbm_ev': None if edges['cbm_ev'] is None else edges['cbm_ev']-reference}
    return levels-reference
