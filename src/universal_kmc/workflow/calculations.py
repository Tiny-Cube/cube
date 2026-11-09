"""Calculation entry registry and dependency-free routing of shared modules.

Register new calculation types here (or through register_calculation). A solver
is a callable(control, preparation_report, *, resume) returning a result dict.
Import expensive solver dependencies inside that callable, never in the parser.
"""
from dataclasses import dataclass
import re
from typing import Callable

EXTRACTION = 'get_dimers_monomers'
QUANTUM_SECTIONS = ('two-mol', 'one-mol', 'reorganization-energy')
MODULE_SECTIONS = (EXTRACTION, *QUANTUM_SECTIONS)
MODULE_ORDER = (EXTRACTION, 'one-mol', 'two-mol', 'reorganization-energy')


@dataclass(frozen=True)
class CalculationType:
    modules: tuple[str, ...]
    solver: Callable | None = None
    unavailable_reason: str | None = None


CALCULATIONS: dict[str, CalculationType] = {}


def register_calculation(name, calculation):
    if not re.fullmatch(r'[a-z][a-z0-9-]*', name) or name in (*MODULE_SECTIONS, 'end', 'get_dimer'):
        raise ValueError('calculation name must be a distinct lowercase section name')
    if name in CALCULATIONS:
        raise ValueError(f'calculation already registered: {name}')
    if not isinstance(calculation, CalculationType):
        raise TypeError('expected a CalculationType')
    if len(set(calculation.modules)) != len(calculation.modules) or not set(calculation.modules) <= set(MODULE_SECTIONS):
        raise ValueError('calculation modules must be unique registered modules')
    if calculation.solver is not None and not callable(calculation.solver):
        raise TypeError('solver must be callable')
    CALCULATIONS[name] = calculation


register_calculation('charge-carrier', CalculationType(MODULE_ORDER))
register_calculation('exciton', CalculationType((EXTRACTION,), unavailable_reason=
    'excited-state parameters, exciton couplings and diffusion solver are not yet connected'))
register_calculation('surface-hopping', CalculationType((EXTRACTION,), unavailable_reason=
    'electronic-state energies, nonadiabatic couplings and surface-hopping dynamics are not yet connected'))


def section_names():
    return (*CALCULATIONS, *MODULE_SECTIONS)


def route_sections(available, entries, only=None, *, module_switches=None):
    """Compose main workflows and independently requested shared modules.

    Submodule do=1 requests execution without a main calculation; do=0 disables
    that submodule. Omitted do retains automatic main-workflow/legacy dispatch.
    --only narrows the selected union and never overrides an explicit do=0.
    """
    switches = module_switches or {}
    requested = {EXTRACTION if s == 'get_dimer' else s for s in (only or ())}
    if not requested <= set(available):
        raise ValueError('--only must name a section present in input')
    entry_filter = requested.intersection(CALCULATIONS)
    module_filter = requested.intersection(MODULE_SECTIONS)
    selected = [name for name in CALCULATIONS if name in entries and (not entry_filter or name in entry_filter)]
    enabled = [name for name in selected if entries[name]['do'] == 1]
    disabled = [name for name in selected if entries[name]['do'] == 0]
    if entries:
        owners = {module: [name for name in enabled if module in CALCULATIONS[name].modules
                           and entries[name].get(module, 1) == 1]
                  for module in MODULE_ORDER}
    else:
        owners = {module: ['legacy'] for module in MODULE_ORDER}
    for module in MODULE_ORDER:
        if switches.get(module) == 0:
            owners[module] = []
        elif switches.get(module) == 1:
            # Module requests are independent of main selection and ownership.
            owners[module] = [name for name in owners[module] if name != 'legacy']+['standalone']
    tasks = [module for module in MODULE_ORDER if module in available and owners[module]
             and (not module_filter or module in module_filter)]
    calculations = enabled if not module_filter else []
    return {'routing': 'calculation_entries' if entries else 'legacy_modules',
            'requested_modules': [module for module in tasks if switches.get(module) == 1],
            'module_switches': {module: switches.get(module) for module in MODULE_ORDER if module in available},
            'enabled_calculations': enabled, 'calculation_tasks': calculations,
            'disabled_tasks': disabled, 'tasks': tasks,
            'module_owners': {module: owners[module] for module in tasks},
            'skipped_modules': {module: ('disabled_by_module_do' if switches.get(module) == 0 else
                                        'no_enabled_calculation_requires_module' if not owners[module]
                                        else 'not_selected_by_only')
                                for module in MODULE_ORDER if module in available and module not in tasks},
            'unsupported_tasks': [name for name in calculations if CALCULATIONS[name].unavailable_reason],
            'pending_solvers': [name for name in calculations if CALCULATIONS[name].solver is None],
            'entry_points': {name: {'do': entries[name]['do'], 'selected': name in selected,
                                   'modules': list(CALCULATIONS[name].modules),
                                   'module_selection': {module: entries[name].get(module, 1)
                                                        for module in CALCULATIONS[name].modules},
                                   'solver_available': CALCULATIONS[name].solver is not None,
                                   'unavailable_reason': CALCULATIONS[name].unavailable_reason}
                             for name in CALCULATIONS if name in entries}}
