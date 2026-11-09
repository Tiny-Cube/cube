"""Calculation gates, shared module ownership and extension dispatch contracts."""
from contextlib import redirect_stdout
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from universal_kmc.workflow.calculations import (CALCULATIONS, CalculationType, register_calculation)
from universal_kmc.workflow.control import load_input
from universal_kmc.workflow.runner import run_workflow

ENTRIES = '$charge-carrier\ndo = {charge}\n$end\n$exciton\ndo = {exciton}\n$end\n$surface-hopping\ndo = {surface}\n$end\n'
# A disabled module must not require its files or validate its method settings.
UNAVAILABLE_MODULES = '''$get_dimers_monomers
tpr = absent.tpr
gro = absent.gro
xtc = absent.xtc
fragment-file = absent.txt
$end
$one-mol
xyz = absent.xyz
qc-engine = not-an-engine
$end
$two-mol
basis-set = intentionally-unused
$end
$reorganization-energy
environment = intentionally-unused
$end
'''


class CalculationEntryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.path = self.root/'input'

    def write(self, charge=0, exciton=0, surface=0, modules=UNAVAILABLE_MODULES):
        self.path.write_text(ENTRIES.format(charge=charge, exciton=exciton, surface=surface)+modules)
        return load_input(self.path)

    def test_all_disabled_ignores_missing_module_files_even_with_only(self):
        control = self.write()
        self.assertEqual(control.plan()['tasks'], [])
        self.assertEqual(set(control.inactive_sections), {'get_dimers_monomers', 'one-mol', 'two-mol', 'reorganization-energy'})
        for only in (None, ['one-mol'], ['get_dimer'], ['charge-carrier', 'one-mol']):
            with self.subTest(only=only), redirect_stdout(io.StringIO()):
                self.assertEqual(run_workflow(control, only=only), 0)
            report = json.loads((self.root/'cube-results/workflow.json').read_text())
            self.assertEqual(report['tasks'], {})
        self.assertFalse((self.root/'mol').exists())

    def test_disabled_module_syntax_still_rejects_typos(self):
        for module in ('$one-mol\nmisspelled = 1\n$end\n', '$one-mol\nxyz = a\nxyz = b\n$end\n'):
            with self.assertRaisesRegex(ValueError, 'unknown or duplicate key'):
                self.write(modules=module)

    def test_exciton_and_surface_do_not_run_ground_state_charge_modules(self):
        # Existing extraction inputs suffice for shared geometry; invalid charge
        # settings remain inactive rather than silently becoming excited states.
        for name in ('a.tpr','a.gro','a.xtc','core.txt'):(self.root/name).write_text('placeholder')
        geometry = '$get_dimers_monomers\ntpr=a.tpr\ngro=a.gro\nxtc=a.xtc\nfragment-file=core.txt\n$end\n'
        others = UNAVAILABLE_MODULES[UNAVAILABLE_MODULES.index('$one-mol'):]
        for enabled in ('exciton','surface'):
            control = self.write(**{enabled:1}, modules=geometry+others)
            self.assertEqual(control.plan()['tasks'], ['get_dimers_monomers'])
            self.assertEqual(set(control.plan()['skipped_modules']), {'one-mol','two-mol','reorganization-energy'})
            with self.assertRaises(NotImplementedError):control.execution_plan()
            self.assertEqual(control.execution_plan(['get_dimer'])['tasks'], ['get_dimers_monomers'])
            self.assertEqual(control.plan(['one-mol'])['tasks'], [])

    def test_new_entry_reuses_shared_modules_once_and_cli_recognizes_it(self):
        xyz = self.root/'h2.xyz';xyz.write_text('2\nH2\nH 0 0 0\nH .74 0 0\n')
        with patch.dict(CALCULATIONS):
            register_calculation('new-calculation', CalculationType(('one-mol',)))
            text = ENTRIES.format(charge=1,exciton=0,surface=0)
            text += '$new-calculation\ndo = 1\n$end\n$one-mol\nxyz=h2.xyz\n$end\n'
            self.path.write_text(text)
            control = load_input(self.path)
            plan = control.plan()
            self.assertEqual(plan['tasks'], ['one-mol'])
            self.assertEqual(plan['module_owners']['one-mol'], ['charge-carrier','new-calculation'])
            self.assertEqual(control.plan(['new-calculation'])['module_owners']['one-mol'], ['new-calculation'])
            from universal_kmc.workflow.cli import main
            with redirect_stdout(io.StringIO()) as output:
                self.assertEqual(main([str(self.path),'--plan','--only','new-calculation']), 0)
            self.assertEqual(json.loads(output.getvalue())['calculation_tasks'], ['new-calculation'])
        self.assertNotIn('new-calculation', CALCULATIONS)

    def test_registered_solver_dispatch_needs_no_core_parser_changes(self):
        calls=[]
        def solver(control, preparation_report, *, resume):
            calls.append((preparation_report['tasks'],resume))
            return {'status':'completed','test_value':7}
        with patch.dict(CALCULATIONS):
            register_calculation('new-calculation', CalculationType((),solver=solver))
            self.path.write_text('$new-calculation\ndo = 1\n$end\n')
            with redirect_stdout(io.StringIO()):
                self.assertEqual(run_workflow(load_input(self.path),resume=True),0)
            self.assertEqual(calls,[({},True)])
            report=json.loads((self.root/'cube-results/workflow.json').read_text())
            self.assertEqual(report['calculations']['new-calculation'],{'status':'completed','test_value':7})

    def test_preparation_does_not_claim_completed_mobility(self):
        control=self.write(charge=1,modules='')
        with redirect_stdout(io.StringIO()):self.assertEqual(run_workflow(control),2)
        report=json.loads((self.root/'cube-results/workflow.json').read_text())
        self.assertEqual(report['status'],'prepared')
        self.assertEqual(report['calculations']['charge-carrier']['reason'],'solver_not_implemented')

    def test_failed_solver_marks_the_entry_failed(self):
        def fail(control, preparation_report, *, resume):
            raise RuntimeError('solver failed')
        with patch.dict(CALCULATIONS):
            register_calculation('new-calculation',CalculationType((),solver=fail))
            self.path.write_text('$new-calculation\ndo=1\n$end\n')
            with self.assertRaisesRegex(RuntimeError,'solver failed'):
                run_workflow(load_input(self.path))
            report=json.loads((self.root/'cube-results/workflow.json').read_text())
            self.assertEqual(report['status'],'failed')
            self.assertEqual(report['calculations']['new-calculation']['status'],'failed')

    def test_only_cannot_borrow_extraction_from_unselected_entry(self):
        for name in ('a.tpr','a.gro','a.xtc','core.txt'):(self.root/name).write_text('placeholder')
        with patch.dict(CALCULATIONS):
            register_calculation('new-calculation',CalculationType(('one-mol',)))
            self.path.write_text(ENTRIES.format(charge=1,exciton=0,surface=0)+
                '$new-calculation\ndo=1\n$end\n$one-mol\n$end\n$get_dimers_monomers\n'
                'tpr=a.tpr\ngro=a.gro\nxtc=a.xtc\nfragment-file=core.txt\n$end\n')
            control=load_input(self.path)
            with self.assertRaisesRegex(ValueError,'owned extraction'):
                control.plan(['new-calculation','one-mol'])

    def test_registry_rejects_duplicate_names_unknown_modules_and_bad_solver(self):
        with patch.dict(CALCULATIONS):
            for name,spec in [('charge-carrier',CalculationType(())),('one-mol',CalculationType(())),
                              ('bad_name',CalculationType(())),('new',CalculationType(('missing',))),
                              ('new',CalculationType(('one-mol','one-mol')))]:
                with self.subTest(name=name),self.assertRaises(ValueError):register_calculation(name,spec)
            with self.assertRaises(TypeError):register_calculation('new',CalculationType((),solver='not callable'))

    def test_entry_configuration_can_move_without_changing_the_plan(self):
        (self.root/'h2.xyz').write_text('2\nH2\nH 0 0 0\nH .74 0 0\n')
        entry='$charge-carrier\ndo=1\n$end\n'
        module='$one-mol\nxyz=h2.xyz\n$end\n'
        self.path.write_text(entry+module)
        first=load_input(self.path).plan()
        self.path.write_text(module+entry)
        second=load_input(self.path).plan()
        self.assertEqual(first['tasks'],second['tasks'])
        self.assertEqual(first['enabled_calculations'],second['enabled_calculations'])
        self.assertEqual(first['module_owners'],second['module_owners'])


if __name__=='__main__':unittest.main()
