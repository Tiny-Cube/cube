"""Main orchestration and independent submodule execution share one dispatcher."""
from contextlib import redirect_stdout
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from universal_kmc.workflow.cli import main
from universal_kmc.workflow.control import load_input
from universal_kmc.workflow.runner import run_workflow

ENTRIES = '$charge-carrier\ndo={charge}\n$end\n$exciton\ndo={exciton}\n$end\n$surface-hopping\ndo=0\n$end\n'
UNUSED = '''$get_dimers_monomers
do=0
tpr=absent.tpr
gro=absent.gro
xtc=absent.xtc
fragment-file=absent.txt
$end
$two-mol
do=0
qc-engine=unused
$end
$reorganization-energy
do=0
environment=unused
$end
'''


class SubmoduleSwitchTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.path=self.root/'input'
        (self.root/'h2.xyz').write_text('2\nH2 synthetic fixture\nH 0 0 0\nH 0 0 .74\n')

    def write(self, module, *, charge=0, exciton=0, others=UNUSED):
        self.path.write_text(ENTRIES.format(charge=charge,exciton=exciton)+module+others)
        return load_input(self.path)

    def test_child_runs_with_all_mains_disabled_and_ignores_unused_files(self):
        control=self.write('$one-mol\ndo=1\nxyz=h2.xyz\n$end\n')
        plan=control.execution_plan()
        self.assertEqual(plan['tasks'],['one-mol'])
        self.assertEqual(plan['requested_modules'],['one-mol'])
        self.assertEqual(plan['calculation_tasks'],[])
        self.assertEqual(plan['module_owners'],{'one-mol':['standalone']})
        self.assertEqual(control.module_switches['one-mol'],1)
        self.assertNotIn('do',control.sections['one-mol'])
        self.assertEqual(set(control.inactive_sections),{'get_dimers_monomers','two-mol','reorganization-energy'})

    @unittest.skipUnless(importlib.util.find_spec('numpy'), 'install NumPy for shared scientific dispatch')
    def test_main_and_standalone_request_share_one_execution(self):
        control=self.write('$one-mol\ndo=1\nxyz=h2.xyz\n$end\n',charge=1)
        self.assertEqual(control.plan()['module_owners']['one-mol'],['charge-carrier','standalone'])
        with patch('universal_kmc.workflow.scf.run_single_point',return_value={'status':'completed'}) as single, redirect_stdout(io.StringIO()):
            self.assertEqual(run_workflow(control),2) # mobility remains unimplemented
        self.assertEqual(single.call_count,1)
        report=json.loads((self.root/'cube-results/workflow.json').read_text())
        self.assertEqual(report['requested_modules'],['one-mol'])

    def test_explicit_zero_disables_even_with_main_and_only(self):
        control=self.write('$one-mol\ndo=0\nxyz=missing.xyz\nqc-engine=unused\n$end\n',charge=1)
        self.assertEqual(control.plan()['tasks'],[])
        self.assertEqual(control.plan(['one-mol'])['tasks'],[])
        self.assertEqual(control.plan()['skipped_modules']['one-mol'],'disabled_by_module_do')

    def test_missing_do_retains_automatic_main_behavior(self):
        module='$one-mol\nxyz=h2.xyz\n$end\n'
        self.assertEqual(self.write(module).plan()['tasks'],[])
        self.assertEqual(self.write(module,charge=1).plan()['tasks'],['one-mol'])

    def test_legacy_controls_respect_module_do_without_main_sections(self):
        self.path.write_text('$one-mol\ndo=1\nxyz=h2.xyz\n$end\n'+UNUSED)
        plan=load_input(self.path).plan()
        self.assertEqual(plan['tasks'],['one-mol'])
        self.assertEqual(plan['requested_modules'],['one-mol'])
        self.assertEqual(plan['module_owners']['one-mol'],['standalone'])

    def test_standalone_module_needs_no_main_ownership(self):
        control=self.write('$one-mol\ndo=1\nxyz=h2.xyz\n$end\n',exciton=1)
        # A full exciton request is still unsupported; selecting the independently
        # requested child does not pretend to run excited-state dynamics.
        with self.assertRaises(NotImplementedError):control.execution_plan()
        self.assertEqual(control.execution_plan(['one-mol'])['tasks'],['one-mol'])
        self.assertEqual(control.plan(['one-mol'])['calculation_tasks'],[])

    def test_required_extraction_is_added_unless_explicitly_disabled(self):
        for name in ('a.tpr','a.gro','a.xtc','core.txt'):(self.root/name).write_text('plan fixture')
        extraction='$get_dimers_monomers\ntpr=a.tpr\ngro=a.gro\nxtc=a.xtc\nfragment-file=core.txt\n$end\n'
        module='$one-mol\ndo=1\n$end\n'
        control=self.write(module,others=extraction)
        self.assertEqual(control.plan()['tasks'],['get_dimers_monomers','one-mol'])
        self.assertEqual(control.plan()['module_owners']['get_dimers_monomers'],['standalone'])
        self.assertEqual(control.plan(['one-mol'])['tasks'],['get_dimers_monomers','one-mol'])
        with self.assertRaisesRegex(ValueError,'do=0.*extracted data'):
            self.write(module,others=extraction.replace('$get_dimers_monomers\n','$get_dimers_monomers\ndo=0\n'))
        # Explicit existing coordinates remove that dependency entirely.
        self.assertEqual(self.write(module.replace('do=1','do=1\nxyz=h2.xyz'),others=UNUSED).plan()['tasks'],['one-mol'])

    def test_bad_do_is_rejected_for_all_children_even_when_mains_disabled(self):
        for name in ('get_dimers_monomers','one-mol','two-mol','reorganization-energy'):
            for value in ('true','2','-1','1.0','01',''):
                with self.subTest(module=name,value=value),self.assertRaisesRegex(ValueError,'do must be 0 or 1'):
                    self.write(f'${name}\ndo={value}\n$end\n',others='')

    def test_cli_plan_shows_standalone_request_without_scientific_execution(self):
        self.write('$one-mol\ndo=1\nxyz=h2.xyz\n$end\n')
        with redirect_stdout(io.StringIO()) as captured:
            self.assertEqual(main([str(self.path),'--plan']),0)
        self.assertEqual(json.loads(captured.getvalue())['requested_modules'],['one-mol'])
        self.assertFalse((self.root/'cube-results').exists())

    @unittest.skipUnless(importlib.util.find_spec('pyscf'), 'install the quantum extra for standalone SCF execution')
    def test_real_standalone_scf_and_resume_with_all_mains_off(self):
        self.write('$one-mol\ndo=1\nxyz=h2.xyz\nmethod=hf\nbasis-set=sto-3g\nthreads=1\n$end\n')
        with redirect_stdout(io.StringIO()):
            self.assertEqual(main([str(self.path)]),0)
            self.assertEqual(main([str(self.path),'--resume']),0)
        report=json.loads((self.root/'cube-results/workflow.json').read_text())
        self.assertEqual(report['status'],'completed')
        self.assertEqual(report['task_order'],['one-mol'])
        self.assertTrue(all(value['status']=='disabled' for value in report['calculations'].values()))
        record=report['tasks']['one-mol']['outputs']['h2']
        self.assertTrue(Path(record['checkpoint']).is_file())
        self.assertTrue(record['reused'])

    @unittest.skipUnless(importlib.util.find_spec('pyscf'), 'install the quantum extra for SCF reuse')
    def test_scientific_output_is_reusable_after_changing_orchestration(self):
        settings='$one-mol\nxyz=h2.xyz\nmethod=hf\nbasis-set=sto-3g\nthreads=1\n$end\n'
        self.write(settings,charge=1)
        with redirect_stdout(io.StringIO()):self.assertEqual(main([str(self.path)]),2)
        self.write(settings.replace('$one-mol\n','$one-mol\ndo=1\n'))
        with redirect_stdout(io.StringIO()):self.assertEqual(main([str(self.path),'--resume']),0)
        report=json.loads((self.root/'cube-results/workflow.json').read_text())
        self.assertTrue(report['tasks']['one-mol']['outputs']['h2']['reused'])


if __name__=='__main__':unittest.main()
