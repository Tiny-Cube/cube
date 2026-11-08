"""Behavioral checks for the standalone viewer's embedded interactive renderer."""
import ast
import json
from pathlib import Path
import shutil
import subprocess
import unittest


@unittest.skipUnless(shutil.which('node'), 'Node.js required to check browser renderer')
class GROViewerRenderingTests(unittest.TestCase):
    def test_full_scene_projection_preview_recovery_and_export(self):
        root = Path(__file__).resolve().parents[1]
        tree = ast.parse((root / 'scripts/view_clusters.py').read_text())
        html = next(ast.literal_eval(node.value) for node in tree.body
                    if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'HTML' for t in node.targets))
        script = html.split('<script>')[1].split('</script>')[0].replace('__INITIAL__', '{"center":1,"cutoff":30}')
        # More than the preview limit exercises bounded sampling and complete restoration.
        count = 13020
        groups = [{'molecule': i + 1, 'label': str(i + 1), 'image': [0, 0, 0],
                   'start': i * 20, 'count': 20, 'minimum_distance': i / 100}
                  for i in range(count // 20)]
        cluster = {'center': 1, 'cutoff': 30, 'coordinates': [[1 + i % 31, i % 27, i % 19] for i in range(count)],
                   'elements': ['H' if i % 3 == 0 else 'C' for i in range(count)],
                   'bonds': [[i, i + 1] for i in range(count - 1) if i // 20 == (i + 1) // 20],
                   'atom_groups': [i // 20 for i in range(count)], 'groups': groups,
                   'central_atoms': 20, 'neighbor_atoms': count - 20, 'neighbor_images': len(groups) - 1,
                   'periodic_images': 0, 'minimum_distance': .01, 'names': ['C'] * count,
                   'serials': list(range(count)), 'global_indices': list(range(count))}
        for size in (120, count):
            sample = {**cluster, 'coordinates': cluster['coordinates'][:size], 'elements': cluster['elements'][:size],
                      'atom_groups': cluster['atom_groups'][:size], 'bonds': [b for b in cluster['bonds'] if b[1] < size],
                      'groups': groups[:size // 20], 'neighbor_atoms': size - 20, 'neighbor_images': size // 20 - 1}
            with self.subTest(atom_count=size):
                result = subprocess.run(['node', str(root / 'tests/viewer_rendering.cjs')],
                                        input=json.dumps({'script': script, 'cluster': sample}),
                                        text=True, capture_output=True, timeout=30)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(json.loads(result.stdout)['checks']['atoms'], size)


if __name__ == '__main__':
    unittest.main()
