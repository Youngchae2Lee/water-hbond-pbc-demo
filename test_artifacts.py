"""Check that the documented entry point produces usable, repeatable artifacts."""
import contextlib
import csv
import io
import json
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path
from unittest.mock import patch

import run_demo


class ArtifactTests(unittest.TestCase):
    def test_end_to_end_generation(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            with patch.object(run_demo, 'ROOT', root), contextlib.redirect_stdout(io.StringIO()):
                run_demo.main()
                first = {p.relative_to(root): p.read_bytes()
                         for p in root.rglob('*') if p.is_file()}
                run_demo.main()
                second = {p.relative_to(root): p.read_bytes()
                          for p in root.rglob('*') if p.is_file()}
            self.assertEqual(first, second)
            self.assertEqual(len(first), 7)
            with (root/'results'/'cutoff_sweep.csv').open(newline='', encoding='utf-8') as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual([int(row['hbond_pair_edges']) for row in rows], [0,2,2,3,3])
            self.assertEqual([int(row['distance_only_edges']) for row in rows], [0,2,2,3,4])
            results = json.loads((root/'results'/'analysis.json').read_text(encoding='utf-8'))
            self.assertEqual(results[3]['degrees'], [1,2,2,1])
            for figure in (root/'results').glob('*.svg'):
                self.assertEqual(ET.parse(figure).getroot().tag, '{http://www.w3.org/2000/svg}svg')


if __name__ == '__main__':
    unittest.main()
