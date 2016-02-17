"""
Tests for trajectory analysis utilities.
"""
import VMD
from mock import sentinel

from pyvmd.analyzer import Analyzer
from pyvmd.molecules import Molecule

from .utils import data, PyvmdTestCase


class TestAnalyzer(PyvmdTestCase):
    """
    Test `Analyzer` class.
    """
    def setUp(self):
        self.mol = Molecule.create()
        self.mol.load(data('water.psf'))
        # Storage for callback data
        self.coords = []
        self.frames = []

    def test_analyze_callback_args(self):
        # Test callback is called with extra arguments
        steps = []

        def callback(step, *args, **kwargs):
            self.assertEqual(step.molecule, self.mol)
            self.assertEqual(args, (sentinel.arg1, sentinel.arg2))
            self.assertEqual(kwargs, {'key': sentinel.value, 'another': sentinel.junk})
            steps.append(step.frame)

        analyzer = Analyzer(self.mol, [data('water.1.dcd')])
        analyzer.add_callback(callback, sentinel.arg1, sentinel.arg2, key=sentinel.value, another=sentinel.junk)
        analyzer.analyze()

        self.assertEqual(steps, range(12))

    def _get_status(self, status):
        # Callback to collect status data
        self.frames.append(status.frame)

    def _get_x(self, status):
        # Callback to collect data
        self.coords.append(VMD.atomsel.atomsel('index 0', molid=status.molecule.molid).get('x')[0])

    def test_analyze(self):
        # Test analyzer works correctly with default parameters
        analyzer = Analyzer(self.mol, [data('water.1.dcd'), data('water.2.dcd')])
        analyzer.add_callback(self._get_status)
        analyzer.add_callback(self._get_x)
        analyzer.analyze()
        result = [-1.4911567, -1.4851371, -1.4858487, -1.4773947, -1.4746015, -1.4673382, -1.4535547, -1.4307435,
                  -1.4120502, -1.3853478, -1.3674825, -1.3421925, -1.3177859, -1.2816998, -1.2579591, -1.2262495,
                  -1.2036057, -1.1834533, -1.174916, -1.1693807, -1.1705244, -1.1722997, -1.1759951, -1.175245]
        self.assertAlmostEqualSeqs(self.coords, result)
        self.assertEqual(self.frames, range(24))

    def test_analyze_params(self):
        # Test load every other frame, all 12 at once
        self.coords = []
        self.frames = []
        analyzer = Analyzer(self.mol, [data('water.1.dcd'), data('water.2.dcd')], step=2, chunk=12)
        analyzer.add_callback(self._get_status)
        analyzer.add_callback(self._get_x)
        analyzer.analyze()
        result = [-1.4911567, -1.4858487, -1.4746015, -1.4535547, -1.4120502, -1.3674825, -1.3177859, -1.2579591,
                  -1.2036057, -1.174916, -1.1705244, -1.1759951]
        self.assertAlmostEqualSeqs(self.coords, result)
        self.assertEqual(self.frames, range(12))
