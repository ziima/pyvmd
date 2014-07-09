"""
Tests for trajectory analysis utilities.
"""
import VMD

from pyvmd.molecules import Molecule
from pyvmd.traj import Loader

from .utils import data, PyvmdTestCase


class TestLoader(PyvmdTestCase):
    """
    Test `Loader` class.
    """
    def setUp(self):
        self.mol = Molecule.create()
        self.mol.load(data('water.psf'))
        # Storage for callback data
        self.coords = []
        self.frames = []

    def _get_status(self, status):
        # Callback to collect status data
        self.frames.append(status.frame)

    def _get_x(self, status):
        # Callback to collect data
        self.coords.append(VMD.atomsel.atomsel('index 0', molid=status.molecule.molid).get('x')[0])

    def test_load(self):
        # Test loader works correctly with default parameters

        loader = Loader(self.mol, [data('water.1.dcd'), data('water.2.dcd')])
        loader.add_callback(self._get_status)
        loader.add_callback(self._get_x)
        loader.run()
        result = [-1.4911567, -1.4851371, -1.4858487, -1.4773947, -1.4746015, -1.4673382, -1.4535547, -1.4307435,
                  -1.4120502, -1.3853478, -1.3674825, -1.3421925, -1.3177859, -1.2816998, -1.2579591, -1.2262495,
                  -1.2036057, -1.1834533, -1.174916, -1.1693807, -1.1705244, -1.1722997, -1.1759951, -1.175245]
        self.assertAlmostEqualSeqs(self.coords, result)
        self.assertEqual(self.frames, range(24))

    def test_load_params(self):
        # Test load every other frame, all 12 at once
        self.coords = []
        self.frames = []
        loader = Loader(self.mol, [data('water.1.dcd'), data('water.2.dcd')], step=2, chunk=12)
        loader.add_callback(self._get_status)
        loader.add_callback(self._get_x)
        loader.run()
        result = [-1.4911567, -1.4858487, -1.4746015, -1.4535547, -1.4120502, -1.3674825, -1.3177859, -1.2579591,
                  -1.2036057, -1.174916, -1.1705244, -1.1759951]
        self.assertAlmostEqualSeqs(self.coords, result)
        self.assertEqual(self.frames, range(12))
