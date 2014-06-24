"""
Tests for trajectory analysis utilities.
"""
import VMD

from pyvmd.traj import Loader

from .utils import data, PyvmdTestCase


X_COORDS = [-1.4911567, -1.4851371, -1.4858487, -1.4773947, -1.4746015, -1.4673382, -1.4535547, -1.4307435, -1.4120502,
            -1.3853478, -1.3674825, -1.3421925, -1.3177859, -1.2816998, -1.2579591, -1.2262495, -1.2036057, -1.1834533,
            -1.174916, -1.1693807, -1.1705244, -1.1722997, -1.1759951, -1.175245]


class TestLoader(PyvmdTestCase):
    """
    Test `Loader` class.
    """
    def test_load(self):
        # Test loader works correctly
        x_coords = []
        frames = []

        # Callback to collect status data
        def _get_status(status):
            frames.append(status.frame)

        # Callback to collect data
        def _get_x(status):
            x_coords.append(VMD.atomsel.atomsel('index 0', molid=status.molecule.molid).get('x')[0])

        loader = Loader(data('water.psf'), [data('water.1.dcd'), data('water.2.dcd')])
        loader.add_callback(_get_status)
        loader.add_callback(_get_x)
        loader.run()
        self.assertAlmostEqualSeqs(x_coords, X_COORDS)
        self.assertEqual(frames, range(24))
