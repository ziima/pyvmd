"""
Tests for trajectory analysis utilities.
"""
import os
import unittest

import VMD

from pyvmd.traj import Loader


def data(filename):
    """
    Return full filename of the test datafile.
    """
    return os.path.join(os.path.dirname(__file__), 'data', filename)


X_COORDS = [-1.4911566972732544, -1.4851371049880981, -1.4858486652374268, -1.477394700050354, -1.4746015071868896,
            -1.46733820438385, -1.4535547494888306, -1.4307434558868408, -1.4120502471923828, -1.385347843170166,
            -1.3674825429916382, -1.342192530632019, -1.3177858591079712, -1.281699776649475, -1.2579591274261475,
            -1.2262494564056396, -1.2036056518554688, -1.1834533214569092, -1.1749160289764404, -1.1693806648254395,
            -1.1705243587493896, -1.1722997426986694, -1.175995111465454, -1.1752450466156006]


class TestLoader(unittest.TestCase):
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
            x_coords.append(VMD.atomsel.atomsel('index 0', molid=status.molecule.id).get('x')[0])

        loader = Loader(data('water.psf'), [data('water.1.dcd'), data('water.2.dcd')])
        loader.add_callback(_get_status)
        loader.add_callback(_get_x)
        loader.run()
        self.assertEqual(x_coords, X_COORDS)
        self.assertEqual(frames, range(24))
