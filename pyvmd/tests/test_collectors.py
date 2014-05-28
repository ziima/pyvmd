"""
Tests for data collectors.
"""
import os
from tempfile import mkstemp
import unittest

from Molecule import Molecule

from pyvmd.traj import Loader
from pyvmd.collectors import RMSDCollector

from .utils import data


class TestRMSDCollector(unittest.TestCase):
    """
    Test RMSD collector.
    """
    def setUp(self):
        dummy, filename = mkstemp(prefix='pyvmd_test_')
        self.tmpfile = filename

    def tearDown(self):
        os.unlink(self.tmpfile)

    def test_rmsd(self):
        # Test RMSD collector returns correct results
        # Set up the reference, collector and loader
        ref = Molecule()
        ref.load(data('water.psf'))
        ref.load(data('water.pdb'))
        rmsd = RMSDCollector(ref, self.tmpfile)
        rmsd.add_selection('all')
        rmsd.add_selection('all and name OH2')
        rmsd.add_selection('all and noh', 'noh')
        loader = Loader(data('water.psf'), [data('water.1.dcd')])
        loader.add_collector(rmsd)
        loader.run()

        # Check the result
        self.assertEqual(open(self.tmpfile).read(), open(data('rmsd.dat')).read())

    def test_invalid_names(self):
        ref = Molecule()
        rmsd = RMSDCollector(ref, self.tmpfile)
        rmsd.add_selection('all', 'unique')
        self.assertRaises(ValueError, rmsd.add_selection, 'noh', 'unique')
