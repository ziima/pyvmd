"""
Tests for data collectors.
"""
from cStringIO import StringIO

from pyvmd.collectors import RMSDCollector
from pyvmd.molecules import Molecule
from pyvmd.traj import Loader

from .utils import data, PyvmdTestCase


class TestRMSDCollector(PyvmdTestCase):
    """
    Test RMSD collector.
    """
    def test_rmsd(self):
        # Test RMSD collector returns correct results
        # Set up the reference, collector and loader
        ref = Molecule.create()
        ref.load(data('water.psf'))
        ref.load(data('water.pdb'))
        mol = Molecule.create()
        mol.load(data('water.psf'))
        rmsd = RMSDCollector(ref)
        rmsd.add_selection('all')
        rmsd.add_selection('all and name OH2')
        rmsd.add_selection('all and noh', 'noh')
        loader = Loader(mol, [data('water.1.dcd')])
        loader.add_collector(rmsd)
        loader.run()

        # Write data to check result
        buf = StringIO()
        rmsd.dataset.write(buf)
        # Check the result
        self.assertEqual(buf.getvalue(), open(data('rmsd.dat')).read())

    def test_invalid_names(self):
        ref = Molecule.create()
        rmsd = RMSDCollector(ref)
        rmsd.add_selection('all', 'unique')
        self.assertRaises(ValueError, rmsd.add_selection, 'noh', 'unique')
