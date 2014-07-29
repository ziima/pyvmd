"""
Tests for data collectors.
"""
from cStringIO import StringIO

from pyvmd.analyzer import Analyzer
from pyvmd.atoms import Selection
from pyvmd.collectors import RMSDCollector
from pyvmd.datasets import DataSet
from pyvmd.molecules import Molecule

from .utils import data, PyvmdTestCase


class TestRMSDCollector(PyvmdTestCase):
    """
    Test RMSD collector.
    """
    def test_rmsd(self):
        # Test RMSD collector returns correct results
        # Set up the reference, collectors, dataset and analyzer
        ref = Molecule.create()
        ref.load(data('water.psf'))
        ref.load(data('water.pdb'))
        mol = Molecule.create()
        mol.load(data('water.psf'))
        dset = DataSet()
        dset.add_collector(RMSDCollector('all', Selection('all', ref)))
        dset.add_collector(RMSDCollector('all and name OH2', Selection('all and name OH2', ref)))
        dset.add_collector(RMSDCollector('all and noh', Selection('all and noh', ref), name='noh'))
        analyzer = Analyzer(mol, [data('water.1.dcd')])
        analyzer.add_dataset(dset)
        analyzer.analyze()

        # Write data to check result
        buf = StringIO()
        dset.write(buf)
        # Check the result
        self.assertEqual(buf.getvalue(), open(data('rmsd.dat')).read())
