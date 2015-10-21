"""
Tests for data collectors.
"""
from cStringIO import StringIO

from pyvmd.analyzer import Analyzer
from pyvmd.atoms import Selection
from pyvmd.collectors import (AngleCollector, DihedralCollector, DistanceCollector, RMSDCollector, XCoordCollector,
                              YCoordCollector, ZCoordCollector)
from pyvmd.datasets import DataSet
from pyvmd.molecules import Molecule

from .utils import data, PyvmdTestCase


class TestCollectors(PyvmdTestCase):
    """
    Test collectors.
    """
    def setUp(self):
        # Prepare molecule
        self.mol = Molecule.create()
        self.mol.load(data('water.psf'))

    def test_coordinate_collectors(self):
        # Test coordinate collector
        dset = DataSet()
        dset.add_collector(XCoordCollector('index 0'))
        dset.add_collector(XCoordCollector('all'))
        dset.add_collector(YCoordCollector('index 0'))
        dset.add_collector(YCoordCollector('all'))
        dset.add_collector(ZCoordCollector('index 0'))
        dset.add_collector(ZCoordCollector('all'))
        analyzer = Analyzer(self.mol, [data('water.1.dcd')])
        analyzer.add_dataset(dset)
        analyzer.analyze()

        # Write data to check result
        buf = StringIO()
        dset.write(buf)
        # Check the result
        self.assertEqual(buf.getvalue(), open(data('coords.dat')).read())

    def test_geometry_collectors(self):
        # Test geometry collectors - distance, angle, dihedral and improper.
        dset = DataSet()
        dset.add_collector(DistanceCollector('index 0', 'index 1'))
        dset.add_collector(DistanceCollector('resid 1', 'all'))
        dset.add_collector(AngleCollector('index 0', 'index 1', 'index 2'))
        dset.add_collector(AngleCollector('resid 1', 'resid 2', 'resid 3'))
        dset.add_collector(DihedralCollector('index 0', 'index 1', 'index 2', 'index 3'))
        dset.add_collector(DihedralCollector('resid 1', 'resid 2', 'resid 3', 'resid 4'))
        analyzer = Analyzer(self.mol, [data('water.1.dcd')])
        analyzer.add_dataset(dset)
        analyzer.analyze()

        # Write data to check result
        buf = StringIO()
        dset.write(buf)
        # Check the result
        self.assertEqual(buf.getvalue(), open(data('geometry.dat')).read())

    def _test_collector_error(self, collector):
        dset = DataSet()
        dset.add_collector(collector)
        analyzer = Analyzer(self.mol, [data('water.1.dcd')])
        analyzer.add_dataset(dset)
        self.assertRaises(ValueError, analyzer.analyze)

    def test_selection_errors(self):
        # Test collectors raises errors on empty selections
        self._test_collector_error(XCoordCollector('none'))
        self._test_collector_error(YCoordCollector('none'))
        self._test_collector_error(ZCoordCollector('none'))
        self._test_collector_error(DistanceCollector('none', 'index 0'))
        self._test_collector_error(DistanceCollector('index 0', 'none'))
        self._test_collector_error(AngleCollector('none', 'index 0', 'index 1'))
        self._test_collector_error(AngleCollector('index 0', 'none', 'index 1'))
        self._test_collector_error(AngleCollector('index 0', 'index 1', 'none'))
        self._test_collector_error(DihedralCollector('none', 'index 0', 'index 1', 'index 2'))
        self._test_collector_error(DihedralCollector('index 0', 'none', 'index 1', 'index 2'))
        self._test_collector_error(DihedralCollector('index 0', 'index 1', 'none', 'index 2'))
        self._test_collector_error(DihedralCollector('index 0', 'index 1', 'index 2', 'none'))

    def test_rmsd_collector(self):
        # Test RMSD collector
        ref = Molecule.create()
        ref.load(data('water.psf'))
        ref.load(data('water.pdb'))
        dset = DataSet()
        dset.add_collector(RMSDCollector('all', Selection('all', ref)))
        dset.add_collector(RMSDCollector('all and name OH2', Selection('all and name OH2', ref)))
        dset.add_collector(RMSDCollector('all and noh', Selection('all and noh', ref), name='noh'))
        analyzer = Analyzer(self.mol, [data('water.1.dcd')])
        analyzer.add_dataset(dset)
        analyzer.analyze()

        # Write data to check result
        buf = StringIO()
        dset.write(buf)
        # Check the result
        self.assertEqual(buf.getvalue(), open(data('rmsd.dat')).read())
