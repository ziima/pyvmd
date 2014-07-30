"""
Tests for measure.
"""
import VMD

from pyvmd.atoms import Atom
from pyvmd.measure import angle, dihedral, distance

from .utils import data, PyvmdTestCase


class TestMeasure(PyvmdTestCase):
    """
    Test measure utilities.
    """
    def setUp(self):
        molid = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        self.molid = molid

    def test_distance(self):
        # Test `distance` function
        a = Atom(0)
        b = Atom(4)
        c = Atom(7)

        self.assertAlmostEqual(distance(a, b), 4.4765141)
        self.assertAlmostEqual(distance(b, a), 4.4765141)
        self.assertAlmostEqual(distance(a, c), 4.0660655)
        self.assertAlmostEqual(distance(c, a), 4.0660655)
        self.assertAlmostEqual(distance(b, c), 4.4653567)
        self.assertAlmostEqual(distance(c, b), 4.4653567)

    def test_angle(self):
        # Test `angle` function
        a = Atom(0)
        b = Atom(4)
        c = Atom(7)

        self.assertAlmostEqual(angle(a, b, c), 54.0939249)
        self.assertAlmostEqual(angle(c, b, a), 54.0939249)
        self.assertAlmostEqual(angle(b, c, a), 63.0930652)
        self.assertAlmostEqual(angle(a, c, b), 63.0930652)
        self.assertAlmostEqual(angle(c, a, b), 62.8130099)
        self.assertAlmostEqual(angle(b, a, c), 62.8130099)

    def test_dihedrals(self):
        # Test `dihedral` function
        a = Atom(0)
        b = Atom(4)
        c = Atom(7)
        d = Atom(10)

        self.assertAlmostEqual(dihedral(a, b, c, d), -80.113001)
        self.assertAlmostEqual(dihedral(d, c, b, a), -80.113001)
        self.assertAlmostEqual(dihedral(a, c, b, d), 80.113001)
        self.assertAlmostEqual(dihedral(d, b, c, a), 80.113001)
