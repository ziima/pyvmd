"""
Tests for analysis.
"""
import VMD

from pyvmd.analysis import hydrogen_bonds, HydrogenBond
from pyvmd.atoms import Atom, Selection

from .utils import data, PyvmdTestCase


class TestAnalysis(PyvmdTestCase):
    """
    Test analysis utilities.
    """
    def setUp(self):
        molid = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        self.molid = molid

    def test_hydrogen_bonds(self):
        # Test `hydrogen_bonds` function
        sel = Selection('noh')

        result = [HydrogenBond(Atom(18), Atom(19), Atom(9)), HydrogenBond(Atom(9), Atom(11), Atom(6)),
                  HydrogenBond(Atom(6), Atom(7), Atom(15))]
        self.assertEqual(list(hydrogen_bonds(sel)), result)
        self.assertEqual(list(hydrogen_bonds(sel, sel)), result)
        result = [HydrogenBond(Atom(18), Atom(19), Atom(9)), HydrogenBond(Atom(9), Atom(11), Atom(6)),
                  HydrogenBond(Atom(6), Atom(7), Atom(9)), HydrogenBond(Atom(6), Atom(7), Atom(15))]
        self.assertEqual(list(hydrogen_bonds(sel, angle=75)), result)
        self.assertEqual(list(hydrogen_bonds(sel, angle=180)), [])
        self.assertEqual(list(hydrogen_bonds(sel, distance=2)), [])

        # If the selections do not share same atoms, check the hydrogen bonds are returned only in correct direction
        sel1 = Selection('index 6')
        sel2 = Selection('index 9')
        self.assertEqual(list(hydrogen_bonds(sel1, sel2, angle=75)), [HydrogenBond(Atom(6), Atom(7), Atom(9))])
