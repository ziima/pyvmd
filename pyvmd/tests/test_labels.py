"""
Tests for labels utilities.
"""
import VMD

from pyvmd.atoms import Atom
from pyvmd.labels import (ANGLE, ANGLE_LABELS, AngleLabel, ATOM, ATOM_LABELS, AtomLabel, BOND, BOND_LABELS, BondLabel,
                          DIHEDRAL, DIHEDRAL_LABELS, DihedralLabel)
from pyvmd.molecules import Molecule, MOLECULES

from .utils import data, PyvmdTestCase


def clear_labels():
    """
    Utility to delete all labels.
    """
    for category in ATOM, BOND, ANGLE, DIHEDRAL:
        for label in VMD.label.listall(category):
            VMD.label.delete(category, label)


class TestLabels(PyvmdTestCase):
    """
    Test labels classes.
    """
    def setUp(self):
        self.mol = Molecule.create()
        self.mol.load(data('water.psf'))
        self.mol.load(data('water.pdb'))
        self.mol.load(data('water.1.dcd'))
        self.mol.frame = 0

        self.other = Molecule.create()
        self.other.load(data('water.psf'))
        self.other.load(data('water.pdb'))

        MOLECULES.top = self.mol

        self.addCleanup(clear_labels)

    def test_atom_label(self):
        # Test atom labels
        a = Atom(0)
        b = Atom(4)
        VMD.label.add(ATOM, (self.mol.molid, ), (0, ))
        VMD.label.add(ATOM, (self.mol.molid, ), (4, ))
        VMD.label.add(ATOM, (self.other.molid, ), (0, ))

        # Check label
        label_1 = AtomLabel(a)
        self.assertEqual(label_1.category, ATOM)
        self.assertEqual(label_1.atoms, (a, ))
        self.assertEqual(label_1.atom, a)
        self.assertEqual(label_1, label_1)

        # Create other instance of the same label
        label_2 = AtomLabel(a)
        self.assertEqual(label_2, label_1)

        # Check other label
        label_3 = AtomLabel(b)
        self.assertEqual(label_3.category, ATOM)
        self.assertEqual(label_3.atoms, (b, ))
        self.assertEqual(label_3.atom, b)
        self.assertNotEqual(label_3, label_1)

        # Check label for other instance of the same atom
        label_4 = AtomLabel(Atom(0))
        self.assertEqual(label_4, label_1)

        # Check label for atom in another frame - it's the same label, because labels are independent on frames
        label_5 = AtomLabel(Atom(0, frame=4))
        self.assertEqual(label_5, label_1)

        # Check label to atom from other molecule is different
        label_6 = AtomLabel(Atom(0, self.other))
        self.assertNotEqual(label_6, label_1)

        # Check ValueError if label doesn't exist
        with self.assertRaises(ValueError):
            AtomLabel(Atom(5))

    def test_bond_label(self):
        # Test bond labels
        a = Atom(0)
        b = Atom(4)
        c = Atom(7)
        VMD.label.add(BOND, (self.mol.molid, self.mol.molid), (0, 4))
        VMD.label.add(BOND, (self.mol.molid, self.mol.molid), (4, 7))
        VMD.label.add(BOND, (self.other.molid, self.other.molid), (0, 4))

        # Check label
        label_1 = BondLabel(a, b)
        self.assertEqual(label_1.category, BOND)
        self.assertEqual(label_1.atoms, (a, b))
        self.assertEqual(label_1, label_1)

        # Create other instance of the same label
        label_2 = BondLabel(a, b)
        self.assertEqual(label_2, label_1)

        # Check other label
        label_3 = BondLabel(b, c)
        self.assertEqual(label_3.category, BOND)
        self.assertEqual(label_3.atoms, (b, c))
        self.assertNotEqual(label_3, label_1)

        # Check label with opposite order of atoms
        label_4 = BondLabel(b, a)
        self.assertEqual(label_4.atoms, (b, a))
        self.assertEqual(label_4, label_1)

        # Create label for atom in another frame - it's the same label, because labels are independent on frames
        label_5 = BondLabel(Atom(0, frame=4), Atom(4, frame=8))
        self.assertEqual(label_5, label_1)

        # Check label to atom from other molecule is different
        label_6 = BondLabel(Atom(0, self.other), Atom(4, self.other))
        self.assertNotEqual(label_6, label_1)

        # Check ValueError if label doesn't exist
        with self.assertRaises(ValueError):
            BondLabel(Atom(5), Atom(6))

    def test_angle_label(self):
        # Test angle labels
        a = Atom(0)
        b = Atom(4)
        c = Atom(7)
        VMD.label.add(ANGLE, (self.mol.molid, self.mol.molid, self.mol.molid), (0, 4, 7))
        VMD.label.add(ANGLE, (self.mol.molid, self.mol.molid, self.mol.molid), (4, 7, 0))
        VMD.label.add(ANGLE, (self.other.molid, self.other.molid, self.other.molid), (0, 4, 7))

        # Check label
        label_1 = AngleLabel(a, b, c)
        self.assertEqual(label_1.category, ANGLE)
        self.assertEqual(label_1.atoms, (a, b, c))
        self.assertEqual(label_1, label_1)

        # Create other instance of the same label
        label_2 = AngleLabel(a, b, c)
        self.assertEqual(label_2, label_1)

        # Check other label
        label_3 = AngleLabel(b, c, a)
        self.assertEqual(label_3.category, ANGLE)
        self.assertEqual(label_3.atoms, (b, c, a))
        self.assertNotEqual(label_3, label_1)

        # Check label with opposite order of atoms
        label_4 = AngleLabel(c, b, a)
        self.assertEqual(label_4.atoms, (c, b, a))
        self.assertEqual(label_4, label_1)

        # Create label for atom in another frame - it's the same label, because labels are independent on frames
        label_5 = AngleLabel(Atom(0, frame=4), Atom(4, frame=8), Atom(7, frame=2))
        self.assertEqual(label_5, label_1)

        # Check label to atom from other molecule is different
        label_6 = AngleLabel(Atom(0, self.other), Atom(4, self.other), Atom(7, self.other))
        self.assertNotEqual(label_6, label_1)

        # Check ValueError if label doesn't exist
        with self.assertRaises(ValueError):
            AngleLabel(Atom(5), Atom(6), Atom(7))

    def test_dihedral_label(self):
        # Test dihedral labels
        a = Atom(0)
        b = Atom(4)
        c = Atom(7)
        d = Atom(10)
        VMD.label.add(DIHEDRAL, (self.mol.molid, self.mol.molid, self.mol.molid, self.mol.molid), (0, 4, 7, 10))
        VMD.label.add(DIHEDRAL, (self.mol.molid, self.mol.molid, self.mol.molid, self.mol.molid), (4, 7, 10, 0))
        VMD.label.add(DIHEDRAL, (self.other.molid, self.other.molid, self.other.molid, self.other.molid), (0, 4, 7, 10))

        # Check label
        label_1 = DihedralLabel(a, b, c, d)
        self.assertEqual(label_1.category, DIHEDRAL)
        self.assertEqual(label_1.atoms, (a, b, c, d))
        self.assertEqual(label_1, label_1)

        # Create other instance of the same label
        label_2 = DihedralLabel(a, b, c, d)
        self.assertEqual(label_2, label_1)

        # Check other label
        label_3 = DihedralLabel(b, c, d, a)
        self.assertEqual(label_3.category, DIHEDRAL)
        self.assertEqual(label_3.atoms, (b, c, d, a))
        self.assertNotEqual(label_3, label_1)

        # Check label with opposite order of atoms
        label_4 = DihedralLabel(d, c, b, a)
        self.assertEqual(label_4.atoms, (d, c, b, a))
        self.assertEqual(label_4, label_1)

        # Create label for atom in another frame - it's the same label, because labels are independent on frames
        label_5 = DihedralLabel(Atom(0, frame=4), Atom(4, frame=8), Atom(7, frame=2), Atom(10, frame=3))
        self.assertEqual(label_5, label_1)

        # Check label to atom from other molecule is different
        label_6 = DihedralLabel(Atom(0, self.other), Atom(4, self.other), Atom(7, self.other), Atom(10, self.other))
        self.assertNotEqual(label_6, label_1)

        # Check ValueError if label doesn't exist
        with self.assertRaises(ValueError):
            DihedralLabel(Atom(5), Atom(6), Atom(7), Atom(8))

    def test_unequal_types(self):
        VMD.label.add(ATOM, (self.mol.molid, ), (0, ))
        atom_label = AtomLabel(Atom(0))
        VMD.label.add(BOND, (self.mol.molid, self.mol.molid), (0, 4))
        bond_label = BondLabel(Atom(0), Atom(4))
        VMD.label.add(ANGLE, (self.mol.molid, self.mol.molid, self.mol.molid), (0, 4, 7))
        angle_label = AngleLabel(Atom(0), Atom(4), Atom(7))
        VMD.label.add(DIHEDRAL, (self.mol.molid, self.mol.molid, self.mol.molid, self.mol.molid), (0, 4, 7, 10))
        dihedral_label = DihedralLabel(Atom(0), Atom(4), Atom(7), Atom(10))

        self.assertNotEqual(atom_label, bond_label)
        self.assertNotEqual(atom_label, angle_label)
        self.assertNotEqual(atom_label, dihedral_label)
        self.assertNotEqual(bond_label, angle_label)
        self.assertNotEqual(bond_label, dihedral_label)
        self.assertNotEqual(angle_label, dihedral_label)

    def test_create(self):
        atom_label = AtomLabel.create(Atom(0))
        self.assertEqual(VMD.label.listall(ATOM),
                         [{'atomid': (0, ), 'molid': (self.mol.molid, ), 'on': 1, 'value': 0.0}])
        self.assertTrue(atom_label.visible)

        bond_label = BondLabel.create(Atom(0), Atom(4))
        data_1 = {'atomid': (0, 4), 'molid': (self.mol.molid, self.mol.molid), 'on': 1, 'value': 4.476513862609863}
        self.assertEqual(VMD.label.listall(BOND), [data_1])
        self.assertTrue(bond_label.visible)

        angle_label = AngleLabel.create(Atom(0), Atom(4), Atom(7))
        data_2 = {'atomid': (0, 4, 7), 'molid': (self.mol.molid, self.mol.molid, self.mol.molid), 'on': 1,
                  'value': 54.093936920166016}
        self.assertEqual(VMD.label.listall(ANGLE), [data_2])
        self.assertTrue(angle_label.visible)

        dihedral_label = DihedralLabel.create(Atom(0), Atom(4), Atom(7), Atom(10))
        data_3 = {'atomid': (0, 4, 7, 10), 'molid': (self.mol.molid, self.mol.molid, self.mol.molid, self.mol.molid),
                  'on': 1, 'value': -80.11302947998047}
        self.assertEqual(VMD.label.listall(DIHEDRAL), [data_3])
        self.assertTrue(dihedral_label.visible)

    def test_delete(self):
        VMD.label.add(ATOM, (self.mol.molid, ), (0, ))
        atom_label = AtomLabel(Atom(0))
        atom_label.delete()
        self.assertEqual(VMD.label.listall(ATOM), [])
        # Label can't be deleted twice
        with self.assertRaises(ValueError):
            atom_label.delete()

        VMD.label.add(BOND, (self.mol.molid, self.mol.molid), (0, 4))
        bond_label = BondLabel(Atom(0), Atom(4))
        bond_label.delete()
        self.assertEqual(VMD.label.listall(BOND), [])
        # Label can't be deleted twice
        with self.assertRaises(ValueError):
            bond_label.delete()
        # Check label values
        with self.assertRaises(ValueError):
            bond_label.distances

        VMD.label.add(ANGLE, (self.mol.molid, self.mol.molid, self.mol.molid), (0, 4, 7))
        angle_label = AngleLabel(Atom(0), Atom(4), Atom(7))
        angle_label.delete()
        self.assertEqual(VMD.label.listall(ANGLE), [])
        # Label can't be deleted twice
        with self.assertRaises(ValueError):
            angle_label.delete()
        # Check label values
        with self.assertRaises(ValueError):
            angle_label.angles

        VMD.label.add(DIHEDRAL, (self.mol.molid, self.mol.molid, self.mol.molid, self.mol.molid), (0, 4, 7, 10))
        dihedral_label = DihedralLabel(Atom(0), Atom(4), Atom(7), Atom(10))
        dihedral_label.delete()
        self.assertEqual(VMD.label.listall(DIHEDRAL), [])
        # Label can't be deleted twice
        with self.assertRaises(ValueError):
            dihedral_label.delete()
        # Check label values
        with self.assertRaises(ValueError):
            dihedral_label.dihedrals

    def test_visibility(self):
        # Test `visible` property
        VMD.label.add(ATOM, (self.mol.molid, ), (0, ))
        label = AtomLabel(Atom(0))
        self.assertEqual(VMD.label.listall(ATOM),
                         [{'atomid': (0, ), 'molid': (self.mol.molid, ), 'on': 1, 'value': 0.0}])
        self.assertTrue(label.visible)

        # Show visible label doesn't change anything
        label.visible = True
        self.assertEqual(VMD.label.listall(ATOM),
                         [{'atomid': (0, ), 'molid': (self.mol.molid, ), 'on': 1, 'value': 0.0}])
        self.assertTrue(label.visible)

        # Hide label
        label.visible = False
        self.assertEqual(VMD.label.listall(ATOM),
                         [{'atomid': (0, ), 'molid': (self.mol.molid, ), 'on': 0, 'value': 0.0}])
        self.assertFalse(label.visible)

        # Hide again doesn't change anything
        label.visible = False
        self.assertEqual(VMD.label.listall(ATOM),
                         [{'atomid': (0, ), 'molid': (self.mol.molid, ), 'on': 0, 'value': 0.0}])
        self.assertFalse(label.visible)

        # Check direct changes to visibility are followed
        VMD.label.show(ATOM, {'atomid': (0, ), 'molid': (self.mol.molid, )})
        self.assertTrue(label.visible)
        VMD.label.hide(ATOM, {'atomid': (0, ), 'molid': (self.mol.molid, )})
        self.assertFalse(label.visible)

        label.delete()
        # Deleted label raises errors
        with self.assertRaises(ValueError):
            label.visible
        with self.assertRaises(ValueError):
            label.visible = True

    def test_values(self):
        VMD.label.add(BOND, (self.mol.molid, self.mol.molid), (0, 1))
        bond_label = BondLabel(Atom(0), Atom(1))
        distances = [1.0100465, 0.9646406, 0.9611206, 0.9662452, 0.9626279, 0.9649956, 0.9611024, 0.9658498, 0.9603688,
                     0.9671795, 0.961661, 0.9632944, 0.9573071]
        self.assertAlmostEqualSeqs(bond_label.distances, distances)

        VMD.label.add(ANGLE, (self.mol.molid, self.mol.molid, self.mol.molid), (0, 1, 2))
        angle_label = AngleLabel(Atom(0), Atom(1), Atom(2))
        angles = [38.9131393, 38.1194229, 38.2041321, 38.341835, 38.3854752, 38.1108742, 37.7457161, 37.9469795,
                  38.6991272, 38.7136955, 38.4580612, 38.6055489, 39.1369514]
        self.assertAlmostEqualSeqs(angle_label.angles, angles)

        VMD.label.add(DIHEDRAL, (self.mol.molid, self.mol.molid, self.mol.molid, self.mol.molid), (0, 1, 2, 3))
        dihedral_label = DihedralLabel(Atom(0), Atom(1), Atom(2), Atom(3))
        dihedrals = [54.5871239, 57.054245, 61.5406532, 68.0926056, 75.0311661, 81.3050842, 86.0759125, 89.2662735,
                     93.1087952, 99.3121338, 110.2715073, 125.2944031, 142.6539612]
        self.assertAlmostEqualSeqs(dihedral_label.dihedrals, dihedrals)


class TestLabelManagers(PyvmdTestCase):
    """
    Test label managers.
    """
    def setUp(self):
        self.addCleanup(clear_labels)

    def test_managers(self):
        molid = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))

        # There are no labels
        self.assertEqual(len(ATOM_LABELS), 0)
        self.assertEqual(len(BOND_LABELS), 0)
        self.assertEqual(len(ANGLE_LABELS), 0)
        self.assertEqual(len(DIHEDRAL_LABELS), 0)

        # Create a few labels directly in VMD
        VMD.label.add(ATOM, (molid, ), (0, ))
        VMD.label.add(ATOM, (molid, ), (1, ))
        VMD.label.add(ATOM, (molid, ), (2, ))
        VMD.label.add(ATOM, (molid, ), (3, ))
        VMD.label.add(BOND, (molid, molid), (2, 3))
        VMD.label.add(BOND, (molid, molid), (3, 4))
        VMD.label.add(BOND, (molid, molid), (4, 5))
        VMD.label.add(ANGLE, (molid, molid, molid), (2, 4, 6))
        VMD.label.add(ANGLE, (molid, molid, molid), (3, 5, 7))
        VMD.label.add(DIHEDRAL, (molid, molid, molid, molid), (8, 5, 3, 2))

        # Check manager's properties
        self.assertEqual(len(ATOM_LABELS), 4)
        self.assertEqual(len(BOND_LABELS), 3)
        self.assertEqual(len(ANGLE_LABELS), 2)
        self.assertEqual(len(DIHEDRAL_LABELS), 1)

        label_1 = AtomLabel(Atom(0))
        self.assertIn(label_1, ATOM_LABELS)
        self.assertNotIn(label_1, BOND_LABELS)
        self.assertNotIn(label_1, ANGLE_LABELS)
        self.assertNotIn(label_1, DIHEDRAL_LABELS)

        self.assertEqual(list(ATOM_LABELS), [label_1, AtomLabel(Atom(1)), AtomLabel(Atom(2)), AtomLabel(Atom(3))])
        self.assertEqual(list(BOND_LABELS),
                         [BondLabel(Atom(2), Atom(3)), BondLabel(Atom(3), Atom(4)), BondLabel(Atom(4), Atom(5))])
        self.assertEqual(list(ANGLE_LABELS),
                         [AngleLabel(Atom(2), Atom(4), Atom(6)), AngleLabel(Atom(3), Atom(5), Atom(7))])
        self.assertEqual(list(DIHEDRAL_LABELS), [DihedralLabel(Atom(8), Atom(5), Atom(3), Atom(2))])
