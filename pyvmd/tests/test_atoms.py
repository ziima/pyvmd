"""
Tests for atom objects.
"""
from numpy import ndarray
import VMD

from pyvmd.atoms import Atom, Chain, Residue, Segment, Selection, NOW
from pyvmd.molecules import Molecule
from .utils import data, PyvmdTestCase


class MoleculePropertyTest(PyvmdTestCase):
    """
    Test `molecule` property in atom objects.
    """
    def _test_molecule_property(self, obj_type, obj_id):
        # Test molecule property works correctly
        molid_1 = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        mol_1 = Molecule(molid_1)
        molid_2 = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        mol_2 = Molecule(molid_2)
        VMD.molecule.set_top(molid_1)

        # If molecule is not defined, top molecule is used
        obj = obj_type(obj_id)
        self.assertEqual(obj.molecule, mol_1)

        # If molecule is defined, it is used
        obj = obj_type(obj_id, mol_2)
        self.assertEqual(obj.molecule, mol_2)

        # Molecule can't be changed
        with self.assertRaises(AttributeError):
            obj.molecule = mol_2

    def test_atom(self):
        self._test_molecule_property(Atom, 0)

    def test_residue(self):
        self._test_molecule_property(Residue, 0)

    def test_chain(self):
        self._test_molecule_property(Chain, 'W')

    def test_segment(self):
        self._test_molecule_property(Segment, 'W1')

    def test_selection(self):
        self._test_molecule_property(Selection, 'index < 10')


class FramePropertyTest(PyvmdTestCase):
    """
    Test `frame` property in atom objects.
    """
    def _test_frame_property(self, obj_type, obj_id, getter):
        # Test frame property works correctly
        molid = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        VMD.molecule.read(molid, 'dcd', data('water.1.dcd'), waitfor=-1)
        mol = Molecule(molid)
        # To check frame works correctly, check the 'x' coordinate of the object's first atom.

        # Create residue linked to the molecule frame
        obj = obj_type(obj_id)
        self.assertEqual(obj.frame, NOW)
        self.assertAlmostEqual(getter(obj), -1.3421925)
        # Change the molecule frame, the object will follow
        mol.frame = 0
        self.assertEqual(obj.frame, NOW)
        self.assertAlmostEqual(getter(obj), -1.493)
        mol.frame = 4
        self.assertEqual(obj.frame, NOW)
        self.assertAlmostEqual(getter(obj), -1.4773947)

        # Define the object's frame
        obj.frame = 9
        self.assertEqual(obj.frame, 9)
        self.assertAlmostEqual(getter(obj), -1.4120502)
        # Change the molecule frame, the object keeps its frame
        mol.frame = 11
        self.assertEqual(obj.frame, 9)
        self.assertAlmostEqual(getter(obj), -1.4120502)

        # Set the object's frame back to the molecule's frame
        obj.frame = NOW
        self.assertEqual(obj.frame, NOW)
        self.assertAlmostEqual(getter(obj), -1.3674825)

    def test_atom(self):
        self._test_frame_property(Atom, 0, lambda obj: obj.x)

    def test_residue(self):
        self._test_frame_property(Residue, 0, lambda obj: list(obj)[0].x)

    def test_chain(self):
        self._test_frame_property(Chain, 'W', lambda obj: list(obj)[0].x)

    def test_segment(self):
        self._test_frame_property(Segment, 'W1', lambda obj: list(obj)[0].x)

    def test_selection(self):
        self._test_frame_property(Selection, 'index < 10', lambda obj: list(obj)[0].x)


class ComparisonTest(PyvmdTestCase):
    """
    Test comparison of atom objects.
    """
    def _test_equality(self, obj_type, obj_id1, obj_id2):
        molid_1 = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        mol_1 = Molecule(molid_1)
        molid_2 = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        mol_2 = Molecule(molid_2)
        VMD.molecule.set_top(molid_1)

        obj1 = obj_type(obj_id1)
        obj2 = obj_type(obj_id1)
        other = obj_type(obj_id2)
        frame = obj_type(obj_id1, frame=0)
        same_mol = obj_type(obj_id1, mol_1)
        other_mol = obj_type(obj_id1, mol_2)

        self.assertEqual(obj1, obj1)
        self.assertEqual(obj1, obj2)
        self.assertEqual(obj1, same_mol)
        self.assertNotEqual(obj1, other)
        self.assertNotEqual(obj2, other)
        self.assertNotEqual(obj1, frame)
        self.assertNotEqual(obj1, other_mol)

    def test_atom(self):
        self._test_equality(Atom, 0, 1)

    def test_residue(self):
        self._test_equality(Residue, 0, 1)

    def test_chain(self):
        self._test_equality(Chain, 'W', 'X')

    def test_segment(self):
        self._test_equality(Segment, 'W1', 'X')

    # Selection object doesn't have equality defined

    def test_different_objects(self):
        molid_1 = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))

        self.assertNotEqual(Atom(0), Residue(0))
        self.assertNotEqual(Chain('X'), Segment('X'))


class HashabilityTest(PyvmdTestCase):
    """
    Test hashability of atoms and residues.
    """
    def _test_hashability(self, obj_type):
        # Test objects are correctly understood by sets
        molid_1 = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        mol_1 = Molecule(molid_1)
        molid_2 = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        mol_2 = Molecule(molid_2)
        VMD.molecule.set_top(molid_1)

        # Same objects
        self.assertEqual({obj_type(0), obj_type(0)}, {obj_type(0)})
        self.assertEqual({obj_type(0, frame=0), obj_type(0, frame=0)}, {obj_type(0, frame=0)})
        self.assertEqual({obj_type(0, mol_1), obj_type(0, mol_1)}, {obj_type(0, mol_1)})
        # Top molecule explicitely and implicitly
        self.assertEqual({obj_type(0), obj_type(0, mol_1)}, {obj_type(0)})
        # Frame differs
        self.assertEqual({obj_type(0), obj_type(0, frame=0)}, {obj_type(0), obj_type(0, frame=0)})
        # Molecule differs
        self.assertEqual({obj_type(0, mol_1), obj_type(0, mol_2)}, {obj_type(0, mol_1), obj_type(0, mol_2)})

    def test_atom(self):
        self._test_hashability(Atom)

    def test_residue(self):
        self._test_hashability(Residue)


class ContainerTest(PyvmdTestCase):
    """
    Test container API for objects which supports it.
    """
    def _test_container(self, obj_type, obj_id, length, atom_ids, out_ids):
        molid_1 = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        mol_1 = Molecule(molid_1)
        molid_2 = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        mol_2 = Molecule(molid_2)
        VMD.molecule.set_top(molid_1)

        obj = obj_type(obj_id)

        # Check __len__
        self.assertEqual(len(obj), length)

        # Check __iter__
        self.assertEqual(list(obj), [Atom(i) for i in atom_ids])

        # Check __contains__
        for atom_id in atom_ids:
            self.assertIn(Atom(atom_id), obj)
        for atom_id in atom_ids:
            # Frame differs
            self.assertNotIn(Atom(atom_id, frame=0), obj)
        for atom_id in atom_ids:
            # Molecule differs
            self.assertNotIn(Atom(atom_id, mol_2), obj)
        for atom_id in out_ids:
            self.assertNotIn(Atom(atom_id), obj)

    def test_residue(self):
        self._test_container(Residue, 0, 3, range(3), range(4, 21))

    def test_chain(self):
        self._test_container(Chain, 'W', 15, range(15), range(16, 21))

    def test_empty_chain(self):
        self._test_container(Chain, 'X', 0, [], range(21))

    def test_segment(self):
        self._test_container(Segment, 'W1', 12, range(12), range(13, 21))

    def test_empty_segment(self):
        self._test_container(Segment, 'X', 0, [], range(21))

    def test_selection(self):
        self._test_container(Selection, 'resid 1 to 3', 9, range(9), range(10, 21))

    def test_empty_selection(self):
        self._test_container(Selection, 'none', 0, [], range(21))


class TestAtom(PyvmdTestCase):
    """
    Test `Atom` class.
    """
    def test_properties(self):
        # Test getters and setters
        molid = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        mol = Molecule(molid)

        atom = Atom(0)
        # Test getters
        self.assertEqual(atom.index, 0)
        self.assertAlmostEqual(atom.x, -1.493)
        self.assertAlmostEqual(atom.y, 1.900)
        self.assertAlmostEqual(atom.z, 1.280)
        self.assertIsInstance(atom.coords, ndarray)
        self.assertAlmostEqualSeqs(list(atom.coords), [-1.493, 1.9, 1.28])
        self.assertEqual(atom.name, 'OH2')
        self.assertEqual(atom.type, 'OT')
        self.assertEqual(atom.element, 'O')
        self.assertAlmostEqual(atom.beta, 0.0)
        self.assertAlmostEqual(atom.occupancy, 1.0)
        self.assertAlmostEqual(atom.mass, 15.9994, places=6)
        self.assertAlmostEqual(atom.charge, -0.834)
        self.assertAlmostEqual(atom.radius, 1.52)
        self.assertEqual(list(atom.bonded), [Atom(1), Atom(2)])
        self.assertEqual(atom.residue, Residue(0))
        self.assertEqual(atom.chain, Chain('W'))
        self.assertEqual(atom.segment, Segment('W1'))

        # Test setters
        with self.assertRaises(AttributeError):
            atom.index = 42
        with self.assertRaises(AttributeError):
            atom.residue = Residue(4)

        sel = VMD.atomsel.atomsel('index 0', molid=molid)

        # Test setters for coordinates
        atom.x = 23.9
        atom.y = -200.45
        atom.z = 0
        #XXX: There are some troubles with rounding in set
        self.assertAlmostEqualSeqs(sel.get('x'), [23.9], places=6)
        self.assertAlmostEqualSeqs(sel.get('y'), [-200.45], places=5)
        self.assertAlmostEqualSeqs(sel.get('z'), [0])
        self.assertAlmostEqual(atom.x, 23.9, places=6)
        self.assertAlmostEqual(atom.y, -200.45, places=5)
        self.assertAlmostEqual(atom.z, 0)
        self.assertAlmostEqualSeqs(list(atom.coords), [23.9, -200.45, 0], places=5)

        # Set complete coordinates
        atom.coords = (-90.56, 42, 17.85)
        #XXX: There are some troubles with rounding in set
        self.assertAlmostEqualSeqs(sel.get('x'), [-90.56], places=5)
        self.assertAlmostEqualSeqs(sel.get('y'), [42])
        self.assertAlmostEqualSeqs(sel.get('z'), [17.85], places=6)
        self.assertAlmostEqual(atom.x, -90.56, places=5)
        self.assertAlmostEqual(atom.y, 42)
        self.assertAlmostEqual(atom.z, 17.85, places=6)
        self.assertAlmostEqualSeqs(list(atom.coords), [-90.56, 42, 17.85], places=5)

        # Test setters for other attributes
        atom.name = 'NEW'
        self.assertEqual(atom.name, 'NEW')
        self.assertEqual(sel.get('name'), ['NEW'])
        atom.type = 'N18'
        self.assertEqual(atom.type, 'N18')
        self.assertEqual(sel.get('type'), ['N18'])
        atom.element = 'Y'
        self.assertEqual(atom.element, 'Y')
        self.assertEqual(sel.get('element'), ['Y'])
        atom.beta = 2.5
        self.assertEqual(atom.beta, 2.5)
        self.assertEqual(sel.get('beta'), [2.5])
        atom.occupancy = -3.8
        self.assertAlmostEqual(atom.occupancy, -3.8)
        self.assertAlmostEqualSeqs(sel.get('occupancy'), [-3.8])
        atom.mass = 42.89
        self.assertAlmostEqual(atom.mass, 42.89, places=5)
        self.assertAlmostEqualSeqs(sel.get('mass'), [42.89], places=5)
        atom.charge = -7.05
        self.assertAlmostEqual(atom.charge, -7.05, places=6)
        self.assertAlmostEqualSeqs(sel.get('charge'), [-7.05], places=6)
        atom.radius = 4.9
        self.assertAlmostEqual(atom.radius, 4.9, places=6)
        self.assertAlmostEqualSeqs(sel.get('radius'), [4.9], places=6)

        atom.chain = Chain('A')
        # Ensure only this atom was moved to the new chain
        self.assertEqual(VMD.atomsel.atomsel('chain A').get('index'), [0])
        self.assertEqual(atom.chain, Chain('A'))

        atom.segment = Segment('S')
        # Ensure only this atom was moved to the new segment
        self.assertEqual(VMD.atomsel.atomsel('segname S').get('index'), [0])
        self.assertEqual(atom.segment, Segment('S'))

    def test_constructors(self):
        # Test various ways of Atom instance creation
        molid1 = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        VMD.molecule.read(molid1, 'dcd', data('water.1.dcd'), waitfor=-1)
        mol1 = Molecule(molid1)
        molid2 = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        mol2 = Molecule(molid2)
        VMD.molecule.set_top(molid2)

        # Top molecule and NOW
        atom = Atom(0)
        self.assertEqual(atom.molecule, mol2)
        self.assertEqual(atom.frame, NOW)
        self.assertAlmostEqual(atom.x, -1.493)

        # Molecule and frame
        atom = Atom(0, molecule=mol1, frame=5)
        self.assertEqual(atom.molecule, mol1)
        self.assertEqual(atom.frame, 5)
        self.assertAlmostEqual(atom.x, -1.4746015)

        # Get atom using selection string
        atom = Atom.pick('resid 2 and name OH2')
        self.assertEqual(atom.molecule, mol2)
        self.assertEqual(atom.frame, NOW)
        self.assertAlmostEqual(atom.x, 0.337)

        # Get atom using selection string, molecule and frame
        atom = Atom.pick('resid 2 and name OH2', molecule=mol1, frame=8)
        self.assertEqual(atom.molecule, mol1)
        self.assertEqual(atom.frame, 8)
        self.assertAlmostEqual(atom.x, 0.3521036)

        # Atom which does not exist
        with self.assertRaises(ValueError):
            Atom(8947)

        # Selection which returns none or too many atoms
        with self.assertRaises(ValueError):
            Atom.pick('all')
        with self.assertRaises(ValueError):
            Atom.pick('none')


class TestResidue(PyvmdTestCase):
    """
    Test `Residue` class.
    """
    def test_properties(self):
        # Test getters and setters
        molid = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        mol = Molecule(molid)

        res = Residue(0)
        # Test getters
        self.assertEqual(res.index, 0)
        self.assertEqual(res.number, 1)
        self.assertEqual(res.name, 'TIP3')

        # Test setters
        with self.assertRaises(AttributeError):
            res.index = 42

        res.number = 42
        sel = VMD.atomsel.atomsel('residue 0', molid=molid)
        self.assertEqual(sel.get('resid'), [42, 42, 42])
        self.assertEqual(res.number, 42)

        res.name = 'WAT'
        self.assertEqual(sel.get('resname'), ['WAT', 'WAT', 'WAT'])
        self.assertEqual(res.name, 'WAT')


class TestChain(PyvmdTestCase):
    """
    Test `Chain` class.
    """
    def test_properties(self):
        # Test getters and setters
        molid = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        mol = Molecule(molid)

        chain = Chain('W')
        # Test getters
        self.assertEqual(chain.name, 'W')

        # Test name setter - when chain is renamed, all atoms still belong to the chain
        chain.name = 'A'
        sel = VMD.atomsel.atomsel('chain A', molid=molid)
        self.assertEqual(sel.get('chain'), ['A'] * 15)
        self.assertEqual(chain.name, 'A')


class TestSegment(PyvmdTestCase):
    """
    Test `Segment` class.
    """
    def test_properties(self):
        # Test getters and setters
        molid = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        mol = Molecule(molid)

        segment = Segment('W1')
        # Test getters
        self.assertEqual(segment.name, 'W1')

        # Test name setter - when segment is renamed, all atoms still belong to the segment
        segment.name = 'A'
        sel = VMD.atomsel.atomsel('segname A', molid=molid)
        self.assertEqual(sel.get('segname'), ['A'] * 12)
        self.assertEqual(segment.name, 'A')


class TestSelection(PyvmdTestCase):
    """
    Test `Selection` class.
    """
    def test_properties(self):
        # Test basic properties
        molid = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        mol = Molecule(molid)

        sel = Selection('resid 1 to 3')
        # Test getters
        self.assertEqual(sel.selection, 'resid 1 to 3')
        # Test setters
        with self.assertRaises(AttributeError):
            sel.selection = 'none'

    def test_selection_updates(self):
        # Test selection updates if frame is changed
        molid = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        VMD.molecule.read(molid, 'dcd', data('water.1.dcd'), waitfor=-1)
        mol = Molecule(molid)

        # Create selection linked to the molecule frame
        sel = Selection('x < -1.4')
        self.assertEqual(list(sel), [Atom(9), Atom(18), Atom(19), Atom(20)])
        # Change the molecule frame, the selection should follow
        mol.frame = 0
        self.assertEqual(list(sel), [Atom(0), Atom(1), Atom(9), Atom(10), Atom(18), Atom(19), Atom(20)])
        mol.frame = 4
        self.assertEqual(list(sel), [Atom(0), Atom(1), Atom(9), Atom(18), Atom(19), Atom(20)])

        # Define the selection's frame
        result = [Atom(0, frame=9), Atom(1, frame=9), Atom(9, frame=9), Atom(18, frame=9), Atom(19, frame=9),
                  Atom(20, frame=9)]
        sel.frame = 9
        self.assertEqual(list(sel), result)
        # Change the molecule frame, the selection keeps its frame
        mol.frame = 11
        self.assertEqual(list(sel), result)

        # Set the selection's frame back to the molecule's frame
        sel.frame = NOW
        self.assertEqual(list(sel), [Atom(9), Atom(18), Atom(19), Atom(20)])

    def test_contacts(self):
        # Test `contacts` method
        VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))

        sel1 = Selection('resid 1 to 3 and noh')
        sel2 = Selection('hydrogen')

        self.assertEqual(list(sel1.contacts(sel2, 1.0)), [])
        self.assertEqual(list(sel1.contacts(sel2, 2.0)), [(Atom(6), Atom(11))])
        self.assertEqual(list(sel2.contacts(sel1, 2.0)), [(Atom(11), Atom(6))])
