"""
Tests for atom objects.
"""
from numpy import ndarray
import VMD

from pyvmd.atoms import Atom, Chain, Residue, Segment, Selection, NOW
from pyvmd.molecules import Molecule
from .utils import data, PyvmdTestCase


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
        self.assertEqual(atom.molecule, mol)
        self.assertEqual(atom.frame, NOW)
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
            atom.molecule = Molecule(molid)
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

    def test_comparison(self):
        # Test atom comparison
        molid_1 = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        mol_1 = Molecule(molid_1)
        molid_2 = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        mol_2 = Molecule(molid_2)
        VMD.molecule.set_top(molid_1)

        atom1 = Atom(0)
        atom2 = Atom(0)
        other = Atom(1)
        frame = Atom(0, frame=0)
        same_mol = Atom(0, mol_1)
        other_mol = Atom(0, mol_2)

        self.assertEqual(atom1, atom1)
        self.assertEqual(atom1, atom2)
        self.assertEqual(atom1, same_mol)
        self.assertNotEqual(atom1, other)
        self.assertNotEqual(atom2, other)
        self.assertNotEqual(atom1, frame)
        self.assertNotEqual(atom1, other_mol)

    def test_hashability(self):
        # Test atoms are correctly understood by sets
        molid_1 = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        mol_1 = Molecule(molid_1)
        molid_2 = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        mol_2 = Molecule(molid_2)
        VMD.molecule.set_top(molid_1)

        self.assertEqual({Atom(0), Atom(0)}, {Atom(0)})
        self.assertEqual({Atom(0, frame=0), Atom(0, frame=0)}, {Atom(0, frame=0)})
        self.assertEqual({Atom(0), Atom(0, frame=0)}, {Atom(0), Atom(0, frame=0)})
        self.assertEqual({Atom(0, mol_1), Atom(0, mol_1)}, {Atom(0, mol_1)})
        self.assertEqual({Atom(0, mol_1), Atom(0, mol_2)}, {Atom(0, mol_1), Atom(0, mol_2)})
        self.assertEqual({Atom(0), Atom(0, mol_1)}, {Atom(0)})

    def test_frame_property(self):
        # Test frame property works correctly
        molid = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        VMD.molecule.read(molid, 'dcd', data('water.1.dcd'), waitfor=-1)
        mol = Molecule(molid)

        # Create atom linked to the molecule frame
        atom = Atom(0)
        self.assertAlmostEqual(atom.x, -1.3421925)
        # Change the molecule frame, the atom will follow
        mol.frame = 0
        self.assertAlmostEqual(atom.x, -1.493)
        mol.frame = 4
        self.assertAlmostEqual(atom.x, -1.4773947)

        # Define the atom's frame
        atom.frame = 9
        self.assertAlmostEqual(atom.x, -1.4120502)
        # Change the molecule frame, the atom keeps its frame
        mol.frame = 11
        self.assertAlmostEqual(atom.x, -1.4120502)

        # Set the atom's frame back to the molecule's frame
        atom.frame = NOW
        self.assertAlmostEqual(atom.x, -1.3674825)

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
        self.assertEqual(res.molecule, mol)
        self.assertEqual(res.frame, NOW)
        self.assertEqual(res.number, 1)
        self.assertEqual(res.name, 'TIP3')

        # Test setters
        with self.assertRaises(AttributeError):
            res.index = 42
        with self.assertRaises(AttributeError):
            res.molecule = Molecule(molid)

        res.number = 42
        sel = VMD.atomsel.atomsel('residue 0', molid=molid)
        self.assertEqual(sel.get('resid'), [42, 42, 42])
        self.assertEqual(res.number, 42)

        res.name = 'WAT'
        self.assertEqual(sel.get('resname'), ['WAT', 'WAT', 'WAT'])
        self.assertEqual(res.name, 'WAT')

        # Test container methods
        self.assertEqual(len(res), 3)
        self.assertEqual(list(res), [Atom(0), Atom(1), Atom(2)])
        self.assertIn(Atom(0), res)
        self.assertIn(Atom(2), res)
        self.assertNotIn(Atom(3), res)
        self.assertNotIn(Atom(15), res)

    def test_comparison(self):
        # Test residue comparison
        molid_1 = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        mol_1 = Molecule(molid_1)
        molid_2 = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        mol_2 = Molecule(molid_2)
        VMD.molecule.set_top(molid_1)

        res1 = Residue(0)
        res2 = Residue(0)
        other = Residue(1)
        frame = Residue(0, frame=0)
        same_mol = Residue(0, mol_1)
        other_mol = Residue(0, mol_2)

        self.assertEqual(res1, res1)
        self.assertEqual(res1, res2)
        self.assertEqual(res1, same_mol)
        self.assertNotEqual(res1, other)
        self.assertNotEqual(res2, other)
        self.assertNotEqual(res1, frame)
        self.assertNotEqual(res1, other_mol)

    def test_hashability(self):
        # Test atoms are correctly understood by sets
        molid_1 = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        mol_1 = Molecule(molid_1)
        molid_2 = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        mol_2 = Molecule(molid_2)
        VMD.molecule.set_top(molid_1)

        self.assertEqual({Residue(0), Residue(0)}, {Residue(0)})
        self.assertEqual({Residue(0, frame=0), Residue(0, frame=0)}, {Residue(0, frame=0)})
        self.assertEqual({Residue(0), Residue(0, frame=0)}, {Residue(0), Residue(0, frame=0)})
        self.assertEqual({Residue(0, mol_1), Residue(0, mol_1)}, {Residue(0, mol_1)})
        self.assertEqual({Residue(0, mol_1), Residue(0, mol_2)}, {Residue(0, mol_1), Residue(0, mol_2)})
        self.assertEqual({Residue(0), Residue(0, mol_1)}, {Residue(0)})

    def test_frame_property(self):
        # Test frame property works correctly
        molid = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        VMD.molecule.read(molid, 'dcd', data('water.1.dcd'), waitfor=-1)
        mol = Molecule(molid)
        # To check frame works correctly, check x coordinate of the residue's first atom.

        # Create residue linked to the molecule frame
        res = Residue(0)
        self.assertAlmostEqual(list(res)[0].x, -1.3421925)
        # Change the molecule frame, the residue will follow
        mol.frame = 0
        self.assertAlmostEqual(list(res)[0].x, -1.493)
        mol.frame = 4
        self.assertAlmostEqual(list(res)[0].x, -1.4773947)

        # Define the residue's frame
        res.frame = 9
        self.assertAlmostEqual(list(res)[0].x, -1.4120502)
        # Change the molecule frame, the residue keeps its frame
        mol.frame = 11
        self.assertAlmostEqual(list(res)[0].x, -1.4120502)

        # Set the residue's frame back to the molecule's frame
        res.frame = NOW
        self.assertAlmostEqual(list(res)[0].x, -1.3674825)


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
        self.assertEqual(chain.molecule, mol)
        self.assertEqual(chain.frame, NOW)

        # Test name setter - when chain is renamed, all atoms still belong to the chain
        chain.name = 'A'
        sel = VMD.atomsel.atomsel('chain A', molid=molid)
        self.assertEqual(sel.get('chain'), ['A'] * 21)
        self.assertEqual(chain.name, 'A')

        # Test container methods
        self.assertEqual(len(chain), 21)
        self.assertEqual(list(chain), [Atom(i) for i in xrange(21)])
        self.assertIn(Atom(0), chain)
        self.assertIn(Atom(2), chain)
        empty_chain = Chain('X')
        self.assertNotIn(Atom(3), empty_chain)
        self.assertNotIn(Atom(15), empty_chain)

    def test_comparison(self):
        # Test residue comparison
        molid_1 = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        mol_1 = Molecule(molid_1)
        molid_2 = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        mol_2 = Molecule(molid_2)
        VMD.molecule.set_top(molid_1)

        chain1 = Chain('W')
        chain2 = Chain('W')
        other = Chain('X')
        frame = Chain('W', frame=0)
        same_mol = Chain('W', mol_1)
        other_mol = Chain('W', mol_2)

        self.assertEqual(chain1, chain1)
        self.assertEqual(chain1, chain2)
        self.assertEqual(chain1, same_mol)
        self.assertNotEqual(chain1, other)
        self.assertNotEqual(chain2, other)
        self.assertNotEqual(chain1, frame)
        self.assertNotEqual(chain1, other_mol)

    def test_frame_property(self):
        # Test frame property works correctly
        molid = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        VMD.molecule.read(molid, 'dcd', data('water.1.dcd'), waitfor=-1)
        mol = Molecule(molid)
        # To check frame works correctly, check x coordinate of the chain's first atom.

        # Create chain linked to the molecule frame
        chain = Chain('W')
        self.assertAlmostEqual(list(chain)[0].x, -1.3421925)
        # Change the molecule frame, the chain will follow
        mol.frame = 0
        self.assertAlmostEqual(list(chain)[0].x, -1.493)
        mol.frame = 4
        self.assertAlmostEqual(list(chain)[0].x, -1.4773947)

        # Define the chain's frame
        chain.frame = 9
        self.assertAlmostEqual(list(chain)[0].x, -1.4120502)
        # Change the molecule frame, the chain keeps its frame
        mol.frame = 11
        self.assertAlmostEqual(list(chain)[0].x, -1.4120502)

        # Set the chain's frame back to the molecule's frame
        chain.frame = NOW
        self.assertAlmostEqual(list(chain)[0].x, -1.3674825)


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
        self.assertEqual(segment.molecule, mol)
        self.assertEqual(segment.frame, NOW)

        # Test name setter - when segment is renamed, all atoms still belong to the segment
        segment.name = 'A'
        sel = VMD.atomsel.atomsel('segname A', molid=molid)
        self.assertEqual(sel.get('segname'), ['A'] * 21)
        self.assertEqual(segment.name, 'A')

        # Test container methods
        self.assertEqual(len(segment), 21)
        self.assertEqual(list(segment), [Atom(i) for i in xrange(21)])
        self.assertIn(Atom(0), segment)
        self.assertIn(Atom(2), segment)
        empty_segment = Segment('X')
        self.assertNotIn(Atom(3), empty_segment)
        self.assertNotIn(Atom(15), empty_segment)

    def test_comparison(self):
        # Test residue comparison
        molid_1 = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        mol_1 = Molecule(molid_1)
        molid_2 = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        mol_2 = Molecule(molid_2)
        VMD.molecule.set_top(molid_1)

        segment1 = Segment('W1')
        segment2 = Segment('W1')
        other = Segment('X')
        frame = Segment('W1', frame=0)
        same_mol = Segment('W1', mol_1)
        other_mol = Segment('W1', mol_2)

        self.assertEqual(segment1, segment1)
        self.assertEqual(segment1, segment2)
        self.assertEqual(segment1, same_mol)
        self.assertNotEqual(segment1, other)
        self.assertNotEqual(segment2, other)
        self.assertNotEqual(segment1, frame)
        self.assertNotEqual(segment1, other_mol)

    def test_frame_property(self):
        # Test frame property works correctly
        molid = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        VMD.molecule.read(molid, 'dcd', data('water.1.dcd'), waitfor=-1)
        mol = Molecule(molid)
        # To check frame works correctly, check x coordinate of the segment's first atom.

        # Create segment linked to the molecule frame
        segment = Segment('W1')
        self.assertAlmostEqual(list(segment)[0].x, -1.3421925)
        # Change the molecule frame, the segment will follow
        mol.frame = 0
        self.assertAlmostEqual(list(segment)[0].x, -1.493)
        mol.frame = 4
        self.assertAlmostEqual(list(segment)[0].x, -1.4773947)

        # Define the segment's frame
        segment.frame = 9
        self.assertAlmostEqual(list(segment)[0].x, -1.4120502)
        # Change the molecule frame, the segment keeps its frame
        mol.frame = 11
        self.assertAlmostEqual(list(segment)[0].x, -1.4120502)

        # Set the segment's frame back to the molecule's frame
        segment.frame = NOW
        self.assertAlmostEqual(list(segment)[0].x, -1.3674825)


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
        self.assertEqual(sel.molecule, mol)
        self.assertEqual(sel.frame, NOW)
        # Test setters
        with self.assertRaises(AttributeError):
            sel.selection = 'none'
        with self.assertRaises(AttributeError):
            sel.molecule = Molecule(molid)

        # Test container methods
        self.assertEqual(len(sel), 9)
        self.assertEqual(list(sel), [Atom(0), Atom(1), Atom(2), Atom(3), Atom(4), Atom(5), Atom(6), Atom(7), Atom(8)])
        self.assertIn(Atom(0), sel)
        self.assertIn(Atom(5), sel)
        self.assertNotIn(Atom(9), sel)
        self.assertNotIn(Atom(15), sel)

        # Test empty selection
        sel = Selection('none')
        self.assertEqual(sel.selection, 'none')
        self.assertEqual(sel.molecule, mol)
        self.assertEqual(sel.frame, NOW)
        # Test container methods
        self.assertEqual(len(sel), 0)
        self.assertEqual(list(sel), [])
        self.assertNotIn(Atom(0), sel)
        self.assertNotIn(Atom(5), sel)

    def test_frame_property(self):
        # Test frame property works correctly
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
