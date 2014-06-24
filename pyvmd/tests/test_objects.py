"""
Tests for objects.
"""
from numpy import ndarray
from Molecule import Molecule as _Molecule
import VMD

from pyvmd.objects import Atom, Molecule, MoleculeManager, Selection, FORMAT_PDB, NOW
from .utils import data, PyvmdTestCase


class TestMolecule(PyvmdTestCase):
    """
    Test `Molecule` class.
    """
    def setUp(self):
        molid = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        VMD.molecule.read(molid, 'dcd', data('water.1.dcd'), waitfor=-1)
        self.molid = molid

    def test_molecule_properties(self):
        # Test basic properties of Molecule
        mol = Molecule(self.molid)

        # Check frame property
        mol.frame = 3
        self.assertEqual(VMD.molecule.get_frame(self.molid), 3)
        self.assertEqual(mol.frame, 3)

        mol.frame = 8
        self.assertEqual(VMD.molecule.get_frame(self.molid), 8)
        self.assertEqual(mol.frame, 8)

        # Check name property
        mol.name = 'My precious'
        self.assertEqual(mol.name, 'My precious')
        self.assertEqual(VMD.molecule.name(self.molid), 'My precious')
        mol.name = 'The Ring'
        self.assertEqual(mol.name, 'The Ring')
        self.assertEqual(VMD.molecule.name(self.molid), 'The Ring')

        # Check molecule property
        self.assertIsInstance(mol.molecule, _Molecule)

        # Check error if molecule does not exists
        self.assertRaises(ValueError, Molecule, 66000)

    def test_molecule_create(self):
        # Test molecule creation
        old_molids = VMD.molecule.listall()
        created = Molecule.create()
        self.assertNotIn(created.molid, old_molids)
        self.assertTrue(VMD.molecule.exists(created.molid))
        self.assertEqual(VMD.molecule.get_filenames(created.molid), [])
        self.assertEqual(VMD.molecule.get_filetypes(created.molid), [])

        # Check create with name
        other = Molecule.create('The One')
        self.assertEqual(other.name, 'The One')
        self.assertTrue(VMD.molecule.exists(other.molid))
        self.assertEqual(VMD.molecule.get_filenames(other.molid), [])
        self.assertEqual(VMD.molecule.get_filetypes(other.molid), [])

    def test_molecule_delete(self):
        # Test molecule deletion
        mol = Molecule(self.molid)
        mol.delete()
        self.assertFalse(VMD.molecule.exists(self.molid))

    def test_molecule_load(self):
        # Test file loading
        mol = Molecule.create()
        mol.load(data('water.psf'))
        self.assertEqual(VMD.molecule.get_filenames(mol.molid), [data('water.psf')])
        self.assertEqual(VMD.molecule.get_filetypes(mol.molid), ['psf'])
        mol.load(data('water.pdb'), FORMAT_PDB)
        self.assertEqual(VMD.molecule.get_filenames(mol.molid), [data('water.psf'), data('water.pdb')])
        self.assertEqual(VMD.molecule.get_filetypes(mol.molid), ['psf', 'pdb'])
        mol.load(data('water.1.dcd'))
        self.assertEqual(VMD.molecule.get_filenames(mol.molid), [data('water.psf'), data('water.pdb'),
                                                                   data('water.1.dcd')])
        self.assertEqual(VMD.molecule.get_filetypes(mol.molid), ['psf', 'pdb', 'dcd'])

        self.assertRaises(ValueError, mol.load, 'no_extension')

    def test_molecule_comparison(self):
        # Test molecule comparison
        mol1 = Molecule(self.molid)
        mol2 = Molecule(self.molid)
        other = Molecule.create()

        self.assertEqual(mol1, mol1)
        self.assertEqual(mol1, mol2)
        self.assertNotEqual(mol1, other)
        self.assertNotEqual(mol2, other)

    def test_frames(self):
        # Test frames wrapper

        # Utility to get x coordinate through the trajectory to check results
        sel = VMD.atomsel.atomsel('index 0', molid=self.molid)

        def _get_x_coord():
            # Utility to get x coordinate through trajectory
            values = []
            for frame in xrange(VMD.molecule.numframes(self.molid)):
                sel.frame = frame
                values.append(sel.get('x')[0])
            return values

        # Check number of frames and iterator
        mol = Molecule(self.molid)
        self.assertEqual(len(mol.frames), 13)
        self.assertEqual(list(mol.frames), range(13))

        # Test frame duplication - duplicate focused frame
        mol.frame = 3
        mol.frames.copy()

        # Check there is one more frame
        self.assertEqual(len(mol.frames), 14)
        coords = [-1.493, -1.4911567, -1.4851371, -1.4858487, -1.4773947, -1.4746015, -1.4673382, -1.4535547,
                  -1.4307435, -1.4120502, -1.3853478, -1.3674825, -1.3421925, -1.4858487]
        self.assertAlmostEqualSeqs(_get_x_coord(), coords)
        # Check molecule is focused to the new frame
        self.assertEqual(mol.frame, 13)

        # Test frame duplication - duplicate defined frame
        mol.frames.copy(4)
        # Check there is one more frame
        self.assertEqual(len(mol.frames), 15)
        coords = [-1.493, -1.4911567, -1.4851371, -1.4858487, -1.4773947, -1.4746015, -1.4673382, -1.4535547,
                  -1.4307435, -1.4120502, -1.3853478, -1.3674825, -1.3421925, -1.4858487, -1.4773947]
        self.assertAlmostEqualSeqs(_get_x_coord(), coords)
        # Molecule is focused to the new frame
        self.assertEqual(mol.frame, 14)

        # Test frame deletion - positive index
        del mol.frames[14]
        self.assertEqual(len(mol.frames), 14)
        coords = [-1.493, -1.4911567, -1.4851371, -1.4858487, -1.4773947, -1.4746015, -1.4673382, -1.4535547,
                  -1.4307435, -1.4120502, -1.3853478, -1.3674825, -1.3421925, -1.4858487]
        self.assertAlmostEqualSeqs(_get_x_coord(), coords)

        # Test frame deletion - negative index
        del mol.frames[-2]
        self.assertEqual(len(mol.frames), 13)
        coords = [-1.493, -1.4911567, -1.4851371, -1.4858487, -1.4773947, -1.4746015, -1.4673382, -1.4535547,
                  -1.4307435, -1.4120502, -1.3853478, -1.3674825, -1.4858487]
        self.assertAlmostEqualSeqs(_get_x_coord(), coords)

        # Check deletion of frame slice
        del mol.frames[2:11:3]
        self.assertEqual(len(mol.frames), 10)
        coords = [-1.493, -1.4911567, -1.4858487, -1.4773947, -1.4673382, -1.4535547, -1.4120502, -1.3853478,
                  -1.3674825, -1.4858487]
        self.assertAlmostEqualSeqs(_get_x_coord(), coords)

        # Invalid key for slice
        with self.assertRaises(TypeError):
            del mol.frames[None]


class TestMoleculeManager(PyvmdTestCase):
    """
    Test `MoleculeManager` class.
    """
    def test_no_molecules(self):
        man = MoleculeManager()

        self.assertEqual(len(man), 0)
        with self.assertRaises(ValueError):
            man[0]
        with self.assertRaises(ValueError):
            man['molecule']
        with self.assertRaises(TypeError):
            man[None]
        with self.assertRaises(ValueError):
            del man[0]
        with self.assertRaises(ValueError):
            del man['molecule']
        self.assertEqual(list(man), [])
        with self.assertRaises(AttributeError):
            man.top

    def test_molecules(self):
        mol1 = Molecule(VMD.molecule.new('Mine'))
        mol2 = Molecule(VMD.molecule.new('Other'))
        VMD.molecule.set_top(mol1.molid)
        man = MoleculeManager()

        self.assertEqual(len(man), 2)
        self.assertEqual(man[mol1.molid], mol1)
        self.assertEqual(man['Mine'], mol1)
        self.assertEqual(man[mol2.molid], mol2)
        self.assertEqual(man['Other'], mol2)
        self.assertEqual(list(man), [mol1, mol2])
        self.assertIn(mol1, man)
        self.assertIn(mol2, man)
        # Check top property
        self.assertEqual(man.top, mol1)
        man.top = mol2
        self.assertEqual(VMD.molecule.get_top(), mol2.molid)
        self.assertEqual(man.top, mol2)

        # Try deletion - using name
        del man['Mine']
        self.assertFalse(VMD.molecule.exists(mol1.molid))
        # Check the other attributes
        self.assertEqual(len(man), 1)
        with self.assertRaises(ValueError):
            man[mol1.molid]
        with self.assertRaises(ValueError):
            man['Mine']
        self.assertEqual(man[mol2.molid], mol2)
        self.assertEqual(man['Other'], mol2)
        self.assertEqual(list(man), [mol2])
        self.assertNotIn(mol1, man)
        self.assertIn(mol2, man)

        # Try deletion - using molid
        del man[mol2.molid]
        self.assertFalse(VMD.molecule.exists(mol2.molid))
        # Check the other attributes
        self.assertEqual(len(man), 0)
        with self.assertRaises(ValueError):
            man[mol2.molid]
        with self.assertRaises(ValueError):
            man['Other']
        self.assertEqual(list(man), [])
        self.assertNotIn(mol1, man)
        self.assertNotIn(mol2, man)

        # Check second deletions raises ValueError
        with self.assertRaises(ValueError):
            del man[mol1.molid]
        with self.assertRaises(ValueError):
            del man[mol2.molid]
        with self.assertRaises(ValueError):
            del man['Mine']
        with self.assertRaises(ValueError):
            del man['Other']

    def test_remote_actions(self):
        # Test manager can cope with molecule changes which happen without its knowledge
        man = MoleculeManager()

        # Create a molecule
        mol1 = Molecule(VMD.molecule.new('Unique'))

        # Manager finds it exists and adds it into name cache
        self.assertEqual(man['Unique'], mol1)

        # Delete the molecule and create other with the same name
        VMD.molecule.delete(mol1.molid)
        mol2 = Molecule(VMD.molecule.new('Unique'))

        # Manager correctly returns the new molecule
        self.assertEqual(man['Unique'], mol2)


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
        self.assertEqual(atom.molecule, mol)
        self.assertEqual(atom.frame, NOW)
        self.assertAlmostEqual(atom.x, -1.493)
        self.assertAlmostEqual(atom.y, 1.900)
        self.assertAlmostEqual(atom.z, 1.280)
        self.assertIsInstance(atom.coords, ndarray)
        self.assertAlmostEqualSeqs(list(atom.coords), [-1.493, 1.9, 1.28])
        self.assertEqual(list(atom.bonded), [Atom(1), Atom(2)])
        # Test setters
        atom.x = 23.9
        atom.y = -200.45
        atom.z = 0
        sel = VMD.atomsel.atomsel('index 0', molid=molid)
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

    def test_atom_comparison(self):
        # Test molecule comparison
        VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        atom1 = Atom(0)
        atom2 = Atom(0)
        other = Atom(1)

        self.assertEqual(atom1, atom1)
        self.assertEqual(atom1, atom2)
        self.assertNotEqual(atom1, other)
        self.assertNotEqual(atom2, other)

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
