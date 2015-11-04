"""
Tests for molecule utilities.
"""
import VMD
from Molecule import Molecule as _Molecule

from pyvmd.molecules import FORMAT_PDB, Molecule, MoleculeManager
from pyvmd.representations import Representation

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
        self.assertEqual(VMD.molecule.get_filenames(mol.molid),
                         [data('water.psf'), data('water.pdb'), data('water.1.dcd')])
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
        with self.assertRaises(ValueError):
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


class TestRepresentationManager(PyvmdTestCase):
    """
    Test `RepresentationManager` class.
    """
    def setUp(self):
        molid = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        self.mol = Molecule(molid)

    def test_property(self):
        # Test Molecule.representations property
        man = self.mol.representations

        self.assertEqual(man.molecule, self.mol)

        with self.assertRaises(AttributeError):
            self.mol.representations = Representation('rep0')

    def test_len(self):
        self.assertEqual(len(self.mol.representations), 1)
        VMD.molrep.addrep(self.mol.molid)
        self.assertEqual(len(self.mol.representations), 2)

    def test_getitem(self):
        VMD.molrep.addrep(self.mol.molid)

        # Test positive indexes
        self.assertEqual(self.mol.representations[0], Representation('rep0'))
        self.assertEqual(self.mol.representations[1], Representation('rep1'))
        with self.assertRaises(IndexError):
            self.mol.representations[2]
        # Check negative indexes
        self.assertEqual(self.mol.representations[-2], Representation('rep0'))
        self.assertEqual(self.mol.representations[-1], Representation('rep1'))
        with self.assertRaises(IndexError):
            self.mol.representations[-3]

        # Test representation names
        self.assertEqual(self.mol.representations['rep0'], Representation('rep0'))
        self.assertEqual(self.mol.representations['rep1'], Representation('rep1'))
        with self.assertRaises(KeyError):
            self.mol.representations['junk']

        # Test slices
        self.assertEqual(self.mol.representations[::-1], [Representation('rep1'), Representation('rep0')])
        self.assertEqual(self.mol.representations[:10], [Representation('rep0'), Representation('rep1')])

        # Test type error
        with self.assertRaises(TypeError):
            self.mol.representations[None]

    def test_iter(self):
        self.assertEqual(list(self.mol.representations), [Representation('rep0')])

        VMD.molrep.addrep(self.mol.molid)

        self.assertEqual(list(self.mol.representations), [Representation('rep0'), Representation('rep1')])
