"""
Tests for objects.
"""
import unittest

from Molecule import Molecule as _Molecule
import VMD

from pyvmd.objects import Molecule, FORMAT_PDB
from .utils import data


class TestMolecule(unittest.TestCase):
    """
    Test `Molecule` class.
    """
    def setUp(self):
        molid = VMD.molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        VMD.molecule.read(molid, 'dcd', data('water.1.dcd'), waitfor=-1)
        self.molid = molid

    def tearDown(self):
        # Delete all molecules
        for molid in VMD.molecule.listall():
            VMD.molecule.delete(molid)

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
            data = []
            for frame in xrange(VMD.molecule.numframes(self.molid)):
                sel.frame = frame
                data.append(sel.get('x')[0])
            return data

        # Check number of frames and iterator
        mol = Molecule(self.molid)
        self.assertEqual(len(mol.frames), 13)
        self.assertEqual(list(mol.frames), range(13))

        # Test frame duplication - duplicate focused frame
        mol.frame = 3
        mol.frames.copy()

        # Check there is one more frame
        self.assertEqual(len(mol.frames), 14)
        coords = [-1.4930000305175781, -1.4911566972732544, -1.4851371049880981, -1.4858486652374268,
                  -1.477394700050354, -1.4746015071868896, -1.46733820438385, -1.4535547494888306, -1.4307434558868408,
                  -1.4120502471923828, -1.385347843170166, -1.3674825429916382, -1.342192530632019, -1.4858486652374268]
        self.assertEqual(_get_x_coord(), coords)
        # Check molecule is focused to the new frame
        self.assertEqual(mol.frame, 13)

        # Test frame duplication - duplicate defined frame
        mol.frames.copy(4)
        # Check there is one more frame
        self.assertEqual(len(mol.frames), 15)
        coords = [-1.4930000305175781, -1.4911566972732544, -1.4851371049880981, -1.4858486652374268,
                  -1.477394700050354, -1.4746015071868896, -1.46733820438385, -1.4535547494888306, -1.4307434558868408,
                  -1.4120502471923828, -1.385347843170166, -1.3674825429916382, -1.342192530632019, -1.4858486652374268,
                  -1.477394700050354]
        self.assertEqual(_get_x_coord(), coords)
        # Molecule is focused to the new frame
        self.assertEqual(mol.frame, 14)

        # Test frame deletion - positive index
        del mol.frames[14]
        self.assertEqual(len(mol.frames), 14)
        coords = [-1.4930000305175781, -1.4911566972732544, -1.4851371049880981, -1.4858486652374268,
                  -1.477394700050354, -1.4746015071868896, -1.46733820438385, -1.4535547494888306, -1.4307434558868408,
                  -1.4120502471923828, -1.385347843170166, -1.3674825429916382, -1.342192530632019, -1.4858486652374268]
        self.assertEqual(_get_x_coord(), coords)

        # Test frame deletion - negative index
        del mol.frames[-2]
        self.assertEqual(len(mol.frames), 13)
        coords = [-1.4930000305175781, -1.4911566972732544, -1.4851371049880981, -1.4858486652374268,
                  -1.477394700050354, -1.4746015071868896, -1.46733820438385, -1.4535547494888306, -1.4307434558868408,
                  -1.4120502471923828, -1.385347843170166, -1.3674825429916382, -1.4858486652374268]
        self.assertEqual(_get_x_coord(), coords)

        # Check deletion of frame slice
        del mol.frames[2:11:3]
        self.assertEqual(len(mol.frames), 10)
        coords = [-1.4930000305175781, -1.4911566972732544, -1.4858486652374268, -1.477394700050354, -1.46733820438385,
                  -1.4535547494888306, -1.4120502471923828, -1.385347843170166, -1.3674825429916382,
                  -1.4858486652374268]
        self.assertEqual(_get_x_coord(), coords)

        # Invalid key for slice
        with self.assertRaises(TypeError):
            del mol.frames[None]
