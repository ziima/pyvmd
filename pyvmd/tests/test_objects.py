"""
Tests for objects.
"""
import unittest

from Molecule import Molecule as _Molecule
import VMD

from pyvmd.objects import Molecule
from .utils import data


class TestMolecule(unittest.TestCase):
    """
    Test `Molecule` class.
    """
    def setUp(self):
        _mol = _Molecule()
        _mol.load(data('water.psf'))
        _mol.load(data('water.pdb'))
        _mol.load(data('water.1.dcd'))
        self.molid = _mol.id

    def tearDown(self):
        VMD.molecule.delete(self.molid)

    def test_molecule(self):
        # Test basic features of Molecule
        mol = Molecule(self.molid)

        # Check active frame
        self.assertEqual(mol.frame, 12)
        # Check change of active frame works
        mol.frame = 3
        self.assertEqual(VMD.molecule.get_frame(self.molid), 3)
        self.assertEqual(mol.frame, 3)

        # Check comparison
        same = Molecule(self.molid)
        self.assertEqual(mol, same)

        _other = _Molecule()
        other = Molecule(_other.id)
        self.assertNotEqual(mol, other)

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
