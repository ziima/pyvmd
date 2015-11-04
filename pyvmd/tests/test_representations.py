"""
Tests for molecule representations.
"""
from VMD import molecule as _molecule, molrep as _molrep

from pyvmd.atoms import Selection
from pyvmd.molecules import Molecule
from pyvmd.representations import (Color, COLOR_COLOR, COLOR_NAME, COLOR_TYPE, DRAW_LINES, DRAW_POINTS, Representation,
                                   Style)

from .utils import data, PyvmdTestCase


class TestRepresentation(PyvmdTestCase):
    """
    Test `Representation` class.
    """
    def setUp(self):
        self.molid = _molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        self.other_molid = _molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        _molecule.set_top(self.molid)

    def test_init(self):
        # Test molecule is used if defined
        rep = Representation('rep0', Molecule(self.other_molid))
        self.assertEqual(rep.molecule, Molecule(self.other_molid))

    def test_init_top(self):
        # Test top molecule is used if molecule isn't defined
        rep = Representation('rep0')
        self.assertEqual(rep.molecule, Molecule(self.molid))

    def test_init_non_existing(self):
        # Test top molecule is used if molecule isn't defined
        with self.assertRaises(ValueError):
            Representation('JUNK')

    def test_comparison(self):
        rep_1 = Representation('rep0')
        rep_2 = Representation('rep0')
        rep_explicit = Representation('rep0', Molecule(self.molid))
        rep_new = Representation.create()
        rep_mol2 = Representation('rep0', Molecule(self.other_molid))

        self.assertEqual(rep_1, rep_1)
        self.assertEqual(rep_1, rep_2)
        self.assertEqual(rep_1, rep_explicit)
        self.assertNotEqual(rep_1, rep_new)
        self.assertNotEqual(rep_1, rep_mol2)

    def test_create(self):
        # Test new representaion is created
        Representation.create(Molecule(self.other_molid))
        self.assertEqual(_molrep.num(self.other_molid), 2)

    def test_create_top(self):
        # Test new representation is created in top molecule
        Representation.create()
        self.assertEqual(_molrep.num(self.molid), 2)

    def test_delete(self):
        rep = Representation('rep0')
        rep.delete()
        self.assertEqual(_molrep.num(self.molid), 0)

    def test_name(self):
        rep = Representation('rep0')
        self.assertEqual(rep.name, 'rep0')
        with self.assertRaises(AttributeError):
            rep.name = 'new_name'

    def test_molecule(self):
        rep = Representation('rep0')
        self.assertEqual(rep.molecule, Molecule(self.molid))
        with self.assertRaises(AttributeError):
            rep.molecule = Molecule(self.molid)

    def test_selection(self):
        rep = Representation('rep0')
        self.assertEqual(rep.selection, Selection('all'))

        rep.selection = Selection('resid 1 to 4')
        self.assertEqual(_molrep.get_selection(self.molid, 0), 'resid 1 to 4')
        self.assertEqual(rep.selection, Selection('resid 1 to 4'))

    def test_visible(self):
        rep = Representation('rep0')
        self.assertTrue(rep.visible)

        rep.visible = False
        self.assertFalse(_molrep.get_visible(self.molid, 0))
        self.assertFalse(rep.visible)

    def test_update_selection(self):
        rep = Representation('rep0')
        self.assertFalse(rep.update_selection)

        rep.update_selection = True
        self.assertTrue(_molrep.get_autoupdate(self.molid, 0))
        self.assertTrue(rep.update_selection)

    def test_update_color(self):
        rep = Representation('rep0')
        self.assertFalse(rep.update_color)

        rep.update_color = True
        self.assertTrue(_molrep.get_colorupdate(self.molid, 0))
        self.assertTrue(rep.update_color)

    def test_style(self):
        rep = Representation('rep0')
        self.assertEqual(rep.style, Style(Representation('rep0')))
        with self.assertRaises(AttributeError):
            rep.style = Style(Representation('rep0'))

    def test_color(self):
        rep = Representation('rep0')
        self.assertEqual(rep.color, Color(Representation('rep0')))
        with self.assertRaises(AttributeError):
            rep.style = Color(Representation('rep0'))


class TestStyle(PyvmdTestCase):
    """
    Test `Style` class.
    """
    def setUp(self):
        self.molid = _molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        self.rep = Representation('rep0')
        self.rep2 = Representation.create()

    def test_comparison(self):
        style_1 = Style(self.rep)
        style_2 = Style(self.rep)
        style_rep2 = Style(self.rep2)

        self.assertEqual(style_1, style_1)
        self.assertEqual(style_1, style_2)
        self.assertNotEqual(style_1, style_rep2)

    def test_method(self):
        style = Style(self.rep)
        self.assertEqual(style.method, DRAW_LINES)

        style.method = DRAW_POINTS
        self.assertEqual(style.method, DRAW_POINTS)
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Points')

    def test_get_parameters(self):
        style = Style(self.rep)
        self.assertEqual(style.get_parameters(), {'size': 1.})

    def test_set_parameters(self):
        style = Style(self.rep)
        style.set_parameters(size=4)
        self.assertEqual(_molrep.get_style(self.molid, 0), 'Lines 4')
        self.assertEqual(style.get_parameters(), {'size': 4})

    def test_set_parameters_no_kwargs(self):
        style = Style(self.rep)
        with self.assertRaises(ValueError):
            style.set_parameters()

    def test_set_parameters_invalid_kwargs(self):
        style = Style(self.rep)
        with self.assertRaises(ValueError):
            style.set_parameters(unknown=78)


class TestColor(PyvmdTestCase):
    """
    Test `Color` class.
    """
    def setUp(self):
        self.molid = _molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        self.rep = Representation('rep0')
        self.rep2 = Representation.create()

    def test_comparison(self):
        color_1 = Color(self.rep)
        color_2 = Color(self.rep)
        color_rep2 = Color(self.rep2)

        self.assertEqual(color_1, color_1)
        self.assertEqual(color_1, color_2)
        self.assertNotEqual(color_1, color_rep2)

    def test_method(self):
        color = Color(self.rep)
        self.assertEqual(color.method, COLOR_NAME)

        color.method = COLOR_TYPE
        self.assertEqual(color.method, COLOR_TYPE)
        self.assertEqual(_molrep.get_color(self.molid, 0), 'Type')

    def test_get_parameters(self):
        color = Color(self.rep)
        self.assertEqual(color.get_parameters(), {})

        color.method = COLOR_COLOR
        self.assertEqual(color.get_parameters(), {'color': 1})

    def test_set_parameters(self):
        color = Color(self.rep)
        color.method = COLOR_COLOR
        color.set_parameters(color=4)
        self.assertEqual(_molrep.get_color(self.molid, 0), 'ColorID 4')
        self.assertEqual(color.get_parameters(), {'color': 4})

    def test_set_parameters_no_kwargs(self):
        color = Color(self.rep)
        with self.assertRaises(ValueError):
            color.set_parameters()

    def test_set_parameters_invalid_kwargs(self):
        color = Color(self.rep)
        with self.assertRaises(ValueError):
            color.set_parameters(unknown=78)
