"""
Tests for coloring methods of molecule representations.
"""
from VMD import molecule as _molecule, molrep as _molrep

from pyvmd.representations import (COLOR_BACKBONE, COLOR_BETA, COLOR_CHAIN, COLOR_CHARGE, COLOR_COLOR,
                                   COLOR_CONFORMATION, COLOR_ELEMENT, COLOR_FRAGMENT, COLOR_INDEX, COLOR_MASS,
                                   COLOR_MOLECULE, COLOR_NAME, COLOR_OCCUPANCY, COLOR_PHYSICAL_TIME, COLOR_POS,
                                   COLOR_POS_X, COLOR_POS_Y, COLOR_POS_Z, COLOR_RESID, COLOR_RESNAME, COLOR_RESTYPE,
                                   COLOR_SEGNAME, COLOR_STRUCTURE, COLOR_THROB, COLOR_TIMESTEP, COLOR_TYPE, COLOR_USER,
                                   COLOR_USER_2, COLOR_USER_3, COLOR_USER_4, COLOR_VELOCITY, COLOR_VOLUME,
                                   Representation)

from .utils import data, PyvmdTestCase


class TestColoringMethods(PyvmdTestCase):
    """
    Test coloring methods are defined correctly.
    """
    def setUp(self):
        self.molid = _molecule.load('psf', data('water.psf'), 'pdb', data('water.pdb'))
        self.rep = Representation('rep0')
        self.color = self.rep.color

    def test_name(self):
        self.color.method = COLOR_NAME
        self.assertEqual(self.color.method, COLOR_NAME)
        self.assertEqual(self.color.get_parameters(), {})
        self.assertEqual(_molrep.get_color(self.molid, 0), 'Name')

    def test_type(self):
        self.color.method = COLOR_TYPE
        self.assertEqual(self.color.method, COLOR_TYPE)
        self.assertEqual(self.color.get_parameters(), {})
        self.assertEqual(_molrep.get_color(self.molid, 0), 'Type')

    def test_element(self):
        self.color.method = COLOR_ELEMENT
        self.assertEqual(self.color.method, COLOR_ELEMENT)
        self.assertEqual(self.color.get_parameters(), {})
        self.assertEqual(_molrep.get_color(self.molid, 0), 'Element')

    def test_resname(self):
        self.color.method = COLOR_RESNAME
        self.assertEqual(self.color.method, COLOR_RESNAME)
        self.assertEqual(self.color.get_parameters(), {})
        self.assertEqual(_molrep.get_color(self.molid, 0), 'ResName')

    def test_restype(self):
        self.color.method = COLOR_RESTYPE
        self.assertEqual(self.color.method, COLOR_RESTYPE)
        self.assertEqual(self.color.get_parameters(), {})
        self.assertEqual(_molrep.get_color(self.molid, 0), 'ResType')

    def test_resid(self):
        self.color.method = COLOR_RESID
        self.assertEqual(self.color.method, COLOR_RESID)
        self.assertEqual(self.color.get_parameters(), {})
        self.assertEqual(_molrep.get_color(self.molid, 0), 'ResID')

    def test_chain(self):
        self.color.method = COLOR_CHAIN
        self.assertEqual(self.color.method, COLOR_CHAIN)
        self.assertEqual(self.color.get_parameters(), {})
        self.assertEqual(_molrep.get_color(self.molid, 0), 'Chain')

    def test_segname(self):
        self.color.method = COLOR_SEGNAME
        self.assertEqual(self.color.method, COLOR_SEGNAME)
        self.assertEqual(self.color.get_parameters(), {})
        self.assertEqual(_molrep.get_color(self.molid, 0), 'SegName')

    def test_conformation(self):
        self.color.method = COLOR_CONFORMATION
        self.assertEqual(self.color.method, COLOR_CONFORMATION)
        self.assertEqual(self.color.get_parameters(), {})
        self.assertEqual(_molrep.get_color(self.molid, 0), 'Conformation')

    def test_molecule(self):
        self.color.method = COLOR_MOLECULE
        self.assertEqual(self.color.method, COLOR_MOLECULE)
        self.assertEqual(self.color.get_parameters(), {})
        self.assertEqual(_molrep.get_color(self.molid, 0), 'Molecule')

    def test_structure(self):
        self.color.method = COLOR_STRUCTURE
        self.assertEqual(self.color.method, COLOR_STRUCTURE)
        self.assertEqual(self.color.get_parameters(), {})
        self.assertEqual(_molrep.get_color(self.molid, 0), 'Structure')

    def test_color(self):
        self.color.method = COLOR_COLOR
        self.assertEqual(self.color.method, COLOR_COLOR)
        self.assertEqual(self.color.get_parameters(), {'color': 1})
        self.assertEqual(_molrep.get_color(self.molid, 0), 'ColorID 1')

        self.color.set_parameters(color=3)
        self.assertEqual(self.color.method, COLOR_COLOR)
        self.assertEqual(self.color.get_parameters(), {'color': 3})
        self.assertEqual(_molrep.get_color(self.molid, 0), 'ColorID 3')

    def test_beta(self):
        self.color.method = COLOR_BETA
        self.assertEqual(self.color.method, COLOR_BETA)
        self.assertEqual(self.color.get_parameters(), {})
        self.assertEqual(_molrep.get_color(self.molid, 0), 'Beta')

    def test_occupancy(self):
        self.color.method = COLOR_OCCUPANCY
        self.assertEqual(self.color.method, COLOR_OCCUPANCY)
        self.assertEqual(self.color.get_parameters(), {})
        self.assertEqual(_molrep.get_color(self.molid, 0), 'Occupancy')

    def test_mass(self):
        self.color.method = COLOR_MASS
        self.assertEqual(self.color.method, COLOR_MASS)
        self.assertEqual(self.color.get_parameters(), {})
        self.assertEqual(_molrep.get_color(self.molid, 0), 'Mass')

    def test_charge(self):
        self.color.method = COLOR_CHARGE
        self.assertEqual(self.color.method, COLOR_CHARGE)
        self.assertEqual(self.color.get_parameters(), {})
        self.assertEqual(_molrep.get_color(self.molid, 0), 'Charge')

    def test_pos(self):
        self.color.method = COLOR_POS
        self.assertEqual(self.color.method, COLOR_POS)
        self.assertEqual(self.color.get_parameters(), {})
        self.assertEqual(_molrep.get_color(self.molid, 0), 'Pos')

    def test_pos_x(self):
        self.color.method = COLOR_POS_X
        self.assertEqual(self.color.method, COLOR_POS_X)
        self.assertEqual(self.color.get_parameters(), {})
        self.assertEqual(_molrep.get_color(self.molid, 0), 'PosX')

    def test_pos_y(self):
        self.color.method = COLOR_POS_Y
        self.assertEqual(self.color.method, COLOR_POS_Y)
        self.assertEqual(self.color.get_parameters(), {})
        self.assertEqual(_molrep.get_color(self.molid, 0), 'PosY')

    def test_pos_z(self):
        self.color.method = COLOR_POS_Z
        self.assertEqual(self.color.method, COLOR_POS_Z)
        self.assertEqual(self.color.get_parameters(), {})
        self.assertEqual(_molrep.get_color(self.molid, 0), 'PosZ')

    def test_user(self):
        self.color.method = COLOR_USER
        self.assertEqual(self.color.method, COLOR_USER)
        self.assertEqual(self.color.get_parameters(), {})
        self.assertEqual(_molrep.get_color(self.molid, 0), 'User')

    def test_user_2(self):
        self.color.method = COLOR_USER_2
        self.assertEqual(self.color.method, COLOR_USER_2)
        self.assertEqual(self.color.get_parameters(), {})
        self.assertEqual(_molrep.get_color(self.molid, 0), 'User2')

    def test_user_3(self):
        self.color.method = COLOR_USER_3
        self.assertEqual(self.color.method, COLOR_USER_3)
        self.assertEqual(self.color.get_parameters(), {})
        self.assertEqual(_molrep.get_color(self.molid, 0), 'User3')

    def test_user_4(self):
        self.color.method = COLOR_USER_4
        self.assertEqual(self.color.method, COLOR_USER_4)
        self.assertEqual(self.color.get_parameters(), {})
        self.assertEqual(_molrep.get_color(self.molid, 0), 'User4')

    def test_fragment(self):
        self.color.method = COLOR_FRAGMENT
        self.assertEqual(self.color.method, COLOR_FRAGMENT)
        self.assertEqual(self.color.get_parameters(), {})
        self.assertEqual(_molrep.get_color(self.molid, 0), 'Fragment')

    def test_index(self):
        self.color.method = COLOR_INDEX
        self.assertEqual(self.color.method, COLOR_INDEX)
        self.assertEqual(self.color.get_parameters(), {})
        self.assertEqual(_molrep.get_color(self.molid, 0), 'Index')

    def test_backbone(self):
        self.color.method = COLOR_BACKBONE
        self.assertEqual(self.color.method, COLOR_BACKBONE)
        self.assertEqual(self.color.get_parameters(), {})
        self.assertEqual(_molrep.get_color(self.molid, 0), 'Backbone')

    def test_throb(self):
        self.color.method = COLOR_THROB
        self.assertEqual(self.color.method, COLOR_THROB)
        self.assertEqual(self.color.get_parameters(), {})
        self.assertEqual(_molrep.get_color(self.molid, 0), 'Throb')

    def test_physical_time(self):
        self.color.method = COLOR_PHYSICAL_TIME
        self.assertEqual(self.color.method, COLOR_PHYSICAL_TIME)
        self.assertEqual(self.color.get_parameters(), {})
        self.assertEqual(_molrep.get_color(self.molid, 0), 'PhysicalTime')

    def test_timestep(self):
        self.color.method = COLOR_TIMESTEP
        self.assertEqual(self.color.method, COLOR_TIMESTEP)
        self.assertEqual(self.color.get_parameters(), {})
        self.assertEqual(_molrep.get_color(self.molid, 0), 'Timestep')

    def test_velocity(self):
        self.color.method = COLOR_VELOCITY
        self.assertEqual(self.color.method, COLOR_VELOCITY)
        self.assertEqual(self.color.get_parameters(), {})
        self.assertEqual(_molrep.get_color(self.molid, 0), 'Velocity')

    def test_volume(self):
        self.color.method = COLOR_VOLUME
        self.assertEqual(self.color.method, COLOR_VOLUME)
        self.assertEqual(self.color.get_parameters(), {'volume_id': 0})
        self.assertEqual(_molrep.get_color(self.molid, 0), 'Volume 0')

        self.color.set_parameters(volume_id=42)
        self.assertEqual(self.color.method, COLOR_VOLUME)
        self.assertEqual(self.color.get_parameters(), {'volume_id': 42})
        self.assertEqual(_molrep.get_color(self.molid, 0), 'Volume 42')
