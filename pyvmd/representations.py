"""
Objects for manipulation with VMD molecule representations.
"""
from collections import namedtuple

from VMD import molrep as _molrep

from pyvmd.atoms import Selection
from pyvmd.molecules import Molecule, MOLECULES

__all__ = ['COLOR_BACKBONE', 'COLOR_BETA', 'COLOR_CHAIN', 'COLOR_CHARGE', 'COLOR_COLOR', 'COLOR_CONFORMATION',
           'COLOR_ELEMENT', 'COLOR_FRAGMENT', 'COLOR_INDEX', 'COLOR_MASS', 'COLOR_MOLECULE', 'COLOR_NAME',
           'COLOR_OCCUPANCY', 'COLOR_PHYSICAL_TIME', 'COLOR_POS', 'COLOR_POS_X', 'COLOR_POS_Y', 'COLOR_POS_Z',
           'COLOR_RESID', 'COLOR_RESNAME', 'COLOR_RESTYPE', 'COLOR_SEGNAME', 'COLOR_STRUCTURE', 'COLOR_THROB',
           'COLOR_TIMESTEP', 'COLOR_TYPE', 'COLOR_USER', 'COLOR_USER_2', 'COLOR_USER_3', 'COLOR_USER_4',
           'COLOR_VELOCITY', 'COLOR_VOLUME', 'COLORING_METHODS', 'DRAW_BEADS', 'DRAW_BONDS', 'DRAW_CARTOON', 'DRAW_CPK',
           'DRAW_DOTTED', 'DRAW_DYNAMIC_BONDS', 'DRAW_FIELD_LINES', 'DRAW_HBONDS', 'DRAW_ISOSURFACE', 'DRAW_LICORICE',
           'DRAW_LINES', 'DRAW_MSMS', 'DRAW_NEW_CARTOON', 'DRAW_NEW_RIBBONS', 'DRAW_ORBITAL', 'DRAW_PAPER_CHAIN',
           'DRAW_POINTS', 'DRAW_POLYHEDRA', 'DRAW_QUICKSURF', 'DRAW_RIBBONS', 'DRAW_SOLVENT', 'DRAW_SURF', 'DRAW_TRACE',
           'DRAW_TUBE', 'DRAW_TWISTER', 'DRAW_VDW', 'DRAW_VOLUME_SLICE', 'DRAWING_METHODS', 'Representation']


def _modrep(representation, **kwargs):
    """
    Utility wrapper over `modrep` to raise exception in case of an error.
    """
    result = _molrep.modrep(representation.molecule.molid, representation.repindex, **kwargs)
    if not result:
        raise ValueError("Error occured while changing molecule representation.")


################################################################################
# Drawing methods
#
# XXX: There is no way to get default parameters for styles from VMD, so we have to define them ourselves.
# name: VMD style keyword
# parameters: List of paramater names. Customized because there is no definition of those in VMD.
# defaults: Dictionary with default values. Values are copied from AtomRep.C.
DrawingMethod = namedtuple('DrawingMethod', ('name', 'parameters', 'defaults'))
# TODO: Even though VMD handles all parameters as floats, some of them actually behaves differently, e.g. as a bool.
# High level interface implemented here should honor the actual behavior of the parameters and interpret float values
# correctly.

# Drawing methods as defined in AtomRep.C
DRAW_LINES = DrawingMethod('Lines', ('size', ), {'size': 1.})
DRAW_BONDS = DrawingMethod('Bonds', ('size', 'resolution'), {'size': 0.3, 'resolution': 10})
DRAW_DYNAMIC_BONDS = DrawingMethod('DynamicBonds', ('cutoff', 'size', 'resolution'),
                                   {'cutoff': 3., 'size': 0.3, 'resolution': 10})
DRAW_HBONDS = DrawingMethod('HBonds', ('cutoff', 'angle_cutoff', 'size'),
                            {'cutoff': 3., 'angle_cutoff': 20., 'size': 1.})
DRAW_POINTS = DrawingMethod('Points', ('size', ), {'size': 1.})
DRAW_VDW = DrawingMethod('VDW', ('size', 'resolution'), {'size': 1., 'resolution': 12})
DRAW_CPK = DrawingMethod('CPK', ('size', 'bond_radius', 'resolution', 'bond_resolution'),
                         {'size': 1., 'bond_radius': 0.3, 'resolution': 12, 'bond_resolution': 10})
DRAW_LICORICE = DrawingMethod('Licorice', ('size', 'resolution', 'bond_resolution'),
                              {'size': 0.3, 'resolution': 10, 'bond_resolution': 10})
DRAW_POLYHEDRA = DrawingMethod('Polyhedra', ('cutoff', ), {'cutoff': 3.})
DRAW_TRACE = DrawingMethod('Trace', ('size', 'resolution'), {'size': 0.3, 'resolution': 10})
DRAW_TUBE = DrawingMethod('Tube', ('size', 'resolution'), {'size': 0.3, 'resolution': 10})
DRAW_RIBBONS = DrawingMethod('Ribbons', ('margin_size', 'resolution', 'size'),
                             {'margin_size': 0.3, 'resolution': 10, 'size': 2.})
# TODO: spline constants
DRAW_NEW_RIBBONS = DrawingMethod('NewRibbons', ('thickness', 'resolution', 'width', 'spline'),
                                 {'thickness': 0.3, 'resolution': 10, 'width': 3., 'spline': 0})
DRAW_CARTOON = DrawingMethod('Cartoon', ('helix_size', 'resolution', 'sheet_size'),
                             {'helix_size': 2.1, 'resolution': 12, 'sheet_size': 5.})
# TODO: spline constants
DRAW_NEW_CARTOON = DrawingMethod('NewCartoon', ('thickness', 'resolution', 'width', 'spline'),
                                 {'thickness': 0.3, 'resolution': 10, 'width': 4.5, 'spline': 0})
DRAW_PAPER_CHAIN = DrawingMethod('PaperChain', ('height', 'max_ring_size'), {'height': 1., 'max_ring_size': 10})
# TODO: start, hide_shared_links constants
DRAW_TWISTER = DrawingMethod(
    'Twister',
    ('start', 'hide_shared_links', 'steps', 'width', 'height', 'max_ring_size', 'linking_distance'),
    {'start': 1, 'hide_shared_links': 0, 'steps': 10, 'width': 0.3, 'height': 0.05, 'max_ring_size': 10,
     'linking_distance': 5},
)
DRAW_QUICKSURF = DrawingMethod('QuickSurf', ('sphere_scale', 'density_isovalue', 'grid_spacing', 'resolution'),
                               {'sphere_scale': 1., 'density_isovalue': 0.5, 'grid_spacing': 1., 'resolution': 0})
# TODO: selection, method constants
DRAW_MSMS = DrawingMethod('MSMS', ('probe_radius', 'density', 'atoms', 'method'),
                          {'probe_radius': 1.5, 'density': 1.5, 'atoms': 0, 'method': 0})
# TODO: method constant
DRAW_SURF = DrawingMethod('Surf', ('probe_radius', 'method'), {'probe_radius': 1.4, 'method': 0})
# TODO: axis, quality constants
DRAW_VOLUME_SLICE = DrawingMethod('VolumeSlice', ('slice', 'volume_id', 'axis', 'quality'),
                                  {'slice': 0.5, 'volume_id': 0, 'axis': 0, 'quality': 2})
# TODO: display, method constants
DRAW_ISOSURFACE = DrawingMethod('Isosurface', ('isovalue', 'volume_id', 'display', 'method', 'step', 'size'),
                                {'isovalue': 0.5, 'volume_id': 0, 'display': 2, 'method': 2, 'step': 1, 'size': 1})
DRAW_FIELD_LINES = DrawingMethod('FieldLines', ('volume_id', 'gradient', 'min_length', 'max_length', 'size'),
                                 {'volume_id': 0, 'gradient': 1.8, 'min_length': 10, 'max_length': 50, 'size': 1})
# TODO: display, method, wavefunction, spin constants
DRAW_ORBITAL = DrawingMethod(
    'Orbital',
    ('isovalue', 'orbital_id', 'display', 'method', 'grid_spacing', 'size', 'wavefunction', 'spin', 'excitation',
     'step'),
    {'isovalue': 0.05, 'orbital_id': 0, 'display': 0, 'method': 0, 'grid_spacing': 0.075, 'size': 1, 'wavefunction': 0,
     'spin': 0, 'excitation': 0, 'step': 1}
)
DRAW_BEADS = DrawingMethod('Beads', ('size', 'resolution'), {'size': 1., 'resolution': 12})
DRAW_DOTTED = DrawingMethod('Dotted', ('size', 'resolution'), {'size': 1., 'resolution': 12})
# TODO: method constants
DRAW_SOLVENT = DrawingMethod('Solvent', ('probe_radius', 'detail', 'method'),
                             {'probe_radius': 0, 'detail': 7, 'method': 1})

DRAWING_METHODS = (DRAW_LINES, DRAW_BONDS, DRAW_DYNAMIC_BONDS, DRAW_HBONDS, DRAW_POINTS, DRAW_VDW, DRAW_CPK,
                   DRAW_LICORICE, DRAW_POLYHEDRA, DRAW_TRACE, DRAW_TUBE, DRAW_RIBBONS, DRAW_NEW_RIBBONS, DRAW_CARTOON,
                   DRAW_NEW_CARTOON, DRAW_PAPER_CHAIN, DRAW_TWISTER, DRAW_QUICKSURF, DRAW_MSMS, DRAW_SURF,
                   DRAW_VOLUME_SLICE, DRAW_ISOSURFACE, DRAW_FIELD_LINES, DRAW_ORBITAL, DRAW_BEADS, DRAW_DOTTED,
                   DRAW_SOLVENT)
DRAWING_METHODS_MAP = {m.name: m for m in DRAWING_METHODS}


def _parse_raw_style(raw_style):
    """
    Parses raw style string and returns drawing method and options.
    """
    assert raw_style
    chunks = raw_style.split()
    method_name = chunks.pop(0)
    if method_name not in DRAWING_METHODS_MAP:
        raise ValueError('Unknown style: %s' % raw_style)
    method = DRAWING_METHODS_MAP[method_name]
    # If parameters are set to default, some of them may be missing from `raw_style`. Return complete set of parameters.
    params = method.defaults.copy()
    params.update({p: float(v) for p, v in zip(method.parameters, chunks)})
    return method, params


class Style(object):
    """
    Represents drawing method and its parameters of a particular representaiton.

    This class is a proxy for molecule represenation style in VMD.
    """
    def __init__(self, representation):
        assert isinstance(representation, Representation)
        self.representation = representation

    def __repr__(self):
        return "<%s: '%r'>" % (type(self).__name__, self.representation)

    def __eq__(self, other):
        return type(self) == type(other) and self.representation == other.representation

    def __ne__(self, other):
        return not self.__eq__(other)

    ############################################################################
    # Style properties
    def _get_method(self):
        raw_style = _molrep.get_style(self.representation.molecule.molid, self.representation.repindex)
        return _parse_raw_style(raw_style)[0]

    def _set_method(self, method):
        # Only assert. There is no function to get the list of available drawing methods, so just send it to VMD.
        assert method in DRAWING_METHODS
        _modrep(self.representation, style=method.name)

    method = property(_get_method, _set_method, doc="Drawing method")

    def get_parameters(self):
        """Returns drawing parameters"""
        raw_style = _molrep.get_style(self.representation.molecule.molid, self.representation.repindex)
        return _parse_raw_style(raw_style)[1]

    def set_parameters(self, **kwargs):
        """Sets drawing parameters"""
        if not kwargs:
            raise ValueError('At least one parameter is required.')
        method = self.method
        if set(kwargs) - set(method.parameters):
            raise ValueError('Unknown parameters: %s' % (set(kwargs) - set(method.parameters)))
        # Get current parameters and update with modifications
        params = self.get_parameters()
        params.update(kwargs)
        # Construct style string
        data = (str(params[p]) for p in method.parameters)
        style_str = '%s %s' % (method.name, ' '.join(data))
        # Set the style in VMD
        _modrep(self.representation, style=style_str)


################################################################################
# Color methods
#
# name: VMD color keyword
# parameters: List of paramater names. Customized because there is no definition of those in VMD.
# defaults: Dictionary with default values.
ColoringMethod = namedtuple('ColoringMethod', ('name', 'parameters', 'defaults'))

# Color methods as defined in AtomColor.C
COLOR_NAME = ColoringMethod("Name", (), {})
COLOR_TYPE = ColoringMethod("Type", (), {})
COLOR_ELEMENT = ColoringMethod("Element", (), {})
COLOR_RESNAME = ColoringMethod("ResName", (), {})
COLOR_RESTYPE = ColoringMethod("ResType", (), {})
COLOR_RESID = ColoringMethod("ResID", (), {})
COLOR_CHAIN = ColoringMethod("Chain", (), {})
COLOR_SEGNAME = ColoringMethod("SegName", (), {})
COLOR_CONFORMATION = ColoringMethod("Conformation", (), {})
COLOR_MOLECULE = ColoringMethod("Molecule", (), {})
COLOR_STRUCTURE = ColoringMethod("Structure", (), {})
COLOR_COLOR = ColoringMethod("ColorID", ('color', ), {'color': 1})
COLOR_BETA = ColoringMethod("Beta", (), {})
COLOR_OCCUPANCY = ColoringMethod("Occupancy", (), {})
COLOR_MASS = ColoringMethod("Mass", (), {})
COLOR_CHARGE = ColoringMethod("Charge", (), {})
COLOR_POS = ColoringMethod("Pos", (), {})
COLOR_POS_X = ColoringMethod("PosX", (), {})
COLOR_POS_Y = ColoringMethod("PosY", (), {})
COLOR_POS_Z = ColoringMethod("PosZ", (), {})
COLOR_USER = ColoringMethod("User", (), {})
COLOR_USER_2 = ColoringMethod("User2", (), {})
COLOR_USER_3 = ColoringMethod("User3", (), {})
COLOR_USER_4 = ColoringMethod("User4", (), {})
COLOR_FRAGMENT = ColoringMethod("Fragment", (), {})
COLOR_INDEX = ColoringMethod("Index", (), {})
COLOR_BACKBONE = ColoringMethod("Backbone", (), {})
COLOR_THROB = ColoringMethod("Throb", (), {})
COLOR_PHYSICAL_TIME = ColoringMethod("PhysicalTime", (), {})
COLOR_TIMESTEP = ColoringMethod("Timestep", (), {})
COLOR_VELOCITY = ColoringMethod("Velocity", (), {})
COLOR_VOLUME = ColoringMethod("Volume", ('volume_id', ), {'volume_id': 0})

COLORING_METHODS = (COLOR_NAME, COLOR_TYPE, COLOR_ELEMENT, COLOR_RESNAME, COLOR_RESTYPE, COLOR_RESID, COLOR_CHAIN,
                    COLOR_SEGNAME, COLOR_CONFORMATION, COLOR_MOLECULE, COLOR_STRUCTURE, COLOR_COLOR, COLOR_BETA,
                    COLOR_OCCUPANCY, COLOR_MASS, COLOR_CHARGE, COLOR_POS, COLOR_POS_X, COLOR_POS_Y, COLOR_POS_Z,
                    COLOR_USER, COLOR_USER_2, COLOR_USER_3, COLOR_USER_4, COLOR_FRAGMENT, COLOR_INDEX, COLOR_BACKBONE,
                    COLOR_THROB, COLOR_PHYSICAL_TIME, COLOR_TIMESTEP, COLOR_VELOCITY, COLOR_VOLUME)
COLORING_METHODS_MAP = {m.name: m for m in COLORING_METHODS}


def _parse_raw_color(raw_color):
    """
    Parses raw color string and returns coloring method and options.
    """
    assert raw_color
    chunks = raw_color.split()
    method_name = chunks.pop(0)
    if method_name not in COLORING_METHODS_MAP:
        raise ValueError('Unknown color: %s' % raw_color)
    method = COLORING_METHODS_MAP[method_name]
    params = {p: int(v) for p, v in zip(method.parameters, chunks)}
    return method, params


class Color(object):
    """
    Represents coloring method and its parameters of a particular representaiton.

    This class is a proxy for molecule represenation color in VMD.
    """
    def __init__(self, representation):
        assert isinstance(representation, Representation)
        self.representation = representation

    def __repr__(self):
        return "<%s: '%r'>" % (type(self).__name__, self.representation)

    def __eq__(self, other):
        return type(self) == type(other) and self.representation == other.representation

    def __ne__(self, other):
        return not self.__eq__(other)

    ############################################################################
    # Color properties
    def _set_color(self, method, parameters):
        # Utility method to format and set the color with all parameters.
        # Construct color string
        if method.parameters:
            data = (str(parameters[p]) for p in method.parameters)
            color_str = '%s %s' % (method.name, ' '.join(data))
        else:
            color_str = method.name
        # Set the style in VMD
        _modrep(self.representation, color=color_str)

    def _get_method(self):
        raw_color = _molrep.get_color(self.representation.molecule.molid, self.representation.repindex)
        return _parse_raw_color(raw_color)[0]

    def _set_method(self, method):
        # Only assert. There is no function to get the list of available coloring methods, so just send it to VMD.
        assert method in COLORING_METHODS
        self._set_color(method, method.defaults)

    method = property(_get_method, _set_method, doc="Coloring method")

    def get_parameters(self):
        """Returns coloring parameters"""
        raw_color = _molrep.get_color(self.representation.molecule.molid, self.representation.repindex)
        return _parse_raw_color(raw_color)[1]

    def set_parameters(self, **kwargs):
        """Sets coloring parameters"""
        if not kwargs:
            raise ValueError('At least one parameter is required.')
        method = self.method
        if set(kwargs) - set(method.parameters):
            raise ValueError('Unknown parameters: %s' % (set(kwargs) - set(method.parameters)))
        # Get current parameters and update with modifications
        params = self.get_parameters()
        params.update(kwargs)
        self._set_color(method, params)


################################################################################
# Representation
#
class Representation(object):
    """
    Representation of molecule graphical representation.

    This class is a proxy for molecule represenation in VMD.
    """
    def __init__(self, name, molecule=None):
        """
        Creates a molecule representation.

        @param name: Representation identifier.
        @type name: str
        @param molecule: Molecule to select from. Top if not provided.
        @type molecule: Molecule or None
        """
        if molecule is None:
            molecule = MOLECULES.top
        else:
            assert isinstance(molecule, Molecule)
        self._molecule = molecule

        assert isinstance(name, str)
        self._name = name
        # If the representation does not exists, complain immediately
        self.repindex

    def __repr__(self):
        return "<%s: '%s' of '%r'>" % (type(self).__name__, self._name, self._molecule)

    def __eq__(self, other):
        return type(self) == type(other) and self._molecule == other.molecule and self._name == other.name

    def __ne__(self, other):
        return not self.__eq__(other)

    @property
    def name(self):
        """Index"""
        return self._name

    @property
    def molecule(self):
        """Molecule"""
        return self._molecule

    @property
    def repindex(self):
        """Utility property to facilitate the VMD calls."""
        repindex = _molrep.repindex(self._molecule.molid, self._name)
        if repindex < 0:
            raise ValueError("Molecule %s does not have a representation '%s'." % (self._molecule, self._name))
        return repindex

    @classmethod
    def create(cls, molecule=None):
        """
        Creates new molecule representation.

        @param molecule: Molecule to select from. Top if not provided.
        @type molecule: Molecule or None
        """
        if molecule is None:
            molecule = MOLECULES.top
        else:
            assert isinstance(molecule, Molecule)
        _molrep.addrep(molecule.molid)
        # XXX: `addrep` has no return value, get the index of the newly created representation from their number.
        index = _molrep.num(molecule.molid) - 1
        name = _molrep.get_repname(molecule.molid, index)

        return cls(name, molecule)

    def delete(self):
        """
        Deletes the molecule representation.
        """
        _molrep.delrep(self._molecule.molid, self.repindex)

    ############################################################################
    # Representation properties
    def _get_selection(self):
        return Selection(_molrep.get_selection(self._molecule.molid, self.repindex), self._molecule)

    def _set_selection(self, selection):
        assert isinstance(selection, Selection)
        assert self._molecule == selection.molecule
        _modrep(self, sel=selection.selection)

    selection = property(_get_selection, _set_selection, doc="Selection")

    def _get_visible(self):
        return _molrep.get_visible(self._molecule.molid, self.repindex)

    def _set_visible(self, value):
        assert isinstance(value, bool)
        _molrep.set_visible(self._molecule.molid, self.repindex, value)

    visible = property(_get_visible, _set_visible, doc="Visibility")

    def _get_update_selection(self):
        return _molrep.get_autoupdate(self._molecule.molid, self.repindex)

    def _set_update_selection(self, value):
        assert isinstance(value, bool)
        _molrep.set_autoupdate(self._molecule.molid, self.repindex, value)

    update_selection = property(_get_update_selection, _set_update_selection, doc="Update selection every frame")

    def _get_update_color(self):
        return _molrep.get_colorupdate(self._molecule.molid, self.repindex)

    def _set_update_color(self, value):
        assert isinstance(value, bool)
        _molrep.set_colorupdate(self._molecule.molid, self.repindex, value)

    update_color = property(_get_update_color, _set_update_color, doc="Update color every frame")

    @property
    def style(self):
        """Graphical style"""
        return Style(self)

    @property
    def color(self):
        """Color"""
        return Color(self)
