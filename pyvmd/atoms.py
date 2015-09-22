"""
Objects for manipulation with VMD atoms.
"""
import itertools

from numpy import array
from atomsel import atomsel as _atomsel
from VMD import molecule as _molecule

from .molecules import Molecule, MOLECULES

__all__ = ['Atom', 'Chain', 'Residue', 'Segment' 'Selection', 'NOW']


# Constant which always references active frame
NOW = -1


class SelectionBase(object):
    """
    Base class for selection-like objects.
    """
    # Large amounts of these objects can be created, slots has some performance benefits.
    __slots__ = ('_molecule', '_frame', '_atomsel')

    def __init__(self, molecule=None, frame=NOW):
        """
        Creates selection-like object.

        @param molecule: Molecule to select from. Top if not provider.
        @type molecule: Molecule or None
        @param frame: Selection frame
        @type frame: Non-negative integer or NOW
        """
        if molecule is None:
            molecule = MOLECULES.top
        else:
            assert isinstance(molecule, Molecule)
        assert frame == NOW or (isinstance(frame, int) and frame >= 0)
        self._molecule = molecule
        self._frame = frame
        self._atomsel = None

    @property
    def molecule(self):
        "Molecule"
        return self._molecule

    def _get_frame(self):
        return self._frame

    def _set_frame(self, frame):
        assert frame == NOW or (isinstance(frame, int) and frame >= 0)
        self.atomsel.frame = frame
        self._frame = frame

    frame = property(_get_frame, _set_frame, doc="Frame")

    @property
    def atomsel(self):
        """
        Returns respective 'VMD.atomsel' instance. Derived class must implement this property.
        """
        raise NotImplementedError

    # Basic getters and setters for selection data
    def _getter(self, name):
        # The getter should be used only for values which are the same through the selection.
        return self.atomsel.get(name)[0]

    def _setter(self, name, value):
        # The setter should be used only for values which are the same through the selection.
        self.atomsel.set(name, value)


class IterableSelectionMixin(object):
    """
    Base class for selection-based objects which represent multiple atoms.
    """
    # Large amounts of these objects can be created, slots has some performance benefits.
    __slots__ = ()

    def __len__(self):
        return len(self.atomsel)

    def __iter__(self):
        for index in self.atomsel:
            yield Atom(index, self._molecule, self._frame)

    def __contains__(self, atom):
        assert isinstance(atom, Atom)
        return self.atomsel[atom.index]


class Selection(IterableSelectionMixin, SelectionBase):
    """
    Selection of atoms.

    This class is a proxy to a selection in VMD.
    Coordinate based selections are automatically updated.
    """
    def __init__(self, selection, molecule=None, frame=NOW):
        """
        Creates selection.

        @param selection: Selection text
        @param molecule: Molecule to select from. Top if not provider.
        @type molecule: Molecule or None
        @param frame: Selection frame
        @type frame: Non-negative integer or NOW
        """
        super(Selection, self).__init__(molecule=molecule, frame=frame)
        self._selection = selection
        # No need to delay creation of the atomsel. This also checks if the selection text makes sense.
        self._atomsel = _atomsel(selection, frame=frame, molid=self._molecule.molid)
        # _update_frame is the frame that have been used to filter coordinate based selection in last update.
        self._update_frame = self._get_active_frame()

    def _get_active_frame(self):
        # Return which frame is used to filter coordinate based selection.
        if self._frame == NOW:
            return self._molecule.frame
        else:
            return self._frame

    def __repr__(self):
        return "<%s: '%s' of '%r' at %d>" % (type(self).__name__, self._selection, self._molecule, self._frame)

    @property
    def selection(self):
        "Selection text"
        return self._selection

    @property
    def atomsel(self):
        """
        Returns respective 'VMD.atomsel' instance.
        """
        active_frame = self._get_active_frame()
        # Selection can be coordinate-based. If update frame and active frame differ, update selection.
        if active_frame != self._update_frame:
            self._atomsel.update()
            self._update_frame = active_frame
        return self._atomsel

    ############################################################################
    # Useful methods
    def contacts(self, other, distance):
        """
        Returns iterator of atom pairs which are closer than distance.
        """
        assert isinstance(other, Selection)
        assert isinstance(distance, (int, float, long)) and distance >= 0
        atoms_self, atoms_other = self.atomsel.contacts(other.atomsel, distance)
        return ((Atom(a, self._molecule, self._frame), Atom(b, other.molecule, other.frame))
                for a, b in itertools.izip(atoms_self, atoms_other))


class StaticSelection(SelectionBase):
    """
    Base class for `Atom` and `Residue`.

    Atoms and residues are static objects in VMD. They can't be created or deleted. They have unique identifier which
    does not change. Atoms can not be moved to a different residue.
    """
    # Large amounts of these objects can be created, slots has some performance benefits.
    __slots__ = ('_index', )

    def __init__(self, index, molecule=None, frame=NOW):
        """
        Creates the object.

        @param index: Index of the selection which can not be changed and provide unique identification.
        @param molecule: Selection's molecule. Top if not provider.
        @type molecule: Molecule or None
        @param frame: Selection's frame
        @type frame: Non-negative integer or NOW
        """
        assert isinstance(index, int) and index >= 0
        super(StaticSelection, self).__init__(molecule=molecule, frame=frame)
        self._index = index

    def __repr__(self):
        return "<%s: %d of '%r' at %d>" % (type(self).__name__, self._index, self._molecule, self._frame)

    def __hash__(self):
        # We can use unique identifier for objects which have one.
        return self._index

    def __eq__(self, other):
        return type(self) == type(other) and self._index == other.index and self._molecule == other.molecule and \
            self._frame == other.frame

    def __ne__(self, other):
        return not self.__eq__(other)

    @property
    def index(self):
        "Index"
        return self._index


def _object_property(keyword, doc):
    """Utility function to map VMD keywords to objects properties."""
    def _getter_wrapper(self):
        return self._getter(keyword)

    def _setter_wrapper(self, value):
        self._setter(keyword, value)

    return property(_getter_wrapper, _setter_wrapper, doc=doc)


class Atom(StaticSelection):
    """
    Atom representation.

    This class is a proxy to a atom in molecule loaded into VMD.
    """
    # Large amounts of these objects can be created, slots has some performance benefits.
    __slots__ = ()

    def __init__(self, index, molecule=None, frame=NOW):
        """
        Creates atom representation.

        @param index: Index of the atom
        @param molecule: Atom's molecule. Top if not provider.
        @type molecule: Molecule or None
        @param frame: Atom's frame
        @type frame: Non-negative integer or NOW
        """
        super(Atom, self).__init__(index, molecule=molecule, frame=frame)
        # Check if index makes sense
        if index >= _molecule.numatoms(self._molecule.molid):
            raise ValueError("Atom %d doesn't exist in '%s' at %s" % (index, self._molecule, frame))

    @classmethod
    def pick(cls, selection, molecule=None, frame=NOW):
        """
        Creates atom from selection text.

        @param selection: Selection text
        @param molecule: Molecule to select from. Top if not defined.
        @type molecule: Molecule
        @param frame: Atom's frame
        @type frame: Non-negative integer or NOW
        """
        if molecule is None:
            molecule = MOLECULES.top
        else:
            assert isinstance(molecule, Molecule)
        assert frame == NOW or (isinstance(frame, int) and frame >= 0)

        sel = _atomsel(selection, frame=frame, molid=molecule.molid)
        if len(sel) != 1:
            raise ValueError("Selection '%s' doesn't define single atom in '%s' at %s" % (selection, molecule, frame))
        self = cls(sel.get('index')[0], molecule, frame)
        return self

    @property
    def atomsel(self):
        """
        Returns respective 'VMD.atomsel' instance.
        """
        if self._atomsel is None:
            self._atomsel = _atomsel('index %d' % self._index, frame=self._frame, molid=self._molecule.molid)
        return self._atomsel

    ############################################################################
    # Atom's data
    # Coordinates
    x = _object_property('x', doc="Coordinate in 'x' dimension.")
    y = _object_property('y', doc="Coordinate in 'y' dimension.")
    z = _object_property('z', doc="Coordinate in 'z' dimension.")

    def _get_coords(self):
        #XXX: This is unintuitive, but very fast. The atom's center is the location of the atom.
        # Apparently getting the coordinates all at once has lower overhead than underlying conputation of the center.
        return array(self.atomsel.center())

    def _set_coords(self, value):
        self.x, self.y, self.z = value

    coords = property(_get_coords, _set_coords, doc="Array of (x, y, z) coordinates.")

    # Other data
    name = _object_property('name', doc="Atom name.")
    type = _object_property('type', doc="Atom type.")
    element = _object_property('element', doc="Atom element.")
    beta = _object_property('beta', doc="Atom beta factor.")
    occupancy = _object_property('occupancy', doc="Atom occupancy.")
    mass = _object_property('mass', doc="Atom mass.")
    charge = _object_property('charge', doc="Atom charge.")
    radius = _object_property('radius', doc="Atom radius.")

    ############################################################################
    # Connections to other objects
    @property
    def bonded(self):
        """
        Returns iterator over Atoms bonded to this one.
        """
        return (Atom(i, self._molecule, self._frame) for i in self.atomsel.bonds[0])

    @property
    def residue(self):
        """
        Returns atom's residue.
        """
        return Residue(self._getter('residue'), self._molecule, self._frame)

    def _get_chain(self):
        return Chain(self._getter('chain'), self._molecule, self._frame)

    def _set_chain(self, chain):
        assert isinstance(chain, Chain)
        self._setter('chain', chain.name)

    chain = property(_get_chain, _set_chain, doc="Atom chain")

    def _get_segment(self):
        return Segment(self._getter('segname'), self._molecule, self._frame)

    def _set_segment(self, segment):
        assert isinstance(segment, Segment)
        self._setter('segname', segment.name)

    segment = property(_get_segment, _set_segment, doc="Atom segment")


class Residue(IterableSelectionMixin, StaticSelection):
    """
    Residue representation.

    This class is a proxy to a residue in molecule loaded into VMD.
    The residue is identified by 'residue' value from VMD.
    """
    # Large amounts of these objects can be created, slots has some performance benefits.
    __slots__ = ()

    #TODO: Check if index makes sense in __init__

    @property
    def atomsel(self):
        """
        Returns respective 'VMD.atomsel' instance.
        """
        if self._atomsel is None:
            self._atomsel = _atomsel('residue %d' % self._index, frame=self._frame, molid=self._molecule.molid)
        return self._atomsel

    ############################################################################
    # Residue's data
    number = _object_property('resid', doc="Residue number.")
    name = _object_property('resname', doc="Residue name.")


class DynamicSelection(IterableSelectionMixin, SelectionBase):
    """
    Base class for `Chain` and `Segment`.

    Unlike residues, chains and segments are not static in VMD. They can be created and deleted. They have identifiers,
    but they can be changed. Atoms can be moved between chains and segments.
    """
    # Large amounts of these objects can be created, slots has some performance benefits.
    __slots__ = ('_name', )

    def __init__(self, name, molecule=None, frame=NOW):
        """
        Creates the object.

        @param name: Name of the selection. Also acts as an identifier.
        @param molecule: Selection's molecule. Top if not provider.
        @type molecule: Molecule or None
        @param frame: Selection's frame
        @type frame: Non-negative integer or NOW
        """
        assert isinstance(name, str)
        super(DynamicSelection, self).__init__(molecule=molecule, frame=frame)
        self._name = name

    def __repr__(self):
        return "<%s: '%s' of '%r' at %d>" % (type(self).__name__, self.name, self._molecule, self._frame)

    def __eq__(self, other):
        return type(self) == type(other) and self.name == other.name and self._molecule == other.molecule and \
            self._frame == other.frame

    def __ne__(self, other):
        return not self.__eq__(other)


class Chain(DynamicSelection):
    """
    Chain representation.

    This class is a proxy to a chain in molecule loaded into VMD.
    The chain is identified by 'chain' value from VMD.
    """
    # Large amounts of these objects can be created, slots has some performance benefits.
    __slots__ = ()

    @property
    def atomsel(self):
        """
        Returns respective 'VMD.atomsel' instance.
        """
        if self._atomsel is None:
            self._atomsel = _atomsel('chain "%s"' % self.name, frame=self._frame, molid=self._molecule.molid)
        return self._atomsel

    def _get_name(self):
        return self._name

    def _set_name(self, value):
        # Rename the segment in all its atoms
        self._setter('chain', value)
        self._name = value

    name = property(_get_name, _set_name, doc="Chain name")


class Segment(DynamicSelection):
    """
    Segment representation.

    This class is a proxy to a segment in molecule loaded into VMD.
    The segment is identified by 'segname' value from VMD.
    """
    # Large amounts of these objects can be created, slots has some performance benefits.
    __slots__ = ()

    @property
    def atomsel(self):
        """
        Returns respective 'VMD.atomsel' instance.
        """
        if self._atomsel is None:
            self._atomsel = _atomsel('segname "%s"' % self.name, frame=self._frame, molid=self._molecule.molid)
        return self._atomsel

    def _get_name(self):
        return self._name

    def _set_name(self, value):
        # Rename the segment in all its atoms
        self._setter('segname', value)
        self._name = value

    name = property(_get_name, _set_name, doc="Segment name")
