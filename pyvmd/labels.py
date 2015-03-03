"""
Interface for VMD labels.
"""
from VMD import label as _label

from .atoms import Atom
from .molecules import Molecule

__all__ = ['AngleLabel', 'AtomLabel', 'BondLabel', 'DihedralLabel', 'ANGLE', 'ANGLE_LABELS', 'ATOM', 'ATOM_LABELS', \
           'BOND', 'BOND_LABELS', 'DIHEDRAL', 'DIHEDRAL_LABELS']


ATOM = _label.ATOM
BOND = _label.BOND
ANGLE = _label.ANGLE
DIHEDRAL = _label.DIHEDRAL
# String representation of category
CATEGORY_NAMES = {
    ATOM: 'ATOM',
    BOND: 'BOND',
    ANGLE: 'ANGLE',
    DIHEDRAL: 'DIHEDRAL',
}


def get_label_data(category, atomids, molids):
    """
    Utility function to get VMD data for label.
    """
    rev_atomids = atomids[::-1]
    rev_molids = molids[::-1]
    for hit in _label.listall(category):
        if hit['atomid'] == atomids and hit['molid'] == molids:
            # Found it, return data
            return hit
        # The lists of atomids and molids may be reversed
        elif hit['atomid'] == rev_atomids and hit['molid'] == rev_molids:
            # Found it, return data
            return hit


class BaseLabel(object):
    """
    Representation of atom label.
    """
    _category = None

    def __init__(self, *atoms):
        self._atoms = tuple(atoms)
        # Cache tuples of atom IDs and molecule IDs for VMD interface
        self._atomids = tuple(a.index for a in self._atoms)
        self._molids = tuple(a.molecule.molid for a in self._atoms)
        if get_label_data(self._category, self._atomids, self._molids) is None:
            raise ValueError("%s for %s doesn't exist" % (type(self).__name__, atoms))

    def __repr__(self):
        return "<%s: %r>" % (type(self).__name__, self._atoms)

    def __eq__(self, other):
        if type(self) != type(other) or self.category != other.category:
            return False
        # Labels are independent on frames, so we can't use atom comparison, we need to check atom and molecule IDs
        self_signature = zip(self._atomids, self._molids)
        other_signature = [(a.index, a.molecule.molid) for a in other.atoms]
        # Labels are equal if atom order is reversed
        return self_signature == other_signature or self_signature == other_signature[::-1]

    def __ne__(self, other):
        return not self.__eq__(other)

    @property
    def category(self):
        """
        Label category
        """
        return self._category

    @property
    def atoms(self):
        """
        Labeled atoms
        """
        return self._atoms

    @classmethod
    def create(cls, *atoms):
        """
        Creates the label if it doesn't exist.
        """
        atomids = tuple(a.index for a in atoms)
        molids = tuple(a.molecule.molid for a in atoms)
        _label.add(cls._category, molids, atomids)
        return cls(*atoms)

    def delete(self):
        """
        Delete the label.
        """
        _label.delete(self._category, {'molid': self._molids, 'atomid': self._atomids})

    def _get_visible(self):
        data = get_label_data(self._category, self._atomids, self._molids)
        if data is None:
            raise ValueError("%s for %s doesn't exist" % (type(self).__name__, self._atoms))
        # Found it, return 'on' key converted to boolean
        return bool(data['on'])

    def _set_visible(self, value):
        if value:
            _label.show(self._category, {'molid': self._molids, 'atomid': self._atomids})
        else:
            _label.hide(self._category, {'molid': self._molids, 'atomid': self._atomids})

    visible = property(_get_visible, _set_visible, doc="Label's visibility")


class AtomLabel(BaseLabel):
    _category = ATOM

    def __init__(self, a):
        assert isinstance(a, Atom)
        super(AtomLabel, self).__init__(a)

    @classmethod
    def create(cls, a):
        return super(AtomLabel, cls).create(a)

    @property
    def atom(self):
        """
        Labeled atom.
        """
        return self._atoms[0]


class BondLabel(BaseLabel):
    _category = BOND

    def __init__(self, a, b):
        assert isinstance(a, Atom)
        assert isinstance(b, Atom)
        super(BondLabel, self).__init__(a, b)

    @property
    def distances(self):
        """
        Returns distances between the atoms throughout the loaded trajectory.
        """
        return _label.getvalues(self._category, {'molid': self._molids, 'atomid': self._atomids})


class AngleLabel(BaseLabel):
    _category = ANGLE

    def __init__(self, a, b, c):
        assert isinstance(a, Atom)
        assert isinstance(b, Atom)
        assert isinstance(c, Atom)
        super(AngleLabel, self).__init__(a, b, c)

    @property
    def angles(self):
        """
        Returns angles between the atoms throughout the loaded trajectory.
        """
        return _label.getvalues(self._category, {'molid': self._molids, 'atomid': self._atomids})


class DihedralLabel(BaseLabel):
    _category = DIHEDRAL

    def __init__(self, a, b, c, d):
        assert isinstance(a, Atom)
        assert isinstance(b, Atom)
        assert isinstance(c, Atom)
        assert isinstance(d, Atom)
        super(DihedralLabel, self).__init__(a, b, c, d)

    @property
    def dihedrals(self):
        """
        Returns dihedral angles between the atoms throughout the loaded trajectory.
        """
        return _label.getvalues(self._category, {'molid': self._molids, 'atomid': self._atomids})


class LabelManager(object):
    """
    Manages all labels in the category.
    """
    def __init__(self, category, label_cls):
        assert category in (ATOM, BOND, ANGLE, DIHEDRAL)
        self.category = category
        self._label_cls = label_cls

    def __repr__(self):
        return "<%s: %s>" % (type(self).__name__, CATEGORY_NAMES[self.category])

    def __len__(self):
        return len(_label.listall(self.category))

    def __iter__(self):
        for label in _label.listall(self.category):
            yield self._label_cls(*(Atom(a, Molecule(m)) for a, m in zip(label['atomid'], label['molid'])))

    def __contains__(self, label):
        # Check for category is fast
        if label.category != self.category:
            return False

        # Get ids for VMD interface
        atomids = tuple(a.index for a in label.atoms)
        molids = tuple(a.molecule.molid for a in label.atoms)
        data = get_label_data(self.category, atomids, molids)
        return data is not None


ATOM_LABELS = LabelManager(ATOM, AtomLabel)
BOND_LABELS = LabelManager(BOND, BondLabel)
ANGLE_LABELS = LabelManager(ANGLE, AngleLabel)
DIHEDRAL_LABELS = LabelManager(DIHEDRAL, DihedralLabel)
