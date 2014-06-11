"""
High level interface to the VMD objects.
"""
import logging

from Molecule import Molecule as _Molecule
import VMD


LOGGER = logging.getLogger(__name__)


class Frames(object):
    """
    Wrapper for molecules' frames.
    """
    def __init__(self, molecule):
        """
        @param molecule: Respective molecule
        @type molecule: Molecule
        """
        # Use molecule instance instead of molid for possible callbacks
        assert isinstance(molecule, Molecule)
        self.molecule = molecule

    def __len__(self):
        return VMD.molecule.numframes(self.molecule.molid)

    def __delitem__(self, key):
        #XXX: For some reason, 'skip' in the 'delframe' function means which frames are left when deleting
        # That is not consistent with python slicing, so we have to avoid that argument
        if isinstance(key, slice):
            start, stop, step = key.indices(len(self))
            # We will delete one by one, so we have to that in reversed order
            frames = reversed(xrange(start, stop, step))
        elif isinstance(key, int):
            if key < 0:
                frames = [len(self) + key]
            else:
                frames = [key]
        else:
            raise TypeError("%s indices must be integers, not %s" % (type(self), type(key)))

        for frame in frames:
            LOGGER.debug("Deleting frame %d", frame)
            VMD.molecule.delframe(self.molecule.molid, beg=frame, end=frame)

    def __iter__(self):
        # Return the iterator over frames
        return iter(xrange(len(self)))

    def copy(self, frame=None):
        """
        Copies frame and moves the molecule to the new frame.
        """
        if frame is None:
            frame = self.molecule.frame
        else:
            assert isinstance(frame, int) and frame >= 0
        VMD.molecule.dupframe(self.molecule.molid, frame)


class Molecule(object):
    """
    Molecule representation.

    This class is a proxy for molecule loaded into VMD.
    """
    def __init__(self, molid):
        """
        Creates a new molecule.

        @param molid: Molecule ID
        """
        assert isinstance(molid, int) and molid >= 0
        if not VMD.molecule.exists(molid):
            raise ValueError("Molecule %d does not exist." % molid)
        self.molid = molid
        self._molecule = None

    def __eq__(self, other):
        return type(self) == type(other) and self.molid == other.molid

    def __ne__(self, other):
        return not self.__eq__(other)

    def molecule(self):
        """
        Returns respective VMD.Molecule instance.
        """
        if self._molecule is None:
            self._molecule = _Molecule(id=self.molid)
        return self._molecule

    def _get_frame(self):
        return VMD.molecule.get_frame(self.molid)

    def _set_frame(self, frame):
        assert isinstance(frame, int) and frame >= 0
        VMD.molecule.set_frame(self.molid, frame)

    frame = property(_get_frame, _set_frame, doc="Molecule's frame")

    @property
    def frames(self):
        """
        Returns frames descriptor.
        """
        return Frames(self)
