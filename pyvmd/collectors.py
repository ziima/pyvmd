"""
Data collectors for trajectory analysis.
"""
import logging

from . import measure
from .atoms import Selection

__all__ = ['AngleCollector', 'Collector', 'DihedralCollector', 'DistanceCollector', 'FrameCollector', 'RMSDCollector',
           'XCoordCollector', 'YCoordCollector', 'ZCoordCollector']


LOGGER = logging.getLogger(__name__)


class Collector(object):
    """
    Base class for collectors. Collects data from trajectory and writes them into file.
    """
    # Format of the column header in the output
    header_fmt = '%10s'
    # Format of the data in the output
    data_fmt = '%10.4f'

    # Counter for automatic name generation.
    auto_name_counter = 0

    def __init__(self, name=None):
        """
        Create the collector.

        @param name: Name of the collector. If not provided, it is generated in form 'data#####'.
        """
        if name is None:
            Collector.auto_name_counter += 1
            name = 'data%05d' % Collector.auto_name_counter
        self.name = name

    def collect(self, step):
        """
        Performs the analysis on the frame.

        Derived class must implement this method.
        """
        raise NotImplementedError


class FrameCollector(Collector):
    """
    Utility collector which collects frame number.
    """
    header_fmt = '%8s'
    data_fmt = '%8d'

    def collect(self, step):
        return step.frame


def _selection_center(selection_text, molecule):
    """
    Utility method to get center of selection.
    """
    sel = Selection(selection_text, molecule)
    if not len(sel):
        raise ValueError("Selection '%s' doesn't match any atoms." % selection_text)
    return measure.center(sel)


class BaseCoordCollector(Collector):
    """
    Base class for collectors of X, Y and Z coordinates.
    """
    def __init__(self, selection, name=None):
        """
        Creates coordinate collector.

        @param selection: Selection text for collector.
        """
        assert isinstance(selection, basestring)
        super(BaseCoordCollector, self).__init__(name)
        self.selection = selection


class XCoordCollector(BaseCoordCollector):
    """
    Collects X coordinate of atom or center of selection.
    """
    def collect(self, step):
        return _selection_center(self.selection, step.molecule)[0]


class YCoordCollector(BaseCoordCollector):
    """
    Collects Y coordinate of atom or center of selection.
    """
    def collect(self, step):
        return _selection_center(self.selection, step.molecule)[1]


class ZCoordCollector(BaseCoordCollector):
    """
    Collects Z coordinate of atom or center of selection.
    """
    def collect(self, step):
        return _selection_center(self.selection, step.molecule)[2]


class DistanceCollector(Collector):
    """
    Collects distance between two atoms or centers of atoms.
    """
    def __init__(self, selection1, selection2, name=None):
        """
        Creates distance collector.

        @param selection1: Selection text for collector.
        @param selection2: Selection text for collector.
        """
        assert isinstance(selection1, basestring)
        assert isinstance(selection2, basestring)
        super(DistanceCollector, self).__init__(name)
        self.selection1 = selection1
        self.selection2 = selection2

    def collect(self, step):
        center1 = _selection_center(self.selection1, step.molecule)
        center2 = _selection_center(self.selection2, step.molecule)
        return measure.coords_distance(center1, center2)


class AngleCollector(Collector):
    """
    Collects angle between three atoms or centers of atoms.
    """
    def __init__(self, selection1, selection2, selection3, name=None):
        """
        Creates distance collector.

        @param selection1: Selection text for collector.
        @param selection2: Selection text for collector.
        @param selection3: Selection text for collector.
        """
        assert isinstance(selection1, basestring)
        assert isinstance(selection2, basestring)
        assert isinstance(selection3, basestring)
        super(AngleCollector, self).__init__(name)
        self.selection1 = selection1
        self.selection2 = selection2
        self.selection3 = selection3

    def collect(self, step):
        center1 = _selection_center(self.selection1, step.molecule)
        center2 = _selection_center(self.selection2, step.molecule)
        center3 = _selection_center(self.selection3, step.molecule)
        return measure.coords_angle(center1, center2, center3)


class DihedralCollector(Collector):
    """
    Collects dihedral angle of four atoms or centers of atoms.
    """
    def __init__(self, selection1, selection2, selection3, selection4, name=None):
        """
        Creates distance collector.

        @param selection1: Selection text for collector.
        @param selection2: Selection text for collector.
        @param selection3: Selection text for collector.
        """
        assert isinstance(selection1, basestring)
        assert isinstance(selection2, basestring)
        assert isinstance(selection3, basestring)
        assert isinstance(selection4, basestring)
        super(DihedralCollector, self).__init__(name)
        self.selection1 = selection1
        self.selection2 = selection2
        self.selection3 = selection3
        self.selection4 = selection4

    def collect(self, step):
        center1 = _selection_center(self.selection1, step.molecule)
        center2 = _selection_center(self.selection2, step.molecule)
        center3 = _selection_center(self.selection3, step.molecule)
        center4 = _selection_center(self.selection4, step.molecule)
        return measure.coords_dihedral(center1, center2, center3, center4)


class RMSDCollector(Collector):
    """
    Collects RMSD data.
    """
    def __init__(self, selection, reference, name=None):
        """
        Creates RMSD collector.

        @param selection: Selection text for RMSD
        @param reference: Reference for RMSD
        @type reference: Selection
        """
        assert isinstance(selection, basestring)
        assert isinstance(reference, Selection)
        super(RMSDCollector, self).__init__(name)
        self.selection = selection
        self.reference = reference

    def collect(self, step):
        # Active frame number of the molecule.
        cur_frame = step.molecule.frame
        # Duplicate the trajectory frame because we will modify the coordinates.
        # This also sets the molecule frame to the duplicated frame.
        step.molecule.frames.copy()
        # Duplicated frame number
        dup_frame = step.molecule.frame

        all_atoms = Selection('all', step.molecule)
        sel = Selection(self.selection, step.molecule)

        # Align coordinates to the reference
        all_atoms.atomsel.move(sel.atomsel.fit(self.reference.atomsel))

        # Measure RMSD
        rmsd = sel.atomsel.rmsd(self.reference.atomsel)

        # Delete the duplicated frame and reset trajectory frame
        del step.molecule.frames[dup_frame]
        step.molecule.frame = cur_frame

        # Return the RMSD value
        return rmsd
