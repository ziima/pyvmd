"""
Data collectors for trajectory analysis.
"""
import logging

import VMD

from .datasets import DataSet
from .molecules import Molecule

__all__ = ['Collector', 'DataCollector', 'RMSDCollector']


LOGGER = logging.getLogger(__name__)


class Collector(object):
    """
    Base class for collectors. Collects data from trajectory and writes them into file.
    """
    def collect(self, status):
        """
        Performs the analysis on the frame.

        Derived class must implement this method.
        """
        raise NotImplementedError


class DataCollector(Collector):
    """
    Base class for collectors which collects data for several different selections.
    """
    def __init__(self):
        """
        Creates selection collector.
        """
        self.dataset = DataSet()
        self.dataset.add_column('frame', '%8d', '%8s')
        # List of selection strings for analysis
        self._selections = []
        # Last used ID for automatically generated names
        self._autoid = 0

    def add_selection(self, selection, name=None):
        """
        Adds selection to analysis.

        @param selection: The selection string.
        @param name: Name of the column in data set. If not provided, it is generated in form 'data#####'.
        """
        if name is None:
            self._autoid += 1
            name = 'data%05d' % self._autoid
        self.dataset.add_column(name)
        self._selections.append(selection)
        LOGGER.debug("Added selection '%s' named '%s'", selection, name)


class RMSDCollector(DataCollector):
    """
    Collects RMSD data.
    """
    def __init__(self, reference):
        """
        Creates RMSD collector.

        @param reference: Reference molecule
        @type reference: Molecule
        """
        assert isinstance(reference, Molecule)
        super(RMSDCollector, self).__init__()
        self.reference = reference

    def collect(self, status):
        # Current frame number
        cur_frame = status.molecule.frame
        # Duplicate the trajectory frame because we will modify the coordinates
        # This also sets the molecule to the duplicated frame
        status.molecule.frames.copy()
        # Duplicated frame number
        dup_frame = status.molecule.frame

        whole = VMD.atomsel.atomsel('all', molid=status.molecule.molid)

        # Collect the data
        data = [status.frame]
        for selection in self._selections:
            ref = VMD.atomsel.atomsel(selection, frame=self.reference.frame, molid=self.reference.molid)
            sel = VMD.atomsel.atomsel(selection, molid=status.molecule.molid)

            # Align coordinates to the reference
            whole.move(sel.fit(ref))

            # Measure RMSD
            data.append(sel.rmsd(ref))

        # Delete the duplicated frame and reset trajectory frame
        del status.molecule.frames[dup_frame]
        status.molecule.frame = cur_frame

        # Store the data
        self.dataset.add_row(data)
