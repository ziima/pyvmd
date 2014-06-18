"""
Data collectors for trajectory analysis.
"""
from collections import namedtuple
import logging

import VMD

from .objects import Molecule


LOGGER = logging.getLogger(__name__)


class Collector(object):
    """
    Base class for collectors. Collects data from trajectory and writes them into file.
    """
    def __init__(self, output):
        """
        Creates collector.

        @param output: Output filename
        """
        self.output = output

    def collect(self, status):
        """
        Performs the analysis on the frame.

        Derived class must implement this method.
        """
        raise NotImplementedError


DataSet = namedtuple('DataSet', ('selection', 'name'))


class DataCollector(Collector):
    """
    Base class for collectors which collects data for several different selections.
    """
    def __init__(self, output):
        """
        Creates selection collector.

        @param output: Output filename
        """
        super(DataCollector, self).__init__(output)
        self._datasets = []
        # Last used ID for automatically generated names
        self._autoid = 0

    def _prepare_output(self):
        """
        Creates or cleans the output file and writes the header.
        """
        with open(self.output, 'w') as out:
            out.write('#Frame  ')
            for dataset in self._datasets:
                out.write(' %10s' % dataset.name)
            out.write('\n')

    def _write_data(self, frame, data):
        """
        Writes data row into output file
        """
        with open(self.output, 'a') as out:
            out.write('%8d' % frame)
            for value in data:
                out.write(' %10.4f' % value)
            out.write('\n')

    def _check_name(self, name):
        for dataset in self._datasets:
            if name == dataset.name:
                raise ValueError("The data set name '%s' is already used." % name)

    def add_selection(self, selection, name=None):
        """
        Adds selection to analysis.

        @param selection: The selection string.
        @param name: Name of the column in output file. If not provided, it is generated in form 'dataXXXXX'.
        """
        if name is None:
            self._autoid += 1
            name = 'data%05d' % self._autoid
        else:
            # Check if name was already used
            self._check_name(name)
        self._datasets.append(DataSet(selection, name))
        LOGGER.debug("Added data set '%s' under name '%s'", selection, name)


class RMSDCollector(DataCollector):
    """
    Collects RMSD data.
    """
    def __init__(self, reference, output):
        """
        Creates RMSD collector.

        @param reference: Reference molecule
        @type reference: Molecule
        @param output: Output filename
        """
        assert isinstance(reference, Molecule)
        super(RMSDCollector, self).__init__(output)
        self.reference = reference

    def collect(self, status):
        if status.frame == 0:
            self._prepare_output()

        data = []
        # Current frame number
        cur_frame = status.molecule.frame
        # Duplicate the trajectory frame because we will modify the coordinates
        # This also sets the molecule to the duplicated frame
        status.molecule.frames.copy()
        # Duplicated frame number
        dup_frame = status.molecule.frame

        whole = VMD.atomsel.atomsel('all', molid=status.molecule.molid)

        for dataset in self._datasets:
            ref = VMD.atomsel.atomsel(dataset.selection, frame=self.reference.frame, molid=self.reference.molid)
            sel = VMD.atomsel.atomsel(dataset.selection, molid=status.molecule.molid)

            # Align coordinates to the reference
            whole.move(sel.fit(ref))

            # Measure RMSD
            data.append(sel.rmsd(ref))

        # Delete the duplicated frame and reset trajectory frame
        del status.molecule.frames[dup_frame]
        status.molecule.frame = cur_frame

        # Write the data
        self._write_data(status.frame, data)
