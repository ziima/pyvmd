"""
Data collectors for trajectory analysis.
"""
from collections import namedtuple
import logging

import VMD


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


RMSDSet = namedtuple('RMSDSet', ('selection', 'name'))


class RMSDCollector(Collector):
    """
    Collects RMSD data.
    """
    def __init__(self, reference, output):
        """
        Create RMSD collector.

        @param reference: Reference molecule
        @type reference: Molecule
        @param output: Output filename
        """
        super(RMSDCollector, self).__init__(output)
        self.reference = reference
        self._rmsdsets = []
        # Last used ID for automatically generated names
        self._autoid = 0

    def _prepare_output(self):
        # Clean/create the output file and write the header
        with open(self.output, 'w') as out:
            out.write('#Frame  ')
            for rmsdset in self._rmsdsets:
                out.write(' %10s' % rmsdset.name)
            out.write('\n')

    def collect(self, status):
        if status.frame == 0:
            self._prepare_output()

        data = []
        # Current frame number
        cur_frame = status.molecule.curFrame()
        # Duplicate the trajectory frame because we will modify the coordinates
        # This also sets the molecule to the duplicated frame
        status.molecule.dupFrame()
        # Duplicated frame number
        dup_frame = status.molecule.curFrame()

        whole = VMD.atomsel.atomsel('all', molid=status.molecule.id)

        for rmsdset in self._rmsdsets:
            ref = VMD.atomsel.atomsel(rmsdset.selection, frame=self.reference.curFrame(), molid=self.reference.id)
            sel = VMD.atomsel.atomsel(rmsdset.selection, molid=status.molecule.id)

            # Align coordinates to the reference
            whole.move(sel.fit(ref))

            # Measure RMSD
            data.append(sel.rmsd(ref))

        # Delete the duplicated frame and reset trajectory frame
        status.molecule.delFrame(first=dup_frame, last=dup_frame)
        status.molecule.setFrame(cur_frame)

        # Write the data
        with open(self.output, 'a') as out:
            out.write('%8d' % status.frame)
            for value in data:
                out.write(' %10.4f' % value)
            out.write('\n')

    def _check_name(self, name):
        for rmsdset in self._rmsdsets:
            if name == rmsdset.name:
                raise ValueError("The data set name '%s' is already used." % name)

    def add_selection(self, selection, name=None):
        """
        Adds selection to RMSD analysis.

        @param selection: The selection for RMSD.
        @param name: Name of the column in output file. If not provided, it is generated in form 'rmsdXXXXX'.
        """
        if name is None:
            self._autoid += 1
            name = 'rmsd%05d' % self._autoid
        else:
            # Check if name was already used
            self._check_name(name)
        LOGGER.debug("Add selection '%s' under name '%s'", selection, name)
        self._rmsdsets.append(RMSDSet(selection, name))
