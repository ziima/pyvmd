"""
Data sets to collect and store data from trajectory analysis.
"""
import logging

import numpy

from .collectors import Collector, FrameCollector

__all__ = ['DataSet']

LOGGER = logging.getLogger(__name__)


class DataSet(object):
    """
    Basic data set. Collects and stores data extracted from trajectory.
    """
    # Number of rows added when increasing size of data array
    step = 1000

    def __init__(self):
        # Array with data
        self._data = None
        # List of collectors
        self.collectors = []
        # Number of rows filled with data
        self._rows = 0

        # Register the frame collector
        self.add_collector(FrameCollector('frame'))

    @property
    def data(self):
        """
        Returns collected data.
        """
        # Return only the collected data, not the pre-allocated space.
        return self._data[:self._rows]

    def add_collector(self, collector):
        """
        Registers collector into this dataset.
        """
        assert isinstance(collector, Collector)
        self.collectors.append(collector)
        LOGGER.debug("Added collector '%s' to dataset '%s'", collector.name, self)

    def collect(self, step):
        """
        Retrieves data from collectors and store them. Callback for Analyzer.
        """
        row = [collector.collect(step) for collector in self.collectors]
        self._add_row(row)

    def _add_row(self, row):
        """
        Stores new row with data.
        """
        num_columns = len(self.collectors)
        num_rows = self._rows + 1

        # Create array if it doesn't exist
        if self._data is None:
            shape = (self.step, num_columns)
            LOGGER.debug("Creating array %s", shape)
            self._data = numpy.empty(shape)
        # Check if the array is large enough. If not, enlarge
        # Get number of rows in data array
        data_rows = self._data.shape[0]
        if num_rows > data_rows:
            # The existing array is not big enough. Create new larger one and copy data.
            shape = (data_rows + self.step, num_columns)
            LOGGER.debug("Enlarging array %s to %s", self._data.shape, shape)
            new_data = numpy.empty(shape)
            new_data[:data_rows] = self._data
            self._data = new_data

        # Store the new data
        self._data[self._rows] = row  # Since arrays use 0-based index, row number equals the old number of rows.
        self._rows = num_rows

    def write(self, output):
        """
        Writes data into output.

        @param output: Filename or file-like object
        """
        if isinstance(output, basestring):
            out = open(output, 'w')
        else:
            out = output

        # Write header
        header_fmt = ' '.join([c.header_fmt for c in self.collectors])
        out.write(header_fmt % tuple(c.name for c in self.collectors))
        out.write('\n')

        # Write data
        formats = [c.data_fmt for c in self.collectors]
        numpy.savetxt(out, self.data, formats)
