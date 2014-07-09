"""
Data sets to collect and store data from trajectory analysis.
"""
from collections import namedtuple
import logging

import numpy

__all__ = ['Column', 'DataSet']

LOGGER = logging.getLogger(__name__)


# Contains meta data for DataSet column
Column = namedtuple('Column', ('name', 'fmt', 'header_fmt'))


class DataSet(object):
    """
    Basic data set. Collects and stores data extracted from trajectory.
    """
    # Number of rows added when increasing size of data array
    step = 1000

    def __init__(self):
        # Array with data
        self._data = None
        # List of columns
        self.columns = []
        # Number of rows filled with data
        self._rows = 0

    @property
    def data(self):
        """
        Returns collected data.
        """
        # Return only the data, not the pre-allocated space.
        return self._data[:self._rows]

    def add_column(self, name, fmt='%10.4f', header_fmt='%10s'):
        """
        Adds column to the list of columns.

        @param name: Unique name of a column
        """
        assert isinstance(name, basestring)
        if name in [c.name for c in self.columns]:
            raise ValueError("The column '%s' already exists." % name)
        column = Column(name, fmt, header_fmt)
        self.columns.append(column)

    def add_row(self, row):
        """
        Adds list of data as a new row.
        """
        num_columns = len(self.columns)
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

        @param output: Filename of file-like object
        """
        if isinstance(output, basestring):
            out = open(output, 'w')
        else:
            out = output

        # Write header
        header_fmt = ' '.join([c.header_fmt for c in self.columns])
        out.write(header_fmt % tuple(c.name for c in self.columns))
        out.write('\n')

        # Write data
        formats = [c.fmt for c in self.columns]
        numpy.savetxt(out, self.data, formats)
