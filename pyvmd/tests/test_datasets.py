"""
Tests for data sets.
"""
import os
from cStringIO import StringIO
from tempfile import mkstemp

import numpy
from mock import Mock

from pyvmd.collectors import Collector
from pyvmd.datasets import DataSet

from .utils import data, PyvmdTestCase


class SimpleTestCollector(Collector):
    """
    Simple collector which returns given data.
    """
    def __init__(self, results, name):
        super(SimpleTestCollector, self).__init__(name)
        self.results = results

    def collect(self, step):
        return self.results.pop(0)


class TestDataSet(PyvmdTestCase):
    """
    Test `DataSet` object.
    """
    def setUp(self):
        dummy, filename = mkstemp(prefix='pyvmd_test_')
        self.tmpfile = filename
        self.addCleanup(lambda: os.unlink(self.tmpfile))

    def test_dataset(self):
        # Test basic dataset workflow - add columns, add data, get data
        dset = DataSet()
        # Register few collectors
        dset.add_collector(SimpleTestCollector([5.0, -0.9, -42.0, -204.54], 'first'))
        dset.add_collector(SimpleTestCollector([3.5, 0.5, 0.01, 0.1], 'second'))
        dset.add_collector(SimpleTestCollector([2.9, 5.9, 8.9, 15.89], 'last'))

        # Collect the data
        dset.collect(Mock(frame=0))
        dset.collect(Mock(frame=1))
        dset.collect(Mock(frame=2))
        dset.collect(Mock(frame=42))

        result = numpy.array(([0, 5.0, 3.5, 2.9],
                              [1, -0.9, 0.5, 5.9],
                              [2, -42.0, 0.01, 8.9],
                              [42, -204.54, 0.1, 15.89]))
        self.assertTrue(numpy.array_equal(dset.data, result))

    def test_data_resize(self):
        # Test adding so many rows, that the data array is resized.
        dset = DataSet()
        dset.step = 2  # Lower the array size step so we can easily enforce resizing

        # Register few collectors
        dset.add_collector(SimpleTestCollector([5.0, -0.9, -42.0, -204.54, -15.4], 'first'))
        dset.add_collector(SimpleTestCollector([3.5, 0.5, 0.01, 0.1, 0.5], 'second'))

        # Collect the data
        dset.collect(Mock(frame=0))
        dset.collect(Mock(frame=1))
        dset.collect(Mock(frame=2))
        dset.collect(Mock(frame=3))
        dset.collect(Mock(frame=4))

        result = numpy.array(([0, 5.0, 3.5],
                              [1, -0.9, 0.5],
                              [2, -42.0, 0.01],
                              [3, -204.54, 0.1],
                              [4, -15.4, 0.5]))
        self.assertTrue(numpy.array_equal(dset.data, result))

    def test_write(self):
        dset = DataSet()

        # Register few collectors
        first = SimpleTestCollector([5.0, -0.9, -42.0, -204.54], 'first')
        # Set custom format for first collector
        first.header_fmt = '%20s'
        first.data_fmt = '% 20d'
        dset.add_collector(first)
        dset.add_collector(SimpleTestCollector([3.5, 0.5, 0.01, 0.1], 'second'))
        dset.add_collector(SimpleTestCollector([2.9, 5.9, 8.9, 15.89], 'last'))

        # Collect the data
        dset.collect(Mock(frame=0))
        dset.collect(Mock(frame=1))
        dset.collect(Mock(frame=2))
        dset.collect(Mock(frame=3))

        # Test filename
        dset.write(self.tmpfile)
        self.assertEqual(open(self.tmpfile).read(), open(data('dataset.dat')).read())

        # Test file-like object
        buf = StringIO()
        dset.write(buf)
        self.assertEqual(buf.getvalue(), open(data('dataset.dat')).read())
