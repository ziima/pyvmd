"""
Tests for data sets.
"""
from cStringIO import StringIO
import os
from tempfile import mkstemp

import numpy

from pyvmd.datasets import Column, DataSet

from .utils import data, PyvmdTestCase


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
        dset.add_column('first')
        dset.add_column('second')
        dset.add_column('last')
        dset.add_row([5.0, 3.5, 2.9])
        dset.add_row([-0.9, 0.5, 5.9])
        dset.add_row([-42.0, 0.01, 8.9])
        dset.add_row([-204.54, 0.1, 15.89])

        result = numpy.array(([5.0, 3.5, 2.9],
                              [-0.9, 0.5, 5.9],
                              [-42.0, 0.01, 8.9],
                              [-204.54, 0.1, 15.89]))
        self.assertTrue(numpy.array_equal(dset.data, result))

    def test_add_column(self):
        # Test `add_column` method deeply
        dset = DataSet()
        dset.add_column('first')
        dset.add_column('custom', '%d', '%30s')
        dset.add_column('another')

        # Check column with the same name can't be added twice
        self.assertRaises(ValueError, dset.add_column, 'first')

        # Check columns
        result = [Column('first', '%10.4f', '%10s'), Column('custom', '%d', '%30s'),
                  Column('another', '%10.4f', '%10s')]
        self.assertEqual(dset.columns, result)

    def test_data_resize(self):
        # Test adding so many rows, that the data array is resized.
        dset = DataSet()
        dset.step = 2
        dset.add_column('first')
        dset.add_column('second')
        dset.add_row([5.0, 3.5])
        dset.add_row([-0.9, 0.5])
        dset.add_row([-42.0, 0.01])
        dset.add_row([-204.54, 0.1])
        dset.add_row([-15.4, 0.5])

        result = numpy.array(([5.0, 3.5],
                              [-0.9, 0.5],
                              [-42.0, 0.01],
                              [-204.54, 0.1],
                              [-15.4, 0.5]))
        self.assertTrue(numpy.array_equal(dset.data, result))

    def test_write(self):
        dset = DataSet()
        dset.add_column('first', '% 20d', '%20s')
        dset.add_column('second')
        dset.add_column('last')
        dset.add_row([5.0, 3.5, 2.9])
        dset.add_row([-0.9, 0.5, 5.9])
        dset.add_row([-42.0, 0.01, 8.9])
        dset.add_row([-204.54, 0.1, 15.89])

        # Test filename
        dset.write(self.tmpfile)
        self.assertEqual(open(self.tmpfile).read(), open(data('dataset.dat')).read())

        # Test file-like object
        buf = StringIO()
        dset.write(buf)
        self.assertEqual(buf.getvalue(), open(data('dataset.dat')).read())
