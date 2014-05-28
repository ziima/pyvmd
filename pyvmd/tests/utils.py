"""
Utilities for tests.
"""
import os


def data(filename):
    """
    Return full filename of the test datafile.
    """
    return os.path.join(os.path.dirname(__file__), 'data', filename)
