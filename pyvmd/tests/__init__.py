"""
Unittests for pyvmd.
"""
import os
import sys
import unittest

if not hasattr(unittest, 'skip'):
    # XXX: Use unittest2 to use python 2.7 unittest features.
    import unittest2 as unittest


def cover_main():
    """
    Run tests with coverage.
    """
    import coverage
    cov = coverage.coverage(branch=True, source=['pyvmd'])
    cov.start()
    unit = unittest.main(exit=False)
    cov.stop()
    cov.save()

    if not hasattr(unit, 'result'):
        # Unittest help exit
        sys.exit(2)
    else:
        sys.exit(not unit.result.wasSuccessful())


if __name__ == '__main__':
    if os.environ.get('COVERAGE'):
        cover_main()
    else:
        unittest.main()
