"""
Unittests for pyvmd.
"""
import os
import sys
import unittest

if not hasattr(unittest, 'skip'):
    #XXX: Use unittest2 to use python 2.7 unittest features.
    import unittest2 as unittest


if __name__ == '__main__':
    # Start coverage if required
    cover = os.environ.get('COVERAGE')
    if cover:
        import coverage
        cov = coverage.coverage(branch=True, source=['pyvmd'])
        cov.start()

    unit = unittest.main(exit=False)

    if cover:
        cov.stop()
        cov.save()

    if not hasattr(unit, 'result'):
        # Unittest help exit
        sys.exit(2)
    else:
        sys.exit(not unit.result.wasSuccessful())
