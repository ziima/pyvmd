# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

from pyvmd import __version__ as version


CLASSIFIERS = [
    'Development Status :: 3 - Alpha',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    'Operating System :: OS Independent',
    'Programming Language :: Python',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Scientific/Engineering :: Chemistry',
    'Topic :: Scientific/Engineering :: Information Analysis',
    'Topic :: Scientific/Engineering :: Visualization',
]


if __name__ == '__main__':
    setup(
        name='PyVMD',
        version=version,
        packages=find_packages(exclude=['*.tests', '*.tests.*']),
        install_requires=['numpy'],
        author="Vlastimil Zíma",
        author_email="vlastimil.zima@gmail.com",
        description="Python tools for Visual Molecular Dynamics",
        url="https://github.com/ziima/pyvmd",
        classifiers=CLASSIFIERS,
        license='GPLv3',
    )
