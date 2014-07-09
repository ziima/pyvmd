pyvmd
=====

Python tools for VMD (Visual Molecular Dynamics)

Pyvmd aims to provide high-level pythonic interface to VMD.

### Requirements ###
 * VMD with python support
 * numpy
 * setuptools - for installation only

### Installation ###
Pyvmd has standard `setup.py` based on setuptools. Install using `python setup.py install`.

### Usage ###
Run `vmd -python` to start VMD with python shell. Then you can use pyvmd as any other python library.

Python scripts for VMD can be run using `vmd -python -e my_script.py`. You can pass arguments to the script using
`-args` option.

### Documentation ###
Documentation with examples can be found in [docs](docs).
Detailed documentation of API can be found using python's builtin `help` function.
