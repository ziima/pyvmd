# Measurements #

Among other utilities pyvmd provides also functions to measure properties of system.
The utilities reside in module `pyvmd.measure`.

## List of functions ##
 * `distance(a, b)` - returns distance between two atoms `a` and `b`.
 * `angle(a, b, c)` - returns angle between three atoms `a`, `b` and `c`.
 * `dihedral(a, b, c, d)` - returns dihedral or improper dihedral angle of four atoms `a`, `b`, `c` and `d`.
 * `center(selection)` - returns coordinates of geometric center of the selection or iterable of atoms.

### Examples ###
```python
from pyvmd import measure
from pyvmd.atoms import Atom, Residue, Selection

# Measure distance between atoms
measure.distance(Atom(0), Atom(42))  #>>> 4.578

# Measure angle between atoms
measure.angle(Atom(0), Atom(42), Atom(1))  #>>> 73.278

# Measure dihedral angle of atoms
measure.dihedral(Atom(0), Atom(42), Atom(1), Atom(2))  #>>> -58.378

# Measure center of selection
measure.center(Selection('all and noh'))  #>>> array([0.25, 0.78, 80.58])
# Measure center of residue
measure.center(Residue(12))  #>>> array([12.89, 57.32, 75.38])
# Measure center of atom iterable
my_atoms = (Atom(i) for i in xrange(10))
measure.center(my_atoms)  #>>> array([8.95, 3.9, -59.3])
# Center of atom is equal to its coordinates
measure.center(Atom(0)) == Atom(0).coords  #>>> array([ True,  True,  True])
```
