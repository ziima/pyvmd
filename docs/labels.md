# Labels #

Labels allow display of specific atoms, bonds, angles and dihedrals in graphical window.
Module [`pyvmd.labels`](../pyvmd/labels.py) provides interface for labels.

## Label objects ##
There are 4 different objects for atoms, bond, angles and dihedrals which are `AtomLabel`, `BondLabel`, `AngleLabel` and `DihedralLabel` respectively.
All of them are proxy object for labels in VMD.

### Label creation ###
New labels are created by `create` class method. Creation of the existing label doesn't raise an error.

### Properties ###
Label objects have following properties
 * `category` - Returns one of constants `ATOM`, `BOND`, `ANGLE`, `DIHEDRAL`.
 * `atoms` - Tuple of atoms used to create the label
 * `visible` - Gets and sets whether label is visible in graphical window.

Specific properties
 * `AtomLabel.atom` - Returns atom which is labeled
 * `BondLabel.distances` - Returns distances between atoms throughout the trajectory
 * `AngleLabel.angles` - Returns angles of the atoms throughout the trajectory
 * `DihedralLabel.dihedrals` - Returns dihedral angles of the atoms throughout the trajectory

### Examples ###
```python
from pyvmd.atoms import Atom
from pyvmd.labels import AtomLabel, BondLabel, AngleLabel, DihedralLabel

# Create new label
atom_label = AtomLabel.create(Atom(0))

# Check label is visible
atom_label.visible  #>>> True

# Hide the label
atom_label.visible = False

# Show the label
atom_label.visible = True

# Delete the label
atom_label.delete()

# Create a bond label
bond_label = BondLabel.create(Atom(0), Atom(1))

# Get distances between atoms
bond_label.distances  #>>> [1.569, 1.678, ...]

# Labels supports equality operators
bond_label == BondLabel(Atom(0), Atom(1))  #>>> True
bond_label != BondLabel(Atom(0), Atom(1))  #>>> False
bond_label == BondLabel.create(Atom(0), Atom(2))  #>>> False
bond_label != BondLabel.create(Atom(0), Atom(2))  #>>> True
```

## Label managers ##
There are 4 managers to manager all label of given category - `ATOM_LABELS`, `BOND_LABELS`, `ANGLE_LABELS`, `DIHEDRAL_LABELS`.
All of them have simple container API.

### Examples ###
```python
from pyvmd.labels import ATOM_LABELS, BOND_LABELS, ANGLE_LABELS, DIHEDRAL_LABELS

# Get number of labels in each category
len(ATOM_LABELS)  #>>> 4
len(BOND_LABELS)  #>>> 3
len(ANGLE_LABELS)  #>>> 1
len(DIHEDRAL_LABELS)  #>>> 0

# Iterate through all labels
for label in ATOM_LABELS:
    do_something()
```
