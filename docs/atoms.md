# Atoms #

The most used objects in VMD are probably atoms and selections. Pyvmd provides interface which emulates chemical
objects, such as atoms, residues, chains and segments.

## Common features ##
All atom based objects has common features. These are:
 * `molecule`
  - When creating an object, molecule can be defined. Default is top molecule.
  - Object has a `molecule` property which returns its molecule.
 * `frame`
  - When creating an object, frame can be defined. Default is molecule's active frame (NOW).
  - Object has a `frame` property which returns its frame or can be used to change the frame.
 * `atomsel`
  - Object has a `atomsel` property which returns respective `VMD.atomsel.atomsel` object if its interface is required.
 * Container methods
  - All objects except `Atom` has basic container features. Function `len()` returns size of the object in atoms,
    `in` returns whether atom belongs to a object and object also works as iterable over its atoms.
 * Comparison
  - All objects have defined equality and inequality operators.
 * Hashability
  - Objects `Atom` and `Residue` are hashable, so they can be used in sets or as keys in dictionaries.

### Examples ###
```python
from pyvmd.atoms import Selection, NOW

# Create selection
sel = Selection('protein and noh')
# Create selection defining molecule and frame
sel = Selection('protein and noh', molecule=Molecule(3), frame=17)

# Get selection's molecule
sel.molecule  #>>> Molecule(3)

# Get selection's frame
sel.frame  #>>> 17
# Set selection's frame
sel.frame = 25
# Set selection's frame to molecule's active frame
sel.frame = NOW

# Get number of atoms in selection
len(sel)  #>>> 2048

# Iterate through atoms in selection
from atom in sel:
    do_something()

# Check if atom is in selection
Atom(789) in sel  #>>> True
```

## Atom ##
Object representing the basic element of a structure.

### Atom properties table ###
| Property  | Type                   | VMD keyword |
|-----------|------------------------|-------------|
| index     | `int`                  | index       |
| x         | `float`                | x           |
| y         | `float`                | y           |
| z         | `float`                | z           |
| coords    | `numpy.array(x, y, z)` |             |
| name      | `str`                  | name        |
| type      | `str`                  | type        |
| element   | `str`                  | element     |
| beta      | `float`                | beta        |
| occupancy | `float`                | occupancy   |
| mass      | `float`                | mass        |
| charge    | `float`                | charge      |
| radius    | `float`                | radius      |
| residue   | `pyvmd.Residue`        |             |
| chain     | `pyvmd.Chain`          |             |
| segment   | `pyvmd.Segment`        |             |

### Examples ###
```python
from pyvmd.atoms import Atom

# Get atom using its index, the `index` value
atom = Atom(12000)

# Find atom using selection
atom = Atom.pick('residue 25 and name CA')
# If selection returns none or more than one atom, ValueError is raised
Atom.pick('none')  # ValueError
Atom.pick('all')  # ValueError

# Get atom's index
atom.index  #>>> 167

# Atom's coordinates
a.x  #>>> 5.3
a.y  #>>> 2.5
a.z  #>>> 17.89
# All three coordinates can be obtained at once in numpy array
a.coords  #>>> array([5.3, 2.5, 17.89])
# Setters are also available
a.x = 34
a.y = -90.5
a.z = 0.0
a.coords = (5.6, -9, 0.003)

# Iterate through bonded atoms
for bonded in atom.bonded:
    do_something()

# Get atom's residue
atom.residue  #>>> Residue(25)

# Get atom's chain
atom.chain  #>>> Chain('C')
# Move atom to a different chain
atom.chain = Chain('D')

# Get atom's segment
atom.segment  #>>> Segment('S')
# Move atom to a different segment
atom.segment = Segment('U')
```

## Residue ##
Object representing single residue.

### Residue properties table ###
| Property | Type                   | VMD keyword |
|----------|------------------------|-------------|
| index    | `int`                  | residue     |
| name     | `str`                  | resname     |
| number   | `int`                  | resid       |

### Examples ###
```python
from pyvmd.atoms import Residue

# Find residue using its index (`residue` value)
res = Residue(42)

# Get index of the residue, the `residue` value.
res.index  #>>> 42
# Get number of the residue, the `resid` value.
res.number  #>>> 157
# Change the number of the residue
res.number = 456
```

## Chains and segments ##
These objects represent single chain and segment, respectively.

### Examples ###
```python
from pyvmd.atoms import Chain, Segment

# Find chain using its name (`chain` value)
chain = Chain('A')

# Get name of the chain, the `chain` value.
chain.name  #>>> 'A'
# Change the name of the chain
chain.name = 'B'

# Find segment using its name (`segname` value)
segment = Segment('S')

# Get name of the segment, the `segname` value.
segment.name  #>>> 'S'
# Change the name of the segment
segment.name = 'U'
```

## Selections ##
Pyvmd also supports generic selections similar to `VMD.atomsel.atomsel`. Selection is automatically updated.

### Examples ###
```python
from pyvmd.atoms import Selection

# Make a selection
sel = Selection('protein and noh')

# Get selection text
sel.selection  #>>> 'protein and noh'

# Iterate through contacts between two selections
for atom1, atom2 in sel.contacts(other_sel, 3.0):
    do_something()

# Selections are automatically updated if frame is set to NOW
sel = Selection('ions within 5 of protein')
len(sel)  #>>> 10
sel.molecule.frame = 17
len(sel)  #>>> 12
```
