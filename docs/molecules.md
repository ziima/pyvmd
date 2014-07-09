# Molecules #

One of the most basic objects in VMD are molecules. Module [`pyvmd.molecules`](../pyvmd/molecules.py) provides interface
for molecule operations.

## Molecule ##
Similar to VMD's `Molecule.Molecule` Molecule object is a proxy for molecule loaded into VMD.

### Loading files ###
Method `Molecule.load` has only one required argument which is a filename. But it has several other optional arguments:
 * `filetype` - Format of the file. By default it is determined from filename extension, but it can be defined manually.
    Preferably using one of `FORMAT_*` constants in `molecule` module. If constant is not available string defining the
    filetype can be used.
 * `start` - First frame to be read. Default is first.
 * `stop` - Last frame to be read. Default is last.
 * `step` - Load only every step'th frame. Default is every frame.
 * `wait` - Whether to wait until the file is completely read. Default is True.
 * `volsets` - Volumetric data sets to be read.

### Examples ###
```python
from pyvmd.molecules import Molecule

# Create Molecule for existing molecule
mol = Molecule(3)  # molid is required

# Create new empty molecule
mol = Molecule.create()
# Load a file.
mol.load('my_structure.psf')
# Load a trajectory.
mol.load('my_structure.dcd')

# Get active frame
mol.frame  #>>> 15
# Set active frame
mol.frame = 1

# Get molecule's name
mol.name  #>>> 'molecule'
# Rename the molecule
mol.name = 'My precious'

# If you need a missing interface, `molecule` property returns instance of
# VMD's `Molecule.Molecule` object.
vmd_mol = mol.molecule

# Delete the molecule
mol.delete()

# Molecules supports equality operators
Molecule(0) == Molecule(0)  #>>> True
Molecule(0) != Molecule(0)  #>>> False
Molecule(0) == Molecule(1)  #>>> False
Molecule(0) != Molecule(1)  #>>> True
```

## Molecule's frames ##
Interface for molecule's frames is grouped in `frames` descriptor.

### Examples ###
```python
# Get number of frames
len(mol.frames)  #>>> 16

# Frames can be used as iterable
for frame in mol.frames:
    do_something()

# Duplicate active frame
mol.frames.copy()
# Duplicate frame 4
mol.frames.copy(4)

# Delete second frame
del mol.frames[1]
# Delete all but last 4 frames
del mol.frames[:-4]
# Delete every third frame
del mol.frames[::3]
```

## Application molecules ##
The interface to manipulate all molecules in application is also present. It is available through
`pyvmd.molecules.MOLECULES` and has a usual container-like interface.

### Examples ###
```python
from pyvmd.molecules import MOLECULES

# Get number of molecules
len(MOLECULES)  #>>> 2

# Iterate through all molecules
for mol in MOLECULES:
    do_something()

# Get top molecule
MOLECULES.top  #>>> Molecule(0)
# Set top molecule
MOLECULES.top = Molecule(4)

# Get molecule by molid
mol = MOLECULES[0]
# Get molecule by name. This may not be reliable if multiple molecules with the same name is present.
mol = MOLECULES['My precious']

# Molecules can be also deleted using indexes
del MOLECULES[0]
del MOLECULES['My precious']

# Check if molecule still exists
mol in MOLECULES  #>>> False
```
