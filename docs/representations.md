# Representations #

Molecule representations specify display options of molecules in graphical window.
Module [`pyvmd.representations`](../pyvmd/representations.py) provides interface for representations.

## Representation objects ##
Object `Representation` contains all details about single representation. All representations are referenced by their name, which is usually in the form `rep<number>`, but unlike index it doesn't change.

### Properties ###
Representation objects have following properties
 * `molecule` - Returns represented molecule
 * `name` - Returns name of the representation
 * `style` - Returns Style object
 * `color` - Returns Color object
 * `selection` - Returns or sets representation selection
 * `visible` - Returns or sets visibility
 * `update_selection` - Returns or sets automatic selection update with change of frame
 * `update_color` - Returns or sets automatic color update with change of frame

### Examples ###
```python
from pyvmd.molecules import Molecule
from pyvmd.representations import Representation

# Get representation from top molecule using its name
rep = Representation('rep4')
# Get representation from particular molecule
rep = Representation('rep7', Molecule(15))

# Get molecule
rep.molecule  #>>> Molecule(15)
# Get name
rep.name  #>>> 'rep7'

# Get selection
rep.selection  #>>> Selection('protein', Molecule(15))
# Change selection
rep.selection = Selection('not water')

# Get visibility
rep.visible  #>>> True
# Hide
rep.visible = False
# Show
rep.visible = True

# Create new representation for top molecule
rep = Representations.create()
# Create new representation for particular molecule
rep = Representations.create(Molecule(15))

# Delete representation
rep.delete()
```

## Style objects ##
Styles in VMD are composed from the drawing method and its parameters. Because of this complexity, styles have their own objects.

All drawing methods which can be available in VMD are defined as constants.
Names of the constants are based on the name of the methods in VMD and have prefix `DRAW_`.
Depending on compile options, not all of them may be available in all VMD binaries.

### Examples ###
```python
from pyvmd.molecules import Molecule
from pyvmd.representations import Representation, DRAW_CPK

rep = Representation('rep4')

# Get drawing method
rep.style.method  #>>> DrawingMethod('Lines', ...)
# Get parameters of the drawing method
rep.style.get_parameters()  #>>> {'size': 1.0}

# Change drawing method
rep.style.method = DRAW_CPK
# List available parameters for the method
rep.style.method.parameters  #>>> ('size', 'bond_radius', 'resolution', 'bond_resolution')
# Modify drawing parameters
rep.style.set_parameters(resolution=50, bond_resolution=50)
```

## Color objects ##
Colors are similar to the styles. They are also composed from the drawing method and its parameters, and thus have API similar to the styles.

All coloring methods which can be available in VMD are defined as constants.
Names of the constants are based on the name of the methods in VMD and have prefix `COLOR_`.

Currently only the color and volume methods have any parameters.

### Examples ###
```python
from pyvmd.molecules import Molecule
from pyvmd.representations import Representation, COLOR_COLOR

rep = Representation('rep4')

# Get coloring method
rep.color.method  #>>> ColoringMethod('Name', ...)
# Get parameters of the coloring method
rep.color.get_parameters()  #>>> {}

# Change coloring method
rep.color.method = COLOR_COLOR
# List available parameters for the method
rep.color.method.parameters  #>>> ('color', )
# Modify coloring parameters
rep.color.set_parameters(color=4)
```

## Molecule representations ##

An interface to manage all molecule representations is also provided.

### Examples ###
```python
from pyvmd.molecules import Molecule

mol = Molecule(15)

# Get number of representations
len(mol.representations)  #>>> 5

# Iterate through representations
for rep in mol.representations:
    do_something(rep)

# Get representation by its name
mol.representations['rep0']  #>>> Representation('rep0', ...)
# Get representation by its index
mol.representations[3]  #>>> Representation('rep3', ...)
# Get representation by negative index
mol.representations[-1]  #>>> Representation('rep6', ...)

# Get slices, e.g. list representations backwards
mol.representations[::-1]  #>>> [Representation('rep6', ...), ...]
```
