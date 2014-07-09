# Trajectory analysis #

Pyvmd attempts to provides tools for trajectory analysis.

## Trajectory loading ##
Pyvmd provides tools for iterative analysis of trajectories bigger than available memory.

You can register any number of callbacks to be run on every snapshot of the trajectory.
Each callback will be provided a `status` object, which contains information about the trajectory.

### Examples ###
```python
from pyvmd.molecule import Molecule
from pyvmd.traj import Loader

def my_callback(status):
    print "The molecule:", status.molecule
    print "Total frame number", status.frame

mol = Molecule(0)  # Molecule for trajectory loading. It will be modified.
loader = Loader(mol, ['foo.dcd', 'bar.dcd'])
# Set the callbacks to be run on every frame
loader.add_callback(my_callback)
# Run the analysis
loader.run()
```

## Collectors ##
Some analysis are fairly common and for those, there are prepared objects which handles them, called collectors.

Currenty, there is only `RMSDCollector` which returns RMSD profiles for several selections. All selections are fitted
to the reference prior to measuring the RMSD.

### Examples ###
```python
from pyvmd.collectors import RMSDCollector
from pyvmd.molecule import Molecule
from pyvmd.traj import Loader

ref = Molecule(0)  # Reference molecule for RMSD
mol = Molecule(1)  # Molecule for trajectory loading. It will be modified.
rmsd = RMSDCollector(ref)
rmsd.add_selection('protein')
rmsd.add_selection('proint and backbone', 'backbone')  # Add selection under name 'backbone'

loader = Loader(mol, ['foo.dcd', 'bar.dcd'])
# Set the callbacks to be run on every frame
loader.add_collector(rmsd)
# Run the analysis
loader.run()

# Dataset contains the data
rmsd.dataset
# Write them into a file
rmsd.dataset.write('rmsd.dat')
```
