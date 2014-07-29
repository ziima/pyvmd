# Trajectory analysis #

Pyvmd attempts to provides rather flexible tools for trajectory analysis.

## Analyzer ##
Analyzer is the basic object in trajectory analysis. It loads trajectory files iteratively, thus allowing analysis for
trajectories bigger than available memory.

You can register any number of callbacks to be run on every snapshot of the trajectory.
Each callback will be provided a `step` object, which contains information about the ongoing analysis.

### Examples ###
```python
from pyvmd.analyzer import Analyzer
from pyvmd.molecule import Molecule

def my_callback(step):
    print "The molecule:", step.molecule
    print "Total frame number", step.frame

mol = Molecule(0)  # Molecule for analysis. It will be modified.
analyzer = Analyzer(mol, ['foo.dcd', 'bar.dcd'])
# Set the callbacks to be run on every frame
analyzer.add_callback(my_callback)
# Run the analysis
analyzer.analyze()
```

## Collectors ##
Some analysis are fairly common and for those, there are prepared objects which handles them, called collectors.

Currenty, there is only `RMSDCollector` which returns RMSD profiles for several selections. All selections are fitted
to the reference prior to measuring the RMSD.

### Examples ###
```python
from pyvmd.analyzer import Analyzer
from pyvmd.collectors import RMSDCollector
from pyvmd.molecule import Molecule

ref = Molecule(0)  # Reference molecule for RMSD
mol = Molecule(1)  # Molecule for analysis. It will be modified.
rmsd = RMSDCollector(ref)
rmsd.add_selection('protein')
rmsd.add_selection('proint and backbone', 'backbone')  # Add selection under name 'backbone'

analyzer = Analyzer(mol, ['foo.dcd', 'bar.dcd'])
# Set the callbacks to be run on every frame
analyzer.add_collector(rmsd)
# Run the analysis
analyzer.analyze()

# Dataset contains the data
rmsd.dataset
# Write them into a file
rmsd.dataset.write('rmsd.dat')
```
