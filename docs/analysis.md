# Trajectory analysis #

Pyvmd attempts to provides rather flexible tools for trajectory analysis.

## Analyzer ##
Analyzer is the basic object in trajectory analysis.
It loads trajectory files iteratively, thus allowing analysis for trajectories bigger than available memory.

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

## Datasets and collectors ##
Some kinds of analysis are fairly common and for those pyvmd provides prepared tools.

Datasets are containers for data collected during the analysis.

Collectors perform the single analysis operations and provide data for datasets.

Currently, there is only `RMSDCollector` which computes RMSD between selection and the reference.
The selection is fitted to the reference prior to measuring the RMSD.

### Examples ###
```python
from pyvmd.analyzer import Analyzer
from pyvmd.atoms import Selection
from pyvmd.collectors import RMSDCollector
from pyvmd.datasets import DataSet
from pyvmd.molecule import Molecule

ref = Molecule(0)  # Reference molecule for RMSD
dset = DataSet()
# Register collector for RMSD between the 'protein' and 'protein' in reference.
dset.add_collector(RMSDCollector('protein', Selection('protein', ref)))
# Register collector for RMSD between the 'backbone' and 'backbone' in reference under name 'backbone'.
dset.add_collector(RMSDCollector('protein', Selection('protein', ref), name='backbone'))

mol = Molecule(1)  # Molecule for analysis. It will be modified.
analyzer = Analyzer(mol, ['foo.dcd', 'bar.dcd'])
# Set the callbacks to be run on every frame
analyzer.add_dataset(dset)
# Run the analysis
analyzer.analyze()

# Get the data in numpy array
dset.data
# Write data into a file
dset.write('rmsd.dat')
# Write data into a file descriptor
import sys
dset.write(sys.stdout)
```
