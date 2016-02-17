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

Callbacks can receive any additional arguments that are passed upon their registration.

```python
def my_callback(step, extra_arg, extra_keyword):
    print step.frame, extra_arg, extra_keyword

analyzer.add_callback(my_callback, extra_arg, extra_keyword=extra_value)
```


## Datasets and collectors ##
Some analyses like RMSD or atom distances are done on regular basis.
For such analyses pyvmd provides tools which extract the data.
There are two basic objects - datasets and collectors.

 * Datasets are containers for data collected during the analysis.
   When the analysis is finished dataset provides data either in numpy array or it can write it to the file-like object.
 * Collectors performs the specified analysis operation at each snapshot of the trajectory and stores the result into dataset.

## List of collectors ##
Every collector takes optinal argument `name`, which is used in the column header if dataset writes data into a file.
If the `name` isn't specified it is generated in form `data12345`.
 * `XCoordCollector(selection, name=None)` - Collects X coordinate of an atom or geometric center of selection.
 * `YCoordCollector(selection, name=None)` - Collects Y coordinate of an atom or geometric center of selection.
 * `ZCoordCollector(selection, name=None)` - Collects Z coordinate of an atom or geometric center of selection.
 * `DistanceCollector(selection1, selection2, name=None)` - Collects distance between two atoms or geometric centers of selections.
 * `AngleCollector(selection1, selection2, selection3, name=None)` -
   Collects angle between three atoms or geometric centers of selections.
 * `DihedralCollector(selection1, selection2, selection3, selection4, name=None)` -
   Collects dihedral of improper dihedral angle of four atoms or geometric centers of selections.
 * `RMSDCollector(selection, reference, name=None)` - Collects RMSD between selection and reference.
   The selection is fitted to the reference prior to measuring the RMSD.

### Examples ###
```python
from pyvmd import collectors
from pyvmd.analyzer import Analyzer
from pyvmd.atoms import Selection
from pyvmd.datasets import DataSet
from pyvmd.molecule import Molecule

dset = DataSet()
# Register collector for distance between protein and LIG residue.
dset.add_collector(collectors.DistanceCollector('protein', 'resname LIG')
ref = Molecule(0)  # Reference molecule for RMSD
# Register collector for protein backbone RMSD named 'backbone'.
dset.add_collector(
    collectors.RMSDCollector('protein and backbone',
                             Selection('protein and backbone', ref),
                             name='backbone')
)

mol = Molecule(1)  # Molecule for analysis. It will be modified.
analyzer = Analyzer(mol, ['foo.dcd', 'bar.dcd'])
# Set the callbacks to be run on every frame
analyzer.add_dataset(dset)
# Run the analysis
analyzer.analyze()

# Get the data in numpy array
dset.data  #>>> array([[0, 12.45, 0.0], ...])
# Write data into a file
dset.write('rmsd.dat')
# Write data into a file descriptor
import sys
dset.write(sys.stdout)
```
