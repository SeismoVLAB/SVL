![SeismoVLAB Logo](../Logo.png)

Pre-Process
===========

The **Pre-Analysis** is a (user friendly) Python interface designed to handle large models in different formats. In this regard, the **Pre-Analysis** allows the user to parse formats from other software such as [ETABS](https://www.csiamerica.com/products/etabs), [SAP2000](https://www.csiamerica.com/products/sap2000), and [VTK](https://vtk.org/). The main **Pre-Analysis** task is orchestrated by two python files:`SeismoVLAB.py` that loads all required modules to generate a model; and `Definitions.py` that contains the data structure (dictionaries) that stores the finite element information. The Python routines also provided in **SVL** handles the degree-of-freedom numbering in `Numberer.py`, the domain partition in `Partition.py`, soil spatial variability in `RandomField.py`, and domain reduction forces in `PlaneWave.py`. Other features can be incorporated to meet the user's need as well.

Installation of **Seismo-VLAB** (Pre-Process) on Linux/MacOS/Windows requires a `python3` environment and the following libraries:

* [Numpy](https://numpy.org/)
* [Scipy](https://www.scipy.org/)
* [Matplotlib](https://matplotlib.org/)
* [JSON](https://www.json.org/json-en.html)

Further information can be obtained at:

* [Official website](http://www.seismovlab.com)
* [GitHub repository](https://github.com/SeismoVLAB/SVL)
* [Official documentation](http://www.seismovlab.com/documentation/index.html)

Starting with the Pre-Process
-----------------------------
The best way to start with the **Pre-Process** is to run some debugging cases on the official [documentation website.](https://github.com/SeismoVLAB/SVL/tree/master/03-Validations)

It is convenient to set the **PYTHONPATH** environment variable to have the **Pre-Process** folder path as:

```Python
export PYTHONPATH="/path/to/SVL/01-Pre_Process"
```

The later will allow to import the the **Pre-Process** python module as follows:

```Python
from Core import SeismoVLAB as SVL
```

Running the Pre-Process
-----------------------
Once the model is created, running the **Pre-Process** is straightforward, just execute:

```Python
python3 'path/to/input/file.py'
```

The previous execution will generate the `json` files required to run the **Run-Process**.

Files Description
=================

* **Core**:
  * `Definitions.py`: Declares the main dictionaries used to store user's input options and model information.
  * `Utilities.py`: Provides with useful functions to print variables, save model, clear variables and more.
  * `Partition.py`: Generates the domain partition using [Metis](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview)
  * `Numberer.py`: This python file assigns the degree of freedom numbering for each Point according to the User's numbering pattern.
  * `Outputs.py`: Writes the **Run-Analysis** input files in *.json format
  * `RandomField.py`: Applies a random field to a background finite element model
  * `PlaneWave.py`: This python routine creates the domain reduction input files for the homogeneous linear elastic half-space case.
  * `SeismoVLAB.py`: Main python file that imports all required modules.
* **Method**
  * `Attach.py`: Functions provided to populate the finite element model
  * `Remove.py`: Functions provided to delete finite element components
  * `Builder.py`: Provides with functions to create simple geometries in 1D, 2D, and 3D
  * `Display.py`: Creates a 3D visualization of the model compatible with [ParaView](https://www.paraview.org/)
  * `Compute.py`: Computes the kinematic constraints for diaphragm, rigid body, and rigid link.
* **Parser**
  * `Formats.py`: Parse a file provided with the format in which is written
  * `ETABS.py`: Parses a [ETABS](https://www.csiamerica.com/products/etabs) input file with .*e2k extension 
  * `SAP2000.py`: Parses a [SAP2000](https://www.csiamerica.com/products/sap2000) input file with .*s2k extension 
