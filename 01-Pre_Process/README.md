![SeismoVLAB Logo](../Logo.png)

Pre-Process
===========

The **Pre-Analysis** is a Python (user friendly) interface designed to handle large models in different formats. In this regard, the **Pre-Analysis** allows the user to parse formats from other software (ETABS, SAP, and GMSH) as well as it provides with several routines allows the user to create a model manually using such interface. The main **Pre-Analysis** task is orchestrated by two python files: <span style="color:blue">SeismoVLAB.py</span> that loads all required modules to generate a model; and <span style="color:blue">Definitions.py</span> that contains the data structure (dictionaries) that stores the finite element information. The Python routines also provided in **SVL** handles the degree-of-freedom numbering in <span style="color:blue">Numberer.py</span>, the domain partition in <span style="color:blue">Partition.py</span>, soil spatial variability in <span style="color:blue">RandomField.py</span>, and domain reduction forces in <span style="color:blue">PlaneWave.py</span>. Other features can be incorporated to meet the user's need as well.

Installation of **Seismo-VLab** (Pre-Process) on Linux/MacOSX/Windows requires a `python3` environment and the following libraries:

* Numpy
* Scipy
* Matplotlib
* VTK
* JSON

Further information can be obatained at:

* Official website: http://www.seismovlab.com
* GitHub repository: https://github.com/SeismoVLAB/SVL
* Official documentation: http://www.seismovlab.com/documentation/index.html

Start with Pre-Process
----------------------
The best way to start with the **Pre-Process** is to run some debugging cases on the official documentation webSite: https://github.com/SeismoVLAB/SVL/tree/master/03-Validations

In these examples, you will note that the Pre-Process' Python module is imported using:

<pre>
from Core import SeismoVLAB as SVL
</pre> 

However, the **PYTHONPATH** needs to be set to the **Pre-Process** Python module address:

<pre>
export PYTHONPATH="/path/to/SeismoVLAB/01-Pre_Process"
</pre>

Running the Pre-Process
-----------------------
Once the model is created, running the **Pre-Process** is straightforward, just execute:

<pre>
python3 'path/to/input/file.py'
</pre>

The previous command will generate the command input line to be used in the **Run-Process**.

Files Description
=================

* **Core**:
  * <span style="color:blue">Definitions.py</span>: Declares the main dictionaries used to store user's input options and model information.
  * <span style="color:blue">Utilities.py</span>: Provides with useful functions to print variables, save model, clear variables and more.
  * <span style="color:blue">Partition.py</span>: Generates the domain partition
  * <span style="color:blue">Numberer.py</span>: This python file assigns the degree of freedom numbering for each Point according to the User's numbering pattern.
  * <span style="color:blue">Outputs.py</span>: Writes the **Run-Analysis** input files in *.svl or *.json format
  * <span style="color:blue">RandomField.py</span>: Applies a random field to a background finite element model
  * <span style="color:blue">PlaneWave.py</span>: This python routine creates the domain reduction input files for the homogeneous linear elastic half-space case.
  * <span style="color:blue">SeismoVLAB.py</span>: Main python file tha imports all required modules.
* **Method**
  * <span style="color:blue">Attach.py</span>: Functions provided to populate the the finite element model
  * <span style="color:blue">Remove.py</span>: Functions provided to delete the finite element model
  * <span style="color:blue">Builder.py</span>: Provides vith functions to create simple geometries in 1D, 2D, and 3D
  * <span style="color:blue">Display.py</span>: Creates a 3D visualization of the model defined
  * <span style="color:blue">Compute.py</span>: Computes the kinematic constraints for diaphrag, rigid body, and rigid link.
* **Parser**
  * <span style="color:blue">Formats.py</span>: Parse a file provided with the format in which is written
  * <span style="color:blue">GMSH.py</span>: Parses a gmsh input file with .*mesh (INRIA) extension 
  * <span style="color:blue">ETABS.py</span>: Parses a ETABS input file with .*e2k extension 
  * <span style="color:blue">SAP2000.py</span>: Parses a SAP input file with .*s2k (INRIA) extension 
