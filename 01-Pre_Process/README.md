![SeismoVLAB Logo](../Logo.png)

Pre-Process
===========

The main objective of the **Pre-Analysis** is to abstract the user from complicated manipulation in creating the input files. Therefore, the **Pre-Analysis** folder is employed to parse a pre-defined user's format into the Run-Analysis format. In Seismo-VLAB the **Pre-Analysis** corresponds to a set of python rutines that are in charge of generating input files, domain (mesh) partition, degree of freedom numbering, spatial variability of soil/structure properties, and parameter identification files. 

Further informatin can be obatained at:

* Official webSite: https://SeismoVLAB.github.io/SVL/
* Official documentation: https://SeismoVLAB.github.io/SVL/01-Pre-Process/docs/index.html

Start with Pre-Process
----------------------
The best way to start with the **Pre-Process** is to run some debugging cases on the official documentation webSite: https://SeismoVLAB.github.io/SVL/03-Validations/

In these examples, you will note that the Pre-Process' input file is defined by fields, for which each field starts and ends with a <span style="color:blue">\*TOKEN</span> and <span style="color:blue">\*END</span> keyword that defines the scope of the variables being defined. Inside each field, variables are defined with a <span style="color:blue">-variable</span> declaration followed by its corresponding instantiation, i.e, name, number, path. Each variable are separated by commas <span style="color:blue">“,”</span> which allows the format to be not only order insensitive, but also case insensitive. However, tags (identifiers) need to be declared at the begining of each line within the field. 


Running the Pre-Process
-----------------------
Once the model is created, running the **Pre-Process** is straightforward, just execute:

<pre>
python3 PreAnalysis.py 'path/to/input/file'
</pre>

The previous command will generate the command input line to be used in the Run-Process.

Files Description
=================

* **PreProcess.py**:
  This is the main python file. 
* **Parser.py**:
  Field from the input file are transformed into python dictionaries. 
* **Metis.py**:
  This python file runs the METIS - Serial Graph Partitioning by George Karypis, and gets the element indeces for each partition.
* **Numbering.py**:
  This python file assigns the degree of freedom numbering for each Point according to the User's numbering pattern.
* **PlaneWave.py**:
  This python routine creates the domain reduction input files for the homogeneous linear elastic half-space case.  
* **SeismoVLab.py**:
  This python file writes the output files for the Pre-Process format. In other words, transform all dictionaries information into text files. 
