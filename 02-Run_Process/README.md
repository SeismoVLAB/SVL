![SeismoVLAB Logo](../Logo.png)

Run-Process
===========

The **Run-Process** is the main core of this software, and it is in charge of generating the finite element model that comes from the Pre-Analysis and perform the desired analysis. This essentially generates the finite element matrices, assembles the contribution of each elements, and solves the resulting systems of equations. The solution generated in such process is then recorded in separated files whose format is discussed in Recorder. In the Run-Analysis component, the software's folder is divided in several sub-folder which contains all the objects needed to create a finite element model and then perform an analysis. This internal folder structure is deliberately chosen so that finding/modifying/adding new classes are simple for any user. These folders at the same time are divided into three categories,

- **Geometry Module**: This module is conformed by Node, Material, Section, Elements, Load, and Mesh classes. On this level the mesh of the finite element model is generated, in which nodes are placed in space, elements are defined through node connectivity, materials are assigned to the elements, and loads are specified on nodes and elements, respectively. The mesh encapsulates all these objects, which is finally used by the assembler object to generate the global mass, damping and stiffness matrices as well as the global force vector. 

- **Solution Module**: This module is conformed by Analysis, Algorithm, Integrator, LinearSystem, and Assembler classes. On this level an incremental analysis is performed for either static or dynamical purpose. The analysis uses the algorithm to take care of how the solution will be evolved for each time step. During each time step, the integrator combines the matrices and vectors using the assembler and gives this information to the solver. The solver finds the solution to the linear system (generated from the integrator) and returns it to the algorithm. The algorithm checks if convergence criteria are met, and continues to the next time step.

- **Input/Output Module**: This module is conformed by Recorder and Parser classes. In this level the Parser reads the domain and analysis input file(s) and generates the Node, Element, Material, Load, Mesh objects to define a finite element problem. It also generates the analysis and how the finite element problem is going to be solved. The Recorders stores the solution at Node or Element obtained in the Solution Module in files specified by the user.

Further informatin can be obatained at:

* Official webSite: https://SeismoVLAB.github.io/SVL/
* Official documentation: https://SeismoVLAB.github.io/SVL/02-Run-Process/docs/index.html

Compiling Run-Analysis
----------------------
Installation of **Seismo-VLab** on Linux/MacOSX requires to download `Eigen C++ library`, `MUMPS Library`, and `Pestc Library`. Also, python3 is needed along with libraries such as numpy, scipy, and matplotlib.

* The **Eigen C++ library** can be downloaded from the website http://eigen.tuxfamily.org/. This package needs to be unzip and its content move (for instnce) to `/usr/include/Eigen`. 
* The **MUMPS library** can be downloaded from the website http://mumps.enseeiht.fr/. This package needs to be unzip and compiled (for instnce) at `/usr/include/Mumps`.
* The **Pestc Library** library can be downloaded at the website https://www.mcs.anl.gov/petsc/. This package needs to be unzip and compiled (for instnce) at `/usr/include/Petsc`.

Assuming the previous libraries are successfully installed, then modify the `Makefile.inc` file such the previous path point to the right libraries:

<pre>
EIGEN_DIR = /usr/include/Eigen
PETSC_DIR = /usr/include/Petsc
MUMPS_DIR = /usr/include/Mumps
</pre>

Also, check that : `MPI_DIR, METIS_DIR, SCOTCH_DIR` have the correct path where these libraries are installed. Make sure libraries such as: *-lscalapack-openmpi -lblacs-openmpi -llapack -lblas* and *-lparmetis -lmetis -lptesmumps -lptscotch -lptscotcherr* are also installed.

Finally, write in terminal:
<pre>
make -s DEBUG=False
</pre>

Some libraries can be easily installed using repositories:
* **LAPACK :** sudo apt-get install liblapack-dev
* **ATLAS  :** sudo apt-get install libatlas-dev libatlas-base-dev
* **BLAS   :** sudo apt-get install libblas-dev libblas-common libblacs-mpi-dev
* **METIS  :** sudo apt-get install metis libmetis-dev libparmetis-dev
* **SCOTCH :** sudo apt-get install libscotch-dev libptscotch-dev libscotchmetis-dev libscotchparmetis-dev

Running the Run-Analysis
------------------------

Once **Seismo-VLab** is compiled, write in a terminal window (Linux/MacOSX) or command prompt (cmd in Windows),

<pre>
./SeismoVLab -help
</pre>

and the following message will be printed:

<pre>
       _____      _                    _    ____          __            
      / ___/___  (_)________ ___  ____| |  / / /   ____  / /_           
      \\__ \\/ _ \\/ / ___/ __ `__ \\/ __ \\ | / / /   / __ `/ __ \\    
     ___/ /  __/ /__  / / / / / / /_/ / |/ / /___/ /_/ / /_/ /          
    /____/\\___/_/____/_/ /_/ /_/\\____/|___/_____/\\__,_/_.___/          
                                                                        
               (Seismo)s (V)irtual (Lab)oratory                         
        Module for Serial and Parallel analysis of seismic              
   wave propagation and soil-structure interaction simulation           
   Copyright (C) 2020, The California Institute of Technology 
                     All Rights Reserved.                               
                                                                        
 Written by:                                         
   Danilo S. Kusanovic (dkusanov@caltech.edu)                           
   Elnaz E. Seylabi    (elnaze@unr.edu)                              
                                                                        
 Supervised by:                                      
   Domniki M. Asimaki  (domniki@caltech.edu)                            
                                                                        
 COMMAND LINE FLAGS TO PROVIDE:                                  
     -dir  : Location of the working directory.                  
     -file : Name of the SeismoVLab input file to be loaded.      
</pre>

Then, we can run a serial/parallel **Seismo-VLab** in the following manner:

* For a single-core execution,
<pre>
./SeismoVLab -dir /path/to/simulation/directory -file model.$.svl
</pre>

* For a multiple-core execution,
<pre>
mpirun -np n ./SeismoVLab -dir /path/to/simulation/directory -file model.$.svl
</pre>

The flag <span style="color:blue">`-np = n`</span> specifies that the number of processors are <span style="color:blue">`n`</span>, this requieres the mesh and simulation files to be partitioned in such number of files. The latter explains the <span style="color:blue">`.$.`</span> token which internally is replaced by the processor number.

Folder Description
==================
* **01-Node**:
  This folder contains the Node and Constraint classes.
* **02-Materials**: 
  This folder contains the Material class separated into Linear and NonLinear.
* **03-Sections**:
  This folder contains the Section class separated into Plain (one material) and Fiber (several materials).
* **04-Elements**:
  This folder contains the Element class. Here, solid and structural elements are defined. The integration class is defined as well.
* **05-Loads**:
  This folder contains the Load and LoadCombo classes.
* **06-Mesh**:
  This folder cotains the Mesh class. Node, Material, Element, Load Containers are defined in this class.
* **07-Assembler**:
    This folder contains the Assembler class. This object is in charge of generating the Mass, Stiffness, damping matrices as well as the force vector. 
* **08-Analysis**:
    This folder contains the Analysis class. Definition of Static and Dynamic analyses are provided in this folder.
* **09-Algorithms**:
  This folder contains the Algorithm class. Definition of Linear (linear) and NewtonRaphson (nonlinear) algorithm are provided in this folder.   
* **10-Integrators**:
  This folder contains the Integrator class. Definition of QuasiStatic, CentralDifference, NewmarkBeta, and Bathe integrators are provided in this folder.
* **11-Solvers**:
  This folder contains the LunearSystem class. Definition of EigenSoler (serial), MumpsSolver (parallel) and PetscSolver (parallel) are provided in this folder.
* **12-Utilities**:
  This folder contains definition of several classes. Damping, Parser, and Recorder are defined in this folder.
