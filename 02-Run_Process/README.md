![SeismoVLAB Logo](../Logo.png)

Run-Process
===========

The **Run-Process** is the main core of this software, and it is in charge of performing the desired analysis. The **Run-Porcess** essentially generates the finite element matrices, assembles the contribution of each elements, and solves the resulting systems of equations. The **Run-Porcess** folder is divided in several sub-folder which contains all the objects needed to create a finite element model and then perform an analysis. This internal folder structure is deliberately chosen so that finding/modifying/adding new classes are simple for any user. These folders at the same time are divided into three categories:

- **Geometry Module**: This module is conformed by `Node`, `Material`, `Section`, `Elements`, `Load`, and `Mesh` classes. On this level the mesh of the finite element model is generated, in which nodes are placed in space, elements are defined through node connectivity, materials are assigned to the elements, and loads are specified on nodes and elements, respectively. The mesh encapsulates all these objects, which is finally used by the `Assembler` object to generate the global mass, damping and stiffness matrices as well as the global force vector. 

- **Solution Module**: This module is conformed by `Analysis`, `Algorithm`, `Integrator`, `LinearSystem`, and `Assembler` classes. On this level an incremental analysis is performed for either static or dynamical purpose. The analysis uses the algorithm to take care of how the solution will be evolved for each time step. During each time step, the integrator combines the matrices and vectors using the assembler and gives this information to the solver. The solver finds the solution to the linear system (generated from the integrator) and returns it to the algorithm. The algorithm checks if convergence criteria are met, and continues to the next time step.

- **Input/Output Module**: This module is conformed by `Recorder` and `Parser` classes. In this level the Parser reads the domain and analysis input file(s) and generates the `Node`, `Element`, `Material`, `Load`, `Mesh` objects to define a finite element problem. It also generates the analysis and how the finite element problem is going to be solved. The Recorders stores the solution at `Node` or `Element` obtained in the Solution Module in files specified by the user.

Further information can be obtained at:

* [Official website](http://www.seismovlab.com)
* [GitHub repository](https://github.com/SeismoVLAB/SVL)
* [Official documentation](http://www.seismovlab.com/documentation/index.html)

Compiling the Run-Analysis
--------------------------
Installation of **Seismo-VLAB** on Linux/MacOSX requires to download `Eigen C++ library`, `MUMPS Library`, and `Pestc Library`. Also, python3 is needed along with libraries such as numpy, scipy, and matplotlib.

* The **Eigen C++ library** can be downloaded from this [website](http://eigen.tuxfamily.org/). This package needs to be unzip and its content move to `/usr/include/eigen`. 
* The **MUMPS library** can be downloaded from this [website](http://mumps.enseeiht.fr/). This package needs to be unzip and compiled at `/usr/include/mumps`.
* The **PETSc Library** can be downloaded at this [website](https://www.mcs.anl.gov/petsc/). This package needs to be unzip and compiled at `/usr/include/petsc`.

Assuming the previous libraries are successfully installed, then modify the `Makefile.mk` file such the previous path point to the right libraries:

```bash
EIGEN_DIR = /usr/include/eigen
PETSC_DIR = /usr/include/petsc
MUMPS_DIR = /usr/include/mumps
```

Also, make sure that libraries such as: *libscalapack-openmpi*, *libblacs-openmpi*, *liblapack*, *libblas*, and *libparmetis*, *libmetis*, *libptscotch*, *libptscotcherr* are also installed.

Finally, write in terminal:
```bash
make -s DEBUG=False
```
A detailed explanation on how to install **SVL** on Windows, MacOS, and Linux can be found in [this link.](http://seismovlab.com/documentation/linkInstallation.html)


Running the Run-Analysis
------------------------

Once **Seismo-VLAB** is compiled, write in a terminal (Linux/MacOSX) or command prompt (cmd in Windows),

```bash
./SeismoVLAB.exe --help
```

and the following message will be printed:

<pre>
                ░                                              
  ███████╗ ░░░░   ░░░░ ░░░░░░░ ░░░░ ██╗   ██╗ ██╗   ░░░  ░░░░  
  ██╔════╝ ░    ░ ░    ░  ░  ░ ░  ░ ██║   ██║ ██║  ░   ░ ░   ░ 
  ███████╗ ░░░  ░ ░░░░ ░  ░  ░ ░███╗██║   ██║ ██║  ░   ░ ░░░░  
  ╚════██║ ░    ░    ░ ░  ░  ░ ░╚═░╝╚██╗ ██╔╝ ██║  ░░░░░ ░   ░ 
  ███████║ ░░░░ ░ ░░░░ ░  ░  ░ ░░░░  ╚████╔╝  ███████╗ ░ ░░░░ 
  ╚══════╝                            ╚═══╝   ╚══════╝           
                                                                        
                  Seismo Virtual Laboratory                         
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
     -file : Name of the SeismoVLAB input file to be loaded.      
</pre>

Then, we can run a serial/parallel **Seismo-VLAB** in the following manner:

* For a single-core execution:
  ```bash
  ./SeismoVLAB.exe -dir '/path/to/Partition/folder' -file 'model.$.json'
  ```

* For a multiple-core execution:
  ```bash
  mpirun -np n ./SeismoVLAB.exe -dir '/path/to/Partition/folder' -file 'model.$.json'
  ```

The flag `-np = n` specifies that the number of processors are `n`, this requires the mesh and simulation files to be partitioned in such number of files. The latter explains the `.$.` token which internally is replaced by the processor number.

Folder Description
==================
* **01-Node**:
  This folder contains the `Node` and `Constraint` classes.
* **02-Materials**: 
  This folder contains the Material class separated into `01-Linear`, `02-NonLinear` and `03-Fiber`.
* **03-Sections**:
  This folder contains the `Section` class separated into `01-Plain` (one material) and `02-Fiber` (several materials).
* **04-Elements**:
  This folder contains the `Element` class. Here, solid and structural elements are defined. The integration class is defined as well.
* **05-Loads**:
  This folder contains the `Load` and `LoadCombo` classes.
* **06-Mesh**:
  This folder cotains the `Mesh` class. Node, Material, Element, Load Containers are stored in this class.
* **07-Assembler**:
    This folder contains the `Assembler` class. This object is in charge of generating the mass, stiffness, damping matrices as well as the force vector. 
* **08-Analysis**:
    This folder contains the `Analysis` class. Definition of `01-Static` and `02-Dynamic` analyses are provided in this folder.
* **09-Algorithms**:
  This folder contains the `Algorithm` class. Definition of `01-Linear` (linear) and `02-Newton` (nonlinear) algorithm are provided in this folder.   
* **10-Integrators**:
  This folder contains the `Integrator` class. Definition of `01-QuasiStatic`, `02-CentralDifference`, `03-Newmark`, and `04-Bathe` integrators are provided in this folder.
* **11-Solvers**:
  This folder contains the `LinearSystem` class. Definition of EigenSolver (serial), MumpsSolver (parallel) and PetscSolver (parallel) are provided in this folder.
* **12-Utilities**:
  This folder contains definition of several classes. `Damping`, `Parser`, and `Recorder` are defined in this folder.
