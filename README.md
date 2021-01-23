![SeismoVLAB Logo](Logo.png)

**Seismo-VLAB** is a simple, fast, and extendable C++ multi-platform finite element software designed to optimize large-scale simulations of dynamic, nonlinear soil-structure interaction (SSI) problems. The software is intented to be used by academics researches in structural and geothecnical field but also for industrials, students, etc.

* Official website: http://www.seismovlab.com
* GitHub repository: https://github.com/SeismoVLAB/SVL
* Official documentation: http://www.seismovlab.com/documentation/index.html

What can Seismo-VLAB do?
------------------------
With **Seismo-VLAB** you can solve:

* Linear and Nonlinear wave propagation problems in shallow crust.
* Linear and Nonlinear soil-structure interaction problems.
* Standard mechanics-based nonlinear structural dynamic problems.

Visit the gallery for examples: http://www.seismovlab.com/gallery.html

Installing Seismo-VLAB
----------------------
Installation of **Seismo-VLab** (Pre-Process) on Linux/MacOSX/Windows requires a `python3` environment and the following libraries:

* Numpy
* Scipy
* Matplotlib
* VTK
* JSON

Installation of **Seismo-VLab** (Run-Process) on Linux/MacOSX/Windows requires to download `Eigen C++ library`, `MUMPS Library`, and `Pestc Library`. Also, python3 is needed along with libraries such as numpy, scipy, and matplotlib.

* The **Eigen C++ library** can be downloaded from the website http://eigen.tuxfamily.org/. This package needs to be unzip and its content move (for instnce) to `/usr/include/Eigen`. 
* The **MUMPS library** can be downloaded from the website http://mumps.enseeiht.fr/. This package needs to be unzip and compiled (for instnce) at `/usr/include/Mumps`.
* The **Pestc Library** library can be downloaded at the website https://www.mcs.anl.gov/petsc/. This package needs to be unzip and compiled (for instnce) at `/usr/include/Petsc`.

Assuming the previous libraries are successfully installed, then modify the `Makefile.inc` file such the previous path point to the right libraries:

```makefile
EIGEN_DIR = /usr/include/Eigen
PETSC_DIR = /usr/include/Petsc
MUMPS_DIR = /usr/include/Mumps
```

Also, check that : `MPI_DIR, METIS_DIR, SCOTCH_DIR` have the correct path where these libraries are installed. Make sure libraries such as: *-lscalapack-openmpi -lblacs-openmpi -llapack -lblas* and *-lparmetis -lmetis -lptesmumps -lptscotch -lptscotcherr* are also installed.

Finally, write in terminal:
```bash
make -s DEBUG=False
```

Some libraries can be easily installed using repositories:
* **LAPACK :** sudo apt-get install liblapack-dev
* **ATLAS  :** sudo apt-get install libatlas-dev libatlas-base-dev
* **BLAS   :** sudo apt-get install libblas-dev libblas-common libblacs-mpi-dev
* **METIS  :** sudo apt-get install metis libmetis-dev libparmetis-dev
* **SCOTCH :** sudo apt-get install libscotch-dev libptscotch-dev libscotchmetis-dev libscotchparmetis-dev

License
=======

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Seismo-VLAB is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Seismo-VLAB is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details http://www.gnu.org/licenses.

<!---
Citation
========
To cite Seismo-VLAB, please use :

Kusanovic Danilo, Seylabi Elnaz, Kottke Albert, and Asimaki Domniki (2020). Seismo-Vlab: A parallel object-oriented platform for reliable nonlinear seismic wave propagation and soil-structure interaction simulation. *Computers and Geotechnics*. [![DOI](https://img.shields.io/badge/DOI-10.1016/j.cma.2009.08.016-green.svg)](https://doi.org/10.1016/j.cma.2009.08.016)

```
@article{Kusanovic2020SeismoVLab,
title   = {Seismo-VLAB: A parallel object-oriented platform for reliable nonlinear seismic wave propagation and soil-structure interaction simulation.},
author  = {Kusanovic Danilo and Seylabi Elnaz and Kottke Albert and Asimaki Domniki},
journal = {To be submitted to Computer Methods in Applied Mechanics and Engineering},
volume  = {},
number  = {},
pages   = {},
year    = {2020},
issn    = {},
doi     = {},
url     = {}
}
```

Kusanovic Danilo, Seylabi Elnaz, and Asimaki Domniki (2021). Seismo-VLAB: A parallel C++ finite element software for structural and soil mechanics. *The Journal of Open Source Software*. [![DOI](https://img.shields.io/badge/DOI-10.1016/j.cma.2009.08.016-green.svg)](https://doi.org/10.1016/j.cma.2009.08.016)

```
@article{Kusanovic2021SeismoVLab,
title   = {Seismo-VLAB: A parallel C++ finite element software for structural and soil mechanics.},
author  = {Kusanovic Danilo and Seylabi Elnaz and Asimaki Domniki},
journal = {To be submitted to SoftwareX},
volume  = {},
number  = {},
pages   = {},
year    = {2021},
issn    = {},
doi     = {},
url     = {}
}
```
--->
