![SeismoVLAB Logo](Logo.png)

**Seismo-VLAB** (a.k.a **SVL**) is a simple, fast, and extendable C++ finite element software designed to optimize meso-scale simulations of linear and nonlinear wave-propagation and soil-structure interaction. **SVL** is intended not only to be used by researchers in structural and geothecnical engineering, but also in industries, laboratories, universities, etc.

* [Official website](http://www.seismovlab.com)
* [GitHub repository](https://github.com/SeismoVLAB/SVL)
* [Official documentation](http://www.seismovlab.com/documentation/index.html)

What can Seismo-VLAB do?
------------------------
With **Seismo-VLAB** you can solve:

* Linear and nonlinear wave propagation problems in shallow crust
* Linear and nonlinear soil-structure interaction problems
* Standard structural-mechanics linear and nonlinear dynamic problems

Visit our gallery to see some [examples](http://www.seismovlab.com/gallery.html) of simulations using **SVL**.

Installing Seismo-VLAB
----------------------
Installation of **Seismo-VLAB** on Linux/MacOS/Windows is perform in two steps: 


* The Pre-Process requires `python3` and the following libraries:
    * [Numpy](https://numpy.org/)
    * [Scipy](https://www.scipy.org/)
    * [Matplotlib](https://matplotlib.org/)
    * [JSON](https://www.json.org/json-en.html)
    
    These librarires can be installed using `pip3` or standard Linux and Mac repositories.

* The Run-Process requires requires to download `Eigen` C++ library, `MUMPS` Library, and `PETSc` Library.
    * The **Eigen C++ library** can be downloaded from this [website](http://eigen.tuxfamily.org/). This package needs to be unzip and its content move to `/usr/include/eigen`. 
    * The **MUMPS library** can be downloaded from this [website](http://mumps.enseeiht.fr/). This package needs to be unzip and compiled at `/usr/include/mumps`.
    * The **Pestc Library** can be downloaded at this [website](https://www.mcs.anl.gov/petsc/). This package needs to be unzip and compiled at `/usr/include/petsc`.
    
    Assuming the previous libraries are successfully installed, then modify the `Makefile.inc` file such the previous path point to the right libraries:
    ```makefile
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

License
=======

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**Seismo-VLAB** is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
**Seismo-VLAB** is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details http://www.gnu.org/licenses.

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
