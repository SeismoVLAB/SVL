//==============================================================================
//
//                    (Seismo)s (V)irtual (Lab)oratory
//             Module for Serial and Parallel Analysis of seismic 
//         wave propagation and soil-structure interaction simulation
//         Copyright (C) 2018, The California Institute of Technology
//                           All Rights Reserved.
//
// Commercial use of this program without express permission of the California
// Institute of Technology, is strictly  prohibited. See  file "COPYRIGHT"  in
// main  directory  for  information on  usage  and  redistribution, and for a
// DISCLAIMER OF ALL WARRANTIES.
//
//==============================================================================
//
// Written by:
//   Danilo S. Kusanovic (dkusanov@caltech.edu)
//   Elnaz E. Seylabi    (elnaze@unr.edu)
//
// Supervised by:
//   Domniki M. Asimaki  (domniki@caltech.edu)
//
// References : 
//   [1]
//
// Description:
///This file sets the global variables to be used during execution. 
//------------------------------------------------------------------------------

#ifndef __DEFINITIONS_HPP__ 
#define __DEFINITIONS_HPP__ 

#include <string>

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      June 9, 2020
/// @version   1.0
/// @file      Definitions.hpp
/// @see       Utilities.hpp and Profiler.hpp

///Define if profiler is active (0: no profile, 1: profile code)
#define PROFILING 0

///Define macro for unused parameter
#define UNUSED(x)

///The processor number.
extern int rank;

///The number of partitions.
extern int size;

///The problem dimension (1D, 2D, 3D). 
extern unsigned int nDimensions;

///The file name to be loaded.
extern std::string fileName;

///The folder path where the file is loaded.
extern std::string filePath;

///The Execution for each simulation.
extern bool FormOfExecution;

///The element mass formulation.
extern bool MassFormulation;

///Maximum memory for lumped storage sparse matrix.
extern unsigned int LumpedStorage;

///Maximum memory for consistent storage sparse matrix.
extern unsigned int ConsistentStorage;

///Total number of free-degree-of-freedom.
extern unsigned int numberOfFreeDofs;

///Total number of total-degree-of-freedom.
extern unsigned int numberOfTotalDofs;

///Total number of constrained-degree-of-freedom.
extern unsigned int numberOfConstrainedDofs;

///Maximum memory for PML 3D storage sparse matrix.
extern unsigned int PMLStorage;

#endif
