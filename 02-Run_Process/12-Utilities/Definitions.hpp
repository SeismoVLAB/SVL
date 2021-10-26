//==============================================================================
//
//                       Seismo Virtual Laboratory
//             Module for Serial and Parallel Analysis of seismic 
//         wave propagation and soil-structure interaction simulation
//         Copyright (C) 2018-2021, The California Institute of Technology
//                         All Rights Reserved.
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
// This file sets the global variables to be used during SeismoVLAB execution. 
//------------------------------------------------------------------------------

#ifndef __DEFINITIONS_HPP__ 
#define __DEFINITIONS_HPP__ 

#include <vector>
#include <string>

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      June 9, 2020
/// @version   1.0
/// @file      Definitions.hpp
/// @see       Utilities.hpp and Profiler.hpp
/// @brief This file sets the global variables to be used during SeismoVLAB execution. 

///Define if profiler is active (0: no profile, 1: profile code)
#define PROFILING 0

///Define macro for unused parameter
#define UNUSED(x)

///Define group of truss elements, i.e, linear and quadratic for ParaView visualization
#define GROUP_ELEMENT_TRUSS 10

///Define group of triangular elements, i.e, linear and quadratic for ParaView visualization
#define GROUP_ELEMENT_TRIA 11

///Define group of quadrilateral elements, i.e, linear and quadratic for ParaView visualization
#define GROUP_ELEMENT_QUAD 12

///Define group of tetrahedra elements, i.e, linear and quadratic for ParaView visualization
#define GROUP_ELEMENT_TETRA 13

///Define group of hexahedra elements, i.e, linear and quadratic for ParaView visualization
#define GROUP_ELEMENT_HEXA 14

///Define group of frame elements, i.e, linear and quadratic for ParaView visualization
#define GROUP_ELEMENT_FRAME 20

///Define group of shell elements, i.e, linear and quadratic for ParaView visualization
#define GROUP_ELEMENT_SHELL 21

///Define group of zero length elements, i.e, linear and quadratic for ParaView visualization
#define GROUP_ELEMENT_ZERO 31

///Define group of perfectly-matched layer elements, i.e, linear and quadratic for ParaView visualization
#define GROUP_ELEMENT_PML 32

///Define group of bouc-wen device elements, i.e, linear and quadratic for ParaView visualization
#define GROUP_ELEMENT_BWEN 41

///Define group of high-damping rubber bearing device elements, i.e, linear and quadratic for ParaView visualization
#define GROUP_ELEMENT_HDRB 42

///Define group of null elements, i.e, linear and quadratic for ParaView visualization
#define GROUP_ELEMENT_NULL 50

///Define ParaView VTK point (1 nodes) element types
#define VTK_POINT 1 

///Define ParaView VTK line (2 nodes) element types
#define VTK_LINEAR_LINE 3

///Define ParaView VTK line (3 nodes) element types
#define VTK_QUADRATIC_LINE 21

///Define ParaView VTK triangular (2 nodes) element types
#define VTK_LINEAR_TRIA 5

///Define ParaView VTK triangular (3 nodes) element types
#define VTK_QUADRATIC_TRIA 22

///Define ParaView VTK quadrilateral (4 nodes) element types
#define VTK_LINEAR_QUAD 9

///Define ParaView VTK quadrilateral (8 nodes) element types
#define VTK_QUADRATIC_QUAD 23

///Define ParaView VTK tetrahedron (4 nodes) element types
#define VTK_LINEAR_TETRA 10

///Define ParaView VTK tetrahedron (10 nodes) element types
#define VTK_QUADRATIC_TETRA 24

///Define ParaView VTK hexahedron (8 nodes) element types
#define VTK_LINEAR_HEXA 12

///Define ParaView VTK hexahedron (20 nodes) element types
#define VTK_QUADRATIC_HEXA 25 

///Define global load pattern constant concentrated load
#define POINTLOAD_CONCENTRATED_CONSTANT 1

///Define global load pattern time dependant concentrated load
#define POINTLOAD_CONCENTRATED_DYNAMIC 2

///Define global load pattern for a node body load
#define POINTLOAD_BODY_CONSTANT 3

///Define global load pattern time dependant for a body load
#define POINTLOAD_BODY_DYNAMIC 4

///Define global load pattern for a constant surface load
#define ELEMENTLOAD_SURFACE_CONSTANT 5

///Define global load pattern for a time dependant surface load
#define ELEMENTLOAD_SURFACE_DYNAMIC 6

///Define global load pattern for a constant body load
#define ELEMENTLOAD_BODY_CONSTANT 7

///Define global load pattern for a time dependant body load
#define ELEMENTLOAD_BODY_DYNAMIC 8

///Define global load pattern time dependant general wave load
#define ELEMENTLOAD_DOMAIN_REDUCTION 9

///Define global load pattern constant/time-dependant support motion load
#define POINTLOAD_SUPPORT_MOTION 10

///The processor number.
extern int rank;

///The number of partitions.
extern int size;

///The problem dimension (1D, 2D, 3D). 
extern unsigned int nDimensions;

///The folder path where the file is loaded.
extern std::string filePath;

///The file name to be loaded.
extern std::vector<std::string> fileName;

///Whether the driver (JSON) file is provided  
extern bool driverFile;

///The element mass formulation.
extern bool MassFormulation;

///The update option for member in Mesh.
extern std::string UpdateOption;

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
