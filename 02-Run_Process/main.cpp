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
// Version      :
// Date         : 10-22-2017
// Revised      : 02-18-2020
// Modified     : Danilo Kusanovic
// Source       : main.cpp
// Keywords     : <Main, Program, Parser, Mesh, Analysis>
// Description  : This is the main file. It creates the mesh object as well as
//                the simulation objects for serial/parallel analysis from a 
//                file using the pre-process module, and it runs the analysis. 
//
// Example of Usage:
//  *SERIAL VERSION
//   ./SeismoVLAB.exe -dir /file/path/ -file file.$.json 
//  *PARALLEL VERSION
//   mpirun -np n ./SeismoVLAB.exe -dir /file/path/ -file file.$.json
//------------------------------------------------------------------------------

#include <memory>
#include <stdexcept> 
#include <petscsys.h>

#include "Parser.hpp"
#include "Utilities.hpp"
#include "Definitions.hpp"
#include "Profiler.hpp"

int main(int argc, char **argv){ 
    bool parsefile = false;
                                           
    //Initialize MPI instance.
    PetscInitialize(&argc, &argv, NULL, NULL);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    //Prints the Seismo-VLAB logo.
    printLogo();

    //Parse command line arguments.    
    CommandLine(argc, argv, parsefile);

    //Initialize the profiler.
    #if PROFILING
    Profiler::Get().BeginSession("Profile");
    #endif

    //Performs the simulation using JSON input file.
    RunFromJSON(parsefile);

    //Finalize the profiler.
    #if PROFILING
    Profiler::Get().EndSession();
    #endif

    //Finalize MPI library.
    PetscFinalize();

    return 0;
}
