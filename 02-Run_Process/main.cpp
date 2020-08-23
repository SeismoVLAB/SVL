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
//   ./SeismoVLab -dir /file/path/ -file file.$.svl 
//  *PARALLEL VERSION
//   mpirun -np n ./SeismoVLab -dir /file/path/ -file file.$.svl
//------------------------------------------------------------------------------

#include <memory>
#include <stdexcept> 
#include <petscsys.h>

#include "Parser.hpp"
#include "Mesh.hpp"
#include "Analysis.hpp"
#include "Utilities.hpp"
#include "Definitions.hpp"
#include "Profiler.hpp"

int main(int argc, char **argv){                                            
    //Auxiliar variables.
    bool Stop;
    std::string File;
    std::string Folder;

    //Initialize MPI instance.
    PetscInitialize(&argc, &argv, NULL, NULL);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    //Prints the SeismoVLab logo.
    printLogo();

    //Parse command line argumets.    
    Stop = CommandLine(argc, argv, Folder, File);

    //Initialize the profiler.
    Profiler::Get().BeginSession("Profile", Folder);

    //Check command line arguments is valid.
    if(!Stop){
        //Creates Subdomain Pointer.
        std::shared_ptr<Mesh> mesh;

        //Creates Analysis Pointer.
        std::unique_ptr<Analysis> analysis;

        //Parse the Mesh/Analysis From User's Files.
        Parser SVL(size, Folder, File);
        Stop = SVL.GetFromFile(mesh, analysis);

        //Starts analysis if mesh/analysis are created.
        if(!Stop) 
            Stop = analysis->Analyze();
    }

    //Finilize the profiler.
    Profiler::Get().EndSession();

    //Finilize MPI library.
    PetscFinalize();

    return 0;
}
