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
///This set of functions gets the command line user's inputs.
//------------------------------------------------------------------------------

#ifndef _UTILITIES_HPP_
#define _UTILITIES_HPP_

#include <ctime>
#include <string>
#include <iostream>
#include <stdlib.h>
#include "Definitions.hpp"

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      December 18, 2018
/// @version   1.0
/// @file      Utilities.hpp
/// @see       main.cpp

///Prints out Command Line help
void 
printHelp(){
    //Command Help. 
    if(rank == 0){
        std::cout << " COMMAND LINE FLAGS TO PROVIDE:                                  \n";
        std::cout << "     -dir  : Location of the working directory.                  \n";
        std::cout << "     -file : Name of the SeismoVLab file to be loaded.           \n";
        std::cout << "                                                                 \n";
        std::cout << " \x1B[33mUse Example \x1B[0m: ./SeismoVLab -dir /path/to/files -file model.$.svl\n";
    }
}

///Prints Seismo-VLAB software header
void 
printLogo(){
    //SeismoVLab Logo. 
    if(rank == 0){
        time_t t = time(0);
        tm *timePtr = localtime(&t);
        unsigned int year = timePtr->tm_year + 1900;

        std::cout << "       _____      _                    _    ____          __            \n";
        std::cout << "      / ___/___  (_)________ ___  ____| |  / / /   ____ _/ /_           \n";
        std::cout << "      \\__ \\/ _ \\/ / ___/ __ `__ \\/ __ \\ | / / /   / __ `/ __ \\    \n";
        std::cout << "     ___/ /  __/ /__  / / / / / / /_/ / |/ / /___/ /_/ / /_/ /          \n";
        std::cout << "    /____/\\___/_/____/_/ /_/ /_/\\____/|___/_____/\\__,_/_.___/        \n";
        std::cout << "                                                                        \n";
        std::cout << "               (Seismo)s (V)irtual (Lab)oratory                         \n";
        std::cout << "        Module for Serial and Parallel Analysis of Seismic              \n";
        std::cout << "   Wave Propagation and Soil-Structure Interaction Simulation           \n";
        std::cout << "   Copyright (C) " << year << ", The California Institute of Technology \n";
        std::cout << "                     All Rights Reserved.                               \n";
        std::cout << "                                                                        \n";
        std::cout << " \033[1;34mWritten by:                                         \033[1;0m\n";
        std::cout << "   Danilo S. Kusanovic (dkusanov@caltech.edu)                           \n";
        std::cout << "   Elnaz E. Seylabi    (elnaze@unr.edu)                                 \n";
        std::cout << "                                                                        \n";
        std::cout << " \033[1;33mSupervised by:                                      \033[1;0m\n";
        std::cout << "   Domniki M. Asimaki  (domniki@caltech.edu)                            \n";
        std::cout << "                                                                        \n";
    }
}

///Parse Command Line Inputs.
bool 
CommandLine(int argc, char **argv, std::string &folder, std::string &file){
    //Auxiliar Variable.
    int iter  = 0; 
    int count = 0;

    //Find help flag is active.
    while(iter < argc){
        //Help Command Line Output is Active.
        if(strcasecmp(argv[iter],"--help") == 0){
            printHelp();
            return true;
        }
        else if(strcasecmp(argv[iter],"--h") == 0){
            printHelp();
            return true;
        }
        iter++;
    }

    //SeismoVLab analysis variables.
    iter = 0;
    while(iter < argc){
        //Name of the Mesh input file.
        if(strcasecmp(argv[iter],"-file") == 0){
            file = std::string(argv[iter + 1]);
            count++;        
        }
        //Location of the Mesh File:
        if(strcasecmp(argv[iter],"-dir") == 0){
            folder = std::string(argv[iter + 1]);    
            count++;    
        }
        iter++;     
    }

    //No Enough Number of Input Arguments.
    if(count != 2){
        if(rank == 0)
            std::cout << "\x1B[31m ERROR: \x1B[0mNOT ENOUGH Command line input arguments. \n";
        printHelp();
        return true;
    }
    else{
        return false;
    }
}

#endif
