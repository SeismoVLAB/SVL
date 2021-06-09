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
// This set of functions gets the command line user's inputs.
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
printHelp(bool &help){
    //Command Help. 
    if(rank == 0){
        help = true;
        std::cout << " COMMAND LINE FLAGS TO PROVIDE:                                  \n";
        std::cout << "     -dir  : Location of the working directory.                  \n";
        std::cout << "     -file : Name of the SeismoVLAB file to be loaded.           \n";
        std::cout << "                                                                 \n";
        std::cout << " \x1B[33mRun \x1B[0m: mpirun -np n ./SeismoVLAB.exe -dir '/path/to/files' -file 'model.$.json'\n";
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

        std::cout << "                                                                       \n";
        std::cout << "                   ░                                                   \n";
        std::cout << "     ███████╗ ░░░░   ░░░░ ░░░░░░░ ░░░░ ██╗   ██╗ ██╗   ░░░  ░░░░       \n";
        std::cout << "     ██╔════╝ ░    ░ ░    ░  ░  ░ ░  ░ ██║   ██║ ██║  ░   ░ ░   ░      \n";
        std::cout << "     ███████╗ ░░░  ░ ░░░░ ░  ░  ░ ░███╗██║   ██║ ██║  ░   ░ ░░░░       \n";
        std::cout << "     ╚════██║ ░    ░    ░ ░  ░  ░ ░╚═░╝╚██╗ ██╔╝ ██║  ░░░░░ ░   ░      \n";
        std::cout << "     ███████║ ░░░░ ░ ░░░░ ░  ░  ░ ░░░░  ╚████╔╝  ███████╗ ░ ░░░░       \n";
        std::cout << "     ╚══════╝                            ╚═══╝   ╚══════╝              \n";
        std::cout << "                    Seismo Virtual Laboratory                          \n";
        std::cout << "          Module for Serial and Parallel Analysis of seismic           \n";
        std::cout << "      wave propagation and soil-structure interaction simulation       \n";
        std::cout << " Copyright (C) " << year << ", The California Institute of Technology  \n";
        std::cout << "                       All Rights Reserved.                            \n";
        std::cout << "                                                                       \n";
        std::cout << " \033[1;34mWritten by:                                        \033[1;0m\n";
        std::cout << "   Danilo S. Kusanovic (dkusanov@caltech.edu)                          \n";
        std::cout << "   Elnaz E. Seylabi    (elnaze@unr.edu)                                \n";
        std::cout << "                                                                       \n";
        std::cout << " \033[1;33mSupervised by:                                     \033[1;0m\n";
        std::cout << "   Domniki M. Asimaki  (domniki@caltech.edu)                           \n";
        std::cout << "                                                                       \n";
    }
}

///Parse Command Line Inputs. 
void
CommandLine(int argc, char **argv, bool &parsefile){
    //Auxiliary Variable.
    int iter  = 0; 
    int count = 0;
    bool help = false;

    //Find help flag is active.
    while(iter < argc){
        //Help Command Line Output is Active.
        if(strcasecmp(argv[iter],"--help") == 0){
            printHelp(help);
        }
        else if(strcasecmp(argv[iter],"--h") == 0){
            printHelp(help);
        }
        iter++;
    }

    //SeismoVLab analysis variables.
    iter = 0;
    while(iter < argc){
        //Name of the Mesh input file.
        if(strcasecmp(argv[iter],"-file") == 0){
            fileName = std::string(argv[iter + 1]);
            count++;        
        }
        //Location of the Mesh File:
        if(strcasecmp(argv[iter],"-dir") == 0){
            filePath = std::string(argv[iter + 1]);    
            count++;    
        }
        iter++;     
    }

    //No Enough Number of Input Arguments.
    if(count == 2){
        parsefile = true;
    }
    else{
        if(!help){
            if(rank == 0)
                std::cout << "\x1B[31m ERROR: \x1B[0mNOT ENOUGH Command line input arguments. \n";
        }
    }
}

#endif
