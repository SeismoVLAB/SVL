#!/usr/bin/python3

import os
import sys
from datetime import date
from Source import Parser
from Source import Numbering
from Source import SeismoVLab

def main(path):
    """
    This is the main file to generate the input files for SeismoVLab 

    Parameters
    ----------
    path : str 
           The location where the user's defines configuration file is located.

    Output
    -------
    file : User['MODEL'].$.mesh
           Finite element mesh partition.
    file : User['MODEL'].$.run
           Finite element analysis partition.
    """
    #Empty Parsing Dictionaries.
    User                   = {}
    Initial  , Support     = {}, {} 
    Mass     , Surface     = {}, {} 
    Point    , Element     = {}, {}
    Material , Section     = {}, {} 
    Diaphragm, Constraint  = {}, {}
    Damping  , Function    = {}, {}
    Load     , Combination = {}, {}
    Recorder , Simulation  = {}, {}

    #PARSE THE INPUT SEISMOVLAB FILE
    try:  
        #Open the provided file
        with open(path, "r") as fileHandler:  
            #Loop over each line in file.
            for lines in fileHandler:
                line = list(filter(None, lines.strip().split())) 
                if line:
                    if line[0].upper() == '*MODEL':
                        Parser.GetModelInformation(fileHandler, User, path)
                    elif line[0].upper() == '*PARTITION':
                        Parser.GetPartitionInformation(fileHandler, User)
                    elif line[0].upper() == '*MATERIAL':
                        Parser.GetMaterialInformation(fileHandler, User, Material)
                    elif line[0].upper() == '*SECTION':
                        Parser.GetSectionInformation(fileHandler, User, Section)
                    elif line[0].upper() == '*POINT':
                        Parser.GetNodesInformation(fileHandler, User, Point)
                    elif line[0].upper() == '*MASS':
                        Parser.GetMassInformation(fileHandler, User, Mass, Point)
                    elif line[0].upper() == '*INITIALSTATE':
                        Parser.GetInitialInformation(fileHandler, User, Initial, Point)
                    elif line[0].upper() == '*SUPPORTMOTION':
                        Parser.GetSupportInformation(fileHandler, User, Support, Point)
                    elif line[0].upper() == '*RESTRAIN':
                        Parser.GetRestrainInformation(fileHandler, User, Point)
                    elif line[0].upper() == '*CONSTRAINT':
                        Parser.GetConstraintInformation(fileHandler, User, Constraint)
                    elif line[0].upper() == '*DIAPHRAGM':
                        Parser.GetDiaphragmInformation(fileHandler, User, Diaphragm)
                    elif line[0].upper() == '*ELEMENT':
                        Parser.GetElementInformation(fileHandler, User, Element)
                    elif line[0].upper() == '*SURFACE':
                        Parser.GetSurfaceInformation(fileHandler, User, Surface)
                    elif line[0].upper() == '*DAMPING':
                        Parser.GetDampingInformation(fileHandler, User, Damping)
                    elif line[0].upper() == '*FUNCTION':
                        Parser.GetFunctionInformation(fileHandler, User, Function) 
                    elif line[0].upper() == '*LOAD':
                        Parser.GetLoadInformation(fileHandler, User, Point, Element, Function, Load)
                    elif line[0].upper() == '*COMBINATION':
                        Parser.GetCombinationInformation(fileHandler, User, Combination)
                    elif line[0].upper() == '*RECORDER':
                        Parser.GetRecorderInformation(fileHandler, User, Point, Element, Recorder)
                    elif line[0].upper() == '*SIMULATION':
                        Parser.GetSimulationInformation(fileHandler, User, Simulation)

            #Corrects Free degree-of-freedom numbering for Constrained Points.
            Parser.ComputeConstraints(User, Point, Constraint, Diaphragm)

    except IOError as e:
        print('\x1B[31m ERROR \x1B[0m: The FILE=' + path + ' could not be opened.')
        print(' Check if both the file path and file name provided are correct.')
        sys.exit(-1)

    #CHECK ALL PROVIDED INFORMATION IS CORRECT.
    Parser.CheckCorrectness(User, Point, Mass, Element, Material, Section, Constraint, Diaphragm, Support, Initial, Damping, Surface, Function, Load, Combination, Recorder, Simulation) 

    #GENERATES DRM INFORMATION IF REQUIRED.
    Parser.GenerateDRMFiles(User, Point, Element, Function, Load)

    #SET NODES DEGREE OF FREEDOM NUMBERING 
    Numbering.SetDegreeOfFreedom(User, Point, Element, Constraint)

    #DIVIDES THE MESH INTO FILES
    SeismoVLab.WriteMeshPartition(User, Point, Mass, Element, Material, Section, Constraint, Damping, Surface, Support, Initial, Function, Load, Combination, Recorder, Simulation)
    
    #COMMAND LINE EXECUTION FOR SEISMOVLAB PROGRAM
    print('\n\033[1;32m SeismoVLab can be executed as:\033[1;0m')
    print(User['EXE'])

if __name__ == "__main__":
    #CLEAR THE TERMINAL COMMAND LINE
    os.system('cls' if os.name == 'nt' else 'clear')

    #GETS THE CURRENT DATE
    today = date.today()
    today = today.strftime("%m/%d/%Y")

    #PRINTS THE SOFTWARE'S HEADER 
    print( "       _____      _                    _    ____          __         " )
    print( "      / ___/___  (_)________ ___  ____| |  / / /   ____ _/ /_        " )
    print( "      \__ \/ _ \/ / ___/ __ `__ \/ __ \ | / / /   / __ `/ __ \       " )
    print( "     ___/ /  __/ /__  / / / / / / /_/ / |/ / /___/ /_/ / /_/ /       " )
    print( "    /____/\___/_/____/_/ /_/ /_/\____/|___/_____/\__,_/_.___/        " )
    print( "                                                                     " )
    print( "               (Seismo)s (V)irtual (Lab)oratory                      " )
    print( "        Module for Serial and Parallel Analysis of seismic           " )
    print( "   wave propagation and soil-structure interaction simulation        " )
    print( " Copyright (C) " + today + ", The California Institute of Technology " )
    print( "                     All Rights Reserved.                            " )
    print( "                                                                     " )
    print( " \033[1;34mWritten by:                                      \033[1;0m" )
    print( "   Danilo S. Kusanovic (dkusanov@caltech.edu)                        " )
    print( "   Elnaz E. Seylabi    (elnaze@unr.edu)                              " )
    print( "                                                                     " )
    print( " \033[1;33mSupervised by:                                   \033[1;0m" )
    print( "   Domniki M. Asimaki  (domniki@caltech.edu)                         " )
    print( "                                                                     " )

    #GETS THE CONFIGURATION INPUT FILE
    if len(sys.argv) == 1:
        print(' Please provide an svl input file and run once more:')
        print('   python3 PreProcess.py /path/to/file/input.svl\n')
    else:
        path = str(sys.argv[1])
        path = path.replace("\"", "")
        main(path)
