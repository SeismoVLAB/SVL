#!/usr/bin/env python3
# -*- coding: Utf-8 -*-

import os
import sys
import copy
import inspect
import numpy as np
from datetime import date
from Core.Definitions import Entities, Options, ConvergeTest, SolverOption

def clc():
    """
    This function clears the screen and prints the SVL header\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Returns
    -------
    None
    """
    printHeader()

def cleanAll():
    """
    Sets all SVL dictionaries to its initial values\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Returns
    -------
    None
    """
    #Sets the Entities to empty dictionaries
    for key in Entities:
        Entities[key] = {}

    #Sets the Options to pre-defined values
    metis = Options['metispath']
    preanalysis = Options['preanalysis']
    runanalysis = Options['runanalysis']
    for key in Options:
        if isinstance(Options[key], str):
            Options[key] = ''
        elif isinstance(Options[key], int):
            Options[key] = 0
        elif isinstance(Options[key], np.ndarray):
            Options[key] = []

    Options['file'       ] = 'SeismoVLAB'
    Options['metispath'  ] = metis
    Options['description'] = '\n'
    Options['allocation' ] = 'NO'
    Options['numbering'  ] = 'Plain'
    Options['massform'   ] = 'consistent'
    Options['nparts'     ] =  1
    Options['preanalysis'] = preanalysis
    Options['runanalysis'] = runanalysis
    setFilePath()

def setFilePath():
    """
    Set the path of the working directory in Options.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Returns
    -------
    None
    """
    Options['path'] = os.path.abspath(os.path.dirname(sys.argv[0]))
    try:
        pythonpath =  os.environ['PYTHONPATH']
        pythonpath =  pythonpath.split('/')
        pythonpath =  pythonpath[:-1]
        pythonpath =  "/".join(pythonpath)
        Options['preanalysis'] = pythonpath + '/01-Pre_Process'
        Options['runanalysis'] = pythonpath + '/02-Run_Process'
        Options['metispath'] =  Options['preanalysis'] + '/Metis/mpmetis'
    except KeyError:
        Options['metispath'] = []
        print('\x1B[31m ERROR \x1B[0m: The \'PYTHONPATH\' variable has not been defined')
    

def printHeader():
    """
    Prints Seismo-VLAB header showing author information, contact 
    information, and Copyright.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Returns
    -------
    None
    """
    os.system('cls' if os.name == 'nt' else 'clear')
    today = date.today()
    today = today.strftime("%m/%d/%Y")

    print( "       _____      _                    _    ____          __         " )
    print(r"      / ___/___  (_)________ ___  ____| |  / / /   ____ _/ /_        " )
    print(r"      \__ \/ _ \/ / ___/ __ `__ \/ __ \ | / / /   / __ `/ __ \       " )
    print( "     ___/ /  __/ /__  / / / / / / /_/ / |/ / /___/ /_/ / /_/ /       " )
    print(r"    /____/\___/_/____/_/ /_/ /_/\____/|___/_____/\__,_/_.___/        " )
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
    
def debugInfo(level):
    """
    Provides information regarding the function called\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    level : int
        Level inside the dictionry to look for information

    Returns
    -------
    info : struct
        Structure containing information of the function called
    """
    callerframerecord = inspect.stack()[level]
    frame = callerframerecord[0]
    info  = inspect.getframeinfo(frame)

    return info

def printFormatted(d, indent=0):
    """
    Format and print the given dictionary\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    d : dict
        Dictionary containing the Entity to be printed
    indent : int
        Number of spaces to leave after a line break

    Returns
    -------
    None
    """
    for key, value in d.items():
       if indent == 0:
           print(" [" + str(key) + "]")
       if isinstance(value, dict):
          printFormatted(value, indent+1)
       else:
          print('  '*(indent+1) + str(key) + ' : ' + str(value))

def printAll(name):
    """
    Prints all members defined in Entities\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    name : str
        String containing the Entity to be printed

    Returns
    -------
    None
    """
    if Entities[name]:
        print("\n List of all " + name + " defined so far:")
        printFormatted(Entities[name])
    else:
        print("\n There are no " + name + " defined so far.\n")

def tryOpenfile(filepath):
    """
    This function attempst to open a file\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    filepath : str
        String containing the path where the file is located

    Returns
    -------
    bool
        Whether the open file operation was successful (False) of failed (True)
    """
    try:
        with open(filepath, "r") as fh:
            fh.close()
            return False
    except IOError:
        return True

def saveAs(filename='output'):
    """
    This function saves in a python file all entities defined so far in command prompt\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2021

    Parameters
    ----------
    filepath : str
        String containing the file name

    Returns
    -------
    None
    """
    #TODO: FINISH THIS PROCESS
    #The SVL output file name
    filepath = Options['path'] + '/' + filename + '.py'
    
    pyFile = open(filepath, "w+")

    #Writes the Entities one by one
    today = date.today()
    today = today.strftime("%m/%d/%Y")

    pyFile.write("%s\n" % "\"\"\"")
    pyFile.write("Description\n")
    pyFile.write("-----------\n")
    pyFile.write("%s\n" % Options['description'])
    pyFile.write("DATE: %s\n" % today)
    pyFile.write("%s\n\n" % "\"\"\"")

    pyFile.write("from Core import SeismoVLAB as SVL\n\n")

    pyFile.write("#User's (pre-defined) options\n")
    pyFile.write("SVL.Options[\'file\'     ] = \'%s\'\n" % Options['file'])
    pyFile.write("SVL.Options[\'format\'   ] = \'%s\'\n" % Options['format'])
    pyFile.write("SVL.Options[\'numbering\'] = \'%s\'\n" % Options['numbering'])
    pyFile.write("SVL.Options[\'metispath\'] = \'%s\'\n" % Options['metispath'])
    pyFile.write("SVL.Options[\'massform\' ] = \'%s\'\n" % Options['massform'])
    pyFile.write("SVL.Options[\'nparts\'   ] = %d\n" % Options['nparts'])
    pyFile.write("SVL.Options[\'dimension\'] = %d\n" % Options['dimension'])
    pyFile.write("\n")

    if Entities['Materials']:
        pyFile.write("#Create Materials\n")
        for mTag in Entities['Materials']:
            material = Entities['Materials'][mTag]
            pyFile.write("SVL.addMaterial(tag=%d, name=\'%s\', attributes=%s)\n" % (mTag, material['name'], str(material['attributes'])))
        pyFile.write("\n")

    if Entities['Sections']:
        pyFile.write("#Create Sections\n")
        for sTag in Entities['Sections']:
            section = Entities['Sections'][sTag]
            pyFile.write("SVL.addSection(tag=%d, name=\'%s\', model=\'%s\', attributes=%s)\n" %(sTag, section['name'], section['model'], section['attributes']))
        pyFile.write("\n")

    if Entities['Nodes']:
        pyFile.write("#Create Nodes\n")
        for nTag in Entities['Nodes']:
            node = Entities['Nodes'][nTag]
            pyFile.write("SVL.addNode(tag=%d, ndof=%d, coords=[%s], freedof=[%s], totaldof=[%s])\n" % (nTag, node['ndof'], ",".join(map(str,node['coords'])), ",".join(map(str,node['freedof'])), ",".join(map(str,node['totaldof']))))
        pyFile.write("\n")

    if Entities['Masses']:
        pyFile.write("#Create Masses\n")
        for nTag in Entities['Masses']:
            dof  = []
            vals = []
            mass = copy.deepcopy(Entities['Masses'][nTag])
            for k,m in enumerate(mass['mass']):
                if m != 0.0:
                    dof.append(k+1)
                    vals.append(m)
            pyFile.write("SVL.addMass(tag=%d, dof=[%s], vals=[%s])\n" % (nTag, ",".join(map(str,dof)), ",".join(map(str,vals))))
        pyFile.write("\n")

    if Entities['Constraints']:
        pyFile.write("#Constraint degree of freedom\n")
        for cTag in Entities['Constraints']:
            constraint = copy.deepcopy(Entities['Constraints'][cTag])
            attributes = {'stag': constraint['stag'], 'sdof': constraint['sdof'], 'mtag': constraint['mtag'], 'mdof': constraint['mdof'], 'factor': constraint['factor']}
            attributes['sdof'] += 1
            attributes['mdof'] = [k+1 for k in attributes['mdof']]
            pyFile.write("SVL.addConstraint(tag=%d, name=\'%s\', attributes=%s)\n" % (cTag, constraint['name'], str(attributes)))
        pyFile.write("\n")

    if Entities['Supports']:
        pyFile.write("#Supports motions\n")
        for cTag in Entities['Supports']:
            support = copy.deepcopy(Entities['Supports'][cTag])
            support['dof'] = [k+1 for k in support['dof']]
            pyFile.write("SVL.addSupportMotion(tag=%d, attributes=%s)\n" % (cTag, str(support)))
        pyFile.write("\n")

    if Entities['Elements']:
        pyFile.write("#Create Elements\n")
        for eTag in Entities['Elements']:
            element = Entities['Elements'][eTag]
            pyFile.write("SVL.addElement(tag=%d, conn=[%s], name=\'%s\', attributes=%s)\n" % (eTag, ",".join(map(str,element['conn'])), element['name'], element['attributes']))
        pyFile.write("\n")

    if Entities['Surfaces']:
        pyFile.write("#Create Surfaces\n")
        for sTag in Entities['Surfaces']:
             pyFile.write("SVL.addSurface(tag=%d, etag=%d, conn=[%s])\n" % (sTag, Entities['Surfaces'][sTag]['etag'], ",".join(map(str,Entities['Surfaces'][sTag]['conn']))))

    if Entities['Dampings']:
        pyFile.write("#Create Dampings\n")
        for dTag in Entities['Dampings']:
            pyFile.write("SVL.addDamping(tag=%d, name=\'%s\', attributes=%s)\n" % (dTag, Entities['Dampings'][dTag]['name'], str(Entities['Dampings'][dTag]['attributes'])))
        pyFile.write("\n")

    if Entities['Functions']:
        pyFile.write("#Create functions\n")
        for fTag in Entities['Functions']:
            pyFile.write("SVL.addFunction(tag=%d, name=\'%s\', attributes=%s)\n" % (fTag, Entities['Functions'][fTag]['name'], str(Entities['Functions'][fTag]['attributes'])))
        pyFile.write("\n")

    if Entities['Loads']:
        pyFile.write("#Create a Loads\n")
        for lTag in Entities['Loads']:
            pyFile.write("SVL.addLoad(tag=%d, name=\'%s\', attributes=%s)\n" % (lTag, Entities['Loads'][lTag]['name'], str(Entities['Loads'][lTag]['attributes'])))
        pyFile.write("\n")

    if Entities['Combinations']:
        pyFile.write("#Create a Combination\n")
        for cTag in Entities['Combinations']:
            pyFile.write("SVL.addCombinationCase(tag=%d, name=\'%s\', attributes=%s)\n" % (cTag, Entities['Combinations'][cTag]['name'], str(Entities['Combinations'][cTag]['attributes'])))
        pyFile.write("\n")

    if Entities['Recorders']:
        pyFile.write("#Create Recorder\n")
        for rTag in Entities['Recorders']:
            pyFile.write("SVL.addRecorder(tag=%d, attributes=%s)\n" % (rTag, str(Entities['Recorders'][rTag])))
        pyFile.write("\n")

    if Entities['Simulations']:
        pyFile.write("#Create Analysis\n")
        for aTag in Entities['Simulations']:
            #Indeces for each Entity
            nTag = Entities['Simulations'][aTag]['attributes']['analysis']
            mTag = Entities['Simulations'][aTag]['attributes']['algorithm']
            iTag = Entities['Simulations'][aTag]['attributes']['integrator']
            sTag = Entities['Simulations'][aTag]['attributes']['solver']

            #copy the algorithm (workspace stays the same)
            algorithm = copy.deepcopy(Entities['Algorithms'][mTag])
            if 'cnvgtest' in algorithm:
                for key in ConvergeTest:
                    if ConvergeTest[key] == algorithm['cnvgtest']:
                        algorithm['cnvgtest'] = key
            #copy the solver (workspace stays the same)
            solver = copy.deepcopy(Entities['Solvers'][sTag])
            if 'update' in solver:
                UPDATE = solver['update']
                solver['update'] = 'OFF'  if UPDATE == 1  else 'ON'
            if 'option' in solver:
                for key in SolverOption:
                    if SolverOption[key] == solver['option']:
                        solver['option'] = key
            pyFile.write("SVL.addAnalysis(tag=%d, attributes=%s)\n" % (nTag, str(Entities['Analyses'][nTag])))
            pyFile.write("SVL.addAlgorithm(tag=%d, attributes=%s)\n" % (mTag, str(algorithm)))
            pyFile.write("SVL.addIntegrator(tag=%d, attributes=%s)\n" % (iTag, str(Entities['Integrators'][iTag])))
            pyFile.write("SVL.addSolver(tag=%d, attributes=%s)\n" % (sTag, str(solver)))
            pyFile.write("SVL.addSimulation(tag=%d, combo=%d, attributes=%s)\n" % (aTag, Entities['Simulations'][aTag]['combo'], str(Entities['Simulations'][aTag]['attributes'])))
    pyFile.close()
