#!/usr/bin/python3

import os
import sys
import subprocess
import numpy as np

def SetMetisInputFile(User, Point, Element, Constraint):
    """
    This function creates the input file for METIS - Serial Graph Partitioning
    by George Karypis.

    Parameters
    ----------
    path : str 
           The location where the output file will be stored
    Point : dict 
           The dictionary containing all point information.
    Element : dict 
           The dictionary containing all the user's relevant information.

    Output
    -------
    Metisfile : Graph.out 
           Element connectivities of the Serial Graph.
    """
    #Total number of elements
    nElems = len(Element)
    User['NELEMS'] = nElems
    
    #Creates the Partition Directory
    dirName = User['FOLDER'] + 'Partition'
    if not os.path.exists(dirName):
        os.mkdir(dirName)

    #Creates a Map between Point Tags
    MetisMap = dict()
    for nTag in Point:
        MetisMap[nTag] = Point[nTag]['TAG'] + 1

    #Assign same Point Tag To Equal Constraints
    #OBJECTIVE: Make a uniform partition.
    #TODO: Be extremely carefull with these lines: 44-48, it gives problem for large models.
    auxMap = MetisMap.copy()
    for cTag in Constraint:
        if Constraint[cTag]['NAME'] == 'EQUAL':
            slave  = Constraint[cTag]['SLAVENODE' ]
            master = Constraint[cTag]['MASTERNODE']
            MetisMap[slave] = auxMap[master[0]]

    #Opens the Partition File 
    Metisfile = open(User['FOLDER'] + 'Partition/Graph.out', "w+")
    Metisfile.write(str(nElems) + '\n')

    #Writes the connectivities of the elements for Metis.
    for eTag in Element:
        connection = Element[eTag]['CONNECTIVITY']
        for conn in connection:
            Metisfile.write("%d " % MetisMap[conn])
        Metisfile.write("\n")

    Metisfile.close()

def GetMetisOutputFile(User):
    """
    This function runs the METIS - Serial Graph Partitioning by George Karypis,
    and gets the element indeces for each partition 

    Parameters
    ----------
    User : dict 
           The dictionary containing all the user's relevant information.

    Returns
    -------
    Partition : list 
           Indeces of the element that belongs to each partition.
    """
    #Gets information of the mesh
    npart = User['NPART']
    path  = User['FOLDER'] + 'Partition/Graph.out'

    #Generates the division according to the number of parts
    if npart > 1:
       #Metis path
       metis = User['METIS']

       #Excecutes METIS - Serial Graph Partitioning.
       cmdline = metis + ' \'' + path + '\' ' + str(npart) 
       subprocess.check_output(cmdline, shell=True)

       #Obtains the partition information
       Partition = np.loadtxt(path + '.epart.' + str(npart), dtype='int', skiprows=0)
    elif npart == 1:
        #Generates a unique partition.
        nElems =  User['NELEMS']

        #Element Partition (Useless)
        open(path + '.epart.' + str(npart), 'a+').close()

        #Node Partition (Useless)
        open(path + '.npart.' + str(npart), 'a+').close()

        #The partition information
        Partition = np.full(nElems, 0, dtype='int')
    else:
        print('\x1B[31m ERROR \x1B[0m: The specified number of partition is not possible.')
        print(' Check that the number of partition are greter than 0.')
        sys.exit(-1)

    return Partition
