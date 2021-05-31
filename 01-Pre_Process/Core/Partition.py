#!/usr/bin/python3
# -*- coding: Utf-8 -*-

import os
import sys
import subprocess
import numpy as np
from Core.Definitions import Entities, Options

def SetMetisInputFile():
    """
    This function creates the input file for METIS - Serial Graph Partitioning
    by George Karypis.\n
    @visit  http://glaros.dtc.umn.edu/gkhome/metis/metis/download\n
    @author George Karypis

    Parameters
    ----------
    Entities : dict 
        The dictionary containing model's information

    Returns
    -------
    file : Graph.out 
        Element connectivity's of the Serial Graph.
    """

    #Total number of elements
    nElems = len(Entities['Elements'])
    
    #Creates the Partition folder if it has not
    dirName  = Options['path'] + '/' + 'Partition'
    if not os.path.exists(dirName):
        os.mkdir(dirName)

    #Creates a Map between Node Tags
    count = 0
    MetisMap = dict()
    for nTag in Entities['Nodes']:
        count += 1
        MetisMap[nTag] = count

    #Assign same Point Tag To Equal Constraints
    #OBJECTIVE: Make a uniform partition (in paricular if PML are used)
    #TODO: Be extremely carefull with these lines: 44-48, it gives problem for large models.
    auxMap = MetisMap.copy()
    for ctag in Entities['Constraints']:
        if Entities['Constraints'][ctag]['name'] == 'EQUAL':
            slave  = Entities['Constraints'][ctag]['stag' ]
            master = Entities['Constraints'][ctag]['mtag']
            MetisMap[slave] = auxMap[master[0]]

    #Opens the Partition File
    MetisPath = dirName + '/' + 'Graph.out'
    MetisFile = open(MetisPath, "w+")
    MetisFile.write(str(nElems) + '\n')

    #Writes the connectivities of the elements for Metis.
    for eTag in Entities['Elements']:
        connection = Entities['Elements'][eTag]['conn']
        for conn in connection:
            MetisFile.write("%d " % MetisMap[conn])
        MetisFile.write("\n")

    MetisFile.close()

def GetMetisOutputFile():
    """
    This function runs the METIS - Serial Graph Partitioning by George Karypis,
    and gets the element indexes for each partition.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Returns
    -------
    Partition : list 
        Indeces of the element that belongs to each partition.
    """
    #Gets information of the mesh
    filePath = Options['path'] + '/' + 'Partition' + '/' + 'Graph.out'

    #Generates the division according to the number of parts
    nparts = Options['nparts']
    if nparts > 1:
       #Executes METIS - Serial Graph Partitioning.
       cmdline = Options['metispath'] + ' \'' + filePath + '\' ' + str(nparts)
       subprocess.check_output(cmdline, shell=True)

       #Obtains the partition information
       Partition = np.loadtxt(filePath + '.epart.' + str(nparts), dtype='int', skiprows=0)
    elif nparts == 1:
        #Generates a unique partition.
        nElems = len(Entities['Elements'])

        #Element Partition (Useless)
        open(filePath + '.epart.' + str(nparts), 'a+').close()

        #Node Partition (Useless)
        open(filePath + '.npart.' + str(nparts), 'a+').close()

        #The partition information
        Partition = np.full(nElems, 0, dtype='int')
    else:
        print('\x1B[31m ERROR \x1B[0m: The requested number of partition is not possible.')
        sys.exit(-1)

    Options['partition'] = Partition