#!/usr/bin/python3
# -*- coding: Utf-8 -*-

import os
import sys
import subprocess
import numpy as np
from Core.Definitions import Entities, Options

def most_frequent(List):
    return max(set(List), key = List.count)

def CheckPartition(Partition):
    """
    This function makes the new partition to be consistent to previous
    one. The idea is too keep the previously clustered elements in the 
    same partition and new elements are re-distributed among partitions. 
    Unfortunately, this naive process may generate discontinuity between
    elements that are clustered.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2021

    Parameters
    ----------
    Partition : array 
        The list/array containing the information

    Returns
    -------
    None
        Updates the 'partition' fields in Options
    """
    nparts = Options['nparts']
    nElems = len(Entities['Elements'])

    if len(Partition) > 0:
        #if elements are added, distribute them in partitions 
        kTag = 0
        elems = dict()
        nodes = {x: set() for x in range(nparts)}
        parts = np.zeros(nElems)

        for m, eTag in enumerate(Entities['Elements']):
            if Options['clustermap'][eTag] != -1:
                parts[m] = Options['clustermap'][eTag]
                nodes[parts[m]].update(Entities['Elements'][eTag]['conn'])
            else:
                elems[eTag] = m
                parts[m] = Partition[kTag]
                Options['clustermap'][eTag] = Partition[kTag]
                kTag += 1 

        #TODO: This loop is very naive and can be improved. Essentially checks
        #the nodes of the element added, and see in which partition fits better
        #Corrects elements in a partition that may have nodes in another
        for eTag in elems:
            count =  0
            pTags = -1
            nTags = Entities['Elements'][eTag]['conn']
            for p in nodes:
                n = len(set(nTags).intersection(nodes[p]))
                if n > 0 and count < n:
                    pTags = p
                    count = n
            if pTags != -1:
                 pid = elems[eTag]
                 parts[pid] = pTags
                 Options['clustermap'][eTag] = pTags

        Options['partition'] = parts
    else:
        #If elements were removed, update partition and clustermap
        ne = len(Options['partition'])
        if ne != nElems:
            kTag = 0
            cluster = dict()
            elparts  = np.zeros(nElems)
            for m, eTag in enumerate(Options['clustermap']):
                if eTag in Entities['Elements']:
                    elparts[kTag] = Options['partition'][m]
                    cluster[eTag] = Options['clustermap'][eTag]
                    kTag += 1 

            Options['partition' ] = elparts
            Options['clustermap'] = cluster

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

    #Total number of (added) elements
    nElems = 0
    for eTag in Entities['Elements']:
        if eTag not in Options['clustermap']:
            nElems += 1   
    
    #Creates the Partition folder if it has not
    dirName  = Options['path'] + '/' + 'Partition'
    if not os.path.exists(dirName):
        os.mkdir(dirName)

    #Creates a Map between Node Tags
    MetisMap = dict()
    for k, nTag in enumerate(Entities['Nodes']):
        MetisMap[nTag] = k+1

    #Assign same Point Tag To Equal Constraints
    #OBJECTIVE: Make a uniform partition (in paricular if PML are used)
    #TODO: Be extremely carefull with these lines: 44-48, it gives problem for large models.
    auxMap = MetisMap.copy()
    for ctag in Entities['Constraints']:
        if Entities['Constraints'][ctag]['name'] == 'EQUAL':
            slave  = Entities['Constraints'][ctag]['stag' ]
            master = Entities['Constraints'][ctag]['mtag']
            MetisMap[slave] = auxMap[master[0]]

    #Creates the Partition File
    MetisPath = dirName + '/' + 'Graph.out'
    if nElems > 0:
        MetisFile = open(MetisPath, "w+")
        MetisFile.write(str(nElems) + '\n')

        #Writes the connectivities of the (added) elements for Metis.
        for eTag in Entities['Elements']:
            if eTag not in Options['clustermap']:
                Options['clustermap'][eTag] = -1
                connection = Entities['Elements'][eTag]['conn']
                for conn in connection:
                    MetisFile.write("%d " % MetisMap[conn])
                MetisFile.write("\n")
        MetisFile.close()
    else:
        #No elements to be partitioned (empty file is created)
        open(MetisPath, 'a+').close()

def GetMetisOutputFile():
    """
    This function runs the METIS - Serial Graph Partitioning by George Karypis,
    and gets the element indexes for each partition.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Returns
    -------
    Partition : list 
        Indexes of the element that belongs to each partition.
    """
    #Gets information of the mesh
    filePath = Options['path'] + '/' + 'Partition' + '/' + 'Graph.out'

    #Generates the division according to the number of parts
    nparts = Options['nparts']

    if nparts > 1:
        #Command line to be excecuted
        cmdline = Options['metispath'] + ' \'' + filePath + '\' ' + str(nparts)

        if os.stat(filePath).st_size != 0:
            #Executes METIS - Serial Graph Partitioning.
            subprocess.check_output(cmdline, shell=True, stderr=subprocess.DEVNULL)

            #Obtains the partition information
            Partition = np.loadtxt(filePath + '.epart.' + str(nparts), dtype='int', skiprows=0)
        else:
            open(filePath + '.epart.' + str(nparts), 'a+').close()
            open(filePath + '.npart.' + str(nparts), 'a+').close()

            #Obtains the partition information
            Partition = []   
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

    #Check partition consistency with previous partition
    CheckPartition(Partition)
