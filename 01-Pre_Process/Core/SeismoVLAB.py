#!/usr/bin/python3
# -*- coding: Utf-8 -*-

import os
import sys
import math
import copy
import numpy as np
from Method.Attach import *
from Method.Remove import *
from Method.Display import *
from Method.Builder import *
from Method.Compute import *
from Parser.Formats import *
from Core.Outputs import *
from Core.Utilities import *
from Core.Numberer import *
from Core.Partition import *
from Core.PlaneWave import *
from Core.Definitions import *

def createFolders():
    """
    This function creates all required folders, these are:
       'Partition' : Strores the domain decomposition and SVL Run-Analysis files
       'Paraview'  : Stores the ParaView VTK animation files for a given load combination
       'Solution'  : Stores the recorded responses for a given load combination
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Returns
    -------
    None
    """
    #Creates the Partition/Paraview/Results folders
    for folder in ['Partition','Paraview','Solution']:
        dirName  = Options['path'] + '/' + folder
        if not os.path.exists(dirName):
            os.mkdir(dirName)

    #Creates the Load combination folders used in Results/Paraview 
    for cTag in Entities['Combinations']:
        for folder in ['Solution','Paraview']:
            dirName  = Options['path'] + '/' + folder + '/' + Entities['Combinations'][cTag]['attributes']['folder']
            if not os.path.exists(dirName):
                os.mkdir(dirName)

def Entities2Processor(matSubdomain, secSubdomain, nodeSubdomain, massSubdomain, conSubdomain, elemSubdomain, surfSubdomain, k):
    """
    This function creates a dictionary that holds all information required to
    be written in the k-th processor  
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    matSubdomain  : list
        The material indexes that belongs to the k-th partition
    secSubdomain  : list
        The section indexes that belongs to the k-th partition
    nodeSubdomain : list
        The node indexes that belongs to the k-th partition
    conSubdomain  : list
        The constraint indexes that belongs to the k-th partition
    elemSubdomain : list
        The material indexes that belongs to the k-th partition
    surfSubdomain : list
        The surface indexes that belongs to the k-th partition
    k : int
        The processor (partition) number

    Returns
    -------
    ToProcessor : dict
        Contains all Entities (information) that belongs to this process
    """
    #Empty dictionary to be written
    ToProcessor = {
        'Global'      : {},
        'Materials'   : {}, 
        'Sections'    : {}, 
        'Nodes'       : {}, 
        'Masses'      : {}, 
        'Supports'    : {}, 
        'Constraints' : {},
        'Surfaces'    : {},
        'Elements'    : {}, 
        'Dampings'    : {}, 
        'Loads'       : {}, 
        'Combinations': {}, 
        'Recorders'   : {},
        'Simulations' : {}
        }
    
    #Global parameters stored in Entities for simulation
    ToProcessor['Global'] = {
        'ndim': Options['dimension'], 
        'ntotal': Options['ntotal'], 
        'nfree': Options['nfree'], 
        'massform': Options['massform'].upper(), 
        'solution': Options['solution'].upper()
    } 

    #Gets the materials for this partition
    for tag in matSubdomain:
        ToProcessor['Materials'][str(tag)] = Entities['Materials'][tag]

    #Gets the sections for this partition
    for tag in secSubdomain:
        ToProcessor['Sections'][str(tag)] = Entities['Sections'][tag]

    #Gets the nodes for this partition
    nTags = sorted(nodeSubdomain)
    for tag in nTags:
        ToProcessor['Nodes'][str(tag)] = Entities['Nodes'][tag]

    #Gets the point masses for this partition
    nTags = sorted(massSubdomain.intersection(nodeSubdomain))
    for tag in nTags:
        ToProcessor['Masses'][str(tag)] = copy.deepcopy(Entities['Masses'][tag])
        massSubdomain.remove(tag)

    #Gets the support motions for this partition
    sTags = set(Entities['Supports'].keys())
    nTags = sorted(sTags.intersection(nodeSubdomain))
    for tag in nTags:
        ToProcessor['Supports'][str(tag)] = Entities['Supports'][tag]

    #Gets the constraints for this partition
    cTags = sorted(conSubdomain, reverse=True)
    for tag in cTags:
        mTags = []
        Constraint = Entities['Constraints'][tag]
        sTags = Entities['Nodes'][Constraint['stag']]['totaldof'][Constraint['sdof']]
        for node, dof in zip(Constraint['mtag'], Constraint['mdof']):
            mTags.append(Entities['Nodes'][node]['freedof'][dof])
        ToProcessor['Constraints'][str(tag)] = {'stag': sTags, 'mtag': mTags, 'factor': Constraint['factor']}

    #Gets the elements for this partition
    Options['nparaview'] = 0
    for tag in elemSubdomain:
        ToProcessor['Elements'][str(tag)] = Entities['Elements'][tag]
        Options['nparaview'] += (len(Entities['Elements'][tag]['conn']) + 1)
    Options['nfeatures'] += Options['nparaview']

    #Gets the surfaces for this partition
    for sTag in surfSubdomain:
        eTag = Entities['Surfaces'][sTag]['etag']
        face = Entities['Surfaces'][sTag]['face']
        ToProcessor['Surfaces'][str(sTag)] = {'element': eTag, 'face': face}

    #Gets the dampings for this partition
    for dTag in Entities['Dampings']:
        eTag = list(set(elemSubdomain).intersection(Entities['Dampings'][dTag]['attributes']['list']))
        if eTag:
             ToProcessor['Dampings'][str(dTag)] = copy.deepcopy(Entities['Dampings'][dTag])
             ToProcessor['Dampings'][str(dTag)]['attributes']['list'] = eTag

    #Gets the loads for this partition
    loadSubdomain = list()
    for lTag in Entities['Loads']:
        name = Entities['Loads'][lTag]['name']
        if name == 'POINTLOAD':
            nTags = sorted(nodeSubdomain.intersection(Entities['Loads'][lTag]['attributes']['list']))
            if nTags:
                fTag  = Entities['Loads'][lTag]['attributes']['fun']
                fname = Entities['Functions'][fTag]['name']
                fdir  = Entities['Functions'][fTag]['attributes']['dir']
                if fname == 'CONSTANT':
                    magnitude = Entities['Functions'][fTag]['attributes']['mag']
                    attributes = {'name': fname, 'mag': magnitude, 'dir': fdir, 'list': nTags}
                    ToProcessor['Loads'][str(lTag)] = {'name': name, 'attributes': attributes}
                elif fname == 'TIMESERIE':
                    filepath = Entities['Functions'][fTag]['attributes']['file']
                    attributes = {'name': fname, 'file': filepath, 'dir': fdir, 'list': nTags}
                    ToProcessor['Loads'][str(lTag)] = {'name': name, 'attributes': attributes}
                loadSubdomain.append(lTag)
                Entities['Loads'][lTag]['attributes']['list'] = list(set(Entities['Loads'][lTag]['attributes']['list']).difference(nTags))
        elif name == 'ELEMENTLOAD':
            #Check the indexes correspond to element or surface
            if Entities['Loads'][lTag]['attributes']['type'] == 'SURFACE':
                eTags = sorted(list(set(surfSubdomain).intersection(Entities['Loads'][lTag]['attributes']['list'])))
            else:
                eTags = sorted(list(set(elemSubdomain).intersection(Entities['Loads'][lTag]['attributes']['list'])))

            if eTags:
                fTag  = Entities['Loads'][lTag]['attributes']['fun']
                fname = Entities['Functions'][fTag]['name']
                lname = Entities['Loads'][lTag]['attributes']['type']
                if fname == 'CONSTANT':
                    fdir  = Entities['Functions'][fTag]['attributes']['dir']
                    magnitude = Entities['Functions'][fTag]['attributes']['mag']
                    attributes = {'name': fname, 'type': lname, 'mag': magnitude, 'dir': fdir, 'list': eTags}
                    ToProcessor['Loads'][str(lTag)] = {'name': name, 'attributes': attributes}
                elif fname == 'TIMESERIE':
                    filepath = Entities['Functions'][fTag]['attributes']['file']
                    if lname == 'GENERALWAVE':
                        attributes = {'name': fname, 'type': lname, 'file': filepath, 'list': eTags}
                        ToProcessor['Loads'][str(lTag)] = {'name': name, 'attributes': attributes}
                    elif lname == 'PLANEWAVE':
                        features = {} #Entities['Functions'][fTag]['features']
                        attributes = {'name': fname, 'type': lname, 'file': filepath, 'features': features, 'list': eTags}
                        ToProcessor['Loads'][str(lTag)] = {'name': name, 'attributes': attributes}
                    else:
                        fdir  = Entities['Functions'][fTag]['attributes']['dir']
                        attributes = {'name': fname, 'type': lname, 'file': filepath, 'dir': fdir, 'list': eTags}
                        ToProcessor['Loads'][str(lTag)] = {'name': name, 'attributes': attributes}
                loadSubdomain.append(lTag)
        elif name == 'SUPPORTMOTION':
            nTags = sorted(nodeSubdomain.intersection(Entities['Loads'][lTag]['attributes']['list']))
            if nTags:
                attributes = {'list': nTags}
                ToProcessor['Loads'][str(lTag)] = {'name': name, 'attributes': attributes}
                loadSubdomain.append(lTag)

    #Gets the load combinations for this partition
    for cTag in Entities['Combinations']:
        name = Entities['Combinations'][cTag]['name']
        lTag = list(set(Entities['Combinations'][cTag]['attributes']['load']).intersection(loadSubdomain))
        if lTag:
            load    = Entities['Combinations'][cTag]['attributes']['load']
            factors = Entities['Combinations'][cTag]['attributes']['factor']
            folder  = Entities['Combinations'][cTag]['attributes']['folder']
            attributes = {'folder': folder, 'load': [], 'factor': []}
            for j in loadSubdomain:
                for i in range(len(load)):
                    if load[i] == j:
                        attributes['load'].append(load[i])
                        attributes['factor'].append(factors[i])
        else:
             attributes = {}
        ToProcessor['Combinations'][str(cTag)] = {'name': name, 'attributes': attributes}

    #Gets the recorders for this partition
    for rTag in Entities['Recorders']:
        #Modifies the Output File to be consistent with Partition
        OUTFILE = Entities['Recorders'][rTag]['file']
        OUTFILE = OUTFILE.replace(".", '.' + str(k) + '.')
        
        if Entities['Recorders'][rTag]['name'] == 'NODE':
            nTags = sorted(nodeSubdomain.intersection(Entities['Recorders'][rTag]['list']))
            if nTags:
                ToProcessor['Recorders'][str(rTag)] = copy.deepcopy(Entities['Recorders'][rTag])
                ToProcessor['Recorders'][str(rTag)]['file'] = OUTFILE
                ToProcessor['Recorders'][str(rTag)]['list'] = nTags.copy()
        elif Entities['Recorders'][rTag]['name'] == 'ELEMENT':
            eTags = sorted(list(set(elemSubdomain).intersection(Entities['Recorders'][rTag]['list'])))
            if eTags:
                ToProcessor['Recorders'][str(rTag)] = copy.deepcopy(Entities['Recorders'][rTag])
                ToProcessor['Recorders'][str(rTag)]['file'] = OUTFILE
                ToProcessor['Recorders'][str(rTag)]['list'] = eTags.copy()
        elif Entities['Recorders'][rTag]['name'] == 'SECTION':
            eTags = sorted(list(set(elemSubdomain).intersection(Entities['Recorders'][rTag]['list'])))
            if eTags:
                ToProcessor['Recorders'][str(rTag)] = copy.deepcopy(Entities['Recorders'][rTag])
                ToProcessor['Recorders'][str(rTag)]['file'] = OUTFILE
                ToProcessor['Recorders'][str(rTag)]['list'] = eTags.copy()
        elif Entities['Recorders'][rTag]['name'] == 'PARAVIEW':
            OUTFILE = Entities['Recorders'][rTag]['file']
            OUTFILE = OUTFILE.split('.')
            OUTFILE = OUTFILE[0] + '_PART' + str(k)
            ToProcessor['Recorders'][str(rTag)] = copy.deepcopy(Entities['Recorders'][rTag])
            ToProcessor['Recorders'][str(rTag)]['file'] = OUTFILE
            ToProcessor['Recorders'][str(rTag)]['features'] = Options['nparaview']

    #Gets the simulation for this partition
    for sTag in Entities['Simulations']:
        ctag = Entities['Simulations'][sTag]['combo']

        tag = Entities['Simulations'][sTag]['attributes']['analysis']
        analysis = Entities['Analyses'][tag]

        tag = Entities['Simulations'][sTag]['attributes']['algorithm']
        algorithm = Entities['Algorithms'][tag]

        tag = Entities['Simulations'][sTag]['attributes']['integrator']
        integrator = Entities['Integrators'][tag]

        tag = Entities['Simulations'][sTag]['attributes']['solver']
        solver = Entities['Solvers'][tag]

        if solver['name'] == 'PETSC':
            solver['d_nz'] = int(Options['d_nz'][k])
            solver['o_nz'] = int(Options['o_nz'][k])

        attributes = {'analysis': analysis, 'algorithm': algorithm, 'integrator': integrator, 'solver': solver}
        ToProcessor['Simulations'][sTag] = {'combo': ctag, 'attributes': attributes} 

    #Removes the empty fields
    keys = list(ToProcessor.keys())
    for key in keys:
        if not ToProcessor[key]:
            del ToProcessor[key]
    return ToProcessor

def createPartitions():
    """
    This function creates the partitions according with the pattern generated
    during the domain decomposition. Basically, goes over the Entities and 
    extract the information\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Returns
    -------
    None
    """
    #Creates the required folder
    createFolders()

    #Creates the domain decomposition input files
    SetMetisInputFile()

    #Reads the generated domain decomposition results
    GetMetisOutputFile()

    eTags = np.zeros(len(Entities['Elements']), dtype=np.uint32)
    for k, tag in enumerate(Entities['Elements']):
        eTags[k] = tag

    massSubdomain = set(Entities['Masses'].keys())

    #Writes the mesh file 
    for k in range(Options['nparts']):
        #Element that belong to this partition
        elemSubdomain = eTags[Options['partition'] == k]

        #Check if the partition has Element
        if len(elemSubdomain) == 0:
            print('\x1B[31m ERROR \x1B[0m: The partition for processor [%s] does not have Element' % k)
            print("\x1B[31m **************** THE PROCESS WILL BE ABORTED ****************\x1B[0m\n")
            sys.exit(-1)

        #Nodes that belong to this partition
        nodeSubdomain = set() 
        for eTag in elemSubdomain:
            connection = Entities['Elements'][eTag]['conn']
            for nTag in connection:
                nodeSubdomain.add(nTag)

        #Materials that belong to this partition
        matSubdomain = set()
        for eTag in elemSubdomain:
            if 'material' in Entities['Elements'][eTag]['attributes']:
                matSubdomain.add(Entities['Elements'][eTag]['attributes']['material']) 

        #Sections that belong to this partition
        secSubdomain = set() 
        for eTag in elemSubdomain:
            if 'section' in Entities['Elements'][eTag]['attributes']:
                sTag = Entities['Elements'][eTag]['attributes']['section']
                if Entities['Sections'][sTag]['model'] == 'PLAIN':
                    mTag = Entities['Sections'][sTag]['attributes']['material']
                    matSubdomain.add(mTag)
                elif  Entities['Sections'][sTag]['model'] == 'FIBER':
                    fTag = list(np.unique(np.array(Entities['Sections'][sTag]['attributes']['fiber'])))
                    matSubdomain.update(fTag)
                secSubdomain.add(sTag)

        #Constraints (Equal, General, Diaphragm) that belong to this partition
        conSubdomain = set() 
        for nTag in nodeSubdomain:
            FreeDofs = Entities['Nodes'][nTag]['freedof']
            for dof in FreeDofs:
                if dof < -1:
                    conSubdomain.add(dof)

        #Constraints information must be contained in this partition
        for cTag in conSubdomain:
            Master = Entities['Constraints'][cTag]['mtag']
            for mNode in Master:
                nodeSubdomain.add(mNode)

        #Surfaces that belong to this partition
        surfSubdomain = set() 
        for sTag in Entities['Surfaces']:
            eTag = Entities['Surfaces'][sTag]['etag']
            if eTag in elemSubdomain:
                surfSubdomain.add(sTag)

        #Sets the Entities that belong to this partition
        ToProcessor = Entities2Processor(matSubdomain,secSubdomain,nodeSubdomain,massSubdomain,conSubdomain,elemSubdomain,surfSubdomain,k)

        #Writes the partition in separated files
        dict2json(ToProcessor, k)
    
    #The generated partition file name (generic) path
    Options['execfile'] = Options['file'] + '.$.json'
    Options['execpath'] = Options['path'] + '/Partition'

    #SeismoVLAB execution command line
    Options['run'] = ' ' + Options['runanalysis'] + '/SeismoVLAB.exe -dir \'' + Options['execpath'] + '\' -file \'' + Options['execfile'] + '\'\n'
    if Options['nparts'] > 1:
        Options['run'] = ' mpirun -np ' + str(Options['nparts']) + Options['run']

    #Cleans generated auxiliary files
    os.remove(Options['execpath'] + '/Graph.out')
    os.remove(Options['execpath'] + '/Graph.out.epart.' + str(nparts))
    os.remove(Options['execpath'] + '/Graph.out.npart.' + str(nparts)) 

def checkWarnings():
    """
    This function checks consistency between the parameters specified by the
    user and the parameters employed in SeismoVLAB. A small report is printed
    showing all possible warnings encounter during the process. If they are 
    found, they are reported as an ALERT to be fixed.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Returns
    -------
    None
    """
    nTags = list(Entities['Nodes'].keys())
    sTags = list(Entities['Sections'].keys())
    mTags = list(Entities['Materials'].keys())
    eTags = list(Entities['Elements'].keys())

    print('\n Checking for warnings:')

    #[1] Check all attributes in NODES are defined in ENTITIES
    if not Entities['Nodes']:
        print("   *** There is no definition of Nodes. ***")
        print("\x1B[31m   *************** THE PROCESS WILL BE ABORTED ***************\x1B[0m\n")
        sys.exit(-1)

    nrestrain = 0
    for nTag in Entities['Nodes']:
        if math.isnan(nTag):
            print("   *** Node[%s] is invalid and should be removed ***" % nTag) 
        if Entities['Nodes'][nTag]['ndof'] == 0:
            print("   *** Node[%s] has ndof=0, fix this or delete it ***" % nTag) 
        ndim = len(Entities['Nodes'][nTag]['coords'])
        if Options['dimension'] != ndim:
            print("   *** Node[%s] coordinate's dimension (=%d) disagrees with Options[\'dimension\'] (=%d) ***" % (nTag,ndim,Options['dimension']))
        for free in Entities['Nodes'][nTag]['freedof']:
            if free == -1:
                nrestrain += 1

    if nrestrain == 0:
        print("   *** There is no restrains applied to Nodes ***")
        print("\x1B[31m   *************** THE PROCESS WILL BE ABORTED ***************\x1B[0m\n")
        sys.exit(-1)

    #[2] Check all attributes in MATERIALS are defined in ENTITIES
    for mTag in Entities['Materials']:
        name = Entities['Materials'][mTag]['name']
        if math.isnan(mTag):
            print("   *** Material[%s] is invalid and should be removed" % mTag) 
        if Entities['Materials'][mTag]['name'] not in SVLclasses['Materials']:
            print("   *** Material[%s] does not have an appropriate class name (%s) ***" % (mTag,Entities['Materials'][mTag]['name']))
        if name not in SVLclasses['Materials']:
            print("   *** Section[%s] does not have an appropriate class name (%s) ***" % (mTag,Entities['Materials'][mTag]['name'])) 
        elif Options['dimension'] not in SVLclasses['Materials'][name]['dim']:
            print("   *** Material[%s] cannot be used in a %sD space ***" % (mTag, Options['dimension']))

    #[3] Check all attributes in SECTIONS are defined in ENTITIES
    for sTag in Entities['Sections']:
        name = Entities['Sections'][sTag]['name']
        if math.isnan(sTag):
            print("   *** Section[%s] is invalid and should be removed" % sTag) 
        if Entities['Sections'][sTag]['name'] not in SVLclasses['Sections']:
            print("   *** Section[%s] does not have an appropriate class name (%s) ***" % (sTag,Entities['Sections'][sTag]['name']))
        if name not in SVLclasses['Sections']:
            print("   *** Section[%s] does not have an appropriate class name (%s) ***" % (sTag,Entities['Sections'][sTag]['name'])) 
        elif Options['dimension'] not in SVLclasses['Sections'][name]['dim']:
            print("   *** Section[%s] cannot be used in a %sD space ***" % (sTag, Options['dimension']))

    #[4] Check all attributes in ELEMENTS are defined in ENTITIES
    if not Entities['Elements']:
        print("   *** There is no definition of Elements. ***")
        print("\x1B[31m   *************** THE PROCESS WILL BE ABORTED ***************\x1B[0m\n")
        sys.exit(-1)

    naux = set(nTags)
    for eTag in Entities['Elements']:
        name = Entities['Elements'][eTag]['name']
        if math.isnan(eTag):
            print("   *** Element[%s] is invalid and should be removed" % nTag) 
        if name not in SVLclasses['Elements']:
            print("   *** Elements[%s] does not have an appropriate class name (%s) ***" % (eTag,Entities['Elements'][eTag]['name'])) 
        elif Options['dimension'] not in SVLclasses['Elements'][name]['dim']:
            print("   *** Element[%s] cannot be used in a %sD space ***" % (eTag, Options['dimension']))
        if 'material' in Entities['Elements'][eTag]['attributes']:
            mtag = Entities['Elements'][eTag]['attributes']['material']
            if mtag not in mTags:
                print("   *** Material[%s] has not been defined in Elements[%s] (%s) ***" % (mtag, eTag,Entities['Elements'][eTag]['name'])) 
        if 'section' in Entities['Elements'][eTag]['attributes']:
            stag = Entities['Elements'][eTag]['attributes']['section']
            if stag not in sTags:
                print("   *** Section[%s] has not been defined in Elements[%s] (%s) ***" % (stag, eTag,Entities['Elements'][eTag]['name']))
        defective = set(Entities['Elements'][eTag]['conn']).difference(naux)
        if defective:
            print('   *** Node[%s] have not been defined in Element[%s]' % (', '.join(str(s) for s in defective), eTag))

    #[5] Check all attributes in SURFACES are defined in ENTITIES
    for sTag in Entities['Surfaces']:
        eTag = Entities['Surfaces'][sTag]['etag']
        if eTag in Entities['Elements']:
            name = Entities['Elements'][eTag]['name']
            surf = Entities['Surfaces'][sTag]['conn']
            conn = np.array(Entities['Elements'][eTag]['conn'])
            Entities['Surfaces'][sTag]['face'] = SurfaceFace(SVLclasses['Elements'][name]['type'], surf, conn)
        else:
            print('   *** Surface[%s] does not belong to Element[%s] ***' % (sTag, eTag))

    #[6] Checks consistency between POINTLOAD/ELEMENTLOAD and other ENTITIES
    for lTag in Entities['Loads']:
        Name = Entities['Loads'][lTag]['name']
        if Name == 'POINTLOAD':
            #Prepare Node indexes for 'ALL' case in list
            if isinstance(Entities['Loads'][lTag]['attributes']['list'], str):
                Entities['Loads'][lTag]['attributes']['list'] = list(Entities['Nodes'].keys())

            #Checks consistency between POINTLOAD and degree-of-freedom in Node
            nTag = Entities['Loads'][lTag]['attributes']['list']
            fTag = Entities['Loads'][lTag]['attributes']['fun']
            nDIR = len(Entities['Functions'][fTag]['attributes']['dir'])
            for n in nTag:
                nDOF = Entities['Nodes'][n]['ndof']
                if nDOF != nDIR:
                    print('   *** Load[%s] (POINTLOAD) with Function[%s] (\'dir\') does not match Node[%s] (\'ndof\') ***' % (lTag,fTag,n))
            
            #Check if TimeSerie file can be opened (local and then global address)
            if 'file' in Entities['Functions'][fTag]['attributes']:
                #Tries the path given by user
                filename = Entities['Functions'][fTag]['attributes']['file']
                if tryOpenfile(filename):
                    #Tries a global path with respect to the main file
                    filepath = Options['path'] + '/' + Entities['Functions'][fTag]['attributes']['file']
                    if tryOpenfile(filepath):
                        lType = Entities['Loads'][lTag]['attributes']['type']
                        print('   *** POINTLOAD (%s) file=\'%s\' in Function[%s] could not be opened ***' % (lType,filename,fTag))
                    else:
                        Entities['Functions'][fTag]['attributes']['file'] = filepath
        elif Name == 'ELEMENTLOAD':
            #Prepare Node indexes for 'ALL' case in list
            if isinstance(Entities['Loads'][lTag]['attributes']['list'], str):
                Entities['Loads'][lTag]['attributes']['list'] = list(Entities['Elements'].keys())

            #Checks consistency between ELEMENTLOAD and model dimension in Options
            fTag = Entities['Loads'][lTag]['attributes']['fun']
            if 'dir' in Entities['Functions'][fTag]['attributes']:
                nDIR = len(Entities['Functions'][fTag]['attributes']['dir'])
                if nDIR != Options['dimension']:
                    print('   *** Load[%s] (ELEMENTLOAD) with Function[%s] (\'dir\') does not match Options (\'dimension\') ***' % (lTag,fTag))

            #Check if TimeSerie file can be opened
            if 'file' in Entities['Functions'][fTag]['attributes']:
                LOAD = Entities['Loads'][lTag]['attributes']['type'].upper()
                filename = Entities['Functions'][fTag]['attributes']['file']
                if LOAD == 'GENERALWAVE':
                    #Gets the DRM Nodes
                    cond = False
                    eTag = Entities['Loads'][lTag]['attributes']['list']
                    nTag = set()
                    for k in eTag:
                        nlist = Entities['Elements'][k]['conn']
                        for n in nlist:
                            nTag.add(n)
                    #Check files can be opened
                    for k in nTag:
                        fname = filename.replace("$", str(k))
                        #Tries the path given by user
                        if tryOpenfile(fname):
                            #Tries a global path with respect to the main file
                            filepath = Options['path'] + '/' + fname
                            if tryOpenfile(filepath):
                                print('   *** ELEMENTLOAD (%s) file=\'%s\' in Function[%s] for Node[%s] could not be opened ***' % (LOAD,fname,fTag,k))
                            else:
                                cond = True
                    if cond:
                        Entities['Functions'][fTag]['attributes']['file'] = Options['path'] + '/' + filename
                elif LOAD == 'PLANEWAVE':
                    #Tries the path given by user
                    if tryOpenfile(filename):
                        #Tries a global path with respect to the main file
                        filepath = Options['path'] + '/' + Entities['Functions'][fTag]['attributes']['file']
                        if tryOpenfile(filepath):
                            print('   *** ELEMENTLOAD (%s) file=\'%s\' in Function[%s] could not be opened ***' % (LOAD,filename,fTag))
                        else:
                            filename = filepath
                            Entities['Functions'][fTag]['attributes']['file'] = filepath
                    #Opens the given file
                    with open(filename, "r") as fileHandler:
                        lines = fileHandler.readlines()
                        nDims = Options['dimension']
                        #Gets the number of DRM nodes
                        line  = list(filter(None, lines[0].strip().split())) 
                        nNode = int(line[2])
                        #Gets the DRM Nodes
                        eTag = Entities['Loads'][lTag]['attributes']['list']
                        allTag = set()
                        for k in eTag:
                            nlist = Entities['Elements'][k]['conn']
                            for n in nlist:
                                allTag.add(n)
                        #Gets the DRM node Tags
                        m = 0
                        nTags = np.empty(shape=(nNode,), dtype=int)
                        for k in range(1, nNode+1):
                            line = list(filter(None, lines[k].strip().split()))
                            nTags[m] = int(line[0])
                            m += 1
                        #Gets the the most distant coordinate
                        xmin = np.empty(shape=(nNode,nDims))
                        for k in range(nNode):
                            xmin[k] = Entities['Nodes'][nTags[k]]['coords']
                        Entities['Functions'][fTag]['attributes']['xmin'] = xmin.min(axis=0)

                        undefined = allTag.difference(nTags)
                        if undefined:
                            print('   *** ELEMENTLOAD (%s) in file=\'%s\' not all DRM nodes have been specified ***' % (LOAD,filename))
                elif LOAD == 'BODY':
                    #Tries the path given by user
                    if tryOpenfile(filename):
                        #Tries a global path with respect to the main file
                        filepath = Options['path'] + '/' + Entities['Functions'][fTag]['attributes']['file']
                        if tryOpenfile(filepath):
                            print('   *** ELEMENTLOAD (%s) file=\'%s\' in Function[%s] could not be opened ***' % (LOAD,filename,fTag))
                        else:
                            Entities['Functions'][fTag]['attributes']['file'] = filepath
                elif LOAD == 'SURFACE':
                    #Tries the path given by user
                    if tryOpenfile(filename):
                        #Tries a global path with respect to the main file
                        filepath = Options['path'] + '/' + Entities['Functions'][fTag]['attributes']['file']
                        if tryOpenfile(filepath):
                            print('   *** ELEMENTLOAD (%s) file=\'%s\' in Function[%s] could not be opened ***' % (LOAD,filename,fTag))
                        else:
                            Entities['Functions'][fTag]['attributes']['file'] = filepath
        elif Name == 'SUPPORTMOTION':
            nTag = Entities['Loads'][lTag]['attributes']['list']
            #Check if TimeSerie file can be opened
            for n in nTag:
                if 'file' in Entities['Supports'][n]:
                    filenames = Entities['Supports'][n]['file']
                    for ftag, filepath in enumerate(filenames):
                        #Tries the path given by user
                        if tryOpenfile(filepath):
                            #Tries a global path with respect to the main file
                            filename = Options['path'] + '/' + filepath
                            if tryOpenfile(filename):
                                print('   *** SUPPORTMOTION file=\'%s\' in Supports[%s] could not be opened ***' % (filepath,ftag))
                            else:
                                Entities['Supports'][n]['file'][ftag] = filename

    #[7] Check all attributes in DAMPING are defined in ENTITIES
    dlist = list(Entities['Elements'].keys())
    if not Entities['Dampings']:
        Entities['Dampings'][1] = {'name': 'FREE', 'attributes': {'list': eTags}}
    else:
        for dTag in Entities['Dampings']:
            if isinstance(Entities['Dampings'][dTag]['attributes']['list'], str):
                if Entities['Dampings'][dTag]['attributes']['list'].upper() == 'ALL':
                    Entities['Dampings'][dTag]['attributes']['list'] = eTags
                else:
                    print("   *** Attribute: 'list'=%d in Dampings[%s] is not recognized ***" % (dTag, Entities['Dampings'][dTag]['attributes']['list']))
            else:
                if Entities['Dampings'][dTag]['name'] == 'RAYLEIGH':
                    if 'am' not in Entities['Dampings'][dTag]['attributes']:
                        Entities['Dampings'][dTag]['attributes']['am'] = 0.0
                    if 'ak' not in Entities['Dampings'][dTag]['attributes']:
                        Entities['Dampings'][dTag]['attributes']['ak'] = 0.0
                elif Entities['Dampings'][dTag]['name'] == 'CAUGHEY':
                    Entities['Dampings'][dTag]['name'] = 'FREE'
                    print("   *** CAUGHEY in Dampings[%s] is not implemented yet ***" % dTag)
                elif Entities['Dampings'][dTag]['name'] == 'CAPPED':
                    Entities['Dampings'][dTag]['name'] = 'FREE'
                    print("   *** CAPPED in Dampings[%s] is not implemented yet ***" % dTag)
            #Finds elements without damping 
            daux  = dlist
            dlist = list(set(daux).difference(Entities['Dampings'][dTag]['attributes']['list']))

        #Creates a new free damping with the elements without damping
        if dlist:
            dTag = 1 + max(list( Entities['Dampings'].keys()))
            Entities['Dampings'][dTag] =  {'name': 'FREE', 'list': dlist}

    #[8] Check all attributes in RECORDERS are defined in ENTITIES
    for rTag in Entities['Recorders']:
        if 'list' in Entities['Recorders'][rTag]:
            if isinstance(Entities['Recorders'][rTag]['list'], str):
                if Entities['Recorders'][rTag]['name'] == 'NODE':
                    Entities['Recorders'][rTag]['list'] = nTags
                elif Entities['Recorders'][rTag]['name'] == 'SECTION':
                    Entities['Recorders'][rTag]['list'] = eTags
                elif Entities['Recorders'][rTag]['name'] == 'ELEMENT':
                    Entities['Recorders'][rTag]['list'] = eTags
        if 'response' in Entities['Recorders'][rTag]:
            if Entities['Recorders'][rTag]['response'] == 'REACTION':
                #TODO: Get all nodes that are fixed
                Entities['Recorders'][rTag]['list'] = [] 
        if Entities['Recorders'][rTag]['name'] == 'PARAVIEW':
            dirName = Options['path'] + '/' + 'Paraview'
            if not os.path.exists(dirName):
                os.mkdir(dirName)

    #[9] Check all attributes in SOLVER are defined in ENTITIES
    for sTag in Entities['Solvers']:
        if Entities['Solvers'][sTag]['name'] == 'PETSC':
            if Options['allocation'] == 'NO':
                print('   *** Solver[%s] uses PETSC (parallel) and memory allocation has not being computed ***' % sTag)
        elif Entities['Solvers'][sTag]['name'] == 'MUMPS':
            if Options['nparts'] == 1:
                print('   *** Solver[%s] uses MUMPS (parallel) for number of partition %d, we recommend using EIGEN (serial) ***' % (sTag,Options['nparts']))
        elif Entities['Solvers'][sTag]['name'] == 'EIGEN':
            if Options['nparts'] != 1:
                print('   *** Solver[%s] uses EIGEN (serial) and number of partition is %d (parallel) ***' % (sTag,Options['nparts']))

    #[10] Check all attributes in SIMULATION are defined in ENTITIES
    for sTag in Entities['Simulations']:
        if Entities['Simulations'][sTag]['attributes']['analysis'] not in Entities['Analyses']:
            print('   *** Simulation[%s] has no defined analysis  ***' % sTag)
        if Entities['Simulations'][sTag]['attributes']['algorithm'] not in Entities['Algorithms']:
            print('   *** Simulation[%s] has no defined algorithm  ***' % sTag)
        if Entities['Simulations'][sTag]['attributes']['integrator'] not in Entities['Integrators']:
            print('   *** Simulation[%s] has no defined integrator  ***' % sTag)
        if Entities['Simulations'][sTag]['attributes']['solver'] not in Entities['Solvers']:
            print('   *** Simulation[%s] has no defined solver ***' % sTag)
        if Entities['Simulations'][sTag]['combo'] not in Entities['Combinations']:
            print('   *** Simulation[%s] has no defined combination ***' % sTag)
           
    print(' Done checking!\n')
    Options['wasChecked'] = True

def CreateRunAnalysisFiles(plot=False):
    """
    This function gathers provided model information and post-process them 
    generating constraints, fiber sections, degree of freedom numbering and
    partitions to be written in Run-Analysis JSON format .\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Returns
    -------
    None
    """
    #Creates the fiber section
    GenerateFiberSection()

    #Enforce diaphragm constraints
    ApplyConstraints()

    #Check if the model is properly done
    if not Options['wasChecked']:
        checkWarnings()

    #Generate DRM input files
    GenerateDRMFiles()

    #Set degree of freedom
    setDegreeOfFreedom(plot)

    #Generate the Entities group
    createPartitions()

    #Prints SVL Run-Analysis execution
    print(Options['run'])


#Functions to be run when SeismoVLAB is imported
printHeader()

setFilePath()
