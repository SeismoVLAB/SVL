#!/usr/bin/env python3
# -*- coding: Utf-8 -*-

import shlex
import copy
import numpy as np
from Core.Utilities import *

def parseETABS(filepath):
    """
    This function parses a *.e2k files from ETABS CSI software. Not all functionalities 
    are available, only coordinates, elements, material and sections are obtained
    during the parsing process.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    filepath: str
        The path were the ETABS file will be read

    Returns
    -------
    mesh : dict
        A dictionary that contains nodes and element with the same structure as Entities
    """
    #The MESH dictionary containing the data is read
    mesh = {'Nodes': {}, 'Materials': {}, 'Sections':{}, 'Elements':{}}

    #Map for line and area identifiers
    mapETABS = {'Material': {}, 'Section': {}}

    #The MESH dictionary containing the data is read
    ETABS = {'STORY': {}, 'DIAPH': {}, 'MATERIAL': {}, 'FRAMESEC': {}, 'SHELLPROP': {}, 'POINT': {}, 'LINE': {}, 'AREA': {}}

    #Open the provided file
    with open(filepath, 'r', errors='replace') as f:
        lines = f.readlines()

    k = 0
    mSet = set()
    sSet = set()
    fwdStory2Tag = {}
    invStory2Tag = {}
    while k < len(lines):
        line = lines[k].strip()
        if line == "$ STORIES - IN SEQUENCE FROM TOP":
            sTag = 0
            while True:
                k += 1
                line = list(shlex.split(lines[k]))
                if line:
                    if line[0].upper() == '$':
                        k -= 1
                        break
                    sTag += 1
                    NAME = line[1]
                    fwdStory2Tag[sTag] = NAME
                    invStory2Tag[NAME] = sTag
                    if NAME == 'BASE':
                        ETABS['STORY'][NAME] = {'ELEV': float(line[3]), 'PASI': {}, 'LASI': {}, 'AASI': {}}
                    else:
                        ETABS['STORY'][NAME] = {'HEIGHT': float(line[3]), 'PASI': {}, 'LASI': {}, 'AASI': {}}
        elif line == "$ DIAPHRAGM NAMES":
            while True:
                k += 1
                line = list(shlex.split(lines[k]))
                if line:
                    if line[0].upper() == '$':
                        k -= 1
                        break
                    NAME = line[1]
                    if NAME not in ETABS['DIAPH']:
                        ETABS['DIAPH'][NAME] = {}
                    #Adds the parameters
                    nparam = (len(line) - 2)//2
                    for j in range(nparam):
                        PARAM = line[2 + j*2]
                        VALUE = line[3 + j*2]
                        ETABS['DIAPH'][NAME][PARAM] = VALUE
        elif line == "$ MATERIAL PROPERTIES":
            matTag = 0
            while True:
                k += 1
                line = list(shlex.split(lines[k]))
                if line:
                    if line[0].upper() == '$':
                        k -= 1
                        break
                    #LINE element material
                    NAME = line[1] + '_1D'
                    if NAME not in ETABS['MATERIAL']:
                        matTag += 1
                        ETABS['MATERIAL'][NAME] = {}
                        mapETABS['Material'][NAME] = matTag
                    nparam = (len(line) - 2)//2
                    for j in range(nparam):
                        PARAM = line[2 + j*2]
                        if PARAM in ['M','E','U','TYPE']:
                            VALUE = line[3 + j*2]
                            ETABS['MATERIAL'][NAME][PARAM] = VALUE
                    #AREA element material
                    NAME = line[1] + '_2D'
                    if NAME not in ETABS['MATERIAL']:
                        matTag += 1
                        ETABS['MATERIAL'][NAME] = {}
                        mapETABS['Material'][NAME] = matTag
                    nparam = (len(line) - 2)//2
                    for j in range(nparam):
                        PARAM = line[2 + j*2]
                        if PARAM in ['M','E','U','TYPE']:
                            VALUE = line[3 + j*2]
                            ETABS['MATERIAL'][NAME][PARAM] = VALUE
                    #VOLUME element material
                    NAME = line[1] + '_2D'
                    if NAME not in ETABS['MATERIAL']:
                        matTag += 1
                        ETABS['MATERIAL'][NAME] = {}
                        mapETABS['Material'][NAME] = matTag
                    nparam = (len(line) - 2)//2
                    for j in range(nparam):
                        PARAM = line[2 + j*2]
                        if PARAM in ['M','E','U','TYPE']:
                            VALUE = line[3 + j*2]
                            ETABS['MATERIAL'][NAME][PARAM] = VALUE
        elif line == "$ FRAME SECTIONS":
            secTag = 0
            while True:
                k += 1
                line = list(shlex.split(lines[k]))
                if line:
                    if line[0].upper() == '$':
                        k -= 1
                        break
                    NAME = line[1]
                    if NAME not in ETABS['FRAMESEC']:
                        secTag += 1
                        ETABS['FRAMESEC'][NAME] = {}
                        mapETABS['Section'][NAME] = secTag
                    #Adds the parameters
                    nparam = (len(line) - 2)//2
                    for j in range(nparam):
                        PARAM = line[2 + j*2]
                        if PARAM in ['MATERIAL','SHAPE','D','B','TF','TW','T']:
                            VALUE = line[3 + j*2]
                            ETABS['FRAMESEC'][NAME][PARAM] = VALUE
        elif line == "$ WALL/SLAB/DECK PROPERTIES":
            secTag = len(mapETABS['Section'])
            while True:
                k += 1
                line = list(shlex.split(lines[k]))
                if line:
                    if line[0].upper() == '$':
                        k -= 1
                        break
                    NAME = line[1]
                    if NAME not in ETABS['SHELLPROP']:
                        secTag += 1
                        ETABS['SHELLPROP'][NAME] = {}
                        mapETABS['Section'][NAME] = secTag
                    #Adds the parameters
                    nparam = (len(line) - 2)//2
                    for j in range(nparam):
                        PARAM = line[2 + j*2]
                        if PARAM in ['MATERIAL','PROPTYPE','TM','TB']:
                            VALUE = line[3 + j*2]
                            ETABS['SHELLPROP'][NAME][PARAM] = VALUE
        elif line.upper() == '$ POINT COORDINATES':
            while True:
                k += 1
                line = list(shlex.split(lines[k]))
                if line:
                    if line[0].upper() == '$':
                        k -= 1
                        break
                    NAME = line[1]
                    if NAME not in ETABS['POINT']:
                        if len(line) == 4:
                            coords = np.array([float(line[2]), float(line[3]), np.nan])
                            ETABS['POINT'][NAME] = {'COORDS': coords}
                        else:
                            coords = np.array([float(line[2]), float(line[3]), float(line[4])])
                            ETABS['POINT'][NAME] = {'COORDS': coords}
                    else:
                        print('\x1B[33m In parseETABS() POINT[%s] was already defined at LINE=%d in parsed file.' %(NAME,k))
        elif line.upper() == '$ LINE CONNECTIVITIES':
            while True:
                k += 1
                line = list(shlex.split(lines[k]))
                if line:
                    if line[0].upper() == '$':
                        k -= 1
                        break
                    NAME = line[1]
                    TYPE = line[2]
                    CONN = [line[3], line[4]]
                    NEXT = [int(line[5]), 0]
                    if NAME not in ETABS['LINE']:
                        ETABS['LINE'][NAME] = {'TYPE': TYPE, 'NPTS': 2, 'CONN': CONN, 'NEXT': NEXT}
                    else:
                        print('\x1B[33m ALERT \x1B[0m: In parseETABS() LINE[%s] was already defined at LINE=%d in parsed file.' %(NAME,k))
        elif line.upper() == '$ AREA CONNECTIVITIES':
            while True:
                k += 1
                line = list(shlex.split(lines[k]))
                if line:
                    if line[0].upper() == '$':
                        k -= 1
                        break
                    NAME = line[1]
                    TYPE = line[2]
                    NPTS = int(line[3])
                    if NAME not in ETABS['AREA']:
                        if NPTS == 3:
                            CONN = [line[4], line[5], line[6]]
                            NEXT = [int(line[7]), int(line[8]), int(line[9])]
                        elif NPTS == 4:
                            CONN = [line[4], line[5], line[6], line[7]]
                            NEXT = [int(line[8]), int(line[9]), int(line[10]), int(line[11])]
                        else:
                            print('\x1B[33m ALERT \x1B[0m: AREA=%s in LINE=%s has more than 4 points.' %(NAME,k))
                    ETABS['AREA'][NAME] = {'TYPE': TYPE, 'NPTS': NPTS, 'CONN': CONN, 'NEXT': NEXT}
        elif line.upper() == '$ POINT ASSIGNS':
            while True:
                k += 1
                line = list(shlex.split(lines[k]))
                if line:
                    if line[0][0].upper() == '$':
                        k -= 1
                        break
                    PTAG = line[1]
                    STAG = line[2]
                    if PTAG not in ETABS['STORY'][STAG]['PASI']:
                        ETABS['STORY'][STAG]['PASI'][PTAG] = {}
                    #Adds the parameters
                    nparam = (len(line) - 3)//2
                    for j in range(nparam):
                        PARAM = line[3 + j*2]
                        if PARAM in ['DIAPH','RESTRAINT']:
                            ETABS['STORY'][STAG]['PASI'][PTAG][PARAM] = line[4 + j*2]
        elif line.upper() == '$ LINE ASSIGNS':
            while True:
                k += 1
                line = list(shlex.split(lines[k]))
                if line:
                    if line[0][0].upper() == '$':
                        k -= 1
                        break
                    LTAG = line[1]
                    STAG = line[2]
                    if LTAG not in ETABS['STORY'][STAG]['LASI']:
                        ETABS['STORY'][STAG]['LASI'][LTAG] = {}
                    #Adds the parameters
                    nparam = (len(line) - 3)//2
                    for j in range(nparam):
                        PARAM = line[3 + j*2]
                        if PARAM in ['SECTION','ANG']:
                            VALUE = line[4 + j*2]
                            ETABS['STORY'][STAG]['LASI'][LTAG][PARAM] = VALUE
        elif line.upper() == '$ AREA ASSIGNS':
            while True:
                k += 1
                line = list(shlex.split(lines[k]))
                if line:
                    if line[0][0].upper() == '$':
                        k -= 1
                        break
                    ATAG = line[1]
                    STAG = line[2]
                    if ATAG not in ETABS['STORY'][STAG]['AASI']:
                        ETABS['STORY'][STAG]['AASI'][ATAG] = {}
                    #Adds the parameters
                    nparam = (len(line) - 3)//2
                    for j in range(nparam):
                        PARAM = line[3 + j*2]
                        if PARAM in ['SECTION','DIAPH']:
                            ETABS['STORY'][STAG]['AASI'][ATAG][PARAM] = line[4 + j*2]
        elif line.upper() == '$ ANALYSIS OPTIONS':
            while True:
                k += 1
                line = list(shlex.split(lines[k]))
                if line:
                    if line[0][0].upper() == '$':
                        k -= 1
                        break
                    ATAG = line[0]
                    if ATAG == 'ACTIVEDOF':
                        DOF = list(filter(None, line[1].strip().split(' ')))
                        ndof = len(DOF)
                        if 'UX' in DOF and 'UY' in DOF and 'UZ' in DOF:
                            ndim = 3
                        if 'UX' in DOF and 'UY' in DOF and 'UZ' not in DOF:
                            ndim = 2
        else:
            k += 1

    #CREATES THE ETABS STORY-POINT-TAG MAPPING
    PointMap = {}
    for sTag in ETABS['STORY']:
        keys = list(ETABS['POINT'].keys())
        PointMap[sTag] = dict.fromkeys(keys, -1)

    tag = 0
    for sTag in ETABS['STORY']:
        for pTag in ETABS['STORY'][sTag]['PASI']:
            if PointMap[sTag][pTag] == -1:
                tag += 1  
                PointMap[sTag][pTag] = tag
        for lTag in ETABS['STORY'][sTag]['LASI']:
            CONN = ETABS['LINE'][lTag]['CONN']
            NEXT = ETABS['LINE'][lTag]['NEXT']
            for k, pTag in enumerate(CONN):
                STORY = fwdStory2Tag[invStory2Tag[sTag] + NEXT[k]]
                if PointMap[STORY][pTag] == -1:
                    tag += 1 
                    PointMap[STORY][pTag] = tag
        for aTag in ETABS['STORY'][sTag]['AASI']:
            CONN = ETABS['AREA'][aTag]['CONN']
            NEXT = ETABS['AREA'][aTag]['NEXT']
            for k, pTag in enumerate(CONN):
                STORY = fwdStory2Tag[invStory2Tag[sTag] + NEXT[k]]
                if PointMap[STORY][pTag] == -1:
                    tag += 1 
                    PointMap[STORY][pTag] = tag

    auxiliar = copy.deepcopy(PointMap)
    for sTag in auxiliar:
        for nTag in auxiliar[sTag]:
            if auxiliar[sTag][nTag] == -1:
                del PointMap[sTag][nTag]

    #COMPUTES THE STORY ELEVATIONS
    nSTORY = len(ETABS['STORY'])
    STORIES = list(range(nSTORY-1, 0, -1))
    for k in STORIES:
        ELEV   = ETABS['STORY'][fwdStory2Tag[k+1]]['ELEV']
        HEIGHT = ETABS['STORY'][fwdStory2Tag[k]]['HEIGHT']
        ETABS['STORY'][fwdStory2Tag[k]]['ELEV'] = ELEV + HEIGHT

    #GENERATES LIST OF MATERIALS
    for NAME in ETABS['MATERIAL']:
        rho = 0.0
        Em  = 0.0
        nu  = 0.0
        #SeismoVLAB material class name
        if NAME[-3:] == '_1D':
            mname = "ELASTIC1DLINEAR"
        elif NAME[-3:] == '_2D':
            mname = "ELASTIC2DPLANESTRESS"
        elif NAME[-3:] == '_3D':
            mname = "ELASTIC3DLINEAR"
        #Construct the material attribute
        mTag = mapETABS['Material'][NAME]
        if 'M' in ETABS['MATERIAL'][NAME]:
            rho = ETABS['MATERIAL'][NAME]['M']
        if 'E' in ETABS['MATERIAL'][NAME]:
            Em = ETABS['MATERIAL'][NAME]['E']
        if 'U' in ETABS['MATERIAL'][NAME]:
            nu = ETABS['MATERIAL'][NAME]['U']
        mesh['Materials'][mTag] = {'name': mname, 'E': Em, 'nu': nu, 'rho': rho}

    #GENERATES LIST OF SECTIONS
    for NAME in ETABS['FRAMESEC']:
        if 'MATERIAL' in ETABS['FRAMESEC'][NAME]:
            MNAME = ETABS['FRAMESEC'][NAME]['MATERIAL'] + '_1D'
            sTag = mapETABS['Section'][NAME]
            mTag = mapETABS['Material'][MNAME]
            mSet.add(mTag)
            
            #Assign section attributes
            attributes = {}
            if ETABS['FRAMESEC'][NAME]['SHAPE'] == "Rectangular":
                SecName = 'LIN' + str(ndim) + 'DRECTANGULAR'
                attributes['h'] = float(ETABS['FRAMESEC'][NAME]['D'])
                attributes['b'] = float(ETABS['FRAMESEC'][NAME]['B'])
            elif ETABS['FRAMESEC'][NAME]['SHAPE'] == "I/Wide Flange":
                SecName = 'LIN' + str(ndim) + 'DWIDEFLANGE'
                attributes['h'] = float(ETABS['FRAMESEC'][NAME]['D'])
                attributes['b'] = float(ETABS['FRAMESEC'][NAME]['B'])
                attributes['tf'] = float(ETABS['FRAMESEC'][NAME]['TF'])
                attributes['tw'] = float(ETABS['FRAMESEC'][NAME]['TW'])
            elif ETABS['FRAMESEC'][NAME]['SHAPE'] == "Box/Tube":
                SecName = 'LIN' + str(ndim) + 'DRECTANGULARTUBE'
                attributes['h'] = float(ETABS['FRAMESEC'][NAME]['D'])
                attributes['b'] = float(ETABS['FRAMESEC'][NAME]['B'])
                attributes['tf'] = float(ETABS['FRAMESEC'][NAME]['TF'])
                attributes['tw'] = float(ETABS['FRAMESEC'][NAME]['TW'])
            elif ETABS['FRAMESEC'][NAME]['SHAPE'] == "CIRCLE":
                SecName = 'LIN' + str(ndim) + 'DCIRCULAR'
                attributes['r'] = float(ETABS['FRAMESEC'][NAME]['D'])/2.0
            elif ETABS['FRAMESEC'][NAME]['SHAPE'] == "Pipe":
                SecName = 'LIN' + str(ndim) + 'DCIRCULARTUBE'
                attributes['re'] = float(ETABS['FRAMESEC'][NAME]['D'])/2.0
                attributes['ri'] = attributes['re'] - float(ETABS['FRAMESEC'][NAME]['T'])

            #Default Section rotation/Insertion point
            attributes['theta'] = 0.0
            attributes['ip'] = 10

            mesh['Sections'][sTag] = {'name': SecName, 'model': 'Plain', 'material': mTag, 'attributes': attributes}

    for NAME in ETABS['SHELLPROP']:
        if 'MATERIAL' in ETABS['SHELLPROP'][NAME]:
            MNAME = ETABS['SHELLPROP'][NAME]['MATERIAL'] + '_2D'
            sTag = mapETABS['Section'][NAME]
            mTag = mapETABS['Material'][MNAME]
            mSet.add(mTag)

            #Assign section attributes
            th = 0.0
            if 'TB' in ETABS['SHELLPROP'][NAME]:
                th = ETABS['SHELLPROP'][NAME]['TB']
            attributes = {'th': th}

            mesh['Sections'][sTag] = {'name': 'LIN3DTHINAREA', 'model': 'Plain', 'material': mTag, 'attributes': attributes}
    
    #CREATE COORDINATE MATRIX
    for sTag in PointMap:
        ELEV = ETABS['STORY'][sTag]['ELEV']
        for pTag in PointMap[sTag]:
            coords = np.copy(ETABS['POINT'][pTag]['COORDS'])
            if np.isnan(coords[2]):
                coords[2] = ELEV
            else:
                coords[2] = ELEV - coords[2]
            nTag = PointMap[sTag][pTag]
            mesh['Nodes'][nTag] = {'ndof': ndof, 'freedof': np.zeros(ndof, dtype=int), 'totaldof': np.zeros(ndof, dtype=int), 'coords': coords}

    #CREATE THE ELEMENT CONNECTIVITY
    eTag = 0
    for sTag in ETABS['STORY']:
        for lTag in ETABS['STORY'][sTag]['LASI']:
            CONN = ETABS['LINE'][lTag]['CONN']
            NEXT = ETABS['LINE'][lTag]['NEXT']
            #The associated SeismoVLAB element class
            NAME = 'LIN' + str(ndim) + 'DFRAME2'
            SECT = ETABS['STORY'][sTag]['LASI'][lTag]['SECTION']
            if SECT != 'NONE':                           
                section = mapETABS['Section'][SECT]
                sSet.add(section)
            else:
                print('\x1B[33m ALERT \x1B[0m: LINE=%s in STORY=%s has SECTION=\'NONE\' assigned.' %(lTag,sTag))
                section = np.nan
            #Constructs the connectivity array of this line element
            conn = []
            for k, nTag in enumerate(CONN):
                STORY = fwdStory2Tag[invStory2Tag[sTag] + NEXT[k]]
                JOINT = PointMap[STORY][nTag]
                conn.append(JOINT)

            eTag += 1
            attributes = {'rule': 'GAUSS', 'np': 3, 'formulation':'Bernoulli','section': section}
            mesh['Elements'][eTag] = {'name': NAME, 'conn': conn, 'attributes': attributes}

        for aTag in ETABS['STORY'][sTag]['AASI']:
            CONN = ETABS['AREA'][aTag]['CONN']
            NEXT = ETABS['AREA'][aTag]['NEXT']
            NPTS = ETABS['AREA'][aTag]['NPTS']
            #The associated SeismoVLAB element class
            NAME = 'LIN' + str(ndim) + 'DSHELL' + str(NPTS)
            nQp  = 7*(NPTS == 3) + 9*(NPTS == 4)
            SECT = ETABS['STORY'][sTag]['AASI'][aTag]['SECTION']
            if SECT != 'NONE':                      
                section = mapETABS['Section'][SECT]
                sSet.add(section)
            else:
                print('\x1B[33m ALERT \x1B[0m: AREA=%s in STORY=%s has SECTION=\'NONE\' assigned.' %(aTag,sTag))
                section = np.nan
            #Constructs the connectivity array of this area element
            conn = []
            for k, nTag in enumerate(CONN):
                STORY = fwdStory2Tag[invStory2Tag[sTag] + NEXT[k]]
                JOINT = PointMap[STORY][nTag]
                conn.append(JOINT)
            
            if NPTS == 3 or NPTS == 4:
                eTag += 1
                attributes = {'rule': 'GAUSS', 'np': nQp, 'section': section}
                mesh['Elements'][eTag] = {'name': NAME, 'conn': conn, 'attributes': attributes}

    #Cleans the unused materials
    mTags = list(mesh['Materials'].keys())
    for m in mTags:
        if m not in mSet:
            del mesh['Materials'][m]

    #Cleans the unused sectionss
    sTags = list(mesh['Sections'].keys())
    for s in sTags:
        if s not in sSet:
            del mesh['Sections'][s]
    
    return mesh