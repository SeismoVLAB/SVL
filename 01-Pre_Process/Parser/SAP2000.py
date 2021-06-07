#!/usr/bin/env python3
# -*- coding: Utf-8 -*-

import shlex
import numpy as np

def parseSAP(filepath):
    """
    This function parses a *.s2k files from SAP2000 CSI software. Not all functionalities 
    are available, only coordinates, elements, links, material and sections are obtained
    during the parsing process.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    filepath: str
        The path were the SAP2000 file will be read

    Returns
    -------
    mesh : dict
        A dictionary that contains nodes and element with the same structure as Entities
    """
    #The MESH dictionary containing the data is read
    mesh = {'Nodes': {}, 'Materials': {}, 'Sections':{}, 'Elements':{}, 'Diaphragms': {}}

    #Map for line and area identifiers
    mapSAP = {'Joint': {}, 'Diaph': {}, 'Material': {}, 'Section': {}, 'Link': {}, 'Frame': {}, 'Area': {}}

    #Open the provided file
    with open(filepath, 'r', errors='replace') as f:
        lines = f.readlines()

    #Loop over the entire file
    k = 0
    eTag = 0
    ndof = 0
    ndim = 0
    mSet = set()
    sSet = set()
    while k < len(lines):
        line = lines[k].strip()
        if line == 'TABLE:  "ACTIVE DEGREES OF FREEDOM"':
            while True:
                k += 1
                line = list(shlex.split(lines[k]))
                if line:
                    if line[0] == 'TABLE:':
                        break
                    UX = list(filter(None, line[0].strip().split('=')))
                    UY = list(filter(None, line[1].strip().split('=')))
                    UZ = list(filter(None, line[2].strip().split('=')))
                    RX = list(filter(None, line[3].strip().split('=')))
                    RY = list(filter(None, line[4].strip().split('=')))
                    RZ = list(filter(None, line[5].strip().split('=')))
                    dof = [UX[1],UY[1],UZ[1],RX[1],RY[1],RZ[1]]

                    #Sets the dimenion of the model
                    if UX[1].upper() == 'YES' and UY[1].upper() == 'YES' and UZ[1].upper() == 'NO':
                        ndim = 2
                    elif UX[1].upper() == 'YES' and UY[1].upper() == 'YES' and UZ[1].upper() == 'YES':
                        ndim = 3

                    #Sets the number of degree of freedom per node
                    for c in dof:
                        if c.upper() == 'YES':
                            ndof += 1
        elif line == 'TABLE:  "MATERIAL PROPERTIES 02 - BASIC MECHANICAL PROPERTIES"':
            mTag = 0
            while True:
                k += 1
                line = list(shlex.split(lines[k]))
                if line:
                    if line[0] == 'TABLE:':
                        break
                    Name = list(filter(None, line[0].strip().split('=')))
                    Mass = list(filter(None, line[2].strip().split('=')))
                    #Standar linear material
                    if  len(line) == 7:
                        E1  = list(filter(None, line[3].strip().split('=')))
                        U12 = list(filter(None, line[5].strip().split('=')))
                        E1  = E1[1]
                        U12 = U12[1]
                    elif  len(line) == 5:
                        E1 = list(filter(None, line[3].strip().split('=')))
                        E1  = E1[1]
                        U12 = '0.0'

                    #Creates a material for possible 1D/2D/3D elements
                    prop = {'E': float(E1), 'nu': float(U12), 'rho': float(Mass[1])}
                    mTag += 1
                    mapSAP['Material'][Name[1]+'_1D'] = mTag
                    mesh['Materials'][mTag] = {'name': 'ELASTIC1DLINEAR', 'attributes': prop}
                    mTag += 1
                    mapSAP['Material'][Name[1]+'_2D'] = mTag
                    mesh['Materials'][mTag] = {'name': 'ELASTIC2DPLANESTRESS', 'attributes': prop}
                    mTag += 1
                    mapSAP['Material'][Name[1]+'_3D'] = mTag
                    mesh['Materials'][mTag] = {'name': 'ELASTIC3DLINEAR', 'attributes': prop}
        elif line == 'TABLE:  "FRAME SECTION PROPERTIES 01 - GENERAL"':
            sTag = 0
            while True:
                k += 1
                line = list(shlex.split(lines[k]))
                if line:
                    if line[0] == 'TABLE:':
                        break

                    #PATCH for stupid break line format in SAP
                    if '_' in line:
                        while line[-1] == '_':
                            line.remove('_')
                            k += 1
                            line += list(filter(None, lines[k].strip().split(' ')))

                    sName = list(filter(None, line[0].strip().split('=')))
                    mName = list(filter(None, line[1].strip().split('=')))
                    Shape = list(filter(None, line[2].strip().split('=')))

                    sName = sName[1]
                    mName = mName[1]+'_1D'
                    sTag += 1
                    mapSAP['Section'][sName] = sTag
                    mTag = mapSAP['Material'][mName]

                    if Shape[1] == 'Rectangular':
                        H = list(filter(None, line[3].strip().split('=')))
                        B = list(filter(None, line[4].strip().split('=')))
                        NAME = 'LIN' + str(ndim) + 'DRECTANGULAR'

                        mSet.add(mTag)
                        prop = {'material': mTag, 'h': float(H[1]),'b': float(B[1]), 'theta': 0.0, 'ip': 10}
                        mesh['Sections'][sTag] = {'name': NAME, 'model': 'PLAIN', 'attributes': prop}
                    elif Shape[1] == '"I/Wide Flange"':
                        H = list(filter(None, line[3].strip().split('=')))
                        B = list(filter(None, line[4].strip().split('=')))
                        tf = list(filter(None, line[5].strip().split('=')))
                        tw = list(filter(None, line[6].strip().split('=')))
                        NAME = 'LIN' + str(ndim) + 'DWIDEFLANGE'

                        mSet.add(mTag)
                        prop = {'material': mTag, 'h': float(H[1]),'b': float(B[1]), 'tf': float(tf[1]),'tw': float(tw[1]), 'theta': 0.0, 'ip': 10}
                        mesh['Sections'][sTag] = {'name': Shape[1].upper(), 'model': 'PLAIN', 'attributes': prop}
        elif line == 'TABLE:  "AREA SECTION PROPERTIES"':
            sTag = len(mapSAP['Section'])
            while True:
                k += 1
                line = list(shlex.split(lines[k]))
                if line:
                    if line[0] == 'TABLE:':
                        break
                    
                    #PATCH for stupid break line format in SAP
                    if '_' in line:
                        while line[-1] == '_':
                            line.remove('_')
                            k += 1
                            line += list(filter(None, lines[k].strip().split(' ')))

                    sName = list(filter(None, line[0].strip().split('=')))
                    mName = list(filter(None, line[1].strip().split('=')))
                    Shape = list(filter(None, line[4].strip().split('=')))
                    thickness = list(filter(None, line[6].strip().split('=')))
                    thickness = float(thickness[1])
                    
                    sName = sName[1]
                    mName = mName[1]+'_2D'
                    sTag += 1
                    mapSAP['Section'][sName] = sTag
                    mTag = mapSAP['Material'][mName]
                    mSet.add(mTag)

                    prop = {'material': mTag, 'th': thickness}
                    mesh['Sections'][sTag] = {'name': 'LIN3DTHINAREA', 'model': 'PLAIN', 'attributes': prop}
        elif line == 'TABLE:  "CONSTRAINT DEFINITIONS - DIAPHRAGM"':
            dTag = 0
            while True:
                k += 1
                line = list(shlex.split(lines[k]))
                if line:
                    if line[0] == 'TABLE:':
                        break
                    Name = list(filter(None, line[0].strip().split('=')))
                    Axis = list(filter(None, line[2].strip().split('=')))

                    dTag += 1
                    sTag = Name[1]
                    mapSAP['Diaph'][sTag] = dTag
                    mesh['Diaphragms'][dTag] = {'tag': np.nan, 'axis': Axis[1].upper(), 'list': []}
        elif line == 'TABLE:  "JOINT COORDINATES"':
            while True:
                k += 1
                line = list(shlex.split(lines[k]))
                if line:
                    if line[0] == 'TABLE:':
                        break
                    Joint   = list(filter(None, line[0].strip().split('=')))
                    GlobalX = list(filter(None, line[7].strip().split('=')))
                    GlobalY = list(filter(None, line[8].strip().split('=')))
                    GlobalZ = list(filter(None, line[9].strip().split('=')))

                    nTag   = int(Joint[1])
                    coords = [float(GlobalX[1]), float(GlobalY[1]), float(GlobalZ[1])]
                    mesh['Nodes'][nTag] = {'ndof': ndof, 'freedof': np.zeros(ndof, dtype=int), 'totaldof': np.zeros(ndof, dtype=int), 'coords': np.array(coords)}
        elif line == 'TABLE:  "CONNECTIVITY - FRAME"':
            while True:
                k += 1
                line = list(shlex.split(lines[k]))
                if line:
                    if line[0] == 'TABLE:':
                        break
                    Frame   = list(filter(None, line[0].strip().split('=')))
                    JointI  = list(filter(None, line[1].strip().split('=')))
                    JointJ  = list(filter(None, line[2].strip().split('=')))

                    #The associated SeismoVLAB element class
                    NAME = 'LIN' + str(ndim) + 'DFRAME2'
                    
                    #Element tag and corresponding mapping                        
                    eTag += 1
                    sTag = int(Frame[1])
                    mapSAP['Frame'][sTag] = eTag
                    conn = np.array([int(JointI[1]), int(JointJ[1])])

                    attributes = {'rule': 'GAUSS', 'np': 3, 'formulation':'Bernoulli','section': np.nan}
                    mesh['Elements'][eTag] = {'name': NAME, 'conn': conn, 'attributes': attributes}
        elif line == 'TABLE:  "CONNECTIVITY - AREA"':
            while True:
                k += 1
                line = list(shlex.split(lines[k]))
                if line:
                    if line[0] == 'TABLE:':
                        break
                    Area = list(filter(None, line[0].strip().split('=')))
                    NumJoints = list(filter(None, line[1].strip().split('=')))

                    nNodes = int(NumJoints[1])
                    if nNodes == 3:
                        Joint1 = list(filter(None, line[2].strip().split('=')))
                        Joint2 = list(filter(None, line[3].strip().split('=')))
                        Joint3 = list(filter(None, line[4].strip().split('=')))
                        conn   = np.array([int(Joint1[1]), int(Joint2[1]), int(Joint3[1])])
                        NAME = 'LIN3DSHELL3'
                        nQp = 7
                    elif nNodes == 4:
                        Joint1 = list(filter(None, line[2].strip().split('=')))
                        Joint2 = list(filter(None, line[3].strip().split('=')))
                        Joint3 = list(filter(None, line[4].strip().split('=')))
                        Joint4 = list(filter(None, line[5].strip().split('=')))
                        conn   = np.array([int(Joint1[1]), int(Joint2[1]), int(Joint3[1]), int(Joint4[1])])
                        NAME = 'LIN3DSHELL4'
                        nQp = 9

                    #Element tag and corresponding mapping                        
                    eTag += 1
                    sTag = int(Area[1])
                    mapSAP['Area'][sTag] = eTag
                    attributes = {'rule': 'GAUSS', 'np': nQp, 'section': np.nan}
                    mesh['Elements'][eTag] = {'name': NAME, 'conn': conn, 'attributes': attributes}
        elif line == 'TABLE:  "CONNECTIVITY - LINK"':
            while True:
                k += 1
                line = list(shlex.split(lines[k]))
                if line:
                    if line[0] == 'TABLE:':
                        break
                    Link   = list(filter(None, line[0].strip().split('=')))
                    JointI = list(filter(None, line[1].strip().split('=')))
                    JointJ = list(filter(None, line[2].strip().split('=')))

                    #Element tag and corresponding mapping                        
                    eTag += 1
                    sTag = int(Link[1])
                    mapSAP['Link'][sTag] = eTag
                    conn = np.array([int(JointI[1]), int(JointJ[1])])

                    attributes = {}
                    mesh['Elements'][eTag] = {'name': '', 'conn': conn, 'attributes': attributes}
        elif line == 'TABLE:  "JOINT CONSTRAINT ASSIGNMENTS"':
            while True:
                k += 1
                line = list(shlex.split(lines[k]))
                if line:
                    if line[0] == 'TABLE:':
                        break
                    Joint = list(filter(None, line[0].strip().split('=')))
                    Type  = list(filter(None, line[2].strip().split('=')))
                    Constraint = list(filter(None, line[1].strip().split('=')))
                    if Type[1] == 'Diaphragm':
                        nTag = int(Joint[1])
                        sTag = Constraint[1]
                        dTag = mapSAP['Diaph'][sTag]
                        mesh['Diaphragms'][dTag]['list'].append(nTag)

            #Asign a new node number for this diaphragm
            dtag = max(list(mesh['Nodes'].keys()))
            for diaph in mesh['Diaphragms']:
                dtag += 1 
                mesh['Diaphragms'][diaph]['tag'] = dtag
        elif line == 'TABLE:  "JOINT RESTRAINT ASSIGNMENTS"':
            while True:
                k += 1
                line = list(shlex.split(lines[k]))
                if line:
                    if line[0] == 'TABLE:':
                        break
                    Joint = list(filter(None, line[0].strip().split('=')))
                    UX = list(filter(None, line[1].strip().split('=')))
                    UY = list(filter(None, line[2].strip().split('=')))
                    UZ = list(filter(None, line[3].strip().split('=')))
                    RX = list(filter(None, line[4].strip().split('=')))
                    RY = list(filter(None, line[5].strip().split('=')))
                    RZ = list(filter(None, line[6].strip().split('=')))

                    #Assign the restrain to a certain dof
                    nTag = int(Joint[1])
                    if ndim == 2:
                        for j, dof in enumerate([UX[1],UY[1],RZ[1]]):
                            if dof.upper() == 'YES':
                                mesh['Nodes'][nTag]['freedof'][j] = -1
                    elif ndim == 3:
                        for j, dof in enumerate([UX[1],UY[1],UZ[1],RX[1],RY[1],RZ[1]]):
                            if dof.upper() == 'YES':
                                mesh['Nodes'][nTag]['freedof'][j] = -1
        elif line == 'TABLE:  "FRAME SECTION ASSIGNMENTS"':
            while True:
                k += 1
                line = list(shlex.split(lines[k]))
                if line:
                    if line[0] == 'TABLE:':
                        break
                    Frame = list(filter(None, line[0].strip().split('=')))
                    sName = list(filter(None, line[3].strip().split('=')))
                    Frame = int(Frame[1])
                    sName = sName[1]

                    sTag = mapSAP['Section'][sName]
                    eTag = mapSAP['Frame'][Frame]
                    sSet.add(sTag)
                    mesh['Elements'][eTag]['attributes']['section'] = sTag
        elif line == 'TABLE:  "AREA SECTION ASSIGNMENTS"':
            while True:
                k += 1
                line = list(shlex.split(lines[k]))
                if line:
                    if line[0] == 'TABLE:':
                        break
                    Area = list(filter(None, line[0].strip().split('=')))
                    sName = list(filter(None, line[1].strip().split('=')))
                    Area = int(Area[1])
                    sName = sName[1]

                    sTag = mapSAP['Section'][sName]
                    eTag = mapSAP['Area'][Area]
                    sSet.add(sTag)
                    mesh['Elements'][eTag]['attributes']['section'] = sTag         
        elif line == 'TABLE:  "LINK PROPERTY ASSIGNMENTS"':
            while True:
                k += 1
                line = list(shlex.split(lines[k]))
                if line:
                    if line[0] == 'TABLE:':
                        break
                    Link  = list(filter(None, line[0].strip().split('=')))
                    sName = list(filter(None, line[1].strip().split('=')))
                    Link  = int(Link[1])
                    sName = sName[1]

                    #Associates the SeismoVLAB link class
                    if sName == "Linear":
                        NAME = 'LINEAR' + str(ndim) + 'DLINK'

                        #Predefined values that need to be changed/obtained
                        attributes = {'Ks': 1.00, 'dir': 1}
                    elif sName == "Rubber Isolator":
                        NAME = 'HDRBYAMAMOTO' + str(ndim) + 'DLINK'

                        #Predefined values that need to be changed/obtained
                        attributes = {'De': 1.00, 'Di': 0.15, 'Hr': 0.165}
                    elif sName == "Plastic (Wen)":
                        NAME = 'UNXBOUCWEN' + str(ndim) + 'DLINK'

                        #Predefined values that need to be changed/obtained
                        attributes = {'alpha': 1.0, 'mu': 2.0, 'eta': 1.0, 'beta': 0.5, 'gamma': 0.5, 'tol': 1E-6, 'fy': 250.0, 'k0': 250.0, 'a1': 0.1, 'a2': 0.0, 'dir': 1}
                    eTag = mapSAP['Link'][Link]
                    mesh['Elements'][eTag]['name'] = NAME
                    mesh['Elements'][eTag]['attributes'] = attributes  
        else:
            k += 1

    #Cleans the unused materials
    mTags = list(mesh['Materials'].keys())
    for m in mTags:
        if m not in mSet:
            del mesh['Materials'][m]

    #Cleans the unused sections
    sTags = list(mesh['Sections'].keys())
    for s in sTags:
        if s not in sSet:
            del mesh['Sections'][s]

    return mesh