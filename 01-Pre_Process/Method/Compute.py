#!/usr/bin/python3
# -*- coding: Utf-8 -*-

import numpy as np
from Core.Definitions import Entities, Options
from Core.Utilities import *

def GetDRMInformation(x0, xl):
    """
    This function finds the DRM elements provided with a rectangular box. The 
    x0 is the coordinate of center of the box, and xl the half-length side in 
    each direction of the box. These sides (planes) should be inside the DRM 
    elements to be identified.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2021

    Parameters
    ----------
    x0: array
        The coordinate of center of the box
    xl: array
        The side half-length in each direction

    Returns
    -------
    Elements  : list
        List with the DRM element Tags 
    Interior  : list
        List with the DRM internal node Tags 
    Exterior  : list
        List with the DRM external node Tags 
    """
    #Identify the DRM elements
    elemDRM = set()
    for eTag in Entities['Elements']:
        count = 0
        nodes = Entities['Elements'][eTag]['conn']
        for nTag in nodes:
            xn = Entities['Nodes'][nTag]['coords']
            cn = np.abs(xn - x0) <= xl
            if cn.all():
                count += 1
        if 0 < count and count < len(nodes):
            elemDRM.add(eTag)

    #Interior/Exterior DRM node lists
    intDRM = set()
    extDRM = set()
    for eTag in elemDRM:
        nodes = Entities['Elements'][eTag]['conn']
        for nTag in nodes:
            xn = Entities['Nodes'][nTag]['coords']
            cn = np.abs(xn - x0) <= xl
            if cn.all():
                intDRM.add(nTag)
            else:
                extDRM.add(nTag)
    #Provides the DRM information in Mesh
    intDRM = list(intDRM)
    extDRM = list(extDRM)
    elemDRM = list(elemDRM)    

    return intDRM, extDRM, elemDRM

def SurfaceFace(name, surf, conn):
    """
    This function associate the face (surface) number to an element depending
    on the connectivity array specified for that face.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    name : str
        The basic element type name
    surf : list
        The connectivity array specified for this face
    conn : list
        The connectivity array of the element

    Returns
    -------
    int
        The associated face number for the given connectivity
    """
    if  name == 'LINE':
        return 1
    elif name == 'TRIA':
        wTag  = list(set(surf).intersection(conn[[0,1]]))
        if len(wTag) == 2:
            return 1
        wTag = list(set(surf).intersection(conn[[1,2]]))
        if len(wTag) == 2:
            return 2
        wTag = list(set(surf).intersection(conn[[2,0]]))
        if len(wTag) == 2:
            return 3
        wTag = list(set(surf).intersection(conn[[0,1,2]]))
        if len(wTag) == 3:
            return 4
    elif name == 'QUAD':
        wTag  = list(set(surf).intersection(conn[[0,1]]))
        if len(wTag) == 2:
            return 1
        wTag = list(set(surf).intersection(conn[[1,2]]))
        if len(wTag) == 2:
            return 2
        wTag = list(set(surf).intersection(conn[[2,3]]))
        if len(wTag) == 2:
            return 3
        wTag = list(set(surf).intersection(conn[[3,0]]))
        if len(wTag) == 2:
            return 4
        wTag = list(set(surf).intersection(conn[[0,1,2,3]]))
        if len(wTag) == 4:
            return 5
    elif name == 'TETRA':
        wTag  = list(set(surf).intersection(conn[[0,2,1]]))
        if len(wTag) == 3:
            return 1
        wTag = list(set(surf).intersection(conn[[0,1,3]]))
        if len(wTag) == 3:
            return 2
        wTag = list(set(surf).intersection(conn[[1,2,3]]))
        if len(wTag) == 3:
            return 3
        wTag = list(set(surf).intersection(conn[[2,0,3]]))
        if len(wTag) == 3:
            return 4
    elif name == 'HEXA':
        wTag = list(set(surf).intersection(conn[[0,1,2,3]]))
        if len(wTag) == 4:
            return 1
        wTag = list(set(surf).intersection(conn[[0,4,5,1]]))
        if len(wTag) == 4:
            return 2
        wTag = list(set(surf).intersection(conn[[1,2,6,5]]))
        if len(wTag) == 4:
            return 3
        wTag = list(set(surf).intersection(conn[[3,7,6,2]]))
        if len(wTag) == 4:
            return 4
        wTag = list(set(surf).intersection(conn[[0,3,7,4]]))
        if len(wTag) == 4:
            return 5
        wTag = list(set(surf).intersection(conn[[4,5,6,7]]))
        if len(wTag) == 4:
            return 6

def GenerateFiberSection(section):
    """
    This function takes information from patches and layers and transform it
    into a fiber section using the JSON Run-Analysis format.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    section : Dictionary with fiber section information
        patch : Area to generate a number of fibers over a geometric cross-section 
        layer : Line to generate a row of fibers along a geometric-arc 

    Returns
    -------
    attributes : dict
        The section transformed into fibers
    """
    zi = list()
    yi = list()
    Ai = list()
    Tags = list()
    if 'patch' in section:
        Patches = section['patch']
        for pTag in Patches:
            patch = Patches[pTag]
            tag = patch['fiber']
            if patch['name'].upper() == 'RECTANGULAR':
                nfibz = patch['nfibz']
                nfiby = patch['nfiby']
                points = np.matrix(patch['coords'])
                DX = (points[1,0] - points[0,0])/nfibz
                DY = (points[1,1] - points[0,1])/nfiby
                for j in range(nfiby):
                    for i in range(nfibz):
                        coords = points[0,:] + np.array([DX*(0.5 + i), DY*(0.5 + j)])
                        zi.append(coords[0,0])
                        yi.append(coords[0,1])
                        Ai.append(abs(DX*DY))
                        Tags.append(tag)
            elif patch['name'].upper() == 'CIRCULAR':
                nfibr = patch['nfibr']
                nfibt = patch['nfibt']
                center = patch['center']
                points = np.matrix(patch['coords'])
                r1 = points[0,0]
                DR = np.abs((points[1,0] - points[0,0]))/nfibr
                DT = np.abs((points[1,1] - points[0,1]))/nfibt
                for j in range(nfibr):
                    for i in range(nfibt):
                        #Center fiber only for Solid Circular section
                        if j == 0 and r1 == 0 and np.abs((points[1,1] - points[0,1])) == 360:
                            zi.append(center[0])
                            yi.append(center[1])
                            Ai.append(np.pi*DR**2)
                            Tags.append(tag)
                            break
                        coords = points[0,:] + np.array([DR*(0.5 + j), DT*(0.5 + i)])
                        zi.append(center[0] - coords[0,0]*np.cos(np.radians(coords[0,1])))
                        yi.append(center[1] + coords[0,0]*np.sin(np.radians(coords[0,1])))
                        Ai.append(0.5*np.radians(DR)*DT*(2.0*r1 + (2*j+1)*DR))
                        Tags.append(tag)
            '''
            elif patch['name'].upper() == 'QUADRILATERAL':
                DX = 2.0/nfibz
                DZ = 2.0/nfiby
                points = np.matrix(patch['coords'])
                for j in range(nfiby):
                    si = 1.0*j + DY/2.0 - 1.0
                    for k in range(nfibz):
                        ri = 1.0*k + DX/2.0 - 1.0
                        H1 = 0.25*(1.0 - ri)*(1.0 - si)
                        H2 = 0.25*(1.0 + ri)*(1.0 - si)
                        H3 = 0.25*(1.0 + ri)*(1.0 + si)
                        H4 = 0.25*(1.0 - ri)*(1.0 + si)
                        Xi = H1*points[0,0] + H2*points[1,0] + H3*points[2,0] + H4*points[3,0]
                        Yi = H1*points[0,1] + H2*points[1,1] + H3*points[2,1] + H4*points[3,1]
                        zi.append(Xi)
                        yi.append(Yi)
                        Ai.append(area)####
                        Tags.append(tag)
            '''
    if 'layer' in section:
        Layers = section['layer']
        for lTag in Layers:
            layer = Layers[lTag]
            tag = layer['fiber']
            nfib = layer['nfib']
            area = layer['area']
            if layer['name'].upper() == 'LINE':
                points = np.matrix(layer['coords'])
                if nfib == 1:
                    xm = 0.5*(points[0,:] + points[1,:])
                    zi.append(xm[0,0])
                    yi.append(xm[0,1])
                    Ai.append(area)
                    Tags.append(tag)
                elif nfib == 2:
                    for k in range(nfib):
                        zi.append(points[k,0])
                        yi.append(points[k,1])
                        Ai.append(area)
                        Tags.append(tag)
                else:
                    D = (points[1,:] - points[0,:])/(nfib - 1.0)
                    for k in range(nfib):
                        coords = points[0,:] + D*k
                        zi.append(coords[0,0])
                        yi.append(coords[0,1])
                        Ai.append(area)
                        Tags.append(tag)
            elif layer['name'].upper() == 'ARCH':
                angles = layer['angle']
                radius = layer['radius']
                center = layer['center']
                if nfib == 1:
                    xm = 0.5*(angles[0] + angles[1])
                    zi.append(center[0] - radius*np.cos(np.radians(xm)))
                    yi.append(center[1] + radius*np.sin(np.radians(xm)))
                    Ai.append(area)
                    Tags.append(tag)
                elif nfib == 2:
                    for k in range(nfib):
                        zi.append(center[0] - radius*np.cos(np.radians(angles[k])))
                        yi.append(center[1] + radius*np.sin(np.radians(angles[k])))
                        Ai.append(area)
                        Tags.append(tag)
                else:
                    D = np.abs((angles[1] - angles[0]))/nfib if abs(angles[1] - angles[0]) == 360 else np.abs((angles[1] - angles[0]))/(nfib - 1)
                    for k in range(nfib):
                        alpha = angles[0] + D*k
                        zi.append(center[0] - radius*np.cos(np.radians(alpha)))
                        yi.append(center[1] + radius*np.sin(np.radians(alpha)))
                        Ai.append(area)
                        Tags.append(tag)
    attributes = {'fiber': Tags, 'zi': zi, 'yi': yi, 'Ai': Ai}
    if 'h' in section:
        attributes['h'] = section['h']
    if 'b' in section:
        attributes['b'] = section['b']
    if 'kappa2' in section:
        attributes['kappa2'] = section['kappa2']
    if 'kappa3' in section:
        attributes['kappa3'] = section['kappa3']
    if 'ip' in section:
        attributes['ip'] = section['ip']
    if 'theta' in section:
        attributes['theta'] = section['theta']

    return attributes

def RigidLinkConstraints():
    """
    This function generates the constraints associated to the applied rigid
    diaphragm.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2021

    Returns
    -------
    None
        Updates the 'Constraints' in Entities computing rigid link constraints
    """
    #Gets the last defined constraint index
    if Entities['Constraints']:
        ind = list(Entities['Constraints'].keys())
        cTag = min(ind)
    else:
        cTag = -1

    #TODO: Implement rigid link constraints
    for rlTag in Entities['RigidLinks']:
        #Rigid link master node and coordinate.
        mTag = Entities['RigidLinks'][rlTag]['tag']
        Center = Entities['Nodes'][mTag]['coords']

        #Rigid link nodes.
        nTags = Entities['RigidLinks'][rlTag]['list']

        for sTag in nTags:
            Coordinates = Entities['Nodes'][sTag]['coords']
            
            if Options['dimension'] == 3:
                dx = Coordinates[0] - Center[0]
                dy = Coordinates[1] - Center[1]
                dz = Coordinates[2] - Center[2]

                #Constraint in direction 1
                cTag -= 1
                Entities['Constraints'][cTag] = {'name': 'RIGIDLINK', 'stag': sTag, 'sdof': 0, 'mtag': [mTag,mTag,mTag], 'mdof': [0,4,5], 'factor': [1.00, dz, -dy]}
                Entities['Nodes'][sTag]['freedof'][0] = cTag

                #Constraint in direction 2
                cTag -= 1
                Entities['Constraints'][cTag] = {'name': 'RIGIDLINK', 'stag': sTag, 'sdof': 1, 'mtag': [mTag,mTag,mTag], 'mdof': [1,3,5], 'factor': [1.00, -dz, dx]}
                Entities['Nodes'][sTag]['freedof'][1] = cTag

                #Constraint in direction 3
                cTag -= 1
                Entities['Constraints'][cTag] = {'name': 'RIGIDLINK', 'stag': sTag, 'sdof': 2, 'mtag': [mTag,mTag,mTag], 'mdof': [2,3,4], 'factor': [1.00, dy, -dx]}
                Entities['Nodes'][sTag]['freedof'][2] = cTag

                if Entities['RigidLinks'][rlTag]['type'] == 'STRUCTURAL':
                    #Constraint in rotation 1
                    cTag -= 1
                    Entities['Constraints'][cTag] = {'name': 'RIGIDLINK', 'stag': sTag, 'sdof': 3, 'mtag': [mTag], 'mdof': [3], 'factor': [1.00]}
                    Entities['Nodes'][sTag]['freedof'][3] = cTag

                    #Constraint in rotation 2
                    cTag -= 1
                    Entities['Constraints'][cTag] = {'name': 'RIGIDLINK', 'stag': sTag, 'sdof': 4, 'mtag': [mTag], 'mdof': [4], 'factor': [1.00]}
                    Entities['Nodes'][sTag]['freedof'][4] = cTag

                    #Constraint in rotation 3
                    cTag -= 1
                    Entities['Constraints'][cTag] = {'name': 'RIGIDLINK', 'stag': sTag, 'sdof': 5, 'mtag': [mTag], 'mdof': [5], 'factor': [1.00]}
                    Entities['Nodes'][sTag]['freedof'][5] = cTag

            elif Options['dimension'] == 2:
                dx = Coordinates[0] - Center[0]
                dy = Coordinates[1] - Center[1]

                #Constraint in direction 1
                cTag -= 1
                Entities['Constraints'][cTag] = {'name': 'RIGIDLINK', 'stag': sTag, 'sdof': 0, 'mtag': [mTag,mTag], 'mdof': [0,2], 'factor': [1.00, -dy]}
                Entities['Nodes'][sTag]['freedof'][0] = cTag

                #Constraint in direction 2
                cTag -= 1
                Entities['Constraints'][cTag] = {'name': 'RIGIDLINK', 'stag': sTag, 'sdof': 1, 'mtag': [mTag,mTag], 'mdof': [1,2], 'factor': [1.00, dx]}
                Entities['Nodes'][sTag]['freedof'][1] = cTag

                if Entities['RigidLinks'][rlTag]['type'] == 'STRUCTURAL':
                    #Constraint in rotation 3
                    cTag -= 1
                    Entities['Constraints'][cTag] = {'name': 'RIGIDLINK', 'stag': sTag, 'sdof': 2, 'mtag': [mTag], 'mdof': [2], 'factor': [1.00]}
                    Entities['Nodes'][sTag]['freedof'][2] = cTag

def DiaphragmConstraints():
    """
    This function generates the constraints associated to the applied rigid
    diaphragm.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Returns
    -------
    None
        Updates the 'Constraints' in Entities computing diaphragm constraints
    """
    #Gets the last defined constraint index
    if Entities['Constraints']:
        ind = list(Entities['Constraints'].keys())
        cTag = min(ind)
    else:
        cTag = -1

    for dTag in Entities['Diaphragms']:
        #Computes the center of Diaphragm.
        count = 0.0
        nTags = Entities['Diaphragms'][dTag]['list']
        DiaphCenter = np.zeros(Options['dimension'], dtype='float')
        for sTag in nTags:
            DiaphCenter += Entities['Nodes'][sTag]['coords']
            count += 1.00
        DiaphCenter = DiaphCenter/count

        #Computes Combinational Factors for Constraints.
        mTag = Entities['Diaphragms'][dTag]['tag']  
        for sTag in nTags:
            Coordinates = Entities['Nodes'][sTag]['coords']

            if Options['dimension'] == 3:
                dx = Coordinates[0] - DiaphCenter[0]
                dy = Coordinates[1] - DiaphCenter[1]
                dz = Coordinates[2] - DiaphCenter[2]

                if Entities['Diaphragms'][dTag]['axis'] == 'Z':
                    sdof = [0, 1, 5]
                    factors = [-dy, dx]
                elif Entities['Diaphragms'][dTag]['axis'] == 'X':
                    sdof = [1, 2, 3]
                    factors = [-dz, dy]
                elif Entities['Diaphragms'][dTag]['axis'] == 'Y':
                    sdof = [0, 2, 4]
                    factors = [dz, -dx]

                #Constraint in direction 1
                cTag -= 1
                Entities['Constraints'][cTag] = {'name': 'DIAPHRAGM', 'stag': sTag, 'sdof': sdof[0], 'mtag': [mTag,mTag], 'mdof': [0,2], 'factor': [1.00, factors[0]]}
                Entities['Nodes'][sTag]['freedof'][sdof[0]] = cTag

                #Constraint in direction 2
                cTag -= 1
                Entities['Constraints'][cTag] = {'name': 'DIAPHRAGM', 'stag': sTag, 'sdof': sdof[1], 'mtag': [mTag,mTag], 'mdof': [1,2], 'factor': [1.00, factors[1]]}
                Entities['Nodes'][sTag]['freedof'][sdof[1]] = cTag

                #Constraint in rotation 3
                if Entities['Nodes'][sTag]['ndof'] == 6:
                    cTag -= 1
                    Entities['Constraints'][cTag] = {'name': 'DIAPHRAGM', 'stag': sTag, 'sdof': sdof[2], 'mtag': [mTag], 'mdof': [2], 'factor': [1.00]}
                    Entities['Nodes'][sTag]['freedof'][sdof[2]] = cTag

            elif Options['dimension'] == 2:
                dx = Coordinates[0] - DiaphCenter[0]
                dy = Coordinates[1] - DiaphCenter[1]

                #Constraint in direction 1
                cTag -= 1
                Entities['Constraints'][cTag] = {'name': 'DIAPHRAGM', 'stag': sTag, 'sdof': 0, 'mtag': [mTag,mTag], 'mdof': [0,2], 'factor': [1.00, -dy]}
                Entities['Nodes'][sTag]['freedof'][0] = cTag

                #Constraint in direction 2
                cTag -= 1
                Entities['Constraints'][cTag] = {'name': 'DIAPHRAGM', 'stag': sTag, 'sdof': 1, 'mtag': [mTag,mTag], 'mdof': [1,2], 'factor': [1.00, dx]}
                Entities['Nodes'][sTag]['freedof'][1] = cTag

                #Constraint in rotation 3
                if Entities['Nodes'][sTag]['ndof'] == 3:
                    cTag -= 1
                    Entities['Constraints'][cTag] = {'name': 'DIAPHRAGM', 'stag': sTag, 'sdof': 2, 'mtag': [mTag], 'mdof': [2], 'factor': [1.00]}
                    Entities['Nodes'][sTag]['freedof'][2] = cTag

        #Creates the Diaphragm Node
        Entities['Nodes'][mTag] = {'ndof': 3, 'freedof': np.zeros(3, dtype=int), 'totaldof': np.zeros(3, dtype=int), 'coords': DiaphCenter}

def RigidBodyConstraints():
    """
    This function generates the constraints associated to a rigid body. A point 'Po'
    must be specified form which the constraints will be generated. 'Po' is the rigid 
    body master node\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Returns
    -------
    None
        Updates the 'Constraints' in Entities computing rigid body constraints
    """
    #Gets the last defined constraint index
    if Entities['Constraints']:
        ind = list(Entities['Constraints'].keys())
        cTag = min(ind)
    else:
        cTag = -1

    for rbTag in Entities['RigidBodies']:
        #Computes the center of Diaphragm.
        nTags = Entities['RigidBodies'][rbTag]['list']
        RotationCenter = Entities['RigidBodies'][rbTag]['center']

        #Computes Combinational Factors for Constraints.
        mTag = Entities['RigidBodies'][rbTag]['tag']
        for sTag in nTags:
            Coordinates = Entities['Nodes'][sTag]['coords']
            
            if Options['dimension'] == 3:
                dx = Coordinates[0] - RotationCenter[0]
                dy = Coordinates[1] - RotationCenter[1]
                dz = Coordinates[2] - RotationCenter[2]

                #Constraint in direction 1
                cTag -= 1
                Entities['Constraints'][cTag] = {'name': 'RIGIDBODY', 'stag': sTag, 'sdof': 0, 'mtag': [mTag,mTag,mTag], 'mdof': [0,4,5], 'factor': [1.00, dz, -dy]}
                Entities['Nodes'][sTag]['freedof'][0] = cTag

                #Constraint in direction 2
                cTag -= 1
                Entities['Constraints'][cTag] = {'name': 'RIGIDBODY', 'stag': sTag, 'sdof': 1, 'mtag': [mTag,mTag,mTag], 'mdof': [1,3,5], 'factor': [1.00, -dz, dx]}
                Entities['Nodes'][sTag]['freedof'][1] = cTag

                #Constraint in direction 3
                cTag -= 1
                Entities['Constraints'][cTag] = {'name': 'RIGIDBODY', 'stag': sTag, 'sdof': 2, 'mtag': [mTag,mTag,mTag], 'mdof': [2,3,4], 'factor': [1.00, dy, -dx]}
                Entities['Nodes'][sTag]['freedof'][2] = cTag

                if Entities['Nodes'][sTag]['ndof'] == 6:
                    #Constraint in rotation 1
                    cTag -= 1
                    Entities['Constraints'][cTag] = {'name': 'RIGIDBODY', 'stag': sTag, 'sdof': 3, 'mtag': [mTag], 'mdof': [3], 'factor': [1.00]}
                    Entities['Nodes'][sTag]['freedof'][3] = cTag

                    #Constraint in rotation 2
                    cTag -= 1
                    Entities['Constraints'][cTag] = {'name': 'RIGIDBODY', 'stag': sTag, 'sdof': 4, 'mtag': [mTag], 'mdof': [4], 'factor': [1.00]}
                    Entities['Nodes'][sTag]['freedof'][4] = cTag

                    #Constraint in rotation 3
                    cTag -= 1
                    Entities['Constraints'][cTag] = {'name': 'RIGIDBODY', 'stag': sTag, 'sdof': 5, 'mtag': [mTag], 'mdof': [5], 'factor': [1.00]}
                    Entities['Nodes'][sTag]['freedof'][5] = cTag

            elif Options['dimension'] == 2:
                dx = Coordinates[0] - RotationCenter[0]
                dy = Coordinates[1] - RotationCenter[1]

                #Constraint in direction 1
                cTag -= 1
                Entities['Constraints'][cTag] = {'name': 'RIGIDBODY', 'stag': sTag, 'sdof': 0, 'mtag': [mTag,mTag], 'mdof': [0,2], 'factor': [1.00, -dy]}
                Entities['Nodes'][sTag]['freedof'][0] = cTag

                #Constraint in direction 2
                cTag -= 1
                Entities['Constraints'][cTag] = {'name': 'RIGIDBODY', 'stag': sTag, 'sdof': 1, 'mtag': [mTag,mTag], 'mdof': [1,2], 'factor': [1.00, dx]}
                Entities['Nodes'][sTag]['freedof'][1] = cTag

                if Entities['Nodes'][sTag]['ndof'] == 3:
                    #Constraint in rotation 3
                    cTag -= 1
                    Entities['Constraints'][cTag] = {'name': 'RIGIDBODY', 'stag': sTag, 'sdof': 2, 'mtag': [mTag], 'mdof': [2], 'factor': [1.00]}
                    Entities['Nodes'][sTag]['freedof'][2] = cTag

        #Creates the Rigid Body Node
        ndof = Entities['RigidBodies'][rbTag]['ndof']
        Entities['Nodes'][mTag] = {'ndof': ndof, 'freedof': np.zeros(ndof, dtype=int), 'totaldof': np.zeros(ndof, dtype=int), 'coords': np.array(RotationCenter)}

def ApplyConstraints():
    """
    This function generates the constraints associated to a rigid body, diaphragm,
    and rigid body. These constraints are appended to Entities['Constraints']\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2021

    Returns
    -------
    None
        Updates the 'Constraints' in Entities
    """
    #Apply Rigid Link Constraints
    RigidLinkConstraints()

    #Apply Diaphragm Constraints
    DiaphragmConstraints()

    #Apply Rigid Body Constraints
    RigidBodyConstraints()
