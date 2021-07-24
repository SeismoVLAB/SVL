#!/usr/bin/env python3
# -*- coding: Utf-8 -*-

import copy
import numpy as np
from Core.Utilities import debugInfo
from Core.Definitions import Options

def reLabel(mesh):
    """
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020
    """
    #Re-Numbers Nodes
    tag = 0
    nMap  = dict()
    nodes = dict()
    for nTag in mesh['Nodes']:
        tag += 1
        nMap[nTag] = tag
        nodes[tag] = mesh['Nodes'][nTag]          

    #Re-Numbers Elements
    tag = 0
    elems = dict()
    for eTag in mesh['Elements']:
        tag += 1
        elems[tag] = mesh['Elements'][eTag]
        elems[tag]['conn'] = [nMap[x] for x in elems[tag]['conn']]

    #Re-Numbers Interface
    if 'Boundary' in mesh:
        bcs = dict()
        for name in mesh['Boundary']:
            nTags = list()
            for nTag in mesh['Boundary'][name]:
                if nTag in mesh['Nodes']:
                    nTags.append(nMap[nTag])
            bcs[name] = nTags

    #Re-Numbers Constraints
    if 'Constraints' in mesh:
        for cTag in mesh['Constraints']:
            stag = mesh['Constraints'][cTag]['stag']
            mesh['Constraints'][cTag]['stag'] = nMap[stag]
            mtag = mesh['Constraints'][cTag]['mtag']
            for k, Tag in enumerate(mtag):
                mtag[k] = nMap[Tag]

    #New Labels
    mesh['Nodes'] = nodes
    mesh['Elements'] = elems
    mesh['Boundary'] = bcs

    return nMap

def makeDomainVolume(options={}):
    """
    Creates a volume domain (in 3D) with specified parameters is option.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    options : dict
        Specify some parameters for constructing the mesh
        'ne'     : (list) The number of elements per dimension, i.e., ne = [nx,ny,nz]
        'ndof'   : (int) number of degree of freedom per node
        'class'  : (str) SeismoVLAB element class name for the generated elements
        'P0'     : (list) The origin point of the volume element
        'P1'     : (list) The corner point (x-side) of the volume element
        'P2'     : (list) The corner point (y-side) of the volume element
        'P3'     : (list) The corner point (z-side) of the volume element
        'elems'  : (str) the type of element during meshing,i.e., TETRA4, TETRA10, HEXA8, HEXA20
        'attributes' : (dict) 
            'rule': (str) Integration rule to be used rule=Gauss, Lobatto
            'np'  : (int) Number of integration point to be used
            'material': (int) the material tag
            'other': ...

    Returns
    -------
    Mesh : dict
        A dictionary that contains nodes and element with the same structure as Entities
    """
    #Defines an empty Mesh dictionary
    Mesh = {'Nodes': {}, 'Elements': {}}

    #Check if attributes is provided
    if not options:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d makeDomainVolume(attributes=?) must be specified.' %(info.filename,info.lineno))
        return Mesh

    if 'attributes' not in options:
        if options['elems'].upper() == 'TETRA4':
            nQp = 4
        elif options['elems'].upper() == 'TETRA10':
            nQp = 7
        elif options['elems'].upper() == 'HEXA8':
            nQp = 8
        elif options['elems'].upper() == 'HEXA20':
            nQp = 27
        attributes = {'rule': 'GAUSS', 'np': nQp, 'material': np.nan}
    else:
        if 'rule' in options['attributes']:
            options['attributes']['rule'] = options['attributes']['rule'].upper()
        attributes = options['attributes']

    #Unpack attribute provided by the user
    ndof = options['ndof']
    name = options['class']

    nx, ny, nz = options['ne'][0], options['ne'][1], options['ne'][2]
    P0 = np.array(options['P0'])
    P1 = np.array(options['P1'])
    P2 = np.array(options['P2'])
    P3 = np.array(options['P3'])

    #Creates the grid increment
    DX = (P1 - P0)/nx
    DY = (P2 - P0)/ny
    DZ = (P3 - P0)/nz

    #Defines the primary grid
    tag = 0
    for k in range(nz+1):
        for j  in range(ny+1):
            for i in range(nx+1):
                tag += 1
                coords = P0 + i*DX + j*DY + k*DZ
                Mesh['Nodes'][tag] = {'ndof': ndof, 'freedof': np.zeros(ndof, dtype=int), 'totaldof': np.zeros(ndof, dtype=int), 'coords': coords}

    if options['elems'].upper() == 'TETRA4':
        #Defines the Elements in Mesh
        tag = 0
        for k  in range(nz):
            for j  in range(ny):
                for i in range(nx):
                    n1 = (nx+1)*(ny+1)*k + j*(nx+1) + i + 1
                    n2 = (nx+1)*(ny+1)*k + j*(nx+1) + i + 2
                    n3 = (nx+1)*(ny+1)*k + (j+1)*(nx+1) + i + 2
                    n4 = (nx+1)*(ny+1)*k + (j+1)*(nx+1) + i + 1
                    n5 = (nx+1)*(ny+1)*(k+1) + j*(nx+1) + i + 1
                    n6 = (nx+1)*(ny+1)*(k+1) + j*(nx+1) + i + 2
                    n7 = (nx+1)*(ny+1)*(k+1) + (j+1)*(nx+1) + i + 2
                    n8 = (nx+1)*(ny+1)*(k+1) + (j+1)*(nx+1) + i + 1

                    tag += 1; Mesh['Elements'][tag] = {'name': name.upper(), 'conn': [n1,n2,n3,n6], 'attributes': attributes}
                    tag += 1; Mesh['Elements'][tag] = {'name': name.upper(), 'conn': [n5,n6,n8,n1], 'attributes': attributes}
                    tag += 1; Mesh['Elements'][tag] = {'name': name.upper(), 'conn': [n6,n7,n8,n3], 'attributes': attributes}
                    tag += 1; Mesh['Elements'][tag] = {'name': name.upper(), 'conn': [n1,n3,n8,n6], 'attributes': attributes}
                    tag += 1; Mesh['Elements'][tag] = {'name': name.upper(), 'conn': [n1,n3,n8,n4], 'attributes': attributes}
    elif options['elems'].upper() == 'TETRA10':
        #Creates the secondary 3D grid
        #TODO: Implement this feature
        pass
    elif options['elems'].upper() == 'HEXA8':   
        #Defines the Elements in Mesh
        tag = 0
        for k  in range(nz):
            for j  in range(ny):
                for i in range(nx):
                    tag += 1
                    n1 = (nx+1)*(ny+1)*k + j*(nx+1) + i + 1
                    n2 = (nx+1)*(ny+1)*k + j*(nx+1) + i + 2
                    n3 = (nx+1)*(ny+1)*k + (j+1)*(nx+1) + i + 2
                    n4 = (nx+1)*(ny+1)*k + (j+1)*(nx+1) + i + 1
                    n5 = (nx+1)*(ny+1)*(k+1) + j*(nx+1) + i + 1
                    n6 = (nx+1)*(ny+1)*(k+1) + j*(nx+1) + i + 2
                    n7 = (nx+1)*(ny+1)*(k+1) + (j+1)*(nx+1) + i + 2
                    n8 = (nx+1)*(ny+1)*(k+1) + (j+1)*(nx+1) + i + 1
                
                    Mesh['Elements'][tag] = {'name': name.upper(), 'conn': [n1, n2, n3, n4, n5, n6, n7, n8], 'attributes': attributes}
    elif options['elems'].upper() == 'HEXA20':
        #Creates the secondary 3D grid
        for k in range(nz+1):
            for j in range(ny+1):
                for i in range(nx):
                    mtag = tag + (nx*(ny+1) + ny*(nx+1))*k + (2*nx + 1)*j + i + 1
                    coords = P0 + 0.5*(2*i+1)*DX + j*DY + k*DZ
                    Mesh['Nodes'][mtag] = {'ndof': ndof, 'freedof': np.zeros(ndof, dtype=int), 'totaldof': np.zeros(ndof, dtype=int), 'coords': coords}
        for k in range(nz+1):
            for j in range(ny):
                for i in range(nx+1):
                    mtag = tag + (nx*(ny+1) + ny*(nx+1))*k + nx + (2*nx + 1)*j + i + 1
                    coords = P0 + i*DX + 0.5*(2*j+1)*DY + k*DZ
                    Mesh['Nodes'][mtag] = {'ndof': ndof, 'freedof': np.zeros(ndof, dtype=int), 'totaldof': np.zeros(ndof, dtype=int), 'coords': coords}
        mtag = tag + (nx*(ny+1) + ny*(nx+1))*(nz+1)
        for k in range(nz):
            for j  in range(ny+1):
                for i in range(nx+1):
                    mtag += 1
                    coords = P0 + i*DX + j*DY + 0.5*(2*k+1)*DZ
                    Mesh['Nodes'][mtag] = {'ndof': ndof, 'freedof': np.zeros(ndof, dtype=int), 'totaldof': np.zeros(ndof, dtype=int), 'coords': coords}

        #Defines the Elements in Mesh
        tag = 0
        for k  in range(nz):
            for j  in range(ny):
                for i in range(nx):
                    tag += 1
                    n1 = (nx+1)*(ny+1)*k + j*(nx+1) + i + 1
                    n2 = (nx+1)*(ny+1)*k + j*(nx+1) + i + 2
                    n3 = (nx+1)*(ny+1)*k + (j+1)*(nx+1) + i + 2
                    n4 = (nx+1)*(ny+1)*k + (j+1)*(nx+1) + i + 1
                    n5 = (nx+1)*(ny+1)*(k+1) + j*(nx+1) + i + 1
                    n6 = (nx+1)*(ny+1)*(k+1) + j*(nx+1) + i + 2
                    n7 = (nx+1)*(ny+1)*(k+1) + (j+1)*(nx+1) + i + 2
                    n8 = (nx+1)*(ny+1)*(k+1) + (j+1)*(nx+1) + i + 1
                    n9 = (nx+1)*(ny+1)*(nz+1) + (nx*(ny+1) + ny*(nx+1))*k + (2*nx + 1)*j + i + 1
                    n10 = (nx+1)*(ny+1)*(nz+1) + (nx*(ny+1) + ny*(nx+1))*k + nx + (2*nx + 1)*j + i + 2
                    n11 = (nx+1)*(ny+1)*(nz+1) + (nx*(ny+1) + ny*(nx+1))*k + (2*nx + 1)*(j+1) + i + 1
                    n12 = (nx+1)*(ny+1)*(nz+1) + (nx*(ny+1) + ny*(nx+1))*k + nx + (2*nx + 1)*j + i + 1
                    n13 = (nx+1)*(ny+1)*(nz+1) + (nx*(ny+1) + ny*(nx+1))*(k+1) + (2*nx + 1)*j + i + 1
                    n14 = (nx+1)*(ny+1)*(nz+1) + (nx*(ny+1) + ny*(nx+1))*(k+1) + nx + (2*nx + 1)*j + i + 2
                    n15 = (nx+1)*(ny+1)*(nz+1) + (nx*(ny+1) + ny*(nx+1))*(k+1) + (2*nx + 1)*(j+1) + i + 1
                    n16 = (nx+1)*(ny+1)*(nz+1) + (nx*(ny+1) + ny*(nx+1))*(k+1) + nx + (2*nx + 1)*j + i + 1
                    n17 = (nx+1)*(ny+1)*(nz+1) + (nx*(ny+1) + ny*(nx+1))*(nz+1) + (nx+1)*(ny+1)*k + j*(nx+1) + i + 1
                    n18 = (nx+1)*(ny+1)*(nz+1) + (nx*(ny+1) + ny*(nx+1))*(nz+1) + (nx+1)*(ny+1)*k + j*(nx+1) + i + 2
                    n19 = (nx+1)*(ny+1)*(nz+1) + (nx*(ny+1) + ny*(nx+1))*(nz+1) + (nx+1)*(ny+1)*k + (j+1)*(nx+1) + i + 2
                    n20 = (nx+1)*(ny+1)*(nz+1) + (nx*(ny+1) + ny*(nx+1))*(nz+1) + (nx+1)*(ny+1)*k + (j+1)*(nx+1) + i + 1
                
                    Mesh['Elements'][tag] = {'name': name.upper(), 'conn': [n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14, n15, n16, n17, n18, n19, n20], 'attributes': attributes}

    #Finds the Boundary Nodes list
    if Options['dimension'] == 3:
        TOLx = DX[0]/1000.0
        TOLy = DY[1]/1000.0
        TOLz = DZ[2]/1000.0
        Find3DBoundaries(Mesh, top=P3[2]-TOLz, bottom=P0[2]+TOLz, left=P0[0]+TOLx, right=P1[0]-TOLx, front=P2[1]-TOLy, back=P0[1]+TOLy)

    return Mesh

def makeDomainArea(options={}):
    """
    Creates an area domain (in 2D or 3D) with specified parameters is option.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    options : dict
        Specify some parameters for constructing the mesh
        'ne'     : (list) The number of elements per dimension, i.e., ne = [nx,ny]
        'ndof'   : (int) number of degree of freedom per node
        'class'  : (str) SeismoVLAB element class name for the generated elements
        'P0'     : (list) The origin point of the area element
        'P1'     : (list) The corner point (x-side) of the area element
        'P2'     : (list) The corner point (y-side) of the area element
        'elems'  : (str) the type of element during meshing,i.e., TRIA3, TRIA6, QUAD4, QUAD8
        'attributes' : (dict) 
            'rule': (str) Integration rule to be used rule=Gauss, Lobatto
            'np'  : (int) Number of integration point to be used
            'material': (int) the material tag
            'other': ...

    Returns
    -------
    Mesh : dict
        A dictionary that contains nodes and element with the same structure as Entities
    """
    #Defines an empty Mesh dictionary
    Mesh = {'Nodes': {}, 'Elements': {}}

    #Check if attributes is provided
    if not options:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d makeDomainArea(attributes=?) must be specified.' %(info.filename,info.lineno))
        return Mesh

    if 'attributes' not in options:
        if options['elems'].upper() == 'TRIA3':
            nQp = 3
        elif options['elems'].upper() == 'TRIA6':
            nQp = 7
        elif options['elems'].upper() == 'QUAD4':
            nQp = 4
        elif options['elems'].upper() == 'QUAD8':
            nQp = 9
        attributes = {'rule': 'Gauss', 'np': nQp, 'material': np.nan}
    else:
        attributes = options['attributes']

    #Unpack attribute provided by the user
    ndof = options['ndof']
    name = options['class']

    #Unpack attribute provided by the user
    nx, ny = options['ne'][0], options['ne'][1]
    P0 = np.array(options['P0'])
    P1 = np.array(options['P1'])
    P2 = np.array(options['P2'])

    #Creates the grid increment
    DX = (P1 - P0)/nx
    DY = (P2 - P0)/ny

    #Defines the primary grid
    tag = 0
    for j in range(ny+1):
        for i in range(nx+1):
            tag += 1
            coords = P0 + DX*i + DY*j
            Mesh['Nodes'][tag] = {'ndof': ndof, 'freedof': np.zeros(ndof, dtype=int), 'totaldof': np.zeros(ndof, dtype=int), 'coords': coords}

    if options['elems'].upper() == 'TRIA3':   
        #Defines the Elements in Mesh
        tag = 0
        for j  in range(ny):
            for i in range(nx):
                n1 = (nx+1)*j+i+1
                n2 = (nx+1)*j+i+2
                n3 = (nx+1)*(j+1)+i+2
                n4 = (nx+1)*(j+1)+i+1

                tag += 1; Mesh['Elements'][tag] = {'name': name.upper(), 'conn': [n1, n2, n4], 'attributes': attributes}
                tag += 1; Mesh['Elements'][tag] = {'name': name.upper(), 'conn': [n2, n3, n4], 'attributes': attributes}
    elif options['elems'].upper() == 'TRIA6':
        #Creates the secondary 2D grid
        for j in range(ny+1):
            for i in range(nx):
                mtag = tag + (2*nx + 1)*j + i + 1
                coords = P0 + 0.5*(2*i+1)*DX + DY*j
                Mesh['Nodes'][mtag] = {'ndof': ndof, 'freedof': np.zeros(ndof, dtype=int), 'totaldof': np.zeros(ndof, dtype=int), 'coords': coords}
        for j in range(ny):
            for i in range(nx+1):
                mtag = tag + nx + (2*nx + 1)*j + i + 1
                coords = P0 + i*DX + 0.5*(2*j+1)*DY
                Mesh['Nodes'][mtag] = {'ndof': ndof, 'freedof': np.zeros(ndof, dtype=int), 'totaldof': np.zeros(ndof, dtype=int), 'coords': coords}
        for j in range(ny):
            for i in range(nx):
                mtag = tag + nx*(ny+1) + (nx+1)*ny + nx*j + i + 1
                coords = P0 + 0.5*(2*i+1)*DX + 0.5*(2*j+1)*DY
                Mesh['Nodes'][mtag] = {'ndof': ndof, 'freedof': np.zeros(ndof, dtype=int), 'totaldof': np.zeros(ndof, dtype=int), 'coords': coords}

        #Defines the Elements in Mesh
        ntags = tag
        tag   = 0
        for j  in range(ny):
            for i in range(nx):
                n1 = (nx+1)*j+i+1
                n2 = (nx+1)*j+i+2
                n3 = (nx+1)*(j+1)+i+2
                n4 = (nx+1)*(j+1)+i+1
                n5 = ntags + (2*nx + 1)*j + i + 1
                n6 = ntags + nx + (2*nx + 1)*j + i + 2
                n7 = ntags + (2*nx + 1)*(j+1) + i + 1
                n8 = ntags + nx + (2*nx + 1)*j + i + 1
                n9 = ntags + nx*(ny+1) + (nx+1)*ny + nx*j + i + 1

                tag += 1; Mesh['Elements'][tag] = {'name': name.upper(), 'conn': [n1, n2, n4, n5, n9, n8], 'attributes': attributes}
                tag += 1; Mesh['Elements'][tag] = {'name': name.upper(), 'conn': [n2, n3, n4, n6, n7, n9], 'attributes': attributes}
    elif options['elems'].upper() == 'QUAD4':  
        #Defines the Elements in Mesh
        tag = 0
        for j  in range(ny):
            for i in range(nx):
                tag += 1
                conn = [(nx+1)*j+i+1, (nx+1)*j+i+2, (nx+1)*(j+1)+i+2, (nx+1)*(j+1)+i+1]
                Mesh['Elements'][tag] = {'name': name.upper(), 'conn': conn, 'attributes': attributes}
    elif options['elems'].upper() == 'QUAD8':
        #Creates the secondary 2D grid
        for j in range(ny+1):
            for i in range(nx):
                mtag = tag + (2*nx + 1)*j + i + 1
                coords = P0 + 0.5*(2*i+1)*DX + DY*j
                Mesh['Nodes'][mtag] = {'ndof': ndof, 'freedof': np.zeros(ndof, dtype=int), 'totaldof': np.zeros(ndof, dtype=int), 'coords': coords}
        for j in range(ny):
            for i in range(nx+1):
                mtag = tag + nx + (2*nx + 1)*j + i + 1
                coords = P0 + i*DX + 0.5*(2*j+1)*DY
                Mesh['Nodes'][mtag] = {'ndof': ndof, 'freedof': np.zeros(ndof, dtype=int), 'totaldof': np.zeros(ndof, dtype=int), 'coords': coords}

        #Defines the Elements in Mesh
        ntags = tag
        tag   = 0
        for j  in range(ny):
            for i in range(nx):
                n1 = (nx+1)*j+i+1
                n2 = (nx+1)*j+i+2
                n3 = (nx+1)*(j+1)+i+2
                n4 = (nx+1)*(j+1)+i+1

                n5 = ntags + (2*nx + 1)*j + i + 1
                n6 = ntags + nx + (2*nx + 1)*j + i + 2
                n7 = ntags + (2*nx + 1)*(j+1) + i + 1
                n8 = ntags + nx + (2*nx + 1)*j + i + 1

                tag += 1
                Mesh['Elements'][tag] = {'name': name.upper(), 'conn': [n1, n2, n3, n4, n5, n6, n7, n8], 'attributes': attributes}
                
    #Finds the Boundary Nodes list
    if Options['dimension'] == 2:
        TOLx = DX[0]/1000.0
        TOLz = DY[1]/1000.0
        Find2DBoundaries(Mesh, top=P2[1]-TOLz, bottom=P0[1]+TOLz, left=P0[0]+TOLx, right=P1[0]-TOLx)

    return Mesh

def makeDomainLine(options={}):
    """
    Creates a line domain (in 1D, 2D, or 3D) with specified parameters is option.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    options : dict
        Specify some parameters for constructing the mesh
        'ne'     : (list) The number of elements per dimension, i.e., ne = [nx]
        'ndof'   : (int) number of degree of freedom per node
        'class'  : (str) SeismoVLAB element class name for the generated elements
        'P0'     : (list) The initial point of the line element, coordinate must be consisten with dimension
        'P1'     : (list) The end point of the line element, coordinate must be consisten with dimension
        'elems'  : (str) the type of element during meshing,i.e.,  LINE2, LINE3
        'attributes' : (dict) 
            'rule': (str) Integration rule to be used rule=Gauss, Lobatto
            'np'  : (int) Number of integration point to be used
            'material': (int) the material tag
            'other': ...

    Returns
    -------
    Mesh : dict
        A dictionary that contains nodes and element with the same structure as Entities
    """
    #Defines an empty Mesh dictionary
    Mesh = {'Nodes': {}, 'Elements': {}}

    #Check if attributes is provided
    if not options:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d makeDomainLine(attributes=?) must be specified.' %(info.filename,info.lineno))
        return Mesh

    if 'attributes' not in options:
        if options['elems'].upper() == 'LINE2':
            nQp = 3
        elif options['elems'].upper() == 'LINE3':
            nQp = 5
        attributes = {'rule': 'Gauss', 'np': nQp, 'material': np.nan}
    else:
        attributes = options['attributes']

    #Unpack attribute provided by the user
    ndof = options['ndof']
    name = options['class']

    #Unpack attribute provided by the user
    nx = options['ne']
    P0 = np.array(options['P0'])
    P1 = np.array(options['P1'])

    #Creates the grid increment
    D = (P1 - P0)/nx

    #Defines the primary grid
    tag = 0
    for i in range(nx+1):
        tag += 1
        coords = P0 + D*i
        Mesh['Nodes'][tag] = {'ndof': ndof, 'freedof': np.zeros(ndof, dtype=int), 'totaldof': np.zeros(ndof, dtype=int), 'coords': coords}

    if options['elems'].upper() == 'LINE2':   
        #Defines the Elements in Mesh
        tag = 0
        for i in range(nx):
            tag += 1 
            Mesh['Elements'][tag] = {'name': name.upper(), 'conn': [i+1, i+2], 'attributes': attributes}
    elif options['elems'].upper() == 'LINE3':
        #Creates the secondary grid
        for i in range(nx):
            tag += 1
            coords = P0 + 0.5*(2*i + 1)*D
            Mesh['Nodes'][tag] = {'ndof': ndof, 'freedof': np.zeros(ndof, dtype=int), 'totaldof': np.zeros(ndof, dtype=int), 'coords': coords}

        #Defines the Elements in Mesh
        tag = 0
        for i in range(nx):
            tag += 1 
            Mesh['Elements'][tag] = {'name': name.upper(), 'conn': [i+1, i+2, nx+2+i], 'attributes': attributes}
    return Mesh

def removeDomain(mesh={}, attributes={}):
    """
    This function removes a domain specifying the geometry\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    mesh : dict
        The mesh that contains 'Nodes' and 'Elements' dictionaries from which
        the deletion will be performed. 
    attributes : dict
        The domain to be removed, for example
            'sides' (list) The length of the sides in x,y, and z
            'center' (list) The center of the rectangular domain

    Returns
    -------
    None
    """
    x0 = np.array(attributes['center'])
    xl = 0.5*np.array(attributes['sides'])

    #Identify the elements to be removed
    rmElems = set()
    cpElems = set()
    for eTag in mesh['Elements']:
        count = 0
        nodes = mesh['Elements'][eTag]['conn']
        nNodes = len(nodes)
        for nTag in nodes:
            xn = mesh['Nodes'][nTag]['coords']
            cn = np.abs(xn - x0) <= xl
            if cn.all():
                count += 1
        if 0 < count and count < nNodes:
            cpElems.add(eTag)
        elif count == nNodes:
            rmElems.add(eTag)

    #List of Nodes to be removed
    removeNodes = set()
    for eTag in rmElems:
        nodes = mesh['Elements'][eTag]['conn']
        for nTag in nodes:
            removeNodes.add(nTag)

    #Identify Nodes of the cut section 
    keepNodes = set()
    for eTag in cpElems:
        nodes = mesh['Elements'][eTag]['conn']
        for nTag in nodes:
            xn = mesh['Nodes'][nTag]['coords']
            cn = np.abs(xn - x0) <= xl
            if cn.all():
                keepNodes.add(nTag)

    rmNodes = removeNodes.difference(keepNodes)

    #Remove the elements and associated nodes
    for nTag in rmNodes:
        del mesh['Nodes'][nTag]
    for eTag in rmElems:
        del mesh['Elements'][eTag]

    #Re-Numbers Nodes and Elements and Interface
    nMap = reLabel(mesh)
    cpNodes = list()
    for k in keepNodes:
        cpNodes.append(nMap[k])

    mesh['Boundary']['Interface'] = cpNodes

    return mesh

def swap(mesh1, mesh2):
    return mesh2, mesh1

def mergeDomain(mesh1={}, mesh2={}, TOL=1E-6):
    """
    This function merge two domain into a third mesh\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    mesh1 : dict
        The mesh that merge will be performed. This is it will be and 
        augmented mesh with information from mesh2 
    mesh2 : dict
        The mesh from which data will be copied. 
    TOL   : double
        The tolerance for which 2 nodes will be consider as close.

    Returns
    -------
    mesh1 : dict
        The mesh with the appended mesh domain
    """
    if 'Interface' in mesh1['Boundary']:
        mesh1, mesh2 = mesh2, mesh1

    nNodes = len(mesh1['Nodes'])
    nElems = len(mesh1['Elements'])

    #Nodes are updated with information in mesh2
    nMap = dict()
    for nTag in mesh2['Nodes']:
        tag = nTag + nNodes
        nMap[nTag] = tag
        mesh1['Nodes'][tag] = mesh2['Nodes'][nTag]

    #Elements are updated with information in mesh2
    for eTag in mesh2['Elements']:
        tag = eTag + nElems
        mesh1['Elements'][tag] = mesh2['Elements'][eTag]
        mesh1['Elements'][tag]['conn'] = [nMap[x] for x in mesh1['Elements'][tag]['conn']]

    #The mesh has a common interface to apply constraints 
    if 'Interface' in mesh2['Boundary']:
        #Gets the Node Tags from main Mesh
        nTags = set()
        for names in mesh1['Boundary']:
            if names in ['bottom','left','right','back','front']:
                nTags.update(mesh1['Boundary'][names])

        iTags = mesh2['Boundary']['Interface']
        iTags = np.array(list(iTags))
        nTags = np.array(list(nTags))
        coordM1 = np.zeros((len(nTags), Options['dimension']))
        coordM2 = np.zeros((len(iTags), Options['dimension']))
        for k, nTag in enumerate(nTags):
            coordM1[k,:] = mesh1['Nodes'][nTag]['coords']
        for k, iTag in enumerate(iTags):
            coordM2[k,:] = mesh2['Nodes'][iTag]['coords']
        
        Interface = dict()
        for k, iTag in enumerate(iTags):
            delta = np.subtract(coordM1, coordM2[k,:])
            cond = np.linalg.norm(delta, axis=1) < TOL
            Interface[nTags[cond][0]] = nMap[iTag]

        #Transfers boundaries from mesh2 to mesh 1
        for names in mesh2['Boundary']:
            if names in ['top','bottom','left','right','back','front']:
                nodes = list()
                for nTag in mesh2['Boundary'][names]:
                    nodes.append(nMap[nTag])
                mesh1['Boundary'][names] = nodes
        
        #Create the constraints to tie both meshes together
        if 'Constraints' not in mesh1:
            mesh1['Constraints'] = dict()
            tag = -1
        else:
            cTags = mesh1['Constraints'].keys()
            tag = min(cTags)

        for mtag in Interface:
            stag = Interface[mtag]
            for k in range(Options['dimension']):
                tag += -1
                mesh1['Constraints'][tag] = {'name': 'EQUAL', 'stag': stag, 'sdof': k, 'mtag': [mtag], 'mdof': [k], 'factor': [1.00]}
                mesh1['Nodes'][stag]['freedof'][k] = tag
    return mesh1

def Find2DBoundaries(mesh, top=float('inf'), bottom=float('inf'), left=float('inf'), right=float('inf')):
    """
    This function is a brute force method to find the nodes that are on the boundary
    assuming a rectangular domain. Note that the inputs top, bottom, left, right must
    specify a value a little bit inside the domain. \n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    mesh : dict
        'Nodes'  : (dict) The Node information using SVL format

    Returns
    -------
    mesh['Boundary']: dict
        'Top' : set
            The list of Nodes on the top boundary, i.e., +Z
        'Bottom' : set
            The list of Nodes on the bottom boundary, i.e., -Z
        'Left' : set
            The list of Nodes on the left boundary, i.e., +X
        'Right' : set
            The list of Nodes on the right boundary, i.e., -X
    """
    #Nodes that belong to boundary 
    Top, Bottom = set(), set()
    Left, Right = set(), set()

    for nTag in mesh['Nodes']:
        xn = mesh['Nodes'][nTag]['coords']

        #Left (+X)/Right (-X) boundaries
        if xn[0] - right > 0.0:
            Right.add(nTag)
        elif left - xn[0] > 0.0:
            Left.add(nTag)
        
        #Top (+Z)/Bottom (-Z) boundaries
        if xn[1] - top > 0.0:
            Top.add(nTag)
        elif bottom - xn[1] > 0.0:
            Bottom.add(nTag)

    #Appends the boundary Nodes information
    mesh['Boundary'] = {'top': list(Top), 'bottom': list(Bottom), 'left': list(Left), 'right': list(Right)}

def Find3DBoundaries(mesh, top=float('inf'), bottom=float('inf'), left=float('inf'), right=float('inf'), front=float('inf'), back=float('inf')):
    """
    This function is a brute force method to find the nodes that are on the boundary
    assuming a rectangular prism domain. Note that the inputs top, bottom, left, right 
    must specify a value a little bit inside the domain. \n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    mesh : dict
        'Nodes'  : (dict) The Node information using SVL format

    Returns
    -------
    mesh['Boundary']: dict
        'Top' : set
            The list of Nodes on the top boundary, i.e., +Z
        'Bottom' : set
            The list of Nodes on the bottom boundary, i.e., -Z
        'Left' : set
            The list of Nodes on the left boundary, i.e., +X
        'Right' : set
            The list of Nodes on the right boundary, i.e., -X
        'Front' : set
            The list of Nodes on the front boundary, i.e., +Y
        'Back' : set
            The list of Nodes on the back boundary, i.e., -Y
    """
    #Nodes that belong to boundary 
    Top, Bottom = set(), set()
    Left, Right = set(), set()
    Back, Front = set(), set()

    for nTag in mesh['Nodes']:
        xn = mesh['Nodes'][nTag]['coords']

        #Left (+X)/Right (-X) boundaries
        if xn[0] - right > 0.0:
            Right.add(nTag)
        elif left - xn[0] > 0.0:
            Left.add(nTag)
        
        #Top (+Z)/Bottom (-Z) boundaries
        if xn[2] - top > 0.0:
            Top.add(nTag)
        elif bottom - xn[2] > 0.0:
            Bottom.add(nTag)

        #Front (+Y)/Back (-Y) boundaries 
        if xn[1] - front > 0.0:
            Front.add(nTag)
        elif back - xn[1] > 0.0:
            Back.add(nTag)

    #Appends the boundary Nodes information
    mesh['Boundary'] = {'top': list(Top), 'bottom': list(Bottom), 'left': list(Left), 'right': list(Right), 'front': list(Front), 'back': list(Back)}

def Coords2Tag(mesh, xp, tol=1E-6):
    """
    This function finds the Node tag associated to the specified coordinates.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2021

    Parameters
    ----------
    mesh : dict
        'Nodes'  : dict
            Dictionary that contains the node information in SVL format
        'Elements'  : dict
            Dictionary that contains the element information in SVL format
    xp: numpy matrix
        The coordinates (matrix) to be find
    tol: float
        Tolerance for which the node is accepted

    Returns
    -------
    cPoints  : list
        List with the node Tags 
    """
    cPoints = []
    for nTag in mesh['Nodes']:
        for x in xp:
            if np.linalg.norm(x - mesh['Nodes'][nTag]['coords']) < tol:
                cPoints.append(nTag)
    return cPoints

def setPMLattributes(mesh, x0, xl):
    if Options['dimension'] == 2:
        for eTag in mesh['Elements']:
            coords = np.full(2,  0.0, dtype=float)
            for nTag in mesh['Elements'][eTag]['conn']:
                coords +=  mesh['Nodes'][nTag]['coords']
            xavg = coords/len( mesh['Elements'][eTag]['conn'])

            attributes = copy.deepcopy( mesh['Elements'][eTag]['attributes'])
            if xavg[0] < (x0[0] - xl[0]):
                if xavg[1] < (x0[1] - xl[1]):
                    attributes['x0'] = [x0[0]-xl[0], x0[1]-xl[1]]
                    attributes['npml'] = [-1.0/np.sqrt(2), -1.0/np.sqrt(2)]
                else:
                    attributes['x0'] = [x0[0]-xl[0], x0[1]-xl[1]/2.0]
                    attributes['npml'] = [-1.0, 0.0]
            elif xavg[0] > (x0[0] + xl[0]):
                if xavg[1] < (x0[1] - xl[1]):
                    attributes['x0'] = [x0[0]+xl[0], x0[1]-xl[1]]
                    attributes['npml'] = [1.0/np.sqrt(2), -1.0/np.sqrt(2)]
                else:
                    attributes['x0'] = [x0[0]+xl[0], x0[1]-xl[1]/2.0]
                    attributes['npml'] = [1.0, 0.0]
            else:
                attributes['x0'] = [x0[0], x0[1]-xl[1]]
                attributes['npml'] = [0.0, -1.0]
            mesh['Elements'][eTag]['attributes'] = attributes
    elif Options['dimension'] == 3:
        for eTag in mesh['Elements']:
            coords = np.full(3,  0.0, dtype=float)
            for nTag in mesh['Elements'][eTag]['conn']:
                coords +=  mesh['Nodes'][nTag]['coords']
            xavg = coords/len(mesh['Elements'][eTag]['conn'])

            attributes = copy.deepcopy(mesh['Elements'][eTag]['attributes'])
            if xavg[0] < (x0[0] - xl[0]):
                if xavg[1] < (x0[1] - xl[1]):
                    if xavg[2] < (x0[2] - xl[2]):
                        attributes['x0'] = [x0[0]-xl[0], x0[1]-xl[1], x0[2]-xl[2]]
                        attributes['npml'] = [-1.0/np.sqrt(3), -1.0/np.sqrt(3), -1.0/np.sqrt(3)]
                    else:
                        attributes['x0'] = [x0[0]-xl[0], x0[1]-xl[1], x0[2]-xl[2]/2.0]
                        attributes['npml'] = [-1.0/np.sqrt(2), -1.0/np.sqrt(2), 0.0]
                elif xavg[1] > (x0[1] + xl[1]):
                    if xavg[2] < (x0[2] - xl[2]):
                        attributes['x0'] = [x0[0]-xl[0], x0[1]+xl[1], x0[2]-xl[2]]
                        attributes['npml'] = [-1.0/np.sqrt(3), 1.0/np.sqrt(3), -1.0/np.sqrt(3)]
                    else:
                        attributes['x0'] = [x0[0]-xl[0], x0[1]+xl[1], x0[2]-xl[2]/2.0]
                        attributes['npml'] = [-1.0/np.sqrt(2), 1.0/np.sqrt(2), 0.0]
                else:
                    if xavg[2] < (x0[2] - xl[2]):
                        attributes['x0'] = [x0[0]-xl[0], x0[1], x0[2]-xl[2]]
                        attributes['npml'] = [-1.0/np.sqrt(2), 0.0, -1.0/np.sqrt(2)]
                    else:
                        attributes['x0'] = [x0[0]-xl[0], x0[1], x0[2]-xl[2]/2.0]
                        attributes['npml'] = [-1.0, 0.0, 0.0]
            elif xavg[0] > (x0[0] + xl[0]):
                if xavg[1] < (x0[1] - xl[1]):
                    if xavg[2] < (x0[2] - xl[2]):
                        attributes['x0'] = [x0[0]+xl[0], x0[1]-xl[1], x0[2]-xl[2]]
                        attributes['npml'] = [1.0/np.sqrt(3), -1.0/np.sqrt(3), -1.0/np.sqrt(3)]
                    else:
                        attributes['x0'] = [x0[0]+xl[0], x0[1]-xl[1], x0[2]-xl[2]/2.0]
                        attributes['npml'] = [1.0/np.sqrt(2), -1.0/np.sqrt(2), 0.0]
                elif xavg[1] > (x0[1] + xl[1]):
                    if xavg[2] < (x0[2] - xl[2]):
                        attributes['x0'] = [x0[0]+xl[0], x0[1]+xl[1], x0[2]-xl[2]]
                        attributes['npml'] = [1.0/np.sqrt(3), 1.0/np.sqrt(3), -1.0/np.sqrt(3)]
                    else:
                        attributes['x0'] = [x0[0]+xl[0], x0[1]+xl[1], x0[2]-xl[2]/2.0]
                        attributes['npml'] = [1.0/np.sqrt(2), 1.0/np.sqrt(2), 0.0]
                else:
                    if xavg[2] < (x0[2] - xl[2]):
                        attributes['x0'] = [x0[0]+xl[0], x0[1], x0[2]-xl[2]]
                        attributes['npml'] = [1.0/np.sqrt(2), 0.0, -1.0/np.sqrt(2)]
                    else:
                        attributes['x0'] = [x0[0]+xl[0], x0[1], x0[2]-xl[2]/2.0]
                        attributes['npml'] = [1.0, 0.0, 0.0]
            else:
                if xavg[1] < (x0[1] - xl[1]):
                    if xavg[2] < (x0[2] - xl[2]):
                        attributes['x0'] = [x0[0], x0[1]-xl[1], x0[2]-xl[2]]
                        attributes['npml'] = [0.0, -1.0/np.sqrt(2), -1.0/np.sqrt(2)]
                    else:
                        attributes['x0'] = [x0[0], x0[1]-xl[1], x0[2]-xl[2]/2.0]
                        attributes['npml'] = [0.0, -1.0, 0.0]
                elif xavg[1] > (x0[1] + xl[1]):
                    if xavg[2] < (x0[2] - xl[2]):
                        attributes['x0'] = [x0[0], x0[1]+xl[1], x0[2]-xl[2]]
                        attributes['npml'] = [0.0, 1.0/np.sqrt(2), -1.0/np.sqrt(2)]
                    else:
                        attributes['x0'] = [x0[0], x0[1]+xl[1], x0[2]-xl[2]/2.0]
                        attributes['npml'] = [0.0, 1.0, 0.0]
                else:
                    attributes['x0'] = [x0[0], x0[1], x0[2]-xl[2]]
                    attributes['npml'] = [0.0, 0.0, -1.0]
            mesh['Elements'][eTag]['attributes'] = attributes

def setPMLDomain(attributes, x0, xl):
    """
    This function creates a PML layer to be attached on a rectangular domain. The PML 
    layer dimension is specified in attributes while x0 and xl are center and sides of 
    the rectangular domain where the soil domain will be placed. Note that the mesh size 
    of the PML and the Soil must coincide at the interface.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2021

    Parameters
    ----------
    attributes : dict
        Same fields required to call makeDomainArea or makeDomainVolume
    x0: array
        The coordinate of center of the box
    xl: array
        The side half-length in each direction

    Returns
    -------
    mesh['Boundary']: dict
        'Top' : set
            The list of Nodes on the top boundary, i.e., +Z
        'Bottom' : set
            The list of Nodes on the bottom boundary, i.e., -Z
        'Left' : set
            The list of Nodes on the left boundary, i.e., +X
        'Right' : set
            The list of Nodes on the right boundary, i.e., -X
        'Front' : set
            The list of Nodes on the front boundary, i.e., +Y
        'Back' : set
            The list of Nodes on the back boundary, i.e., -Y
        'Interface'  : list
            List with the Nodes Tags on the interface cut 
    """
    #Create a Mesh 
    if Options['dimension'] == 2:
        mesh = makeDomainArea(options=attributes)
    elif Options['dimension'] == 3:
        mesh = makeDomainVolume(options=attributes)

    #Identify the elements to be removed
    rmElems = set()
    cpElems = set()
    for eTag in mesh['Elements']:
        count = 0
        nodes = mesh['Elements'][eTag]['conn']
        nNodes = len(nodes)
        for nTag in nodes:
            xn = mesh['Nodes'][nTag]['coords']
            cn = np.abs(xn - x0) <= xl
            if cn.all():
                count += 1
        if 0 < count and count < nNodes:
            cpElems.add(eTag)
        elif count == nNodes:
            rmElems.add(eTag)

    #List of Nodes to be removed
    removeNodes = set()
    for eTag in rmElems:
        nodes = mesh['Elements'][eTag]['conn']
        for nTag in nodes:
            removeNodes.add(nTag)

    #Identify Nodes of the cut section 
    keepNodes = set()
    for eTag in cpElems:
        nodes = mesh['Elements'][eTag]['conn']
        for nTag in nodes:
            xn = mesh['Nodes'][nTag]['coords']
            cn = np.abs(xn - x0) <= xl
            if cn.all():
                keepNodes.add(nTag)

    rmNodes = removeNodes.difference(keepNodes)

    #Remove the elements and associated nodes
    for nTag in rmNodes:
        del mesh['Nodes'][nTag]
    for eTag in rmElems:
        del mesh['Elements'][eTag]

    #Re-Numbers Nodes and Elements and Interface
    nMap = reLabel(mesh)
    cpNodes = list()
    for k in keepNodes:
        cpNodes.append(nMap[k])

    setPMLattributes(mesh, x0, xl)

    #Finds the Boundary Nodes list
    if Options['dimension'] == 2:
        TOLx = (attributes['P1'][0] - attributes['P0'][0])/attributes['ne'][0]/1000.0
        TOLz = (attributes['P2'][1] - attributes['P0'][1])/attributes['ne'][1]/1000.0
        Find2DBoundaries(mesh, top=attributes['P2'][1]-TOLz, bottom=attributes['P0'][1]+TOLz, left=attributes['P0'][0]+TOLx, right=attributes['P1'][0]-TOLx)
    elif Options['dimension'] == 3:
        TOLx = (attributes['P1'][0] - attributes['P0'][0])/attributes['ne'][0]/1000.0
        TOLy = (attributes['P2'][1] - attributes['P0'][1])/attributes['ne'][1]/1000.0
        TOLz = (attributes['P3'][2] - attributes['P0'][2])/attributes['ne'][2]/1000.0
        Find3DBoundaries(mesh, top=attributes['P3'][2]-TOLz, bottom=attributes['P0'][2]+TOLz, left=attributes['P0'][0]+TOLx, right=attributes['P1'][0]-TOLx, front=attributes['P2'][1]-TOLy, back=attributes['P0'][1]+TOLy)

    mesh['Boundary']['Interface'] = cpNodes

    return mesh

def setDRMDomain(mesh, x0, xl):
    """
    This function finds the DRM elements provided with a rectangular box. The 
    x0 is the coordinate of center of the box, and xl the half-length side in 
    each direction of the box. These sides (planes) should be inside the DRM 
    elements to be identified.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2021

    Parameters
    ----------
    mesh : dict
        'Nodes'  : dict
            Dictionary that contains the node information in SVL format
        'Elements'  : dict
            Dictionary that contains the element information in SVL format
    x0: array
        The coordinate of center of the box
    xl: array
        The side half-length in each direction

    Returns
    -------
    DRM : dict
        Elements  : list
            List with the DRM element Tags 
        Interior  : list
            List with the DRM internal node Tags 
        Exterior  : list
            List with the DRM external node Tags 
    """
    #Identify the DRM elements
    elemDRM = set()
    for eTag in mesh['Elements']:
        count = 0
        nodes = mesh['Elements'][eTag]['conn']
        for nTag in nodes:
            xn = mesh['Nodes'][nTag]['coords']
            cn = np.abs(xn - x0) <= xl
            if cn.all():
                count += 1
        if 0 < count and count < len(nodes):
            elemDRM.add(eTag)

    #Interior/Exterior DRM node lists
    intDRM = set()
    extDRM = set()
    for eTag in elemDRM:
        nodes = mesh['Elements'][eTag]['conn']
        for nTag in nodes:
            xn = mesh['Nodes'][nTag]['coords']
            cn = np.abs(xn - x0) <= xl
            if cn.all():
                intDRM.add(nTag)
            else:
                extDRM.add(nTag)

    #Provides the DRM information in Mesh
    if 'DRM' not in mesh:
         mesh['DRM'] = dict()

    mesh['DRM']['Interior'] = list(intDRM)
    mesh['DRM']['Exterior'] = list(extDRM)
    mesh['DRM']['Elements'] = list(elemDRM)

def setRestrains(mesh, dof=[], bc=[]):
    """
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020
    """
    for names in bc:
        if names in ['top','bottom','left','right','back','front']:
            for nTag in mesh['Boundary'][names]:
                for n in dof:
                    mesh['Nodes'][nTag]['freedof'][n-1] = -1

def setLysmerDomain(mesh, bc=[], mat=[]):
    """
    This is an approximation of the Lysmer absorbent condition, since it double or triple counts the dashpots at the boundary interseptions
    """
    ndof = Options['dimension']
    nNodes = max(mesh['Nodes'].keys())
    nElems = max(mesh['Elements'].keys())
    
    #Material indexes
    if ndof == 2:
        ind = [[0,1],[1,0]]
    elif ndof == 3:
        ind = [[0,1,2],[2,0,1],[0,2,1]]
    
    #Gets how many nodes belongs to one or more interfaces
    bcTags = dict()
    for names in bc:
        nTags = mesh['Boundary'][names]
        for nTag in nTags:
            if nTag not in bcTags:
                bcTags[nTag] = [names]
            else:
                bcTags[nTag].append(names)

    nTag = nNodes
    eTag = nElems
    for n in bcTags:
        #Creates the new Node
        nTag += 1
        free  = np.full(ndof, -1, dtype=int)
        total = np.full(ndof,  0, dtype=int)
        mesh['Nodes'][nTag] = {'ndof': ndof, 'freedof': free, 'totaldof': total, 'coords': mesh['Nodes'][n]['coords']}

        for bc in bcTags[n]:
            if bc == 'top' or bc == 'bottom':
                for k in range(ndof):
                    eTag += 1
                    mesh['Elements'][eTag] = {'name': 'ZEROLENGTH1D', 'conn': [n, nTag], 'attributes': {'material': mat[ind[0][k]], 'dir': k}}
            elif bc == 'left' or bc == 'right':
                for k in range(ndof):
                    eTag += 1
                    mesh['Elements'][eTag] = {'name': 'ZEROLENGTH1D', 'conn': [n, nTag], 'attributes': {'material': mat[ind[1][k]], 'dir': k}}
            elif bc == 'back' or bc == 'front':
                for k in range(ndof):
                    eTag += 1
                    mesh['Elements'][eTag] = {'name': 'ZEROLENGTH1D', 'conn': [n, nTag], 'attributes': {'material': mat[ind[2][k]], 'dir': k}}
