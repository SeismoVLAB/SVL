#!/usr/bin/env python3
# -*- coding: Utf-8 -*-

import numpy as np
from Core.Utilities import debugInfo
from Core.Definitions import Options

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
                print(tag-1, ' ', [n1, n2, n4, n5, n9, n8])
                print(tag  , ' ', [n2, n3, n4, n6, n7, n9])
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

def removeDomain(mesh={}, atributes={}):
    """
    This function removes a domain specifying the geometry\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    mesh : dict
        The mesh that contains 'Nodes' and 'Elements' dictionaries from which
        the deletion will be performed. 
    atributes : dict
        The domain to be removed, for example
        'name' = Rectangular
            'sides' (list) The length of the sides in x,y, and z
            'center' (list) The center of the rectangular domain
        'name' = Circular
            'radius' (float) The radius of the circle/sphere
            'center' (list) The center of the circular domain

    Returns
    -------
    None
    """
    name = atributes['name'].upper()
    
    s1 = set()
    s2 = set()
    eElems = set()
    if name == 'RECTANGULAR':
        x0 = np.array(atributes['center'])
        xl = np.array(atributes['sides'])
        nd = x0.size
        for eTag in mesh['Elements']:
            nodes = set()
            xm = np.zeros(nd, dtype=float)
            for nTag in mesh['Elements'][eTag]['conn']:
                nodes.add(nTag)
                xm += mesh['Nodes'][nTag]['coords']
            xm = xm/len(nodes)

            #The element is inside/outside the remove domain
            cond = np.abs(xm - x0) < xl/2.0
            if cond.all():
                eElems.add(eTag)
                s1.update(nodes)
            else:
                s2.update(nodes)
    elif name == 'CIRCULAR':
        x0 = np.array(atributes['center'])
        r  = np.array(atributes['radius'])
        nd = x0.size
        for eTag in mesh['Elements']:
            nodes = set()
            xm = np.zeros(nd, dtype=float)
            for nTag in mesh['Elements'][eTag]['conn']:
                nodes.add(nTag)
                xm += mesh['Nodes'][nTag]['coords']
            xm = xm/len(nodes)
            
            #The element is inside/outside the remove domain
            if np.linalg.norm(xm - x0) < r:
                eElems.add(eTag)
                s1.update(nodes)
            else:
                s2.update(nodes)
    else:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d removeDomain(optiatributesons[\'name\']=%s) is not recognized.' %(info.filename,info.lineno,name))

    #Nodes (degenerate) to be removed
    eNodes = s1.difference(s2)

    #Remove the elements and associated nodes
    for nTag in eNodes:
        del mesh['Nodes'][nTag]
    for eTag in eElems:
        del mesh['Elements'][eTag]
    return mesh

def mergeDomain(mesh1={}, mesh2={}):
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

    Returns
    -------
    None
    """
    nNodes = len(mesh1['Nodes'])
    nElems = len(mesh1['Elements'])

    #Nodes are updated with information in mesh2
    for nTag in mesh2['Nodes']:
        tag = nTag + nNodes
        mesh1['Nodes'][tag] = mesh2['Nodes'][nTag] 

    #Elements are updated with information in mesh2
    for eTag in mesh2['Elements']:
        tag = eTag + nElems
        mesh1['Elements'][tag] = mesh2['Elements'][eTag]
        mesh1['Elements'][tag]['conn'] = [x+nNodes for x in mesh1['Elements'][tag]['conn']]

def setRestraints(Mesh={}, options={}):
    """
    Apply constraints to the requested faces.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    options : dict
        Specify some parameters for constructing the mesh
        'ndof'   : (int) number of degree of freedom per node
        'elems'  : (str) the type of element during meshing,i.e., LINE, LINE3, QUAD4, QUAD8, HEXA8, HEXA20 etc 
        'Borders': (dict) 
            'name' : (str)'restrain' or 'constraint'  
            'faces': (list) face number to be restrained 
                     LINE   : 1:left, 2:right
                     AREA   : 1:left, 2:bottom, 3: right, 4:top
                     VOLUME : 1:left (-X), 2:bottom (-Y), 3: right (+X), 4:top (+Y), 5: bottom (-Z), 6:top (+Z)
            'dof'  : (list) Degree of freedom to restrain, i.e.,  0: Free, 1: Restrain

    Returns
    -------
    None
    """
    #Checks that inputs are not empty
    if not Mesh:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d setRestraints(Mesh=?) must be specified.' %(info.filename,info.lineno))
        return False
    if not options:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d setRestraints(options=?) must be specified.' %(info.filename,info.lineno))
        return False
    if 'Borders' not in options:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d the options[\'Borders\'] must be specified.' %(info.filename,info.lineno))
        return False

    #Unpack attribute provided by the user
    name = options['Borders']['name']
    dof  = options['Borders']['dof']

    if name.upper() != 'RESTRAIN':
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d the options[\'Borders\'][\'name\']=restrain must be specified.' %(info.filename,info.lineno))
        return False

    #Gets the nodes to be restrained
    nTag = boundaryNodes(options)

    for n in nTag:
        for k, d in enumerate(dof):
            if d == 1:
                Mesh['Nodes'][n]['freedof'][k] = -1
    return True

def boundaryNodes(options):
    """
    This auxiliar function gets the face node for a simple domain\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    options : dict
        Specify some parameters for finding the boundary nodes
        'elems'  : (str) the type of element during meshing,i.e., LINE, LINE3, QUAD4, QUAD8, HEXA8, HEXA20 etc 
        'Borders': (dict) 
            'name' : (str)'restrain' or 'constraint'  
            'faces': (list) face number to be restrained
                     LINE   : 1:left, 2:right
                     AREA   : 1:left, 2:bottom, 3: right, 4:top
                     VOLUME : 1:left (-X), 2:bottom (-Y), 3: right (+X), 4:top (+Y), 5: bottom (-Z), 6:top (+Z)
            'dof'  : (list) Degree of freedom to restrain, i.e.,  0: Free, 1: Restrain
    """
    #TODO: Implement for TETRA10, HEXA20
    #Faces to identify nodes
    elems = options['elems']
    faces = options['Borders']['faces']

    #The type of domain for which nodes are detected
    if elems.upper() == 'LINE2' or elems.upper() == 'LINE3':
        nTag = []
        nx = options['ne']
        for f in faces:
            if f == 1:
                nTag.append(1)
            elif f == 2:
                nTag.append(nx+1)
    elif elems.upper() == 'TRIA3' or elems.upper() == 'QUAD4':
        nTag = set()
        nx = options['ne'][0]
        ny = options['ne'][1]
        for f in faces:
            if f == 1:
                nTag.update(range(1, (nx+2)*ny, nx+1))
            elif f == 2:
                nTag.update(range(1, nx+2))
            elif f == 3:
                nTag.update(range(nx+1, (nx+2)*(ny+1), nx+1))
            elif f == 4:
                nTag.update(range((nx+1)*ny+1, (nx+1)*(ny+1)+1))
    elif elems.upper() == 'TRIA6' or elems.upper() == 'QUAD8':
        nTag = set()
        nx = options['ne'][0]
        ny = options['ne'][1]
        for f in faces:
            if f == 1:
                nTag.update(range(1, (nx+2)*ny, nx+1))
                nTag.update(range((nx+1)*(ny+1)+nx+1, (nx+1)*(ny+1)+nx*(ny+1)+(nx+1)*ny, 2*nx+1))
            elif f == 2:
                nTag.update(range(1, nx+2))
                nTag.update(range((nx+1)*(ny+1)+1, (nx+1)*(ny+2)))
            elif f == 3:
                nTag.update(range(nx+1, (nx+2)*(ny+1), nx+1))
                nTag.update(range((nx+1)*(ny+1)+2*nx+1, (nx+1)*(ny+1)+nx*(ny+1)+(nx+1)*ny, 2*nx+1))
            elif f == 4:
                nTag.update(range((nx+1)*ny+1, (nx+1)*(ny+1)+1))
                nTag.update(range((nx+1)*(ny+1)+(2*nx+1)*ny+1, (nx+1)*(ny+1)+nx*(ny+1)+(nx+1)*ny+1))
    elif elems.upper() == 'TETRA4' or elems.upper() == 'HEXA8':
        nTag = set()
        nx = options['ne'][0]
        ny = options['ne'][1]
        nz = options['ne'][2]
        for f in faces:
            if f == 1:
                nTag.update(range(1,(nx+1)*(ny+1)*(nz+1),nx+1))
            elif f == 2:
                for k in range(nz+1):
                    nTag.update(range(1+(nx+1)*(ny+1)*k, 1+(nx+1)*((ny+1)*k + 1)))
            elif f == 3:
                nTag.update(range(nx+1,(nx+1)*(ny+1)*(nz+1)+1,nx+1))
            elif f == 4:
                for k in range(nz+1):
                    nTag.update(range(1+(nx+1)*(ny + (ny+1)*k), 1+(nx+1)*(ny + 1 + (ny+1)*k)))
            elif f == 5:
                nTag.update(range(1,(nx+1)*(ny+1)+1))
            elif f == 6:
                nTag.update(range((nx+1)*(ny+1)*nz+1,(nx+1)*(ny+1)*(nz+1)+1))
    return nTag

def setPMLDomain():
    """
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020
    """
    #TODO: Implement this feature
    mesh = {}
    return mesh

def setDRMDomain():
    """
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020
    """
    #TODO: Implement this feature
    mesh = {}
    return mesh

def setConstraints():
    """
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020
    """
    #TODO: Implement this feature
    pass
