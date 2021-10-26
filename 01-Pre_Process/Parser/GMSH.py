#!/usr/bin/env python3
# -*- coding: Utf-8 -*-

import shlex
import copy
import numpy as np
from Core.Utilities import *

def parseGMSH(filepath):
    """
    This function parses a *.mesh (INRIA) files from GMSH software. Not all 
    functionalities are available, only coordinates, elements, boundaries and 
    interfaces can be specified using "Physical Tags." The properties must
    be associated to this value during the parsing process. The user must use 
    this information in order to create the model; therefore, Entities such as 
    Nodes, Materials, Elements, Loads etc must be generated in the python script.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2021

    Parameters
    ----------
    filepath: str
        The path were the GMSH (.mesh) file will be read

    Returns
    -------
    mesh : dict
        A dictionary that contains nodes and element information. This structure
        DOES NOT coincide with the Entities structure
    """
    #The MESH dictionary containing the data is read
    mesh = {'Vertices': {}, 'Edges': {}, 'Triangles':{}, 'Quadrilaterals':{}, 'Tetrahedra':{}, 'Hexahedra':{}}

    #Open the provided file
    with open(filepath, 'r', errors='replace') as f:
        lines = f.readlines()

    #PARSE THE DATA INTO ARRAYS
    k = 0
    while k < len(lines):
        line = lines[k].strip()
        if line == "Vertices":
            k += 1
            ne = int(shlex.split(lines[k])[0])
            Vertices = np.zeros((ne, 3), dtype=np.float32)

            for j in range(ne):
                k += 1
                vals = list(shlex.split(lines[k]))
                Vertices[j,:] = np.array([float(vals[0]), float(vals[1]), float(vals[2])])
            mesh['Vertices'] = Vertices;
            k += 1
        elif line == "Edges":
            k += 1
            ne = int(shlex.split(lines[k])[0])
            Edges = np.zeros((ne, 3), dtype=np.int32)

            for j in range(ne):
                k += 1
                vals = list(shlex.split(lines[k]))
                Edges[j,:] = np.array([int(vals[0]), int(vals[1]), int(vals[2])])
            mesh['Edges']['Data'] = Edges;
            k += 1
        elif line == "Triangles":
            k += 1
            ne = int(shlex.split(lines[k])[0])
            Tria = np.zeros((ne, 4), dtype=np.int32)

            for j in range(ne):
                k += 1
                vals = list(shlex.split(lines[k]))
                Tria[j,:] = np.array([int(vals[0]), int(vals[1]), int(vals[2]), int(vals[3])])
            mesh['Triangles']['Data'] = Tria;
            k += 1
        elif line == "Quadrilaterals":
            k += 1
            ne = int(shlex.split(lines[k])[0])
            Quads = np.zeros((ne, 5), dtype=np.int32)

            for j in range(ne):
                k += 1
                vals = list(shlex.split(lines[k]))
                Quads[j,:] = np.array([int(vals[0]), int(vals[1]), int(vals[2]), int(vals[3]), int(vals[4])])
            mesh['Quadrilaterals']['Data'] = Quads;
            k += 1
        elif line == "Tetrahedra":
            k += 1
            ne = int(shlex.split(lines[k])[0])
            Tetra = np.zeros((ne, 5), dtype=np.int32)

            for j in range(ne):
                k += 1
                vals = list(shlex.split(lines[k]))
                Tetra[j,:] = np.array([int(vals[0]), int(vals[1]), int(vals[2]), int(vals[3]), int(vals[4])])
            mesh['Tetrahedra']['Data'] = Tetra;
            k += 1
        elif line == "Hexahedra":
            k += 1
            ne = int(shlex.split(lines[k])[0])
            Hexas = np.zeros((ne, 9), dtype=np.int32)

            for j in range(ne):
                k += 1
                vals = list(shlex.split(lines[k]))
                Hexas[j,:] = np.array([int(vals[0]), int(vals[1]), int(vals[2]), int(vals[3]), int(vals[4]), int(vals[5]), int(vals[6]), int(vals[7]), int(vals[8])])
            mesh['Hexahedra']['Data'] = Hexas;
            k += 1
        else:
            k += 1

    #PROCESS INDEXES ACCORDING TO PHYSICAL GROUP
    for option in mesh:
        if option == 'Edges':
            if 'Data' in mesh['Edges']:
                indexes = np.unique(mesh['Edges']['Data'][:,2])

                ElemEdges = {k: set() for k in indexes}
                NodeEdges = {k: set() for k in indexes}
                for k, ed in enumerate(mesh['Edges']['Data']):
                    NodeEdges[ed[2]].add(int(ed[0]))
                    NodeEdges[ed[2]].add(int(ed[1]))
                    ElemEdges[ed[2]].add(k)
                mesh['Edges']['Node'] = NodeEdges
                mesh['Edges']['Elem'] = ElemEdges
        elif option == 'Triangles':
            if 'Data' in mesh['Triangles']:
                indexes = np.unique(mesh['Triangles']['Data'][:,3])

                ElemTrias = {k: set() for k in indexes}
                NodeTrias = {k: set() for k in indexes}
                for k, ed in enumerate(mesh['Triangles']['Data']):
                    NodeTrias[ed[3]].add(int(ed[0]))
                    NodeTrias[ed[3]].add(int(ed[1]))
                    NodeTrias[ed[3]].add(int(ed[2]))
                    ElemTrias[ed[3]].add(k)
                mesh['Triangles']['Node'] = NodeTrias
                mesh['Triangles']['Elem'] = ElemTrias
        elif option == 'Quadrilaterals':
            if 'Data' in mesh['Quadrilaterals']:
                indexes = np.unique(mesh['Quadrilaterals']['Data'][:,4])

                ElemQuads = {k: set() for k in indexes}
                NodeQuads = {k: set() for k in indexes}
                for k, ed in enumerate(mesh['Quadrilaterals']['Data']):
                    NodeQuads[ed[4]].add(int(ed[0]))
                    NodeQuads[ed[4]].add(int(ed[1]))
                    NodeQuads[ed[4]].add(int(ed[2]))
                    NodeQuads[ed[4]].add(int(ed[3]))
                    ElemQuads[ed[4]].add(k)
                mesh['Quadrilaterals']['Node'] = NodeQuads
                mesh['Quadrilaterals']['Elem'] = ElemQuads
        elif option == 'Tetrahedra':
            if 'Data' in mesh['Tetrahedra']:
                indexes = np.unique(mesh['Tetrahedra']['Data'][:,4])

                ElemTetras = {k: set() for k in indexes}
                NodeTetras = {k: set() for k in indexes}
                for k, ed in enumerate(mesh['Tetrahedra']['Data']):
                    NodeTetras[ed[4]].add(int(ed[0]))
                    NodeTetras[ed[4]].add(int(ed[1]))
                    NodeTetras[ed[4]].add(int(ed[2]))
                    NodeTetras[ed[4]].add(int(ed[3]))
                    ElemTetras[ed[4]].add(k)
                    mesh['Tetrahedra']['Node'] = NodeTetras
                    mesh['Tetrahedra']['Elem'] = ElemTetras
        elif option == 'Hexahedra':
            if 'Data' in mesh['Hexahedra']:
                indexes = np.unique(mesh['Hexahedra']['Data'][:,8])

                ElemHexas = {k: set() for k in indexes}
                NodeHexas = {k: set() for k in indexes}
                for k, ed in enumerate(mesh['Hexahedra']['Data']):
                    NodeHexas[ed[8]].add(int(ed[0]))
                    NodeHexas[ed[8]].add(int(ed[1]))
                    NodeHexas[ed[8]].add(int(ed[2]))
                    NodeHexas[ed[8]].add(int(ed[3]))
                    NodeHexas[ed[8]].add(int(ed[4]))
                    NodeHexas[ed[8]].add(int(ed[5]))
                    NodeHexas[ed[8]].add(int(ed[6]))
                    NodeHexas[ed[8]].add(int(ed[7]))
                    ElemHexas[ed[8]].add(k)
                mesh['Hexahedra']['Node'] = NodeHexas
                mesh['Hexahedra']['Elem'] = ElemHexas

    return mesh
