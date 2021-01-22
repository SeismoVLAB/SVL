#!/usr/bin/env python3
# -*- coding: Utf-8 -*-

import numpy as np

def parseGMSH(filepath):
    """
    This function parses a *.mesh file generated from GMSH software. Not all functionalities 
    are available, only coordinates, and elements are obtained during the parsing process.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    filepath: str
        The path were the GMSH file will be read

    Returns
    -------
    mesh : dict
        A dictionary that contains nodes and element with the same structure as Entities
    """
    #The MESH dictionary containing the data is read
    mesh = {'Nodes': {}, 'Materials': {}, 'Sections':{}, 'Elements':{}}

    #Open the provided file
    with open(filepath, 'r', errors='replace') as f:
        lines = f.readlines()

    #TODO:
    return mesh