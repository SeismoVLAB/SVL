#!/usr/bin/env python3
# -*- coding: Utf-8 -*-

import shlex
import numpy as np

def parseABAQUS(filepath):
    """
    This function parses a *.inp file generated from ABAQUS software. Not all functionalities 
    are available, only coordinates, and elements are obtained during the parsing process.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    filepath: str
        The path were the ABAQUS file will be read

    Returns
    -------
    mesh : dict
        A dictionary that contains nodes and element with the same structure as Entities
    """
    #TODO: Implement ABAQUS parser 
    print('\x1B[31m ERROR \x1B[0m: The ABAQUS parser is not implemented yet')

    #The MESH dictionary containing the data is read
    mesh = {'Nodes': {}, 'Materials': {}, 'Sections':{}, 'Elements':{}}

    #Open the provided file
    with open(filepath, 'r', errors='replace') as f:
        lines = f.readlines()

    #Loop over the entire file
    k = 0
    while k < len(lines):
        line = list(shlex.split(lines[k]))
        k += 1

    return mesh