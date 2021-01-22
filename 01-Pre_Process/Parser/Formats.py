#!/usr/bin/python3
# -*- coding: Utf-8 -*-

import json
from Parser.GMSH import parseGMSH
from Parser.ETABS import parseETABS
from Parser.SAP2000 import parseSAP
from Core.Utilities import debugInfo

def parseFile(filepath=None, fileformat=None):
    """
    Parse the input file with a specific style/format\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    filepath : str
        The full-path where the files is located
    format : str
        The file format or layout style (style=SAP, ETABS, gmsh, other)

    Returns
    -------
    bool
        Whether the file was successful (True) of failed (False) during parsing
    """
    if fileformat.upper() == 'JSON':
        mesh = parseJSON(filepath)
    elif fileformat.upper() == 'GMSH':
        mesh = parseGMSH(filepath)
    elif fileformat.upper() == 'ETABS':
        mesh = parseETABS(filepath)
    elif fileformat.upper() == 'SAP':
        mesh = parseSAP(filepath)
    else:
        mesh = {}
        info = debugInfo(2)
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d the format=%s is not implemented yet.' %(info.filename,info.lineno,fileformat))
    return mesh

def parseJSON(filepath=None):
    """
    Parse the input file with a specific style/format\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    filepath : str
        The full-path where the files is located

    Returns
    -------
    bool
        Whether the file was successful (True) of failed (False) during parsing
    """
    with open(filepath, 'r') as myfile:
        JSONdata = myfile.read()
    mesh = json.loads(JSONdata)

    return mesh