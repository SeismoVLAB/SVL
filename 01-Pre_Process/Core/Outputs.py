#!/usr/bin/python3
# -*- coding: Utf-8 -*-

import sys
import json
import numpy as np
from json import JSONEncoder
from Core.Definitions import Entities, Options

class NumpyArrayEncoder(JSONEncoder):
    """
    This simple class allow to serialize numpy multi-dimensional
    arrays when using JSON python module. Essentially it transform
    a numpy array into a json list.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020
    """
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(JSONEncoder, self).default(obj)
        return JSONEncoder.default(self, obj)

def dict2json(d, np):
    """
    This function writes a dictionary using JSON format. The serialization
    used was taken from: https://pynative.com/python-serialize-numpy-ndarray-into-json/\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    d : dict
        Dictionary to be transformed into json format
    np : int
        The processor number

    Returns
    -------
    None
    """
    #The JSON output file name
    filename = Options['path'] + '/' + 'Partition' + '/' + Options['file'] + '.' + str(np) + '.json'

    #Serializing json
    JSONdata = json.dumps(d, cls=NumpyArrayEncoder, indent=4)

    #Writes the file in json format
    with open(filename, "w") as outfile: 
        outfile.write(JSONdata)
