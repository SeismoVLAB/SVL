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

def dict2svl(d, np):
    """
    This function writes a dictionary using SVL old format.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    d : dict
        Dictionary to be transformed into SVL old format
    np : int
        The processor number

    Returns
    -------
    None
    """
    #The SVL output file name
    filepath = Options['path'] + '/' + 'Partition' + '/' + Options['file'] + '.' + str(np) + '.svl'
    
    SeismoVLABfile = open(filepath, "w+")

    #[1] Global model/simulation parameters
    SeismoVLABfile.write("MODEL %d %d %d %s\n" % (Options['dimension'], Options['ntotal'], Options['nfree'], Options['massform'].upper()))

    #[2] Materials
    for mTag in d['Materials']:
        #Unpack fields for simpler writing
        Material = d['Materials'][mTag]
        attributes = Material['attributes']
        NAME = Material['name'].upper()

        SeismoVLABfile.write("MATERIAL %s %s " %(mTag, NAME))

        #The possible Materials in SeismoVLAB
        if NAME == "ELASTIC1DLINEAR":
            SeismoVLABfile.write("%E %1.5f %1.5f\n" % (attributes['E'], attributes['nu'], attributes['rho']))
        elif NAME == "HERTZIAN1DLINEAR":
            SeismoVLABfile.write("%E %E %E %E\n" % (attributes['k1'], attributes['k2'], attributes['k3'], attributes['rho']))
        elif NAME == "VISCOUS1DLINEAR":
            SeismoVLABfile.write("%E\n" % attributes['eta'])
        elif NAME == "ELASTIC2DPLANESTRAIN":
            SeismoVLABfile.write("%E %1.5f %1.5f\n" % (attributes['E'], attributes['nu'], attributes['rho']))
        elif NAME == "ELASTIC2DPLANESTRESS":
            SeismoVLABfile.write("%E %1.5f %1.5f\n" % (attributes['E'], attributes['nu'], attributes['rho']))
        elif NAME == "ELASTIC3DLINEAR":
            SeismoVLABfile.write("%E %1.5f %1.5f\n" % (attributes['E'], attributes['nu'], attributes['rho']))
        elif NAME == "PLASTIC1DJ2":
            SeismoVLABfile.write("%E %1.5f %1.5f %E %E %E\n" % (attributes['E'], attributes['nu'], attributes['rho'], attributes['k'], attributes['h'], attributes['Sy']))
        elif NAME == "PLASTICPLANESTRAINJ2":
            SeismoVLABfile.write("%E %1.5f %1.5f %E %1.5f %E\n" % (attributes['K'], attributes['G'], attributes['rho'], attributes['h'], attributes['beta'], attributes['Sy']))
        elif NAME == "PLASTICPLANESTRAINBA":
            SeismoVLABfile.write("%E %E %1.5f %1.5f %1.5f %1.5f %E %1.5f\n" % (attributes['K'], attributes['G'], attributes['rho'], attributes['h0'], attributes['h'], attributes['m'], attributes['Su'], attributes['beta']))
        elif NAME == "PLASTIC3DJ2":
            SeismoVLABfile.write("%E %1.5f %1.5f %1.5f %1.5f %E\n" % (attributes['K'], attributes['G'], attributes['rho'], attributes['h'], attributes['beta'], attributes['Sy']))
        elif NAME == "PLASTIC3DBA":
            SeismoVLABfile.write("%E %E %1.5f %1.5f %1.5f %1.5f %E %1.5f\n" % (attributes['K'], attributes['G'], attributes['rho'], attributes['h0'], attributes['h'], attributes['m'], attributes['Su'], attributes['beta']))
        else:
            print('\x1B[31m ERROR \x1B[0m: In dict2svl() MATERIAL=' + NAME + ' could not be recognized.')
            sys.exit(-1)

    #[3] Sections 
    for sTag in d['Sections']:
        #Unpack fields for simpler writing
        Section = d['Sections'][sTag]
        attributes = Section['attributes']
        NAME = Section['name'].upper()

        SeismoVLABfile.write("SECTION %s %s " % (sTag, Section['model']))

        if Section['model'] == 'PLAIN':
            SeismoVLABfile.write(str(attributes['material']) + ' ' + NAME)

            #Loop Over Section's Prarmeters values
            if NAME == "LIN2DRECTANGULAR" or NAME == "LIN3DRECTANGULAR":
                SeismoVLABfile.write(" %1.5f %1.5f %1.5f %d\n" % (attributes['h'], attributes['b'], attributes['theta'], attributes['ip']))
            elif NAME == "LIN2DCIRCULAR" or NAME == "LIN3DCIRCULAR":
                SeismoVLABfile.write(" %1.5f %1.5f %d\n" % (attributes['r'], attributes['theta'], attributes['ip']))
            elif NAME == "LIN2DANGLE" or NAME == "LIN2DCHANNEL" or NAME == "LIN2DTEE" or NAME == "LIN2DWIDEFLANGE":
                SeismoVLABfile.write(" %1.5f %1.5f %1.5f %1.5f %1.5f %d\n" % (attributes['h'], attributes['b'], attributes['tw'], attributes['tf'], attributes['theta'], attributes['ip']))
            elif NAME == "LIN3DANGLE" or NAME == "LIN3DCHANNEL" or NAME == "LIN3DTEE" or NAME == "LIN3DWIDEFLANGE":
                SeismoVLABfile.write(" %1.5f %1.5f %1.5f %1.5f %1.5f %d\n" % (attributes['h'], attributes['b'], attributes['tw'], attributes['tf'], attributes['theta'], attributes['ip']))
            elif NAME == "LIN2DRECTANGULARTUBE" or NAME == "LIN3DRECTANGULARTUBE":
                SeismoVLABfile.write(" %1.5f %1.5f %1.5f %1.5f %1.5f %d\n" % (attributes['h'], attributes['b'], attributes['tw'], attributes['tf'], attributes['theta'], attributes['ip']))
            elif NAME == "LIN2DCIRCULARTUBE" or NAME == "LIN3DCIRCULARTUBE":
                SeismoVLABfile.write(" %1.5f %1.5f %1.5f %d\n" % (attributes['re'], attributes['ri'], attributes['theta'], attributes['ip']))
            elif NAME == "LIN3DTHINAREA":
                SeismoVLABfile.write(" %1.5f\n" % attributes['th'])
            else:
                print('\x1B[31m ERROR \x1B[0m: In dict2svl() PLAIN SECTION=' + NAME + ' could not be recognized.')
                sys.exit(-1)

        elif Section['model'] == 'FIBER':
            print('\x1B[33m ALERT \x1B[0m: In Section[%d] \'model\'=%s has not been implemented yet.' % (sTag, Section['model']))

    #[4] Nodes
    for nTag in d['Nodes']:
        #Unpack dictionary for simpler writing
        Node = d['Nodes'][nTag]

        SeismoVLABfile.write("NODE %s %d " % (nTag, Node['ndof']))

        #Loop Over Node's Total-Degree-of-Freedom
        total =  Node['totaldof']
        for dof in total:
            SeismoVLABfile.write("%d " % dof)

        #Loop Over Node's Free-Degree-of-Freedom
        free =  Node['freedof']
        for dof in free:
            SeismoVLABfile.write("%d " % dof)

        #Loop Over Node's Coordinates
        coordinates =  Node['coords']
        for coord in coordinates:
            SeismoVLABfile.write("%s " % coord)
        
        SeismoVLABfile.write('\n')

    #[5] Mass
    for nTag in d['Masses']:
        #Unpack dictionary for simpler writing
        Mass = d['Masses'][nTag]

        SeismoVLABfile.write("MASS %s %d" % (nTag, Mass['ndof']))
        for mass in Mass['mass']:
            SeismoVLABfile.write(" %f" % mass)
        SeismoVLABfile.write('\n')

    #[6] Support Motion
    for nTag in d['Supports']:
        #Unpack dictionary for simpler writing
        Support = d['Supports'][nTag]
        for dof in range(len(Support['dof'])):
            SeismoVLABfile.write("SUPPORTMOTION %s %s " % (nTag, Support['type']))
            if(Support['type'] == 'CONSTANT'):
                SeismoVLABfile.write("%E %d" % (Support['value'][dof], Support['dof'][dof]))
            elif(Support['type'] == 'TIMESERIE'):
                SeismoVLABfile.write("%s %d" % (Support['file'][dof], Support['dof'][dof]))
            SeismoVLABfile.write('\n')
    
    #[7] Constraints
    for cTag in d['Constraints']:
        #Unpack dictionary for simpler writing
        Constraint = d['Constraints'][cTag]
        factor = Constraint['factor']
        nCombinations = len(factor)       
        SeismoVLABfile.write("CONSTRAINT %s %d %d" %(cTag, Constraint['stag'], nCombinations))

        for k in range(nCombinations):
            SeismoVLABfile.write(" %d %f" % (Constraint['mtag'][k], factor[k]))
        SeismoVLABfile.write('\n')

    #[8] Elements
    nParaview = 0 
    for eTag in d['Elements']:
        #Unpack dictionary for simpler writing
        Element = d['Elements'][eTag]
        NAME = Element['name'].upper()
        conn = Element['conn']
        attributes = Element['attributes']

        #Writes a text line with the element information
        SeismoVLABfile.write("ELEMENT %s %s " %(eTag, NAME))

        for node in conn:
            SeismoVLABfile.write("%d " % node)

        if NAME == 'ZEROLENGTH1D':
            nParaview += 3
            SeismoVLABfile.write("%d %d %d\n" %(attributes['material'], Options['dimension'], attributes['dir']))
        elif  NAME == 'LIN2DTRUSS2':
            nParaview += 3
            SeismoVLABfile.write("%d %1.5f\n" %(attributes['material'], attributes['area']))
        elif  NAME == 'KIN2DTRUSS2':
            nParaview += 3
            SeismoVLABfile.write("%d %1.5f\n" %(attributes['material'], attributes['area']))
        elif NAME == 'LIN3DTRUSS2':
            nParaview += 3
            SeismoVLABfile.write("%d %1.5f\n" %(attributes['material'], attributes['area']))
        elif NAME == 'KIN3DTRUSS2':
            nParaview += 3
            SeismoVLABfile.write("%d %1.5f\n" %(attributes['material'], attributes['area']))
        elif NAME == 'LIN2DTRUSS3':
            nParaview += 4
            SeismoVLABfile.write("%d %1.5f %s %d\n" %(attributes['material'], attributes['area'], attributes['rule'], attributes['np']))
        elif NAME == 'LIN3DTRUSS3':
            nParaview += 4
            SeismoVLABfile.write("%d %1.5f %s %d\n" %(attributes['material'], attributes['area'], attributes['rule'], attributes['np']))
        elif NAME == 'LIN2DQUAD4':
            nParaview += 5
            SeismoVLABfile.write("%d %1.5f %s %d\n" %(attributes['material'], attributes['th'], attributes['rule'], attributes['np']))
        elif NAME == 'LIN2DQUAD8':
            nParaview += 9
            SeismoVLABfile.write("%d %1.5f %s %d\n" %(attributes['material'], attributes['th'], attributes['rule'], attributes['np']))
        elif NAME == 'PML2DQUAD4':
            nParaview += 5
            X0  = attributes['x0']
            DIR = attributes['npml']
            SeismoVLABfile.write("%d %1.5f %d %1.5f %E %1.5f %1.5f %1.5f %1.5f %s %d\n" %(attributes['material'], attributes['th'], attributes['n'], attributes['L'], attributes['R'], X0[0], X0[1], DIR[0], DIR[1], attributes['rule'], attributes['np']))
        elif NAME == 'PML2DQUAD8':
            nParaview += 9
            X0  = attributes['x0']
            DIR = attributes['npml']
            SeismoVLABfile.write("%d %1.5f %d %1.5f %E %1.5f %1.5f %1.5f %1.5f %s %d\n" %(attributes['material'], attributes['th'], attributes['n'], attributes['L'], attributes['R'], X0[0], X0[1], DIR[0], DIR[1], attributes['rule'], attributes['np']))
        elif  NAME == 'LIN2DFRAME2':
            nParaview += 3
            SeismoVLABfile.write("%d %s %s %d\n" %(attributes['section'], attributes['formulation'], attributes['rule'], attributes['np']))
        elif  NAME == 'KIN2DFRAME2':
            nParaview += 3
            SeismoVLABfile.write("%d\n" % attributes['section'])
        elif  NAME == 'LIN3DFRAME2':
            nParaview += 3
            SeismoVLABfile.write("%d %s %s %d\n" %(attributes['section'], attributes['formulation'], attributes['rule'], attributes['np']))
        elif NAME == 'LIN3DSHELL4':
            nParaview += 5
            SeismoVLABfile.write("%d %s %d\n" %(attributes['section'], attributes['rule'], attributes['np']))
        elif NAME == 'LIN3DHEXA8':
            nParaview += 9
            SeismoVLABfile.write("%d %s %d\n" %(attributes['material'], attributes['rule'], attributes['np']))
        elif NAME == 'KIN3DHEXA8':
            nParaview += 9
            SeismoVLABfile.write("%d %s %d\n" %(attributes['material'], attributes['rule'], attributes['np']))
        elif NAME == 'PML3DHEXA8':
            nParaview += 9
            X0  = attributes['x0']
            DIR = attributes['npml']
            SeismoVLABfile.write("%d %d %1.5f %E %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %s %d\n" %(attributes['material'], attributes['n'], attributes['L'], attributes['R'], X0[0], X0[1], X0[2], DIR[0], DIR[1], DIR[2], attributes['rule'], attributes['np']))
        elif NAME == 'LIN3DHEXA20':
            nParaview += 21
            SeismoVLABfile.write("%d %s %d\n" %(attributes['material'], attributes['rule'], attributes['np']))
        elif  NAME == 'UNXBOUCWEN2DLINK':
            nParaview += 3
            SeismoVLABfile.write("%1.5f %1.5f %1.5f %1.5f %1.5f %E %E %E %1.5f %1.5f %d %d\n" %(attributes['alpha'], attributes['mu'], attributes['eta'], attributes['beta'], attributes['gamma'], attributes['tol'], attributes['Fy'], attributes['k0'], attributes['a1'], attributes['a2'], attributes['dim'], attributes['dir']))
        elif  NAME == 'UNXBOUCWEN3DLINK':
            nParaview += 3
            SeismoVLABfile.write("%1.5f %1.5f %1.5f %1.5f %1.5f %E %E %E %1.5f %1.5f %d %d\n" %(attributes['alpha'], attributes['mu'], attributes['eta'], attributes['beta'], attributes['gamma'], attributes['tol'], attributes['Fy'], attributes['k0'], attributes['a1'], attributes['a2'], attributes['dim'], attributes['dir']))
        elif  NAME == 'HDRBYAMAMOTO2DLINK':
            nParaview += 3
            SeismoVLABfile.write("%1.5f %1.5f %1.5f %d\n" %(attributes['De'], attributes['Di'], attributes['Hr'], attributes['dim']))
        elif  NAME == 'HDRBYAMAMOTO3DLINK':
            nParaview += 3
            SeismoVLABfile.write("%1.5f %1.5f %1.5f %d\n" %(attributes['De'], attributes['Di'], attributes['Hr'], attributes['dim']))
        elif NAME == 'EQLIN2DQUAD4':
            nParaview +=5
            SeismoVLABfile.write("%d %1.5f %s %d %s %1.5f %1.5f %1.5f \n" %(attributes['material'], attributes['th'], attributes['rule'], attributes['np'], attributes['type'], attributes['zref'], attributes['cf1'], attributes['cf2']))
        elif NAME == 'TIEQLIN2DQUAD4':
            nParaview +=5
            SeismoVLABfile.write("%d %1.5f %s %d %s %1.5f %1.5f %1.5f %1.10f \n" %(attributes['material'], attributes['th'], attributes['rule'], attributes['np'], attributes['type'], attributes['zref'], attributes['cf1'], attributes['cf2'], attributes['eref']))
        elif  NAME == 'NULL2DFRAME2':
            nParaview += 3
            SeismoVLABfile.write("\n")
        elif  NAME == 'NULL3DFRAME2':
            nParaview += 3
            SeismoVLABfile.write("\n")
        else:
            print('\x1B[31m ERROR \x1B[0m: In dict2svl() the ELEMENT=' + NAME + ' could not be recognized.')
            sys.exit(-1)

    Options['nparaview']  = nParaview
    Options['nfeatures'] += nParaview

    #[9] Damping
    for dTag in d['Dampings']:
        #Unpack dictionary for simpler writing
        Damping = d['Dampings'][dTag]
        name = Damping['name']
        attributes = Damping['attributes']

        SeismoVLABfile.write("DAMPING %s %s" %(dTag, name))
        if name == 'RAYLEIGH':
            SeismoVLABfile.write(" %E %E" %(attributes['am'],attributes['ak']))
        elif name == 'CAUGHEY':
            pass
        elif name == 'CAPPED':
            pass

        SeismoVLABfile.write(" %d" % len(attributes['list']))
        for k in attributes['list']:
            SeismoVLABfile.write(" %d" % k)
        SeismoVLABfile.write('\n')

    #[10] Load
    for lTag in d['Loads']:
        #Unpack dictionary for simpler writing
        Load = d['Loads'][lTag]
        NAME = Load['name']
        LIST = Load['attributes']['list']
            
        #Point load case
        if NAME == 'POINTLOAD':
            fName = Load['attributes']['name']

            SeismoVLABfile.write("%s %s %s " % (NAME, lTag, fName))
            #The Load Condition
            DIR = Load['attributes']['dir'] 
            if 'mag' in Load['attributes']: 
                SeismoVLABfile.write("%E %d" % (Load['attributes']['mag'], len(DIR)))
            elif 'file' in Load['attributes']:
                path = Load['attributes']['file']
                SeismoVLABfile.write("%s %d" % (path, len(DIR)))
            #The Load direction
            for k in DIR:
                SeismoVLABfile.write(" %E" % k)
            #Points that are loaded
            SeismoVLABfile.write(" %d" % len(LIST))
            for k in LIST:
                SeismoVLABfile.write(" %d" % k)
            SeismoVLABfile.write("\n")

        #Element Load case
        elif NAME == 'ELEMENTLOAD':
            fName = Load['attributes']['name']
            lType = Load['attributes']['type']

            SeismoVLABfile.write("%s %s %s " % (NAME, lTag, fName))
            #The Load Condition
            if lType == 'SURFACE':
                DIR = Load['attributes']['dir']
                #The Load Condition
                if 'mag' in Load['attributes']:
                    SeismoVLABfile.write("%s %E %d" % (lType, Load['attributes']['mag'], len(DIR)))
                    #The Load direction
                    for k in DIR:
                        SeismoVLABfile.write(" %1.5f" % k)
                    #Elements that are loaded
                    SeismoVLABfile.write(" %d" % len(LIST))
                    for k in LIST:
                        SeismoVLABfile.write(" %d" % k)
                    #Surfaces that are loaded
                    for k in LIST:
                        SeismoVLABfile.write(" %d" % Entities['Surfaces'][k]['face'])
                    SeismoVLABfile.write("\n")
                else:
                    print(' \x1B[31m ERROR \x1B[0m: In dict2svl() SURFACE (Dynamic) ELEMENTLOAD is not Implemented yet')
                    sys.exit(-1)
            elif lType == 'BODY':
                DIR = Load['attributes']['dir']
                #The Load Condition
                if 'mag' in Load['attributes']:
                    SeismoVLABfile.write("%s %E %d" % (lType, Load['attributes']['mag'], len(DIR)))
                elif 'file' in Load['attributes']:
                    path = Load['attributes']['file']
                    SeismoVLABfile.write("%s %s %d" % (lType, path, len(DIR)))
                #The Load direction
                for k in DIR:
                    SeismoVLABfile.write(" %E" % k)
                #Elements that are loaded
                SeismoVLABfile.write(" %d" % len(LIST))
                for k in LIST:
                    SeismoVLABfile.write(" %d" % k)
                SeismoVLABfile.write("\n")
            elif lType == 'PLANEWAVE':
                path = Load['attributes']['file']
                SeismoVLABfile.write("GENERALWAVE %s %d" % (path, len(LIST)))
                #Elements that are employed for DRM genral forces
                for k in LIST:
                    SeismoVLABfile.write(" %d" % k)
                SeismoVLABfile.write("\n")
            elif lType == 'GENERALWAVE':
                path = Load['attributes']['file']
                SeismoVLABfile.write("%s %s %d" % (lType, path, len(LIST)))
                #Elements that are employed for DRM general forces
                for k in LIST:
                    SeismoVLABfile.write(" %d" % k)
                SeismoVLABfile.write("\n")

        #Support motion case
        if NAME == 'SUPPORTMOTION':
            SeismoVLABfile.write("POINTLOAD %s %s %d" % (lTag, NAME, len(LIST)))

            #Points that have motion
            for k in LIST:
                SeismoVLABfile.write(" %d" % k)
            SeismoVLABfile.write("\n")

    #[11] Load Combinations
    for cTag in d['Combinations']:
        #Unpack dictionary for simpler writing
        Combination = d['Combinations'][cTag]
        name = Combination['name']
        attributes = Combination['attributes']

        SeismoVLABfile.write("COMBINATION %s %s" % (cTag, name))

        #Gets Load that are in this combination
        if attributes:
            load  = attributes['load']
            factors = attributes['factor']
            SeismoVLABfile.write(" %d" % len(load))

            for l, f in zip(load, factors):
                SeismoVLABfile.write(" %1.5f %d" % (f, l))         
            SeismoVLABfile.write("\n")
        else:
            SeismoVLABfile.write(" -1\n")

    #[12] Recorder
    for rTag in d['Recorders']:
        #Unpack dictionary for simpler writing
        Recorder = d['Recorders'][rTag]
        OUTFILE = Recorder['file']
        NSAMPLE = Recorder['nsamp']
        NDIGITS = Recorder['ndps']

        #Writes the Recorder's according to Point or Element
        if Recorder['name'] == 'NODE':
            nTags = Recorder['list']
            if nTags:
                SeismoVLABfile.write("RECORDER %s %s %d %d %s %d" % (Recorder['name'], OUTFILE, NSAMPLE, NDIGITS, Recorder['resp'], len(nTags)))
                for n in nTags:
                    SeismoVLABfile.write(" %d" % n)
                SeismoVLABfile.write("\n")
        elif Recorder['name'] == 'ELEMENT':
            eTags = Recorder['list']
            if eTags:
                SeismoVLABfile.write("RECORDER %s %s %d %d %s %d" % (Recorder['name'], OUTFILE, NSAMPLE, NDIGITS, Recorder['resp'], len(eTags)))
                for e in eTags:
                    SeismoVLABfile.write(" %d" % e)
                SeismoVLABfile.write("\n")
        elif Recorder['name'] == 'SECTION':
            eTags = Recorder['list']
            if eTags:
                pTags = Recorder['pos']
                SeismoVLABfile.write("RECORDER %s %s %d %d %s %d" % (Recorder['name'], OUTFILE, NSAMPLE, NDIGITS, Recorder['resp'], len(pTags)))
                for p in pTags:
                    SeismoVLABfile.write(" %E" % p)
                SeismoVLABfile.write(" %d" % len(eTags))
                for e in eTags:
                    SeismoVLABfile.write(" %d" % e)
                SeismoVLABfile.write("\n")
        elif Recorder['name'] == 'PARAVIEW':
            SeismoVLABfile.write("RECORDER %s %s %d %d %d\n" % (Recorder['name'], OUTFILE, NSAMPLE, NDIGITS, Options['nparaview']))

    #[13] Simulation
    for sTag in d['Simulations']:
        #Unpack the simulation
        Simulation = d['Simulations'][sTag]

        #Gets the analysis input parameters
        NAME = Simulation['attributes']['analysis']['name']
        NSTEPS = Simulation['attributes']['analysis']['nt']
        SeismoVLABfile.write("ANALYSIS %s %s %d " % (Simulation['combo'], NAME, NSTEPS))

        #Gets the algorithm input parameters
        NSTEPS  = Simulation['attributes']['algorithm']['nstep']
        CNVGTOL = Simulation['attributes']['algorithm']['cnvgtol']
        CNVGTEST = Simulation['attributes']['algorithm']['cnvgtest']
        SeismoVLABfile.write("ALGORITHM %s %E %d %d " % (Simulation['attributes']['algorithm']['name'], CNVGTOL, NSTEPS, CNVGTEST))
        
        #Gets the integrator input parameters
        TIMESTEP = Simulation['attributes']['integrator']['dt']
        MASSTOL  = Simulation['attributes']['integrator']['mtol']
        STIFFTOL = Simulation['attributes']['integrator']['ktol']
        FORCETOL = Simulation['attributes']['integrator']['ftol']
        SeismoVLABfile.write("INTEGRATOR %s %E %E %E %E " % (Simulation['attributes']['integrator']['name'], MASSTOL, STIFFTOL, FORCETOL, TIMESTEP))

        #Sets the solver input parameters
        NAME = Simulation['attributes']['solver']['name']
        if NAME == 'EIGEN':
            SeismoVLABfile.write("SOLVER %s %d\n" % (NAME, Simulation['attributes']['solver']['update']))
        elif NAME == 'MUMPS':
            SeismoVLABfile.write("SOLVER %s %d %d\n" % (NAME, Simulation['attributes']['solver']['option'], Simulation['attributes']['solver']['update']))
        elif NAME == 'PETSC':
            SeismoVLABfile.write("SOLVER %s %d %E %d %d\n" % (NAME, Simulation['attributes']['solver']['option'], Simulation['attributes']['solver']['tol'], Options['d_nz'][np], Options['o_nz'][np]))
    SeismoVLABfile.close()
