#!/usr/bin/python3

import os
import sys
import numpy as np
from Source import PlaneWave

def GenerateDRMFiles(User, Point, Element, Function, Load):
    """
    """
    for lTag in Load:
        if 'FUNCTION' in Load[lTag]:
            fTag = Load[lTag]['FUNCTION']
            if Function[fTag]['NAME'] == 'PLANEWAVE':
                #Gets the Domain Reduction Method Information.
                nodes, conditions, time, disp, vels, accel, dt, option = PlaneWave.ParseDRMFile(Function[fTag])

                #Computes the input displacement, velocities and accelerations.
                disp, vels, accel = PlaneWave.ComputeField(disp, vels, accel, dt, option)

                #Computes and Generates the Domain Reduction files for each node. 
                PlaneWave.Driver(User, Point, Function[fTag], time, disp, vels, accel, nodes, conditions)

            if Function[fTag]['NAME'] == 'GENERALWAVE':
                if 'UNDEFINED' in Function[fTag]:
                    filepath = Function[fTag]['PATH']
                    if Function[fTag]['NC'] == 6:
                        values = "0.000 0.000 0.000 0.000 0.000 0.000\n" * Function[fTag]['NT']
                    if Function[fTag]['NC'] == 9:
                        values = "0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000\n" * Function[fTag]['NT']

                    for k in Function[fTag]['UNDEFINED']:
                        filename = filepath.replace("$", str(k))

                        #Domain Reduction file format.
                        DRMfile = open(filename, "w+")
                        DRMfile.write("%d %d 0\n" % (Function[fTag]['NT'], Function[fTag]['NC']))
                        DRMfile.write("%s" % values)
                        DRMfile.close()

def CheckCorrectness(User, Point, Mass, Element, Material, Section, Constraint, Diaphragm, Support, Initial, Damping, Surface, Function, Load, Combination, Recorder, Simulation):
    """
    This function checks consistency between the parameters specified by the
    user and the parameters employed in SeismoSSI. 

    Parameters
    ----------
    User : dict 
           The dictionary containing all the user's information.

    Point : dict
           The dictionary containing all point information.

    Constraint : dict
           The dictionary containing all constraints information.

    Diaphragm : dict
           The dictionary containing all diaphragm information.

    Material : dict
           The dictionary containing all material information.

    Section : dict
           The dictionary containing all section information.

    Element : dict
           The dictionary containing all element information.

    Function : dict
           The dictionary containing all Loading Patterns.

    Load : dict
           The dictionary containing all Point/Element Loads.

    Combination : dict
           The dictionary containing how Loads are combined.

    Recorder : dict
           The dictionary containing responses to be stored.

    Simulation : dict
           The dictionary containing the simulations to be performed.

    Returns
    -------
    Checks constistency and modifies the 'ALL' fields in Recorder, and Load.
    """
    #Point information was not specified.
    if User['POINT'] == 'NO':
        print('\x1B[31m ERROR \x1B[0m: There are no POINT provided. Elements cannot be generated.')
        sys.exit(-1)

    #Restrain information was not specified.
    if User['RESTRAIN'] == 'NO':
        print('\x1B[31m ERROR \x1B[0m: There are no RESTRAIN provided. Model is ill-conditioned.')
        sys.exit(-1)

    #Material information was not specified.
    #if User['MATERIAL'] == 'NO':
        #print('\x1B[31m ERROR \x1B[0m: There are no MATERIAL provided. Elements cannot be generated.')
        #sys.exit(-1)

    #Element information was not specified.
    if User['ELEMENT'] == 'NO':
        print('\x1B[31m ERROR \x1B[0m: There are no ELEMENT provided. Mesh cannot be created.')
        sys.exit(-1)

    #Section information was not specified.
    if User['DIAPHRAGM'] == 'NO':
        print('\x1B[33m ALERT \x1B[0m: There is no DIAPHRAGM specified.')

    #Section information was not specified.
    if User['SECTION'] == 'NO':
        print('\x1B[33m ALERT \x1B[0m: There is no SECTION specified. Check there are no structural elements.')

    #Recorder information was not specified.
    if User['RECORDER'] == 'NO':
        print('\x1B[33m ALERT \x1B[0m: There is no RECORDER specified. No solution will be stored.')

    #Simulation information was not specified.
    if User['SIMULATION'] == 'NO':
        print('\x1B[33m ALERT \x1B[0m: There is no SIMULATION specified. No simulation will be run.')

    #Damping information was not specified.
    if User['DAMPING'] == 'NO':
        print('\x1B[33m ALERT \x1B[0m: There is no DAMPING\'s file specified. Free damping model will be assumed.')
        damping = {'NAME': 'FREE', 'LIST': list(Element.keys()) }
        Damping[1] = damping
    else:
        #Checks Consistency between Damping Model and Input Data
        for dTag in Damping:
            if isinstance(Damping[dTag]['LIST'], str):
                Damping[dTag]['LIST'] = list(Element.keys())

            if Damping[dTag]['NAME'].upper() == 'RAYLEIGH':
                if 'ALPHA' not in Damping[dTag]:
                    Damping[dTag]['ALPHA'] = 0.0
                if 'BETA' not in Damping[dTag]:
                    Damping[dTag]['BETA'] = 0.0
            elif Damping[dTag]['NAME'].upper() == 'CAUGHEY':
                print('\x1B[31m ERROR \x1B[0m: CAUGHEY DAMPING has not been implemented yet.')
                sys.exit(-1)
            elif Damping[dTag]['NAME'].upper() == 'CAPPED':
                print('\x1B[31m ERROR \x1B[0m: CAPPED DAMPING has not been implemented yet.')
                sys.exit(-1)

        #Finds elements without damping 
        eTags = list(Element.keys())
        for dTag in Damping:
            dlist = eTags
            eTags = list(set(dlist).difference(Damping[dTag]['LIST']))

        if eTags:
            dTags = max(list(Damping.keys())) + 1
            Damping[dTags] = {'NAME': 'Free', 'LIST': eTags}

    #Consistency in Link
    for eTag in Element:
        LinkName = Element[eTag]['NAME'].upper()
        if LinkName == 'UNXBOUCWEN2DLINK' or LinkName == 'UNXBOUCWEN3DLINK' or LinkName == 'HDRBYAMAMOTO2DLINK' or LinkName == 'HDRBYAMAMOTO3DLINK':
            conn = Element[eTag]['CONNECTIVITY']
            nDofi = Point[conn[0]]['NDOF']
            nDofj = Point[conn[1]]['NDOF']
            if nDofi != nDofj:
                print('\x1B[31m ERROR \x1B[0m: The number of degree-of-freedom in Point[' + str(conn[0]) + '] does not match Restrain[' + str(conn[1]) + '].')
                sys.exit(-1)            
            Element[eTag]['DIMENSION'] = nDofi

    #Checks consistency between SURFACE/ELEMENT and Assign face
    if Surface:
        condition = False
        for sTag in Surface:
            if Element[sTag]:
                name = Element[sTag]['NAME'].upper()
                surface    = Surface[sTag]['CONNECTIVITY']
                connection = np.array(Element[sTag]['CONNECTIVITY'])

                if  name == 'LIN2DTRUSS2'  or name == 'KIN2DTRUSS2' or name == 'LIN3DTRUSS2' or name == 'KIN3DTRUSS2' or name == 'LIN2DTRUSS3' or name == 'LIN3DTRUSS3':
                    Surface[sTag]['FACE'] = 1
                elif name == 'LIN2DFRAME2' or name == 'KIN2DFRAME2' or name == 'LIN3DFRAME2' or name == 'KIN2DFRAME2':
                    Surface[sTag]['FACE'] = 1
                elif name == 'LIN2DQUAD4'  or name == 'KIN2DQUAD4'  or name == 'LIN2DQUAD8'  or name == 'KIN2DQUAD8':
                    wTag =list(set(surface).intersection(connection[[0,1]]))
                    if len(wTag) == 2:
                        Surface[sTag]['FACE'] = 1
                    wTag = list(set(surface).intersection(connection[[1,2]]))
                    if len(wTag) == 2:
                        Surface[sTag]['FACE'] = 2
                    wTag = list(set(surface).intersection(connection[[2,3]]))
                    if len(wTag) == 2:
                        Surface[sTag]['FACE'] = 3
                    wTag = list(set(surface).intersection(connection[[3,0]]))
                    if len(wTag) == 2:
                        Surface[sTag]['FACE'] = 4
                elif name == 'LIN3DSHELL4':
                    wTag = list(set(surface).intersection(connection[[0,1,2,3]]))
                    if len(wTag) == 4:
                        Surface[sTag]['FACE'] = 5
                    else:
                        print('\x1B[31mERROR \x1B[0m: The ELEMENT=' + Element[k]['NAME'].upper() + ' only supports normal surface load. Connectivity of surface must be four.')
                        condition = True 
                elif name == 'LIN3DHEXA8' or name == 'KIN3DHEXA8' or name == 'LIN3DHEXA20':
                    wTag = list(set(surface).intersection(connection[[0,1,2,3]]))
                    if len(wTag) == 4:
                        Surface[sTag]['FACE'] = 1
                    wTag = list(set(surface).intersection(connection[[0,4,5,1]]))
                    if len(wTag) == 4:
                        Surface[sTag]['FACE'] = 2
                    wTag = list(set(surface).intersection(connection[[1,2,6,5]]))
                    if len(wTag) == 4:
                        Surface[sTag]['FACE'] = 3
                    wTag = list(set(surface).intersection(connection[[3,7,6,2]]))
                    if len(wTag) == 4:
                        Surface[sTag]['FACE'] = 4
                    wTag = list(set(surface).intersection(connection[[0,3,7,4]]))
                    if len(wTag) == 4:
                        Surface[sTag]['FACE'] = 5
                    wTag = list(set(surface).intersection(connection[[4,5,6,7]]))
                    if len(wTag) == 4:
                        Surface[sTag]['FACE'] = 6
                else:
                    print('\x1B[31mERROR \x1B[0m: The ELEMENT=' + Element[k]['NAME'].upper() + ' does not support surface load.')    
            else:
                print('\x1B[31m ERROR \x1B[0m: The SURFACE [' + str(sTag) + '] does not belong to an ELEMENT.')
                condition = True
        if condition:
            sys.exit(-1)

    #Checks consistency between POINTLOAD/ELEMENTLOAD and Point Degree-of-freedom 
    for lTag in Load:
        Name = Load[lTag]['NAME'].upper()
        if Name == 'POINTLOAD':
            nTag = Load[lTag]['LIST']
            fTag = Load[lTag]['FUNCTION']
            nDIR = len(Function[fTag]['DIRECTION'])
            for n in nTag:
                nDOF = Point[n]['NDOF']
                if nDOF != nDIR:
                    print('\x1B[31m ERROR \x1B[0m: The number of degree-of-freedom in LOAD[' + str(lTag) + '] for Point[' + str(n) + '] does not match Function[' + str(fTag) + '] Direction.')
                    sys.exit(-1)
        elif Name == 'ELEMENTLOAD':
            nTag = Load[lTag]['LIST']
            fTag = Load[lTag]['FUNCTION']
            if 'DIRECTION' in Load[lTag]:
                nDIR = len(Function[fTag]['DIRECTION'])
                if nDIR != User['DIMENSION']:
                    print('\x1B[31m ERROR \x1B[0m: The dimension of LOAD[' + str(lTag) + '] DIRECTION for ELEMENT LIST does not match MODEL DIMENSION.')
                    sys.exit(-1)

    #Check the GeneralWave Inputs are properly assigned
    for lTag in Load:
        if 'TYPE' in Load[lTag]:
            if Load[lTag]['TYPE'] == 'GENERALWAVE':
                eTag = Load[lTag]['LIST']
                fTag = Load[lTag]['FUNCTION']
                path = Function[fTag]['PATH']

                #Gets the DRM Nodes
                nTag = set()
                for k in eTag:
                    nlist = Element[k]['CONNECTIVITY']
                    for n in nlist:
                        nTag.add(n)

                #Checks if the file can be opened.
                nt = 0
                nc = 0
                nopen = 0
                undefine = set()
                for k in nTag:
                    filename = path.replace("$",str(k))

                    try:
                        with open(filename, "r") as fh:
                            nopen = k
                            fh.close()
                    except IOError as e:
                        undefine.add(k)

                #Open file and extract information
                filename = path.replace("$", str(nopen))
                with open(filename, "r") as fh:
                    lines = fh.readlines()

                    #Parse the first line.
                    line = list(filter(None, lines[0].strip().split()))
                    nt, nc = int(line[0]), int(line[1])
                    fh.close()

                #Some DRM nodes are undefined, they are set to be zero.
                if undefine:
                    Function[fTag]['NC'] = nc
                    Function[fTag]['NT'] = nt
                    Function[fTag]['UNDEFINED'] = undefine
                    print('\x1B[33m ALERT \x1B[0m: In GENERALWAVE there are DRM nodes which have no field specified.')

    #Prepare Point/Element Indeces for 'ALL' case
    for lTag in Load:
        if isinstance(Load[lTag]['LIST'], str):
            if Load[lTag]['NAME'] == 'POINTLOAD':
                Load[lTag]['LIST'] = list(Point.keys())
            elif Load[lTag]['NAME'] == 'ELEMENTLOAD':
                Load[lTag]['LIST'] = list(Element.keys())


    #Prepare Point/Element Indeces for 'ALL' case
    for lTag in Load:
        if isinstance(Load[lTag]['LIST'], str):
            if Load[lTag]['NAME'] == 'POINTLOAD':
                Load[lTag]['LIST'] = list(Point.keys())
            elif Load[lTag]['NAME'] == 'ELEMENTLOAD':
                Load[lTag]['LIST'] = list(Element.keys())

    #Signal Time shift for PLANEWAVE function
    for lTag in Load:
        if 'FUNCTION' in Load[lTag]:
            fTag = Load[lTag]['FUNCTION']
            if Function[fTag]['NAME'] == 'PLANEWAVE':
                with open(Function[fTag]['PATH'], "r") as fileHandler:
                    #Gets the number of DRM nodes
                    lines  = fileHandler.readlines()
                    line   = list(filter(None, lines[0].strip().split())) 
                    nDims  = User['DIMENSION']
                    nNodes = int(line[2])
                    nc     = 6*(nDims == 2) + 9*(nDims == 3)  
                    nt     = int(line[0])

                    #Gets ALL DRM nodes
                    eTag = Load[lTag]['LIST']
                    allTag = set()
                    for k in eTag:
                        nlist = Element[k]['CONNECTIVITY']
                        for n in nlist:
                            allTag.add(n)

                    #Gets the DRM node Tags
                    nTags = np.empty(shape=(nNodes,), dtype=int)

                    m = 0
                    for k in range(1, nNodes+1):
                        line = list(filter(None, lines[k].strip().split()))
                        nTags[m] = int(line[0])
                        m += 1

                    undefined = allTag.difference(nTags)
                    if undefined:
                        Function[fTag]['NC'] = nc
                        Function[fTag]['NT'] = nt
                        Function[fTag]['UNDEFINED'] = undefined
                        print('\x1B[33m ALERT \x1B[0m: In PLANEWAVE there are DRM nodes which have no field specified.')
                        
                    #Gets the the most distant coordinate
                    xmin = np.empty(shape=(nNodes,nDims))
                    for k in range(nNodes):
                        xmin[k] = Point[nTags[k]]['COORDINATES']

                    Function[fTag]['XMIN'] = xmin.min(axis=0)
                    fileHandler.close()

    #Checks consistency 'ALL' case and 'Paraview'
    for rTag in Recorder:
        #Prepare Point/Element Indeces for 'ALL' case
        if 'LIST' in Recorder[rTag]:
            if isinstance(Recorder[rTag]['LIST'], str):
                if Recorder[rTag]['NAME'].upper() == 'NODE':
                    Recorder[rTag]['LIST'] = list(Point.keys())
                elif Recorder[rTag]['NAME'].upper() == 'ELEMENT':
                    Recorder[rTag]['LIST'] = list(Element.keys())
                elif Recorder[rTag]['NAME'].upper() == 'SECTION':
                    Recorder[rTag]['LIST'] = list(Element.keys())

        if 'RESPONSE' in Recorder[rTag]:
            if Recorder[rTag]['RESPONSE'] == 'REACTION':
                newTags = []
                oldTags = Recorder[rTag]['LIST'];
                for n in oldTags:
                    if Point[n]['FIXED']:
                        newTags.append(n)
                Recorder[rTag]['LIST'] = newTags;

        if Recorder[rTag]['NAME'].upper() == 'PARAVIEW':
            dirName = User['FOLDER'] + 'Paraview/'
            if not os.path.exists(dirName):
                os.mkdir(dirName)

    #Check consistency between Serial and Parallel Analysis
    for sTag in Simulation:
        if Simulation[sTag]['SOLVER']['NAME'].upper() == 'EIGEN' and User['NPART'] != 1:
            Simulation[sTag]['SOLVER']['NAME'] = 'MUMPS'
            print('\x1B[33m ALERT \x1B[0m: A partition number greater than one requires a parallel solver.')

def GetSimulationInformation(fileHandler, User, Simulation):
    """
    This function creates the Simulation which handles how and what solution
    for Points/Elements will be stored. 

    Parameters
    ----------
    User : dict 
           The dictionary containing all the user's relevant information.

    Returns
    -------
    Simulation : dict
           The dictionary containing the simulations to be performed.
    """
    #Parse the partition input.
    while True:
        line = list(filter(None, fileHandler.readline().strip().split(',')))
        if line:
            #Checks end of Model's scope.
            if line[0].upper() == '*END':
                break

            #Auxiliar parsing dictionaries.
            analysis   = {'NAME': '', 'NSTEP': 1}
            algorithm  = {'NAME': '', 'NSTEP': 1, 'CNVGTOL': 1E-3, 'CNVGTEST': 'UNBALANCEFORCE'}
            integrator = {'NAME': '', 'TIMESTEP': 0.0, 'MTOL': 1E-12, 'KTOL': 1E-12 , 'FTOL': 1E-12 }
            solver     = {'NAME': '', 'TYPE': '', 'OPTION': ''}
            
            #Parse Simulation input's flag.
            sTag = int(line[0])
            for item in line:
                token = list(filter(None, item.strip().split()))

                if token[0].upper() == '-ANALYSIS':
                    analysis['NAME'] = token[1].upper()
                    for k in range(2, len(token)):
                        if token[k].upper() == '-NSTEP':
                            analysis['NSTEP'] = int(token[k+1])
                elif token[0].upper() == '-ALGORITHM':
                    algorithm['NAME'] = token[1].upper()
                    for k in range(2, len(token)):
                        if token[k].upper() == '-NSTEP':
                            algorithm['NSTEP'] = int(token[k+1])
                        elif token[k].upper() == '-CNVGTOL':
                            algorithm['CNVGTOL'] = float(token[k+1])
                        elif token[k].upper() == '-CNVGTEST':
                            algorithm['CNVGTEST'] = token[k+1].upper()
                elif token[0].upper() == '-INTEGRATOR':
                    integrator['NAME'] = token[1].upper()
                    for k in range(2, len(token)):
                        if token[k].upper() == '-DT':
                            integrator['TIMESTEP'] = float(token[k+1])
                        elif token[k].upper() == '-MTOL':
                            integrator['MTOL'] = float(token[k+1])
                        elif token[k].upper() == '-KTOL':
                            integrator['KTOL'] = float(token[k+1])
                        elif token[k].upper() == '-FTOL':
                            integrator['FTOL'] = float(token[k+1])
                elif token[0].upper() == '-SOLVER':
                    solver['NAME'] = token[1].upper()
                    for k in range(2, len(token)):
                        if token[k].upper() == '-TYPE':
                            solver['TYPE'] = token[k+1].upper()
                        elif token[k].upper() == '-UPDATE':
                            solver['UPDATE'] = token[k+1].upper()
                        elif token[k].upper() == '-TOL':
                            solver['TOLERANCE'] = float(token[k+1])
                        if token[k].upper() == '-OPTION':
                            if token[k+1].upper() == 'SPD':
                                solver['OPTION'] = 0
                            elif token[k+1].upper() == 'SYM':
                                solver['OPTION'] = 1
                            elif token[k+1].upper() == 'USYM':
                                solver['OPTION'] = 2
                            elif token[k+1].upper() == 'KSPCG':
                                solver['OPTION'] = 0
                                User['PETSC'] = 'YES'
                            elif token[k+1].upper() == 'KSPBCGS':
                                solver['OPTION'] = 1
                                User['PETSC'] = 'YES'
                            elif token[k+1].upper() == 'KSPCGS':
                                solver['OPTION'] = 2
                                User['PETSC'] = 'YES'
                            elif token[k+1].upper() == 'KSPBICG':
                                solver['OPTION'] = 3
                                User['PETSC'] = 'YES'

            #The Simulation dictionary.
            Simulation[sTag] = {'ANALYSIS': analysis, 'ALGORITHM': algorithm, 'INTEGRATOR': integrator, 'SOLVER': solver}

    #Simulation file was provided
    User['SIMULATION'] = 'YES'

def GetRecorderInformation(fileHandler, User, Point, Element, Recorder):
    """
    This function creates the Recorder which handles how and what solution
    for Points/Elements will be stored. 

    Parameters
    ----------
    User : dict 
           The dictionary containing all the user's relevant information.

    Returns
    -------
    Recorder : dict
           The dictionary containing response to be stored.
    """
    #Parse the partition input.
    while True:
        line = list(filter(None, fileHandler.readline().strip().split(',')))
        if line:
            #Checks end of Model's scope.
            if line[0].upper() == '*END':
                break

            #Parse Recorder input's flag.
            recorder = {'NSAMPLE': 1, 'PRECISION': 8}
            rTag = int(line[0])
            for item in line:
                token = list(filter(None, item.strip().split()))
                if token[0].upper() == '-NAME':
                    recorder['NAME'] = token[1].upper()
                elif token[0].upper() == '-FILE':
                    recorder['FILE'] = token[1]
                elif token[0].upper() == '-NDGTS':
                    recorder['PRECISION'] = int(token[1])
                elif token[0].upper() == '-NSAMP':
                    recorder['NSAMPLE'] = int(token[1])
                elif token[0].upper() == '-FILE':
                    recorder['FILE'] = token[1]
                elif token[0].upper() == '-RESP':
                    recorder['RESPONSE'] = token[1]
                elif token[0].upper() == '-POS':
                    Tags = []
                    for k in range(1, len(token)):
                        Tags.append(float(token[k]))
                    recorder['POSITION'] = Tags
                elif token[0].upper() == '-LIST':
                    if token[1].upper() == 'ALL':
                        recorder['LIST'] = 'ALL'
                    else:
                        Tags = []
                        for k in range(1, len(token)):
                            Tags.append(int(token[k]))
                        recorder['LIST'] = Tags

                Recorder[rTag] = recorder

    #Recorder file was provided
    User['RECORDER'] = 'YES'

def GetCombinationInformation(fileHandler, User, Combination):
    """
    This function creates a Combination of Load to be applied on the finite 
    element mesh. Combination factors specifies how these load are combined

    Parameters
    ----------
    User : dict 
           The dictionary containing all the user's relevant information.

    Returns
    -------
    Combination : dict
           The dictionary containing how Loads are combined.
    """
    #Parse the partition input.
    while True:
        line = list(filter(None, fileHandler.readline().strip().split(',')))
        if line:
            #Checks end of Model's scope.
            if line[0].upper() == '*END':
                break

            #Parse Combination input's flag.
            combination = {}
            cTag = int(line[0])
            for item in line:
                token = list(filter(None, item.strip().split()))
                if token[0].upper() == '-NAME':
                    combination['NAME'] = token[1]
                elif token[0].upper() == '-LOAD':
                    load = []
                    for k in range(1, len(token)):
                        load.append(int(token[k]))
                    combination['LOAD'] = load
                elif token[0].upper() == '-FACTOR':
                    factors = []
                    for k in range(1, len(token)):
                        factors.append(float(token[k]))
                    combination['FACTOR'] = factors

                Combination[cTag] = combination

    #Combination file was provided
    User['COMBINATION'] = 'YES'

def GetLoadInformation(fileHandler, User, Point, Element, Function, Load):
    """
    This function creates a Load dictionary that contains the Points or Nodes 
    for which this load is going to act on.

    Parameters
    ----------
    User : dict 
           The dictionary containing all the user's relevant information.

    Returns
    -------
    Load : dict
           The dictionary containing all Point/Element Loads.
    """
    #Parse the partition input.
    while True:
        line = list(filter(None, fileHandler.readline().strip().split(',')))
        if line:
            #Checks end of Model's scope.
            if line[0].upper() == '*END':
                break

            #Parse Load input's flag.
            load = {}
            lTag = int(line[0])
            for item in line:
                token = list(filter(None, item.strip().split()))
                if token[0].upper() == '-NAME':
                    load['NAME'] = token[1].upper()
                elif token[0].upper() == '-TYPE':
                    load['TYPE'] = token[1].upper()
                elif token[0].upper() == '-FUN':
                    load['FUNCTION'] = int(token[1])
                elif token[0].upper() == '-LIST':
                    if token[1].upper() == 'ALL':
                        load['LIST'] = 'ALL'
                    else:
                        Tags = []
                        for k in range(1, len(token)):
                            Tags.append(int(token[k]))
                        load['LIST'] = Tags

                Load[lTag] = load

    #Load file was provided
    User['LOAD'] = 'YES'

def GetFunctionInformation(fileHandler, User, Function):
    """
    This routine creates a Loading Function Pattern dictionary that contains 
    the type of load to be assigned to Points or Elements.

    Parameters
    ----------
    User : dict 
           The dictionary containing all the user's relevant information.

    Returns
    -------
    Function : dict
           The dictionary containing all Loading Patterns.
    """
    #Parse the function input.
    while True:
        line = list(filter(None, fileHandler.readline().strip().split(',')))
        if line:
            #Checks end of Model's scope.
            if line[0].upper() == '*END':
                break

            #Parse Function input's flag.
            function = {}
            fTag = int(line[0])
            for item in line:
                token = list(filter(None, item.strip().split()))
                if token[0].upper() == '-NAME':
                    function['NAME'] = token[1].upper()
                elif token[0].upper() == '-TYPE':
                    function['TYPE'] = token[1].upper()
                elif token[0].upper() == '-PATH':
                    path = ""
                    if len(token) > 2:
                        path = token[1] + " " + token[2]
                        for k in range(3, len(token)):
                            path = path + " " + token[k]
                    else:
                        path = token[1]
                    function['PATH'] = path
                elif token[0].upper() == '-MAG':
                    function['MAGNITUDE'] = float(token[1])
                elif token[0].upper() == '-THETA':
                    function['THETA'] = float(token[1])
                elif token[0].upper() == '-PHI':
                    function['PHI'] = float(token[1])
                elif token[0].upper() == '-VS':
                    function['VS'] = float(token[1])
                elif token[0].upper() == '-NU':
                    function['NU'] = float(token[1])
                elif token[0].upper() == '-X0':
                    x0 = []
                    for k in range(1, len(token)):
                        x0.append(float(token[k]))
                    function['X0'] = x0
                elif token[0].upper() == '-DIR':
                    direction = []
                    for k in range(1, len(token)):
                        direction.append(float(token[k]))
                    function['DIRECTION'] = direction
                Function[fTag] = function

    #Constructs the file's global path.
    for fTag in Function:
        if 'PATH' in Function[fTag]:
            Function[fTag]['PATH'] = User['FOLDER'] + Function[fTag]['PATH']

            #Checks if the file can be opened.
            try:
                with open(Function[fTag]['PATH'], "r") as fh:
                    fh.close()
            except IOError as e:
                if Function[fTag]['NAME'] != 'GENERALWAVE':
                    print('\x1B[31m ERROR \x1B[0m: In *FUNCTION the FILE=' + Function[fTag]['PATH'] + ' could not be opened.')
                    print(' Check if both the file path and file name provided are correct.')
                    sys.exit(-1)

    #Function file was provided
    User['FUNCTION'] = 'YES'

def GetDampingInformation(fileHandler, User, Damping):
    """
    This function creates a Damping dictionary that will be applied to some
    elements of the mesh. The Damping models currently supported are 'Free', 
    and 'Rayleigh'.

    Parameters
    ----------
    User : dict 
           The dictionary containing all the user's relevant information.
    Element : dict
           The dictionary containing all element information.

    Returns
    -------
    Damping : dict
           The dictionary containing all damping models.
    """
    #Parse the partition input.
    while True:
        line = list(filter(None, fileHandler.readline().strip().split(',')))
        if line:
            #Checks end of Model's scope.
            if line[0].upper() == '*END':
                break

            #Parse Damping input's flag.
            damping = {}
            dTag = int(line[0])
            for item in line:
                token = list(filter(None, item.strip().split()))
                if token[0].upper() == '-NAME':
                    damping['NAME'] = token[1]
                if token[0].upper() == '-A0':
                    damping['ALPHA'] = float(token[1])
                if token[0].upper() == '-A1':
                    damping['BETA'] = float(token[1])
                elif token[0].upper() == '-LIST':
                    if token[1].upper() == 'ALL':
                        damping['LIST'] = 'ALL'
                    else:
                        elem = []
                        for k in range(1, len(token)):
                            elem.append(int(token[k]))
                        damping['LIST'] = elem

            Damping[dTag] = damping

    #Damping file was provided
    User['DAMPING'] = 'YES'

def GetSurfaceInformation(fileHandler, User, Surface):
    """
    """
    while True:
        line = list(filter(None, fileHandler.readline().strip().split(',')))
        if line:
            #Checks end of Model's scope.
            if line[0].upper() == '*END':
                break

            #Parse Function input's flag.
            surface = {}
            sTag = int(line[0])
            for item in line:
                token = list(filter(None, item.strip().split()))
                if token[0].upper() == '-CONN':
                    connectivity = []
                    for k in range(1, len(token)):
                        connectivity.append(int(token[k]))
                    surface['CONNECTIVITY'] = connectivity
                elif token[0].upper() == '-MAG':
                    magnitudes = []
                    for k in range(1, len(token)):
                        magnitudes.append(float(token[k]))
                    surface['MAGNITUDE'] = magnitudes
                Surface[sTag] = surface

    #Element file was provided
    User['SURFACE'] = 'YES'

def GetElementInformation(fileHandler, User, Element):
    """
    This function creates a Elements dictionary that contains the element 
    identifier, name of the element, connectivity array, material, section, 
    and other relevant element information.

    Parameters
    ----------
    User : dict 
           The dictionary containing all the user's relevant information.

    Returns
    -------
    Element : dict
           The dictionary containing all element information.
    """
    #Parse the partition input.
    while True:
        line = list(filter(None, fileHandler.readline().strip().split(',')))
        if line:
            #Checks end of Model's scope.
            if line[0].upper() == '*END':
                break

            #Parse Function input's flag.
            element = {}
            eTag = int(line[0])
            for item in line:
                token = list(filter(None, item.strip().split()))
                if token[0].upper() == '-NAME':
                    element['NAME'] = token[1]
                elif token[0].upper() == '-CONN':
                    connectivity = []
                    for k in range(1, len(token)):
                        connectivity.append(int(token[k]))
                    element['CONNECTIVITY'] = connectivity
                elif token[0].upper() == '-MAT':
                    element['MATERIAL'] = int(token[1])
                elif token[0].upper() == '-SEC':
                    element['SECTION'] = int(token[1])
                elif token[0].upper() == '-NP':
                    element['NPOINTS'] = int(token[1])
                elif token[0].upper() == '-RULE':
                    element['QUADRATURE'] = token[1]
                elif token[0].upper() == '-TH':
                    element['THICKNESS'] = float(token[1])
                elif token[0].upper() == '-A':
                    element['AREA'] = float(token[1])
                elif token[0].upper() == '-L':
                    element['L'] = float(token[1])
                elif token[0].upper() == '-R':
                    element['R'] = float(token[1])
                elif token[0].upper() == '-N':
                    element['DEGREE'] = float(token[1])
                elif token[0].upper() == '-DIR':
                    element['DIRECTION'] = int(token[1]) - 1
                elif token[0].upper() == '-NPML':
                    npml = []
                    for k in range(1, len(token)):
                        npml.append(float(token[k]))
                    element['NORMAL'] = npml
                elif token[0].upper() == '-X0':
                    x0 = []
                    for k in range(1, len(token)):
                        x0.append(float(token[k]))
                    element['X0'] = x0
                elif token[0].upper() == '-TYPE':
                    element['TYPE'] = token[1].upper()
                elif token[0].upper() == '-CF1':
                    element['CF1'] = float(token[1])
                elif token[0].upper() == '-CF2':
                    element['CF2'] = float(token[1])
                elif token[0].upper() == '-ZREF':
                    element['ZREF'] = float(token[1])
                elif token[0].upper() == '-EREF':
                    element['EREF'] = float(token[1])
                elif token[0].upper() == '-FORM':
                    element['FORMULATION'] = token[1]
                elif token[0].upper() == '-DE':
                    element['DE'] = float(token[1])
                elif token[0].upper() == '-DI':
                    element['DI'] = float(token[1])
                elif token[0].upper() == '-HR':
                    element['HR'] = float(token[1])
                elif token[0].upper() == '-ALPHA':
                    element['ALPHA'] = float(token[1])
                elif token[0].upper() == '-MU':
                    element['MU'] = float(token[1])
                elif token[0].upper() == '-ETA':
                    element['ETA'] = float(token[1])
                elif token[0].upper() == '-BETA':
                    element['BETA'] = float(token[1])
                elif token[0].upper() == '-GAMMA':
                    element['GAMMA'] = float(token[1])
                elif token[0].upper() == '-TOL':
                    element['TOL'] = float(token[1])
                elif token[0].upper() == '-FY':
                    element['FY'] = float(token[1])
                elif token[0].upper() == '-K0':
                    element['K0'] = float(token[1])
                elif token[0].upper() == '-ALPHA1':
                    element['ALPHA1'] = float(token[1])
                elif token[0].upper() == '-ALPHA2':
                    element['ALPHA2'] = float(token[1])

            Element[eTag] = element

    #Element file was provided
    User['ELEMENT'] = 'YES'

def GetSectionInformation(fileHandler, User, Section):
    """
    This function creates a Section dictionary.

    Parameters
    ----------
    User : dict 
           The dictionary containing all the user's relevant information.

    Returns
    -------
    Section : dict
           The dictionary containing all section information.
    """
    #Parse the partition input.
    while True:
        line = list(filter(None, fileHandler.readline().strip().split(',')))
        if line:
            #Checks end of Model's scope.
            if line[0] == '*END':
                break

            #Parse Restrain input's flag.
            section = {}
            sTag = int(line[0])
            for item in line:
                token = list(filter(None, item.strip().split()))
                if token[0].upper() == '-NAME':
                    section['NAME'] = token[1]
                elif token[0].upper() == '-TYPE':
                    section['TYPE'] = token[1]
                elif token[0].upper() == '-PROP':
                    section['PROPERTY'] = token[1]
                elif token[0].upper() == '-MAT':
                    section['MATERIAL'] = int(token[1])
                elif token[0].upper() == '-H':
                    section['h'] = float(token[1])
                elif token[0].upper() == '-B':
                    section['b'] = float(token[1])
                elif token[0].upper() == '-TW':
                    section['tw'] = float(token[1])
                elif token[0].upper() == '-TF':
                    section['tf'] = float(token[1])
                elif token[0].upper() == '-RE':
                    section['re'] = float(token[1])
                elif token[0].upper() == '-RI':
                    section['ri'] = float(token[1])
                elif token[0].upper() == '-R':
                    section['r'] = float(token[1])
                elif token[0].upper() == '-T':
                    section['t'] = float(token[1])
                elif token[0].upper() == '-THETA':
                    section['theta'] = float(token[1])
                elif token[0].upper() == '-IP':
                    section['ip'] = int(token[1])

            Section[sTag] = section

    #Section file was provided
    User['SECTION'] = 'YES'

def GetMaterialInformation(fileHandler, User, Material):
    """
    This function creates a Material dictionary.

    Parameters
    ----------
    User : dict 
           The dictionary containing all the user's relevant information.

    Returns
    -------
    Material : dict
           The dictionary containing all material information.
    """
    #Parse the partition input.
    while True:
        line = list(filter(None, fileHandler.readline().strip().split(',')))
        if line:
            #Checks end of Model's scope.
            if line[0] == '*END':
                break

            #Parse Restrain input's flag.
            material = {'NAME': '', 'RHO': 0.000}
            mTag = int(line[0])
            for item in line:
                token = list(filter(None, item.strip().split()))
                if token[0].upper() == '-NAME':
                    material['NAME'] = token[1]
                elif token[0].upper() == '-E':
                    material['E'] = float(token[1])
                elif token[0].upper() == '-NU':
                    material['NU'] = float(token[1])
                elif token[0].upper() == '-K1':
                    material['K1'] = float(token[1])
                elif token[0].upper() == '-K2':
                    material['K2'] = float(token[1])
                elif token[0].upper() == '-K3':
                    material['K3'] = float(token[1])
                elif token[0].upper() == '-K':
                    material['K'] = float(token[1])
                elif token[0].upper() == '-G':
                    material['G'] = float(token[1])
                elif token[0].upper() == '-ETA':
                    material['ETA'] = float(token[1])
                elif token[0].upper() == '-RHO':
                    material['RHO'] = float(token[1])
                elif token[0] == '-H':
                    material['H'] = float(token[1])
                elif token[0].upper() == '-SY':
                    material['SY'] = float(token[1])
                elif token[0].upper() == '-BETA':
                    material['BETA'] = float(token[1])
                elif token[0].upper() == '-SU':
                    material['SU'] = float(token[1])
                elif token[0].upper() == '-H0':
                    material['H0'] = float(token[1])
                elif token[0] == '-h':
                    material['h'] = float(token[1])
                elif token[0] == '-m':
                    material['m'] = float(token[1])

            Material[mTag] = material

    #Material file was provided.
    User['MATERIAL'] = 'YES'

def ComputeConstraints(User, Point, Constraint, Diaphragm):
    """
    This function generates the constraints for the complete model.

    Parameters
    ----------
    User : dict 
           The dictionary containing all the user's relevant information.

    Output
    -------
    Adds Constraint and Point from Diaphragm.
    """
    #Corrects Free degree-of-freedom numbering for Constrained Points
    for cTag in Constraint:
        nTag = Constraint[cTag]['SLAVENODE']
        sDir = Constraint[cTag]['SLAVEDOF' ]
        Point[nTag]['FREE'][sDir] = cTag

    #Corrects Free degree-of-freedom number for Diaphragm Points
    gTag = User['NPOINTS']
    cTag = User['NCONSTRAINT']
    for dTag in Diaphragm:
        #Computes the center of Diaphragm.
        count = 0.0
        nTag  = Diaphragm[dTag]['NODE']
        DiaphCenter = np.zeros(User['DIMENSION'], dtype='float')
        for node in nTag:
            DiaphCenter += Point[node]['COORDINATES']
            count += 1.00

        DiaphCenter = DiaphCenter/count

        #Computes Combinational Factors for Constraints.
        if Diaphragm[dTag]['AXIS'] == 'Z':
            for node in nTag:
                free        = Point[node]['FREE']
                Coordinates = Point[node]['COORDINATES']
                CombXCoeff  = DiaphCenter[1] - Coordinates[1] 
                CombYCoeff  = Coordinates[0] - DiaphCenter[0]

                #Constraint in direction X
                cTag -= 1
                free[0] = cTag
                constraint = {'NAME': 'DIAPHRAGM', 'SLAVENODE': node, 'SLAVEDOF': 0, 'MASTERNODE': [dTag, dTag], 'MASTERDOF': [0, 2], 'FACTORS': [1.00, CombXCoeff]}
                Constraint[cTag] = constraint

                #Constraint in direction Y
                cTag -= 1
                free[1] = cTag
                constraint = {'NAME': 'DIAPHRAGM', 'SLAVENODE': node, 'SLAVEDOF': 1, 'MASTERNODE': [dTag, dTag], 'MASTERDOF': [1, 2], 'FACTORS': [1.00, CombYCoeff]}
                Constraint[cTag] = constraint

                #Constraint in 3D for Rotation about Z-axis.
                if User['DIMENSION'] == 3 and Point[node]['NDOF'] == 6:
                    cTag -= 1
                    free[5] = cTag
                    constraint = {'NAME': 'DIAPHRAGM', 'SLAVENODE': node, 'SLAVEDOF': 5, 'MASTERNODE': [dTag], 'MASTERDOF': [2], 'FACTORS': [1.00]}
                    Constraint[cTag] = constraint

                #Constraint in 2D for Rotation about Z-axis.
                elif User['DIMENSION'] == 2 and Point[node]['NDOF'] == 3:
                    cTag -= 1
                    free[2] = cTag
                    constraint = {'NAME': 'DIAPHRAGM', 'SLAVENODE': node, 'SLAVEDOF': 2, 'MASTERNODE': [dTag], 'MASTERDOF': [2], 'FACTORS': [1.00]}
                    Constraint[cTag] = constraint

                #Corrects Free degree-of-freedom number.
                Point[node]['FREE'] = free

        elif Diaphragm[dTag]['NODE'] == 'Y':
            for node in nTag:
                Coordinates = Point[node]['COORDINATES']
                CombXCoeff[k] = DiaphCenter[2] - Coordinates[2] 
                CombYCoeff[k] = Coordinates[0] - DiaphCenter[0]
                #TODO: 

        elif Diaphragm[dTag]['NODE'] == 'X':
            for node in nTag:
                Coordinates = Point[node]['COORDINATES']
                CombXCoeff[k] = DiaphCenter[2] - Coordinates[2] 
                CombYCoeff[k] = Coordinates[1] - DiaphCenter[1]
                #TODO: 

        #Creates Diaphragm Nodes
        Node = {'TAG': gTag, 'NDOF': 3, 'COORDINATES': DiaphCenter, 'FREE': np.zeros(3, dtype='int'), 'TOTAL': np.zeros(3, dtype='int')}
        Point[dTag] = Node
        gTag += 1

    #Updated number of Points in model
    User['NPOINTS'] = len(Point)

def GetDiaphragmInformation(fileHandler, User, Diaphragm):
    """
    This function creates a Diaphragm dictionary that contains the nodes 
    that belongs to it, and the axis of action.

    Parameters
    ----------
    User : dict 
           The dictionary containing all the user's relevant information.

    Returns
    -------
    Diaphragm : dict
           The dictionary containing all diaphragm information.
    """
    #Parse the partition input.
    while True:
        line = list(filter(None, fileHandler.readline().strip().split(',')))
        if line:
            #Checks end of Model's scope.
            if line[0].upper() == '*END':
                break

            #Parse Constraint input's flag.
            diaphragm = {}
            dTag = int(line[0])
            for item in line:
                token = list(filter(None, item.strip().split()))
                if token[0].upper() == '-AXIS':
                    diaphragm['AXIS'] = token[1].upper()
                elif token[0].upper() == '-LIST':
                    node = []
                    for k in range(1, len(token)):
                        node.append(int(token[k]))
                    diaphragm['NODE'] = node

            Diaphragm[dTag] = diaphragm

    #Diaphragm file was provided.
    User['DIAPHRAGM'] = 'YES'

def GetConstraintInformation(fileHandler, User, Constraint):
    """
    This function creates a Constraint dictionary.

    Parameters
    ----------
    User : dict 
           The dictionary containing all the user's relevant information.

    Returns
    -------
    Constraint : dict
           The dictionary containing all constraints information.
    """
    #Empty Constraint Dictionary
    gTag = -1

    #Parse the Constraint input.
    while True:
        line = list(filter(None, fileHandler.readline().strip().split(',')))
        if line:
            #Checks end of Model's scope.
            if line[0].upper() == '*END':
                break

            #Parse Constraint input's flag.
            gTag -= 1
            constraint = {}
            cTag = int(line[0])
            for item in line:
                token = list(filter(None, item.strip().split()))
                if token[0].upper() == '-NAME':
                    constraint['NAME'] = token[1].upper()
                elif token[0].upper() == '-SNODE':
                    constraint['SLAVENODE'] = int(token[1])
                elif token[0].upper() == '-SDOF':
                    constraint['SLAVEDOF'] = int(token[1]) - 1
                elif token[0].upper() == '-MNODE':
                    nodes  = []
                    for k in range(1, len(token)):
                        nodes.append(int(token[k]))
                    constraint['MASTERNODE'] = nodes
                elif token[0].upper() == '-MDOF':
                    dofs  = []
                    for k in range(1, len(token)):
                        dofs.append(int(token[k]) - 1)
                    constraint['MASTERDOF'] = dofs
                elif token[0].upper() == '-FACTOR':
                    factors  = []
                    for k in range(1, len(token)):
                        factors.append(float(token[k]))
                    constraint['FACTORS'] = factors

            Constraint[gTag] = constraint

    #Provides factor coefficient to EQUAL Constraint
    for cTag in Constraint:
        if Constraint[cTag]['NAME'] == 'EQUAL':
            Constraint[cTag]['FACTORS'] = [1.000]

    #Total number of imposed constraint 
    User['NCONSTRAINT'] = gTag

    #Constraint file was provided.
    User['CONSTRAINT'] = 'YES'

def GetRestrainInformation(fileHandler, User, Point):
    """
    This function creates a Restrain dictionary that contains the point identifier, 
    and the degree of freedom to be restrained.

    Parameters
    ----------
    User : dict 
           The dictionary containing all the user's relevant information.

    Output
    -------
    Modifies the Point['FREE'] field
    """
    #Parse the partition input.
    while True:
        line = list(filter(None, fileHandler.readline().strip().split(',')))
        if line:
            #Checks end of Model's scope.
            if line[0].upper() == '*END':
                break

            #Parse Restrain input's flag.
            nTag = int(line[0])
            for item in line:
                token = list(filter(None, item.strip().split()))
                if token[0].upper() == '-NDOF':
                    nDOF = int(token[1])
                elif token[0].upper() == '-DOF':
                    dof  = []
                    for k in range(1, len(token)):
                        dof.append(int(token[k]))

            if Point[nTag]['NDOF'] != len(dof):
                print('\x1B[31m ERROR \x1B[0m: The number of degree-of-freedom in Point[' + str(nTag) + '] does not match Restrain[' + str(nTag) + '].')
                sys.exit(-1)

            restrain = np.zeros(nDOF, dtype='int')
            for k in range(nDOF):
                if dof[k] == 1:
                    restrain[k] = -1

            Point[nTag]['FREE' ] = restrain
            Point[nTag]['FIXED'] = True

    #Restrain file was provided.
    User['RESTRAIN'] = 'YES'

def GetMassInformation(fileHandler, User, Mass, Point):
    """
    This function adds a Mass Point to the specific point.

    Parameters
    ----------
    User : dict 
           The dictionary containing all the user's relevant information.

    Returns
    -------
    Point : dict
           The dictionary containing all point information.
    """
    #Parse the partition input.
    while True:
        line = list(filter(None, fileHandler.readline().strip().split(',')))
        if line:
            #Checks end of Model's scope.
            if line[0].upper() == '*END':
                break

            #Parse Model input's flag.
            nTag = int(line[0])
            for item in line:
                token = list(filter(None, item.strip().split()))
                if token[0].upper() == '-NDOF':
                    nDOF = int(token[1])
                elif token[0].upper() == '-MASS':
                    mass  = []
                    for k in range(1, len(token)):
                        mass.append(float(token[k]))

            if Point[nTag]['NDOF'] != len(mass):
                print('\x1B[31m ERROR \x1B[0m: The number of degree-of-freedom in Point[' + str(nTag) + '] does not match Mass[' + str(nTag) + '].')
                sys.exit(-1)

            Mass[nTag] = {'NDOF': nDOF, 'MASS': mass}

    #Total number of points in the mesh
    User['NMASS'] = len(Mass)

    #Point file was provided.
    User['MASS'] = 'YES'

def GetSupportInformation(fileHandler, User, Support, Point):
    """
    This function adds a Initial condition to the specific point.

    Parameters
    ----------
    User : dict 
           The dictionary containing all the user's relevant information.

    Returns
    -------
    Support : dict
           The dictionary containing all specified support motion information.
    """
    #Parse the partition input.
    while True:
        line = list(filter(None, fileHandler.readline().strip().split(',')))
        if line:
            #Checks end of Model's scope.
            if line[0].upper() == '*END':
                break

            #Parse Model input's flag.
            file  = []
            value = []
            DOF   = []
            nTag  = int(line[0])
            for item in line:
                token = list(filter(None, item.strip().split()))
                if token[0].upper() == '-DOF':
                    DOF.append(int(token[1]) - 1)
                elif token[0].upper() == '-TYPE':
                    cond = token[1].upper()
                elif token[0].upper() == '-PATH':
                    file.append(User['FOLDER'] + token[1])
                elif token[0].upper() == '-VAL':
                    value.append(float(token[1]))

            #Checks consistency between the DOF and number of DOF in point
            if DOF[0] > Point[nTag]['NDOF']:
                print('\x1B[31m ERROR \x1B[0m: The number of degree-of-freedom in Point[' + str(nTag) + '] does not match SupportMotion[' + str(nTag) + '].')
                sys.exit(-1)

            #Extend the list if the point exists already
            if nTag in Support:
                Support[nTag]['DOF'].extend(DOF)
                if 'VALUE' in Support[nTag]:
                    Support[nTag]['VALUE'].extend(value)
                if 'PATH' in Support[nTag]:
                    Support[nTag]['PATH' ].extend(file)
            else:
                Support[nTag] = {'DOF': DOF, 'TYPE': cond, 'VALUE': value, 'PATH': file}

        #Checks if the file can be opened.
        if(Support[nTag]['TYPE'] == 'DYNAMIC'):
            try:
                with open(Support[nTag]['PATH'][-1], "r") as fh:
                    fh.close()
            except IOError as e:
                print('\x1B[31m ERROR \x1B[0m: In *SUPPORTMOTION the FILE=' + Support[nTag]['PATH'] + ' could not be opened.')
                print(' Check if both the file path and file name provided are correct.')
                sys.exit(-1)

    #Total number of points in the mesh
    User['NSUPPORT'] = len(Support)

    #Point file was provided.
    User['SUPPORT'] = 'YES'

def GetInitialInformation(fileHandler, User, Initial, Point):
    """
    This function adds a Initial condition to the specific point.

    Parameters
    ----------
    User : dict 
           The dictionary containing all the user's relevant information.

    Returns
    -------
    Initial : dict
           The dictionary containing all initial condition information.
    """
    #Parse the partition input.
    while True:
        line = list(filter(None, fileHandler.readline().strip().split(',')))
        if line:
            #Checks end of Model's scope.
            if line[0].upper() == '*END':
                break

            #Parse Model input's flag.
            nTag = int(line[0])
            for item in line:
                token = list(filter(None, item.strip().split()))
                if token[0].upper() == '-NDOF':
                    nDOF = int(token[1])
                elif token[0].upper() == '-TYPE':
                    if token[1].upper() == 'DISP':
                        cond = 1
                    elif token[1].upper() == 'VEL':
                        cond = 2
                    elif token[1].upper() == 'ACCEL':
                        cond = 3
                    else:
                        print('\x1B[31m ERROR \x1B[0m: The initial state specied in Point[' + str(nTag) + '] is not recognized.')
                elif token[0].upper() == '-VALS':
                    initial  = []
                    for k in range(1, len(token)):
                        initial.append(float(token[k]))

            if Point[nTag]['NDOF'] != len(initial):
                print('\x1B[31m ERROR \x1B[0m: The number of degree-of-freedom in Point[' + str(nTag) + '] does not match InitialCondition[' + str(nTag) + '].')
                sys.exit(-1)

            Initial[nTag] = {'NDOF': nDOF, 'TYPE': cond, 'VALUES': initial}

    #Total number of points in the mesh
    User['NINITIAL'] = len(Initial)

    #Point file was provided.
    User['INITIAL'] = 'YES'

def GetNodesInformation(fileHandler, User, Point):
    """
    This function creates a Points dictionary that contains the point identifier, 
    number of degree of freedom, and the space coordinates.

    Parameters
    ----------
    User : dict 
           The dictionary containing all the user's relevant information.

    Returns
    -------
    Point : dict
           The dictionary containing all point information.
    """
    #Inverse Mapping counter
    gTag  = 0

    #Parse the partition input.
    while True:
        line = list(filter(None, fileHandler.readline().strip().split(',')))
        if line:
            #Checks end of Model's scope.
            if line[0].upper() == '*END':
                break

            #Parse Model input's flag.
            nTag = int(line[0])
            for item in line:
                token = list(filter(None, item.strip().split()))
                if token[0].upper() == '-NDOF':
                    nDOF = int(token[1])
                elif token[0].upper() == '-POS':
                    if len(token) == 3:
                        coordinates = [float(token[1]), float(token[2])]
                    elif len(token) == 4:
                        coordinates = [float(token[1]), float(token[2]), float(token[3])]

            Point[nTag] = {'TAG': gTag, 'NDOF': nDOF, 'COORDINATES': coordinates, 'FREE': np.zeros(nDOF, dtype='int'), 'TOTAL': np.zeros(nDOF, dtype='int'), 'FIXED': False, 'DEFECTIVE': True}
            gTag += 1

    #Total number of points in the mesh
    User['NPOINTS'] = len(Point)

    #Point file was provided.
    User['POINT'] = 'YES'

def GetPartitionInformation(fileHandler, User):
    """
    This function adds to User dictionary information related to partition
    and ordering scheme. This information is used to generate the SeismoSSI 
    input files as well as Solution and Paraview directories.

    Parameters
    ----------
    path : str 
           The file location to be parsed

    Returns
    -------
    User : dict
           The dictionary containing all relevant information.
    """
    #Parse the partition input.
    while True:
        line = list(filter(None, fileHandler.readline().strip().split(',')))
        if line:
            #Checks end of Partition's scope.
            if line[0].upper() == '*END':
                break

            #Parse Model's Partition requirements.
            for item in line:
                token = list(filter(None, item.strip().split()))
                if token[0].upper() == '-METIS':
                    User['METIS'] = token[1]
                elif token[0].upper() == '-NP':
                    User['NPART'] = int(token[1])
                elif token[0].upper() == '-ORDERING':
                    User['ORDERING'] = token[1].upper()

def GetModelInformation(fileHandler, User, path):
    """
    This function creates a User dictionary that contains all user's 
    pre-difined information. This information is used to generate the
    SeismoSSI input files as well as Solution and Paraview directories.

    Parameters
    ----------
    path : str 
           The file location to be parsed

    Returns
    -------
    User : dict
           The dictionary containing all relevant information.
    """
    #Gets the Working Directory Path
    folder = path
    folder = folder.split('/')
    folder[-1] = ''
    folder = '/'.join(folder)

    #User's Dictionary and Pre-Define Values
    User['FILE'       ] = path
    User['FOLDER'     ] = folder
    User['MODEL'      ] = 'SeismoVLab'
    User['NPART'      ] = 1
    User['NPARAVIEW'  ] = 0 
    User['NFEATURES'  ] = 0
    User['ORDERING'   ] = 'PLAIN'
    User['MASSFORM'   ] = 'CONSISTENT'
    User['POINT'      ] = 'NO'
    User['MATERIAL'   ] = 'NO'
    User['MASS'       ] = 'NO'
    User['INITIAL'    ] = 'NO'
    User['SUPPORT'    ] = 'NO'
    User['SECTION'    ] = 'NO'
    User['RESTRAIN'   ] = 'NO'
    User['CONSTRAINT' ] = 'NO' 
    User['DIAPHRAGM'  ] = 'NO'
    User['DAMPING'    ] = 'NO'
    User['ELEMENT'    ] = 'NO'
    User['FUNCTION'   ] = 'NO'
    User['LOAD'       ] = 'NO'
    User['COMBINATION'] = 'NO'
    User['RECORDER'   ] = 'NO'
    User['SIMULATION' ] = 'NO'
    User['PETSC'      ] = 'NO'
    User['NPOINTS'    ] = 0
    User['NMASS'      ] = 0
    User['NINITIAL'   ] = 0
    User['NCONSTRAINT'] = 0 
    User['NFREEDOF'   ] = 0
    User['NTOTALDOF'  ] = 0

    #Parse the model input.
    while True:
        line = list(filter(None, fileHandler.readline().strip().split(',')))
        if line:
            #Checks end of Model's scope.
            if line[0].upper() == '*END':
                break

            #Parse Model's general information.
            for item in line:
                token = list(filter(None, item.strip().split()))
                if token[0].upper() == '-NAME':
                    User['MODEL'] = token[1]
                elif token[0].upper() == '-NDIM':
                    User['DIMENSION'] = int(token[1])
                elif token[0].upper() == '-MASS':
                    User['MASSFORM'] = token[1].upper()

    #CREATES THE SOLUTION DIRECTORY
    dirName = User['FOLDER'] + 'Solution/'
    if not os.path.exists(dirName):
        os.mkdir(dirName) 

    #CREATES THE PARAVIEW DIRECTORY
    dirName = User['FOLDER'] + 'Paraview/'
    if not os.path.exists(dirName):
        os.mkdir(dirName)
