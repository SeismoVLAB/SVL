#!/usr/bin/python3

import os
import sys
import numpy as np
from Metis import Metis as mt

def WriteSimulationPartition(SeismoVLabfile, User, Simulation, k):
    """
    This function appends the Simulation specification to run the
    analysis for a particular load combination pattern. 

    Parameters
    ----------
    SeismoVLabfile : file 
           The partition file where material's information is appended.
    Simulation : dict 
           The dictionary containing all simulation information

    Output
    -------
    None
    """
    for sTag in Simulation:
        #GETS THE ANALYSIS INPUT PARAMETERS
        NSTEPS = Simulation[sTag]['ANALYSIS']['NSTEP']
        SeismoVLabfile.write("ANALYSIS %d %s %d " % (sTag, Simulation[sTag]['ANALYSIS']['NAME'], NSTEPS))

        #GETS THE ALGORITHM INPUT PARAMETERS
        NSTEPS  = Simulation[sTag]['ALGORITHM']['NSTEP']
        CNVGTOL = Simulation[sTag]['ALGORITHM']['CNVGTOL']
        if Simulation[sTag]['ALGORITHM']['CNVGTEST'].upper() == 'UNBALANCEFORCE':
            CNVGTEST = 1
        elif Simulation[sTag]['ALGORITHM']['CNVGTEST'].upper() == 'INCREMENTALDISPLACEMENT':
            CNVGTEST = 2
        elif Simulation[sTag]['ALGORITHM']['CNVGTEST'].upper() == 'RELATIVEUNBALANCEFORCE':
            CNVGTEST = 3
        elif Simulation[sTag]['ALGORITHM']['CNVGTEST'].upper() == 'RELATIVEINCREMENTALDISPLACEMENT':
            CNVGTEST = 4
        SeismoVLabfile.write("ALGORITHM %s %E %d %d " % (Simulation[sTag]['ALGORITHM']['NAME'], CNVGTOL, NSTEPS, CNVGTEST))
        
        #GETS THE INETGRATOR INPUT PARAMETERS
        TIMESTEP = Simulation[sTag]['INTEGRATOR']['TIMESTEP']
        MASSTOL  = Simulation[sTag]['INTEGRATOR']['MTOL']
        STIFFTOL = Simulation[sTag]['INTEGRATOR']['KTOL']
        FORCETOL = Simulation[sTag]['INTEGRATOR']['FTOL']
        SeismoVLabfile.write("INTEGRATOR %s %E %E %E %E " % (Simulation[sTag]['INTEGRATOR']['NAME'], MASSTOL, STIFFTOL, FORCETOL, TIMESTEP))

        #SETS THE SOLVER INPUT PARAMETERS
        if Simulation[sTag]['SOLVER']['NAME'] == 'EIGEN':
            if Simulation[sTag]['SOLVER']['UPDATE'] == 'OFF':
                SeismoVLabfile.write("SOLVER %s 1\n" % Simulation[sTag]['SOLVER']['NAME'])
            else:
                SeismoVLabfile.write("SOLVER %s 0\n" % Simulation[sTag]['SOLVER']['NAME'])
        elif Simulation[sTag]['SOLVER']['NAME'] == 'MUMPS':
            if Simulation[sTag]['SOLVER']['UPDATE'] == 'OFF':
                SeismoVLabfile.write("SOLVER %s %d 1\n" % (Simulation[sTag]['SOLVER']['NAME'], Simulation[sTag]['SOLVER']['OPTION']))
            else:
                SeismoVLabfile.write("SOLVER %s %d 0\n" % (Simulation[sTag]['SOLVER']['NAME'], Simulation[sTag]['SOLVER']['OPTION']))
        elif Simulation[sTag]['SOLVER']['NAME'] == 'PETSC':
            SeismoVLabfile.write("SOLVER %s %d %E %d %d\n" % (Simulation[sTag]['SOLVER']['NAME'], Simulation[sTag]['SOLVER']['OPTION'], Simulation[sTag]['SOLVER']['TOLERANCE'], User['D_NZ'][k], User['O_NZ'][k]))

def WriteRecordersPartition(SeismoVLabfile, User, Recorder, nodeSubdomain, elemSubdomain, k):
    """
    This function appends the Recorder needed to store information for
    a particular combination load pattern. 

    Parameters
    ----------
    SeismoVLabfile : file 
           The partition file where material's information is appended.
    User : dict 
           The dictionary containing all the user's relevant information.
    Recorder : dict 
           The dictionary containing all recording information
    k : int
           The number of the partition.

    Output
    -------
    None
    """
    for rTag in Recorder:
        OUTFILE = Recorder[rTag]['FILE']
        NSAMPLE = Recorder[rTag]['NSAMPLE']
        NDIGITS = Recorder[rTag]['PRECISION']

        #Modifies the Output File to be consistent with Partition
        if OUTFILE.find('.') != -1:
            OUTFILE = OUTFILE.replace(".", '.' + str(k) + '.')
        else:
            OUTFILE = OUTFILE + '.' + str(k)
        #Writes the Recorder's according to Point or Element
        if Recorder[rTag]['NAME'].upper() == 'NODE':
            nTags  = sorted(nodeSubdomain.intersection(Recorder[rTag]['LIST']))
            if nTags:
                SeismoVLabfile.write("RECORDER %s %s %d %d %s %d" % (Recorder[rTag]['NAME'], OUTFILE, NSAMPLE, NDIGITS, Recorder[rTag]['RESPONSE'], len(nTags)))
                for n in nTags:
                    SeismoVLabfile.write(" %d" % n)
                SeismoVLabfile.write("\n")
        elif Recorder[rTag]['NAME'].upper() == 'ELEMENT':
            eTags  = sorted(list(set(elemSubdomain).intersection(Recorder[rTag]['LIST'])))
            if eTags:
                SeismoVLabfile.write("RECORDER %s %s %d %d %s %d" % (Recorder[rTag]['NAME'], OUTFILE, NSAMPLE, NDIGITS, Recorder[rTag]['RESPONSE'], len(eTags)))
                for e in eTags:
                    SeismoVLabfile.write(" %d" % e)
                SeismoVLabfile.write("\n")
        elif Recorder[rTag]['NAME'].upper() == 'SECTION':
            eTags = sorted(list(set(elemSubdomain).intersection(Recorder[rTag]['LIST'])))
            if eTags:
                pTags = Recorder[rTag]['POSITION']
                SeismoVLabfile.write("RECORDER %s %s %d %d %s %d" % (Recorder[rTag]['NAME'], OUTFILE, NSAMPLE, NDIGITS, Recorder[rTag]['RESPONSE'], len(pTags)))
                for p in pTags:
                    SeismoVLabfile.write(" %E" % p)
                SeismoVLabfile.write(" %d" % len(eTags))
                for e in eTags:
                    SeismoVLabfile.write(" %d" % e)
                SeismoVLabfile.write("\n")
        elif Recorder[rTag]['NAME'].upper() == 'PARAVIEW':
            FILENAME = Recorder[rTag]['FILE']
            FILENAME = FILENAME.split('.')
            FILENAME = FILENAME[0] + '_PART' + str(k)
            SeismoVLabfile.write("RECORDER %s %s %d %d %d\n" % (Recorder[rTag]['NAME'], FILENAME, NSAMPLE, NDIGITS, User['NPARAVIEW']))

def WriteLoadCombosPartition(SeismoVLabfile, loadSubdomain, Combination):
    """
    This function appends the Load Combination that are used to perform 
    the desired analysis. 

    Parameters
    ----------
    SeismoVLabfile : file 
           The partition file where material's information is appended.
    Combination : dict 
           The dictionary containing all combination to be performed

    Output
    -------
    None
    """
    for cTag in Combination:
        SeismoVLabfile.write("COMBINATION %d %s" % (cTag, Combination[cTag]['NAME']))
        #Gets Load that are in this combination
        lTag = list(set(Combination[cTag]['LOAD']).intersection(loadSubdomain))
        if lTag:
            lCombos  = Combination[cTag]['LOAD']
            lFactors = Combination[cTag]['FACTOR']
            SeismoVLabfile.write(" %d" % len(lTag))
            for j in loadSubdomain:
                for i in range(len(lCombos)):
                    if lCombos[i] == j:
                        SeismoVLabfile.write(" %1.5f %d" % (lFactors[i], lCombos[i]))
            SeismoVLabfile.write("\n")
        else:
            SeismoVLabfile.write(" -1\n")

def WriteLoadsPartition(SeismoVLabfile, Load, Function, Surface, nodeSubdomain, elemSubdomain):
    """
    This function appends the Loads that are used in a combination 
    state to perform the analysis. 

    Parameters
    ----------
    SeismoVLabfile : file 
           The partition file where material's information is appended.
    Load : dict 
           The dictionary containing all load information
    Material : dict
           The dictionary containing all material information.

    Output
    -------
    None
    """
    #Loads that belong to this partition
    LoadTags = list()
    for lTag in Load:
        lName = Load[lTag]['NAME']
            
        #PARSE POINT LOAD CASE
        if lName == 'POINTLOAD':
            fTag  = Load[lTag]['FUNCTION']
            fName = Function[fTag]['NAME']
            fType = Function[fTag]['TYPE']
            Tags  = sorted(nodeSubdomain.intersection(Load[lTag]['LIST']))
            #There are Loads in this partition
            if Tags:
                SeismoVLabfile.write("%s %d %s " % (lName, lTag, fType))
                #The Load Condition
                DIR = Function[fTag]['DIRECTION']
                if fType == 'STATIC':
                    SeismoVLabfile.write("%E %d" % (Function[fTag]['MAGNITUDE'], len(DIR)))
                elif fType == 'DYNAMIC':
                    path = Function[fTag]['PATH']
                    path = path.replace(" ", "~")
                    SeismoVLabfile.write("%s %d" % (path, len(DIR)))
                #The Load direction
                for k in DIR:
                    SeismoVLabfile.write(" %1.5f" % k)
                #Points that are loaded
                SeismoVLabfile.write(" %d" % len(Tags))
                for k in Tags:
                    SeismoVLabfile.write(" %d" % k)
                SeismoVLabfile.write("\n")

                #Guarantees Point Load is considered in only one partition
                LoadTags.append(lTag)
                Load[lTag]['LIST'] = list(set(Load[lTag]['LIST']).difference(Tags))

        #PARSE POINT LOAD CASE
        elif lName == 'ELEMENTLOAD':
            fTag  = Load[lTag]['FUNCTION']
            fName = Function[fTag]['NAME']
            fType = Function[fTag]['TYPE']
            Tags  = sorted(list(set(elemSubdomain).intersection(Load[lTag]['LIST'])))
            #There are Loads in this partition
            if Tags:
                SeismoVLabfile.write("%s %d %s " % (lName, lTag, fType))
                #The Load Condition
                lType = Load[lTag]['TYPE']
                if lType == 'SURFACE':
                    DIR = Function[fTag]['DIRECTION']
                    #The Load Condition
                    if fType == 'STATIC':
                        SeismoVLabfile.write("%s %E %d" % (lType, Function[fTag]['MAGNITUDE'], len(DIR)))
                        #The Load direction
                        for k in DIR:
                            SeismoVLabfile.write(" %1.5f" % k)
                        #Elements that are loaded
                        SeismoVLabfile.write(" %d" % len(Tags))
                        for k in Tags:
                            SeismoVLabfile.write(" %d" % k)
                        #Surfaces that are loaded
                        for k in Tags:
                            SeismoVLabfile.write(" %d" % Surface[k]['FACE'])
                        SeismoVLabfile.write("\n")
                    else:
                        print(' \x1B[31mERROR \x1B[0m:SURFACE (Dynamic) ELEMENTLOAD is not Implemented yet')
                        sys.exit(-1)
                elif lType == 'BODY':
                    DIR = Function[fTag]['DIRECTION']
                    #The Load Condition
                    if fType == 'STATIC':
                        SeismoVLabfile.write("%s %E %d" % (lType, Function[fTag]['MAGNITUDE'], len(DIR)))
                    elif fType == 'DYNAMIC':
                        path = Function[fTag]['PATH']
                        path = path.replace(" ", "~")
                        SeismoVLabfile.write("%s %s %d" % (lType, path, len(DIR)))
                    #The Load direction
                    for k in DIR:
                        SeismoVLabfile.write(" %1.5f" % k)
                    #Elements that are loaded
                    SeismoVLabfile.write(" %d" % len(Tags))
                    for k in Tags:
                        SeismoVLabfile.write(" %d" % k)
                    SeismoVLabfile.write("\n")
                elif lType == 'PLANEWAVE':
                    path = Function[fTag]['PATH']
                    path = path.replace(" ", "~")
                    SeismoVLabfile.write("GENERALWAVE %s %d" % (path, len(Tags)))
                    #Elements that are employed for DRM genral forces
                    for k in Tags:
                        SeismoVLabfile.write(" %d" % k)
                    SeismoVLabfile.write("\n")
                elif lType == 'GENERALWAVE':
                    path = Function[fTag]['PATH']
                    path = path.replace(" ", "~")
                    SeismoVLabfile.write("%s %s %d" % (lType, path, len(Tags)))
                    #Elements that are employed for DRM genral forces
                    for k in Tags:
                        SeismoVLabfile.write(" %d" % k)
                    SeismoVLabfile.write("\n")

                LoadTags.append(lTag)
        #PARSE SUPPORT MOTION CASE
        if lName == 'SUPPORTMOTION':
            Tags    = sorted(nodeSubdomain.intersection(Load[lTag]['LIST']))
            #There are Loads in this partition
            if Tags:
                SeismoVLabfile.write("POINTLOAD %d %s %d" % (lTag, lName, len(Tags)))

                #Points that have motion
                for k in Tags:
                    SeismoVLabfile.write(" %d" % k)
                SeismoVLabfile.write("\n")
                LoadTags.append(lTag)

    return LoadTags

def WriteMaterialPartition(SeismoVLabfile, Tags, Material):
    """
    This function appends the Material that are used to construct either 
    solid/structural Element to the SeismoVLab mesh input files. 

    Parameters
    ----------
    SeismoVLabfile : file 
           The partition file where material's information is appended.
    Tags : array 
           The Material indeces that belong to this partition
    Material : dict
           The dictionary containing all material information.

    Output
    -------
    None
    """
    #Loop over the materials in partition  
    for k in Tags:
        SeismoVLabfile.write("MATERIAL %s %s " %(k, Material[k]['NAME'].upper()))

        #Loop Over the possible Materials in SeismoVLab
        if Material[k]['NAME'].upper() == "ELASTIC1DLINEAR":
            SeismoVLabfile.write("%E %1.5f %1.5f\n" % (Material[k]['E'], Material[k]['NU'], Material[k]['RHO']))
        elif Material[k]['NAME'].upper() == "HERTZIAN1DLINEAR":
            SeismoVLabfile.write("%E %E %E %E\n" % (Material[k]['K1'], Material[k]['K2'], Material[k]['K3'], Material[k]['RHO']))
        elif Material[k]['NAME'].upper() == "VISCOUS1DLINEAR":
            SeismoVLabfile.write("%E\n" % Material[k]['ETA'])
        elif Material[k]['NAME'].upper() == "ELASTIC2DPLANESTRAIN":
            SeismoVLabfile.write("%E %1.5f %1.5f\n" % (Material[k]['E'], Material[k]['NU'], Material[k]['RHO']))
        elif Material[k]['NAME'].upper() == "ELASTIC2DPLANESTRESS":
            SeismoVLabfile.write("%E %1.5f %1.5f\n" % (Material[k]['E'], Material[k]['NU'], Material[k]['RHO']))
        elif Material[k]['NAME'].upper() == "ELASTIC3DLINEAR":
            SeismoVLabfile.write("%E %1.5f %1.5f\n" % (Material[k]['E'], Material[k]['NU'], Material[k]['RHO']))
        elif Material[k]['NAME'].upper() == "PLASTIC1DJ2":
            SeismoVLabfile.write("%E %1.5f %1.5f %1.5f %1.5f %E\n" % (Material[k]['E'], Material[k]['NU'], Material[k]['RHO'], Material[k]['K'], Material[k]['H'], Material[k]['SY']))
        elif Material[k]['NAME'].upper() == "PLASTICPLANESTRAINJ2":
            SeismoVLabfile.write("%E %1.5f %1.5f %1.5f %1.5f %E\n" % (Material[k]['K'], Material[k]['G'], Material[k]['RHO'], Material[k]['H'], Material[k]['BETA'], Material[k]['SY']))
        elif Material[k]['NAME'].upper() == "PLASTICPLANESTRAINBA":
            SeismoVLabfile.write("%E %E %1.5f %1.5f %1.5f %1.5f %E %1.5f\n" % (Material[k]['K'], Material[k]['G'], Material[k]['RHO'], Material[k]['H0'], Material[k]['h'], Material[k]['m'], Material[k]['SU'], Material[k]['BETA']))
        elif Material[k]['NAME'].upper() == "PLASTIC3DJ2":
            SeismoVLabfile.write("%E %1.5f %1.5f %1.5f %1.5f %E\n" % (Material[k]['K'], Material[k]['G'], Material[k]['RHO'], Material[k]['H'], Material[k]['BETA'], Material[k]['SY']))
        elif Material[k]['NAME'].upper() == "PLASTIC3DBA":
            SeismoVLabfile.write("%E %E %1.5f %1.5f %1.5f %1.5f %E %1.5f\n" % (Material[k]['K'], Material[k]['G'], Material[k]['RHO'], Material[k]['H0'], Material[k]['h'], Material[k]['m'], Material[k]['SU'], Material[k]['BETA']))
        else:
            print('\x1B[31mERROR \x1B[0m: The MATERIAL=' + Material[k]['NAME'].upper() + ' could not be recognized.')
            sys.exit(-1)
    
def WriteSectionsPartition(SeismoVLabfile, Tags, Section):
    """
    This function appends the Section that are used to construct structural
    Element to the SeismoVLab mesh input files. 

    Parameters
    ----------
    SeismoVLabfile : file 
           The partition file where section's information is appended.
    Tags : array 
           The Section indeces that belong to this partition
    Section : dict
           The dictionary containing all section information.

    Output
    -------
    None
    """
    #Loop over the sections in partition  
    for k in Tags:
        SeismoVLabfile.write("SECTION %s %d %s " % (Section[k]['TYPE'], k, Section[k]['PROPERTY']))

        if Section[k]['TYPE'].upper() == 'PLAIN':
            SeismoVLabfile.write(str(Section[k]['MATERIAL']) + ' ' + Section[k]['NAME'].upper())

            if 'theta' not in Section[k]:
                Section[k]['theta'] = 0.0
            if 'ip' not in Section[k]:
                Section[k]['ip'] = 10

            #Loop Over Section's Prarmeters values
            if Section[k]['NAME'].upper() == "RECTANGULAR":
                SeismoVLabfile.write(" %1.5f %1.5f %1.5f %d\n" % (Section[k]['h'], Section[k]['b'], Section[k]['theta'], Section[k]['ip']))
            elif Section[k]['NAME'].upper() == "CIRCULAR":
                SeismoVLabfile.write(" %1.5f %1.5f %d\n" % (Section[k]['r'], Section[k]['theta'], Section[k]['ip']))
            elif Section[k]['NAME'].upper() == "ANGLE" or Section[k]['NAME'].upper() == "CHANNEL" or Section[k]['NAME'].upper() == "TEE" or Section[k]['NAME'].upper() == "WIDEFLANGE":
                SeismoVLabfile.write(" %1.5f %1.5f %1.5f %1.5f %1.5f %d\n" % (Section[k]['h'], Section[k]['b'], Section[k]['tw'], Section[k]['tf'], Section[k]['theta'], Section[k]['ip']))
            elif Section[k]['NAME'].upper() == "RECTANGULARTUBE":
                SeismoVLabfile.write(" %1.5f %1.5f %1.5f %d\n" % (Section[k]['h'], Section[k]['b'], Section[k]['tw'], Section[k]['tf'], Section[k]['theta'], Section[k]['ip']))
            elif Section[k]['NAME'].upper() == "CIRCULARTUBE":
                SeismoVLabfile.write(" %1.5f %1.5f %1.5f %d\n" % (Section[k]['re'], Section[k]['ri'], Section[k]['theta'], Section[k]['ip']))
            elif Section[k]['NAME'].upper() == "THINAREA":
                SeismoVLabfile.write(" %1.5f\n" % Section[k]['t'])
            else:
                print('\x1B[31mERROR \x1B[0m: The PLAIN SECTION=' + Section[k]['NAME'].upper() + ' could not be recognized.')
                sys.exit(-1)

        elif Section[k]['TYPE'] == 'FIBER':
            print('TODO: FIBER SECTION')

def WriteMassPartition(SeismoVLabfile, Tags, Mass):
    """
    This function appends the Masses in the finite element mesh to the 
    SeismoVLab mesh input files. 

    Parameters
    ----------
    SeismoVLabfile : file 
           The partition file where point's information is appended.
    Tags : array 
           The Point indeces that belong to this partition
    Mass : dict
           The dictionary containing all point mass information.

    Output
    -------
    None
    """
    mTags = set(Mass.keys())
    nTags = sorted(mTags.intersection(Tags))
    if nTags:
        for node in nTags:
            SeismoVLabfile.write("MASS %d %d" % (node, Mass[node]['NDOF']))
            for mass in Mass[node]['MASS']:
                SeismoVLabfile.write(" %f" % mass)
            SeismoVLabfile.write('\n')

        for node in nTags:
            del Mass[node]

def WriteSupportPartition(SeismoVLabfile, Tags, Support):
    """
    This function appends the Support motions in the finite element mesh to 
    the SeismoVLab mesh input files. 

    Parameters
    ----------
    SeismoVLabfile : file 
           The partition file where point's information is appended.
    Tags : array 
           The Point indeces that belong to this partition
    Support : dict
           The dictionary containing all specified support motions.

    Output
    -------
    None
    """
    mTags = set(Support.keys())
    nTags = sorted(mTags.intersection(Tags))
    if nTags:
        for node in nTags:
            for k in range(len(Support[node]['DOF'])):
                SeismoVLabfile.write("SUPPORTMOTION %d %s " % (node, Support[node]['TYPE']))
                if(Support[node]['TYPE'] == 'STATIC'):
                    SeismoVLabfile.write("%E %d" % (Support[node]['VALUE'][k], Support[node]['DOF'][k]))
                elif(Support[node]['TYPE'] == 'DYNAMIC'):
                    SeismoVLabfile.write("%s %d" % (Support[node]['PATH'][k], Support[node]['DOF'][k]))
                SeismoVLabfile.write('\n')

def WriteInitialPartition(SeismoVLabfile, Tags, Initial):
    """
    This function appends the Initial conditions in the finite element mesh to 
    the SeismoVLab mesh input files. 

    Parameters
    ----------
    SeismoVLabfile : file 
           The partition file where point's information is appended.
    Tags : array 
           The Point indeces that belong to this partition
    Initial : dict
           The dictionary containing all initial condition information.

    Output
    -------
    None
    """
    mTags = set(Initial.keys())
    nTags = sorted(mTags.intersection(Tags))
    if nTags:
        for node in nTags:
            SeismoVLabfile.write("INITIALSTATE %d %d %d" % (node, Initial[node]['TYPE'], Initial[node]['NDOF']))
            for initial in Initial[node]['VALUES']:
                SeismoVLabfile.write(" %f" % initial)
            SeismoVLabfile.write('\n')

def WriteNodePartition(SeismoVLabfile, Tags, Point):
    """
    This function appends the Point in the finite element mesh to the 
    SeismoVLab mesh input files. 

    Parameters
    ----------
    SeismoVLabfile : file 
           The partition file where point's information is appended.
    Tags : array 
           The Point indeces that belong to this partition
    Point : dict
           The dictionary containing all point information.

    Output
    -------
    None
    """
    nTags = sorted(Tags) 
    #Loop over the Point in partition  
    for k in nTags:
        SeismoVLabfile.write("NODE %d %d " % (k, Point[k]['NDOF']))

        #Loop Over Node's Total-Degree-of-Freedom
        total = Point[k]['TOTAL']
        for dof in total:
            SeismoVLabfile.write("%d " % dof)

        #Loop Over Node's Free-Degree-of-Freedom
        free = Point[k]['FREE']
        for dof in free:
            SeismoVLabfile.write("%d " % dof)

        #Loop Over Node's Coordinates
        coordinates = Point[k]['COORDINATES']
        for coord in coordinates:
            SeismoVLabfile.write("%s " % coord)

        SeismoVLabfile.write('\n')

def WriteConstraintsPartition(SeismoVLabfile, Tags, Constraint, Point):
    """
    This function appends the Constraints that are applied to the 
    degree-of-freedom for a particular node. Constraints are 
    represented by negative degree-of freedom numbering. 

    Parameters
    ----------
    SeismoVLabfile : file 
           The partition file where consraints's information is appended.
    Tags : array 
           The Consraints indeces that belong to this partition
    Constraint : dict
           The dictionary containing all consraints information.

    Output
    -------
    None
    """
    cTags = sorted(Tags)
    for k in cTags:
        #Slave Point Information
        SlaveNode = Constraint[k]['SLAVENODE']
        SlaveDOF  = Constraint[k]['SLAVEDOF' ]

        #Master Point information
        MasterNode = Constraint[k]['MASTERNODE']
        MasterDOF  = Constraint[k]['MASTERDOF' ]
        CombCoeffs = Constraint[k]['FACTORS']

        nCombinations = len(CombCoeffs)
        SeismoVLabfile.write("CONSTRAINT %d %d %d" %(k, Point[SlaveNode]['TOTAL'][SlaveDOF], nCombinations))

        for i in range(nCombinations):
            node = MasterNode[i]
            dof  = MasterDOF[i]
            SeismoVLabfile.write(" %d %f" % (Point[node]['FREE'][dof], CombCoeffs[i]))

        SeismoVLabfile.write('\n')

def WriteDampingsPartition(SeismoVLabfile, Damping, eTags, User):
    """
    This function appends the Damping models applied to element in the 
    SeismoVLab mesh input files. 

    Parameters
    ----------
    SeismoVLabfile : file 
           The partition file where element's information is appended.
    Element : dict
           The dictionary containing all damping models.
    eTags : array 
           The Element indeces that belong to this partition
    Element : dict
           The dictionary containing all element information.

    Output
    -------
    None
    """
    for dTag in Damping:
        eList  = list(set(eTags).intersection(Damping[dTag]['LIST']))

        if eList:
            Name = Damping[dTag]['NAME']
            SeismoVLabfile.write("DAMPING %d %s" %(dTag, Name))
            if Name.upper() == 'RAYLEIGH':
                SeismoVLabfile.write(" %1.15f %1.15f" %(Damping[dTag]['ALPHA'], Damping[dTag]['BETA']))
            elif Name.upper() == 'CAUGHEY':
                pass
            elif Name.upper() == 'CAPPED':
                pass

            SeismoVLabfile.write(" %d" % len(eList))
            for k in eList:
                SeismoVLabfile.write(" %d" % k)
            SeismoVLabfile.write('\n')

def WriteElementPartition(SeismoVLabfile, Tags, Element, User):
    """
    This function appends the Element in the finite element mesh to the 
    SeismoVLab mesh input files. 

    Parameters
    ----------
    SeismoVLabfile : file 
           The partition file where element's information is appended.
    Tags : array 
           The Element indeces that belong to this partition
    Element : dict
           The dictionary containing all element information.

    Output
    -------
    None
    """
    #Loop over the Element in partition
    nParaview = 0 
    DIMS = User['DIMENSION']
    eTag = sorted(Tags)
    for k in eTag:
        name = Element[k]['NAME'].upper()
        connection = Element[k]['CONNECTIVITY']

        #Writes a text line with the element information
        SeismoVLabfile.write("ELEMENT %s %s " %(k, name))

        for node in connection:
            SeismoVLabfile.write("%d " % node)

        if name == 'ZEROLENGTH1D':
            nParaview += 3
            DIR = Element[k]['DIRECTION']
            SeismoVLabfile.write("%d %d %d\n" %(Element[k]['MATERIAL'], DIMS, DIR))
        elif  name == 'LIN2DTRUSS2':
            nParaview += 3
            SeismoVLabfile.write("%d %1.5f\n" %(Element[k]['MATERIAL'], Element[k]['AREA']))
        elif  name == 'KIN2DTRUSS2':
            nParaview += 3
            SeismoVLabfile.write("%d %1.5f\n" %(Element[k]['MATERIAL'], Element[k]['AREA']))
        elif name == 'LIN3DTRUSS2':
            nParaview += 3
            SeismoVLabfile.write("%d %1.5f\n" %(Element[k]['MATERIAL'], Element[k]['AREA']))
        elif name == 'KIN3DTRUSS2':
            nParaview += 3
            SeismoVLabfile.write("%d %1.5f\n" %(Element[k]['MATERIAL'], Element[k]['AREA']))
        elif name == 'LIN2DTRUSS3':
            nParaview += 4
            if 'NPOINTS' not in Element[k]:
                Element[k]['NPOINTS'] = 3
            if Element[k]['NPOINTS'] < 3:
                Element[k]['NPOINTS'] = 3
            SeismoVLabfile.write("%d %1.5f %d\n" %(Element[k]['MATERIAL'], Element[k]['AREA'], Element[k]['NPOINTS']))
        elif name == 'LIN3DTRUSS3':
            nParaview += 4
            if 'NPOINTS' not in Element[k]:
                Element[k]['NPOINTS'] = 3
            if Element[k]['NPOINTS'] < 3:
                Element[k]['NPOINTS'] = 3
            SeismoVLabfile.write("%d %1.5f %d\n" %(Element[k]['MATERIAL'], Element[k]['AREA'], Element[k]['NPOINTS']))
        elif name == 'LIN2DQUAD4':
            nParaview += 5
            if 'NPOINTS' not in Element[k]:
                Element[k]['NPOINTS'] = 4
            if Element[k]['NPOINTS'] < 4:
                Element[k]['NPOINTS'] = 4
            SeismoVLabfile.write("%d %1.5f %d\n" %(Element[k]['MATERIAL'], Element[k]['THICKNESS'], Element[k]['NPOINTS']))
        elif name == 'LIN2DQUAD8':
            nParaview += 9
            if 'NPOINTS' not in Element[k]:
                Element[k]['NPOINTS'] = 9
            if Element[k]['NPOINTS'] < 9:
                Element[k]['NPOINTS'] = 9
            SeismoVLabfile.write("%d %1.5f %d\n" %(Element[k]['MATERIAL'], Element[k]['THICKNESS'], Element[k]['NPOINTS']))
        elif name == 'PML2DQUAD4':
            nParaview += 5
            if 'NPOINTS' not in Element[k]:
                Element[k]['NPOINTS'] = 4
            if Element[k]['NPOINTS'] < 4:
                Element[k]['NPOINTS'] = 4
            X0  = Element[k]['X0']
            DIR = Element[k]['NORMAL']
            SeismoVLabfile.write("%d %1.5f %d %1.5f %E %1.5f %1.5f %1.5f %1.5f %d\n" %(Element[k]['MATERIAL'], Element[k]['THICKNESS'], Element[k]['DEGREE'], Element[k]['L'], Element[k]['R'], X0[0], X0[1], DIR[0], DIR[1], Element[k]['NPOINTS']))
        elif name == 'PML2DQUAD8':
            nParaview += 9
            if 'NPOINTS' not in Element[k]:
                Element[k]['NPOINTS'] = 9
            if Element[k]['NPOINTS'] < 9:
                Element[k]['NPOINTS'] = 9
            X0  = Element[k]['X0']
            DIR = Element[k]['NORMAL']
            SeismoVLabfile.write("%d %1.5f %d %1.5f %E %1.5f %1.5f %1.5f %1.5f %d\n" %(Element[k]['MATERIAL'], Element[k]['THICKNESS'], Element[k]['DEGREE'], Element[k]['L'], Element[k]['R'], X0[0], X0[1], DIR[0], DIR[1], Element[k]['NPOINTS']))
        elif  name == 'LIN2DFRAME2':
            nParaview += 3
            if 'NPOINTS' not in Element[k]:
                Element[k]['NPOINTS'] = 3
            if Element[k]['NPOINTS'] < 3:
                Element[k]['NPOINTS'] = 3
            if 'FORMULATION' not in Element[k]:
                Element[k]['FORMULATION'] = "Bernoulli"
            SeismoVLabfile.write("%d %s %d\n" %(Element[k]['SECTION'], Element[k]['FORMULATION'], Element[k]['NPOINTS']))
        elif  name == 'KIN2DFRAME2':
            nParaview += 3
            SeismoVLabfile.write("%d\n" % Element[k]['SECTION'])
        elif  name == 'LIN3DFRAME2':
            nParaview += 3
            if 'NPOINTS' not in Element[k]:
                Element[k]['NPOINTS'] = 3
            if Element[k]['NPOINTS'] < 3:
                Element[k]['NPOINTS'] = 3
            if 'FORMULATION' not in Element[k]:
                Element[k]['FORMULATION'] = "Bernoulli"
            SeismoVLabfile.write("%d %s %d\n" %(Element[k]['SECTION'], Element[k]['FORMULATION'], Element[k]['NPOINTS']))
        elif name == 'LIN3DSHELL4':
            nParaview += 5
            if 'NPOINTS' not in Element[k]:
                Element[k]['NPOINTS'] = 9
            if Element[k]['NPOINTS'] < 9:
                Element[k]['NPOINTS'] = 9
            SeismoVLabfile.write("%d %d\n" %(Element[k]['SECTION'], Element[k]['NPOINTS']))
        elif name == 'LIN3DHEXA8':
            nParaview += 9
            if 'NPOINTS' not in Element[k]:
                Element[k]['NPOINTS'] = 8
            if Element[k]['NPOINTS'] < 8:
                Element[k]['NPOINTS'] = 8
            SeismoVLabfile.write("%d %d\n" %(Element[k]['MATERIAL'], Element[k]['NPOINTS']))
        elif name == 'KIN3DHEXA8':
            nParaview += 9
            if 'NPOINTS' not in Element[k]:
                Element[k]['NPOINTS'] = 8
            if Element[k]['NPOINTS'] < 8:
                Element[k]['NPOINTS'] = 8
            SeismoVLabfile.write("%d %d\n" %(Element[k]['MATERIAL'], Element[k]['NPOINTS']))
        elif name == 'LIN3DHEXA20':
            nParaview += 21
            if 'NPOINTS' not in Element[k]:
                Element[k]['NPOINTS'] = 27
            if Element[k]['NPOINTS'] < 27:
                Element[k]['NPOINTS'] = 27
            SeismoVLabfile.write("%d %d\n" %(Element[k]['MATERIAL'], Element[k]['NPOINTS']))
        elif  name == 'UNXBOUCWEN2DLINK':
            nParaview += 3
            SeismoVLabfile.write("%1.5f %1.5f %1.5f %1.5f %1.5f %E %E %E %1.5f %1.5f %d %d\n" %(Element[k]['ALPHA'], Element[k]['MU'], Element[k]['ETA'], Element[k]['BETA'], Element[k]['GAMMA'], Element[k]['TOL'], Element[k]['FY'], Element[k]['K0'], Element[k]['ALPHA1'], Element[k]['ALPHA2'], Element[k]['DIMENSION'], Element[k]['DIRECTION']))
        elif  name == 'UNXBOUCWEN3DLINK':
            nParaview += 3
            SeismoVLabfile.write("%1.5f %1.5f %1.5f %1.5f %1.5f %E %E %E %1.5f %1.5f %d %d\n" %(Element[k]['ALPHA'], Element[k]['MU'], Element[k]['ETA'], Element[k]['BETA'], Element[k]['GAMMA'], Element[k]['TOL'], Element[k]['FY'], Element[k]['K0'], Element[k]['ALPHA1'], Element[k]['ALPHA2'], Element[k]['DIMENSION'], Element[k]['DIRECTION']))
        elif  name == 'HDRBYAMAMOTO2DLINK':
            nParaview += 3
            SeismoVLabfile.write("%1.5f %1.5f %1.5f %d\n" %(Element[k]['DE'], Element[k]['DI'], Element[k]['HR'], Element[k]['DIMENSION']))
        elif  name == 'HDRBYAMAMOTO3DLINK':
            nParaview += 3
            SeismoVLabfile.write("%1.5f %1.5f %1.5f %d\n" %(Element[k]['DE'], Element[k]['DI'], Element[k]['HR'], Element[k]['DIMENSION']))
        elif name == 'EQLIN2DQUAD4':
            nParaview +=5
            if 'NPOINTS' not in Element[k]:
                Element[k]['NPOINTS'] = 4
            if Element[k]['NPOINTS'] < 4:
                Element[k]['NPOINTS'] = 4
            SeismoVLabfile.write("%d %1.5f %d %s %1.5f %1.5f %1.5f \n" %(Element[k]['MATERIAL'], Element[k]['THICKNESS'], Element[k]['NPOINTS'], Element[k]['TYPE'], Element[k]['ZREF'], Element[k]['CF1'], Element[k]['CF2']))
        elif name == 'TIEQLIN2DQUAD4':
            nParaview +=5
            if 'NPOINTS' not in Element[k]:
                Element[k]['NPOINTS'] = 4
            if Element[k]['NPOINTS'] < 4:
                Element[k]['NPOINTS'] = 4
            SeismoVLabfile.write("%d %1.5f %d %s %1.5f %1.5f %1.5f %1.10f \n" %(Element[k]['MATERIAL'], Element[k]['THICKNESS'], Element[k]['NPOINTS'], Element[k]['TYPE'], Element[k]['ZREF'], Element[k]['CF1'], Element[k]['CF2'], Element[k]['EREF']))
        elif  name == 'NULL2DFRAME2':
            nParaview += 3
            SeismoVLabfile.write("\n")
        elif  name == 'NULL3DFRAME2':
            nParaview += 3
            SeismoVLabfile.write("\n")
        else:
            print('\x1B[31mERROR \x1B[0m: The ELEMENT=' + Element[k]['NAME'].upper() + ' could not be recognized.')
            sys.exit(-1)

    User['NPARAVIEW']  = nParaview
    User['NFEATURES'] += nParaview

def WriteParaviewFile(User, Point, Element, Material, Section):
    """
    This function writes the Paraview mesh file. This file displays
    properties such as Partition, Sections and Materials.

    Parameters
    ----------
    User : dict 
           The dictionary containing all the user's relevant information.
    Point : dict
           The dictionary containing all point information.
    Element : dict
           The dictionary containing all element information.
    Material : dict
           The dictionary containing all material information.
    Section : dict
           The dictionary containing all section information.

    Output
    -------
    file : User['MODEL'].$.vtk
           Finite element mesh model.
    """
    #WRITES THE PARTITIONED MESH FILE
    dirName = User['FOLDER'] + 'Paraview/'
    path    = dirName + User['MODEL'] + '.vtk'

    #Loads the Mesh Partition
    Partition = mt.GetMetisOutputFile(User)

    #Gets the Paraview Point Tag
    ParaviewMap = dict()
    for nTag in Point:
        ParaviewMap[Point[nTag]['TAG']] = nTag 

    #Opens the Partition File 
    Paraviewfile = open(path, "w+")
    Paraviewfile.write("# vtk DataFile Version 4.0\n")
    Paraviewfile.write("SeismoVLab: FINITE ELEMENT MODEL\n")
    Paraviewfile.write("ASCII\n")
    Paraviewfile.write("DATASET UNSTRUCTURED_GRID\n\n")

    #Point Coordinates.
    Paraviewfile.write("POINTS %d float\n" % len(Point))
    if User['DIMENSION'] == 2:
        for nTag in ParaviewMap:
            Coordinates = Point[ParaviewMap[nTag]]['COORDINATES']
            Paraviewfile.write("%1.5f %1.5f 0.00000\n" % (Coordinates[0], Coordinates[1]))

    if User['DIMENSION'] == 3:
        for nTag in ParaviewMap:
            Coordinates = Point[ParaviewMap[nTag]]['COORDINATES']
            Paraviewfile.write("%1.5f %1.5f %1.5f\n" % (Coordinates[0], Coordinates[1], Coordinates[2]))
    Paraviewfile.write("\n")
    
    #Element Connectivities.
    Paraviewfile.write("CELLS %d %d\n" % (len(Element), User['NFEATURES']))
    for eTag in Element:
        connection = Element[eTag]['CONNECTIVITY']
        Paraviewfile.write("%d" % len(connection))
        for nTag in connection:
            Paraviewfile.write(" %d" % Point[nTag]['TAG'])
        Paraviewfile.write("\n")
    Paraviewfile.write("\n")

    #Element Rendering Type.
    Paraviewfile.write("CELL_TYPES %d\n" % len(Element))
    for eTag in Element:
        if Element[eTag]['NAME'].upper() == 'ZEROLENGTH1D':
            Paraviewfile.write("3\n")
        elif Element[eTag]['NAME'].upper() == 'LIN2DTRUSS2':
            Paraviewfile.write("3\n")
        elif Element[eTag]['NAME'].upper() == 'KIN2DTRUSS2':
            Paraviewfile.write("3\n")
        elif Element[eTag]['NAME'].upper() == 'LIN3DTRUSS2':
            Paraviewfile.write("3\n")
        elif Element[eTag]['NAME'].upper() == 'KIN3DTRUSS2':
            Paraviewfile.write("3\n")
        elif Element[eTag]['NAME'].upper() == 'LIN3DTRUSS3':
            Paraviewfile.write("21\n")
        elif Element[eTag]['NAME'].upper() == 'LIN2DQUAD4':
            Paraviewfile.write("9\n")
        elif Element[eTag]['NAME'].upper() == 'LIN2DQUAD8':
            Paraviewfile.write("23\n")
        elif Element[eTag]['NAME'].upper() == 'PML2DQUAD4':
            Paraviewfile.write("9\n")
        elif Element[eTag]['NAME'].upper() == 'PML2DQUAD8':
            Paraviewfile.write("23\n")
        elif Element[eTag]['NAME'].upper() == 'LIN2DFRAME2':
            Paraviewfile.write("3\n")
        elif Element[eTag]['NAME'].upper() == 'KIN2DFRAME2':
            Paraviewfile.write("3\n")
        elif Element[eTag]['NAME'].upper() == 'LIN3DFRAME2':
            Paraviewfile.write("3\n")
        elif Element[eTag]['NAME'].upper() == 'LIN3DSHELL4':
            Paraviewfile.write("9\n")
        elif Element[eTag]['NAME'].upper() == 'LIN3DHEXA8':
            Paraviewfile.write("12\n")
        elif Element[eTag]['NAME'].upper() == 'LIN3DHEXA20':
            Paraviewfile.write("25\n")
        elif Element[eTag]['NAME'].upper() == 'PML3DHEXA8':
            Paraviewfile.write("12\n")
        elif Element[eTag]['NAME'].upper() == 'PML3DHEXA20':
            Paraviewfile.write("25\n")
        elif Element[eTag]['NAME'].upper() == 'UNXBOUCWEN2DLINK':
            Paraviewfile.write("3\n")
        elif Element[eTag]['NAME'].upper() == 'UNXBOUCWEN3DLINK':
            Paraviewfile.write("3\n")
        elif Element[eTag]['NAME'].upper() == 'HDRBYAMAMOTO2DLINK':
            Paraviewfile.write("3\n")
        elif Element[eTag]['NAME'].upper() == 'HDRBYAMAMOTO3DLINK':
            Paraviewfile.write("3\n")
        elif Element[eTag]['NAME'].upper() == 'EQLIN2DQUAD4':
            Paraviewfile.write("9\n")
        elif Element[eTag]['NAME'].upper() == 'TIEQLIN2DQUAD4':
            Paraviewfile.write("9\n")
        elif Element[eTag]['NAME'].upper() == 'NULL2DFRAME2':
            Paraviewfile.write("3\n")
        elif Element[eTag]['NAME'].upper() == 'NULL3DFRAME2':
            Paraviewfile.write("3\n")
    Paraviewfile.write("\n")

    #Scalar Attributes applied to Element.
    Paraviewfile.write("CELL_DATA %d\n" % len(Element))
    Paraviewfile.write("FIELD attributes 3\n")

    #Mesh Partition
    Paraviewfile.write("Partition 1 %d int\n" % len(Element))
    for ind in Partition:
        Paraviewfile.write("%d\n" % ind)
    Paraviewfile.write("\n")

    #Element's Material
    Paraviewfile.write("Materials 1 %d int\n" % len(Element))
    for eTag in Element:
        if 'MATERIAL' in Element[eTag]:
            Paraviewfile.write("%d\n" % Element[eTag]['MATERIAL'])
        elif 'SECTION' in Element[eTag]:
            sTag = Element[eTag]['SECTION']
            mTag = Section[sTag]['MATERIAL']
            Paraviewfile.write("%d\n" % mTag)
        else:
            Paraviewfile.write("-1\n")
    Paraviewfile.write("\n")

    #Element's Section
    Paraviewfile.write("Sections 1 %d int\n" % len(Element))
    for eTag in Element:
        if 'SECTION' in Element[eTag]:
            Paraviewfile.write("%d\n" % Element[eTag]['SECTION'])
        else:
            Paraviewfile.write("-1\n")
    Paraviewfile.write("\n")

    Paraviewfile.close()    

def WriteMeshPartition(User, Point, Mass, Element, Material, Section, Constraint, Damping, Surface, Support, Initial, Function, Load, Combination, Recorder, Simulation):
    """
    This function writes the SeismoVLab (*.svl) file. It first reads the METIS 
    Graph output files, then select the domain partition elements, and finally 
    writes the partition files (name.$.svl) in SeismoVLab format.

    Parameters
    ----------
    User : dict 
           The dictionary containing all the user's relevant information.
    Point : dict
           The dictionary containing all point information.
    Mass  : dict
           The dictionary containing all point mass information.
    Element : dict
           The dictionary containing all element information.
    Material : dict
           The dictionary containing all material information.
    Section : dict
           The dictionary containing all section information.
    Constraint : dict
           The dictionary containing all constraints information.
    Damping : dict
           The dictionary containing all damping models.
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

    Output
    -------
    file : User['MODEL'].$.mesh
           Finite element mesh partition.
    file : User['MODEL'].$.run
           Finite element analysis partition.
    """
    #Creates the METIS Graph Input Files.
    mt.SetMetisInputFile(User, Point, Element, Constraint)

    #Computes the mesh partition using METIS Graph Partitioning.
    Partition = mt.GetMetisOutputFile(User)

    #Gets the Element Tag
    count = 0
    eTags = np.zeros(len(Partition), dtype=np.uint32)
    for k in Element:
        eTags[count] = k
        count += 1
        
    #Writes the Partitioned Mesh Domain.
    for k in range(User['NPART']):
        #Element that belong to this partition
        elemSubdomain = eTags[Partition == k]

        #Check if the partition has Element
        if len(elemSubdomain) == 0:
            print('\x1B[31mERROR \x1B[0m: The partition [' + str(k) + '] does not have Element.')
            sys.exit(-1)

        #Nodes that belong to this partition
        nodeSubdomain = set() 
        for eTag in elemSubdomain:
            connection = Element[eTag]['CONNECTIVITY']
            for nTag in connection:
                nodeSubdomain.add(nTag) 

        #Materials that belong to this partition
        matSubdomain = set()
        for eTag in elemSubdomain:
            if 'MATERIAL' in Element[eTag]:
                matSubdomain.add(Element[eTag]['MATERIAL']) 

        #Sections that belong to this partition
        secSubdomain = set() 
        for eTag in elemSubdomain:
            if 'SECTION' in Element[eTag]:
                sTag = Element[eTag]['SECTION']
                mTag = Section[sTag]['MATERIAL']
                matSubdomain.add(mTag) 
                secSubdomain.add(sTag) 

        #Constraints (Equal, General, Diaphragm) that belong to this partition
        conSubdomain = set() 
        for nTag in nodeSubdomain:
            FreeDofs = Point[nTag]['FREE']
            for dof in FreeDofs:
                if dof < -1:
                    conSubdomain.add(dof)

        #Node Constraints information must be contained in this partition
        for cTag in conSubdomain:
            Master = Constraint[cTag]['MASTERNODE']
            for mNode in Master:
                nodeSubdomain.add(mNode)

        #########
        #print("Elems : ", len(elemSubdomain),", Nodes : ", len(nodeSubdomain), ", Mat : ", len(matSubdomain), ", Secs : ", len(secSubdomain), ", constrs : ", len(conSubdomain))
        #########

        #WRITES THE PARTITIONED MESH FILE
        path  = User['FOLDER'] + 'Partition/' + User['MODEL'] + '.' + str(k) + '.svl'
        
        #Opens the Partition File 
        SeismoVLabfile = open(path, "w+")
        SeismoVLabfile.write("MODEL %d %d %d %s\n" % (User['DIMENSION'], User['NTOTALDOF'], User['NFREEDOF'], User['MASSFORM'].upper()))

        WriteMaterialPartition(SeismoVLabfile, matSubdomain, Material)
        WriteSectionsPartition(SeismoVLabfile, secSubdomain, Section)
        WriteNodePartition(SeismoVLabfile, nodeSubdomain, Point)
        WriteMassPartition(SeismoVLabfile, nodeSubdomain, Mass)
        WriteInitialPartition(SeismoVLabfile, nodeSubdomain, Initial)
        WriteSupportPartition(SeismoVLabfile, nodeSubdomain, Support)
        WriteConstraintsPartition(SeismoVLabfile, conSubdomain, Constraint, Point)
        WriteElementPartition(SeismoVLabfile, elemSubdomain, Element, User)
        WriteDampingsPartition(SeismoVLabfile, Damping, elemSubdomain, User)

        loadSubdomain = WriteLoadsPartition(SeismoVLabfile, Load, Function, Surface, nodeSubdomain, elemSubdomain)

        WriteLoadCombosPartition(SeismoVLabfile, loadSubdomain, Combination)
        WriteRecordersPartition(SeismoVLabfile, User, Recorder, nodeSubdomain, elemSubdomain, k)
        WriteSimulationPartition(SeismoVLabfile, User, Simulation, k)
        SeismoVLabfile.close()

    #THE GENERATED PARTITION NAME
    User['SVL' ] = User['MODEL'] + '.$.svl'

    #WRITES THE MESH IN PARAVIEW FORMAT
    WriteParaviewFile(User, Point, Element, Material, Section)

    #GENERATES SEISMOVLAB EXECUTION COMMAND LINE 
    file = User['SVL'].replace(" ", "~")
    path = User['FOLDER'].replace(" ", "~")
    if User['NPART'] == 1:
        User['EXE'] = ' ./SeismoVLab.exe -dir ' + path + 'Partition' + ' -file ' + file + '\n'
    else:
        User['EXE'] = ' mpirun -np ' + str(User['NPART']) + ' ./SeismoVLab.exe -dir ' + path + 'Partition' + ' -file ' + file + '\n'

    #CLEANS METIS GENERATED FILES
    os.remove(User['FOLDER'] + 'Partition/Graph.out')
    os.remove(User['FOLDER'] + 'Partition/Graph.out.epart.'+ str(User['NPART']))
    os.remove(User['FOLDER'] + 'Partition/Graph.out.npart.'+ str(User['NPART']))
