#!/usr/bin/python3
# -*- coding: Utf-8 -*-

import os
import numpy as np
from scipy import integrate
from Core.Utilities import *
from Core.Definitions import *

def GetDerivative(f, dx):
    """
    This function computes the numerical derivative of a time serie\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2021

    Parameters
    ----------
    f  : array
        Vector (function) that contains the time serie values
    dx : float 
        The time step for the given time serie 

    Returns
    -------
    dfdx : array 
        The derivative of the given function
    """
    dfdx = np.gradient(f, dx)
    return dfdx

def GetIntegration(f, dx):
    """
    This function computes the numerical integration of a time serie\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2021

    Parameters
    ----------
    f  : array
        Vector (function) that contains the time serie values
    dx : float 
        The time step for the given time serie 

    Returns
    -------
    Ix : array 
        The integrated values of the given function
    """
    n  = len(f)
    x  = np.linspace(0.0, (n-1)*dx, n)
    Ix = integrate.cumtrapz(f, x, initial=0)
    return Ix

def GenerateTimeSeries(disp, vels, accel, dt, option):
    """
    This function computes the missing time serie for the displacement, 
    velocity and acceleration depending on the option provied by the user\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2021

    Parameters
    ----------
    disp   : array
        The provided displacement time serie
    vels   : array
        The provided velocity time serie
    accel  : array
        The provided acceleration time serie
    dt     : float
        The time step for the given time serie 
    option : str
        User's time serie data, option=DISP, VEL, ACCEL, or ALL

    Returns
    -------
    None
    """
    if  option.upper() == 'DISP':
        vels  = GetDerivative(disp, dt)
        accel = GetDerivative(vels, dt)
    elif option.upper() == 'VEL':
        disp  = GetIntegration(vels, dt)
        accel = GetDerivative (vels, dt)
    elif option.upper() == 'ACCEL':
        vels  = GetIntegration(accel, dt)
        accel = GetIntegration(vels, dt)

def ParseDRMFile(Function):
    """
    This function parses the DRM file for plane-wave provided in the Entities['Functions'].  
    The routine reads the displacement, velocity, or accelertion input signal depending on 
    the option=ALL,DISP,VEL,ACCEL provied (at the header) by the user.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2021

    Parameters
    ----------
    Function  : dict
        Contains the information (attributes) to generate the plane-wave

    Returns
    -------
    nodes  : array
        The DRM nodes list provided by the user
    conds  : array
        The DRM node condition (0: interior, 1: exterior)
    time   : array
        The time vector
    disp   : array
        The displacement time serie (generated or provided)
    vels   : array
        The velocity time serie (generated or provided)
    accel  : array
        The acceleration time serie (generated or provided)
    dt     : float
        The time step for the given time serie 
    option : str
        User's time serie data, option=DISP, VEL, ACCEL, or ALL
    """
    #Open the provided file
    path = Function['attributes']['file']

    with open(path, "r") as fileHandler:  
        #Loop over each line in file.
        lines = fileHandler.readlines()

        #Parse the first line.
        line = list(filter(None, lines[0].strip().split()))
        nt, dt, nn, option = int(line[0]), float(line[1]), int(line[2]), line[3]

        time = np.zeros(nt)
        disp = np.zeros(nt)
        vels = np.zeros(nt)
        accel = np.zeros(nt)
        nodes = np.zeros(nn, dtype=int)
        conds = np.zeros(nn, dtype=int)

        #Parse DRM node information.
        for m,k in enumerate(range(1, nn+1)):
            line = list(filter(None, lines[k].strip().split()))
            nodes[m], conds[m] = int(line[0]), int(line[1])

        #Parse DRM input signal information.
        if option.upper() == 'ALL':
            for m,k in enumerate(range(nn+1, nn+nt+1)):
                line = list(filter(None, lines[k].strip().split()))
                time[m], disp[m], vels[m], accel[m] = float(line[0]), float(line[1]), float(line[2]), float(line[3])
        elif option.upper() == 'DISP':
            for m,k in enumerate(range(nn+1, nn+nt+1)):
                line = list(filter(None, lines[k].strip().split()))
                time[m], disp[m] = float(line[0]), float(line[1])
        elif option.upper() == 'VEL':
            for m,k in enumerate(range(nn+1, nn+nt+1)):
                line = list(filter(None, lines[k].strip().split()))
                time[m], vels[m] = float(line[0]), float(line[1])
        elif option.upper() == 'ACCEL':
            for m,k in enumerate(range(nn+1, nn+nt+1)):
                line = list(filter(None, lines[k].strip().split()))
                time[m], accel[m] = float(line[0]), float(line[1])
    return nodes, conds, time, disp, vels, accel, dt, option

def WriteDRMFile(filepath, filename, fTag, Disp, Vels, Accel, nt, nc, n, option):
    """
    This function writes the nodal DRM information to a *.drm file.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2021

    Parameters
    ----------
    filepath  : str
        The path where the DRM file will be written
    filename  : str
        The DRM file name
    fTag  : int
        The function tag in Entities
    Disp   : array
        The displacement time serie (generated or provided)
    Vels   : array
        The velocity time serie (generated or provided)
    Accel  : array
        The acceleration time serie (generated or provided)
    nt  : int
        The number of time steps in the time serie
    nc  : int
        The number of DRM columns to be written (2D: 6, 3D: 9)
    n  : int
        The DRM Node tag
    option : str
        If the DRM node is interior (0) or exterior (1)

    Returns
    -------
    None
    """
    #The output file.
    path = filepath + "/" + filename + "-" + str(fTag) + "." + str(n) + ".drm"

    #Domain Reduction file format.
    DRMfile = open(path, "w+")
    DRMfile.write("%d %d %d\n" % (nt, nc, option))

    if nc == 6:
        for k in range(nt):
            DRMfile.write("%E %E %E %E %E %E\n" % (Disp[k,0], Disp[k,1], Vels[k,0], Vels[k,1], Accel[k,0], Accel[k,1]))
    elif nc == 9:
        for k in range(nt):
            DRMfile.write("%E %E %E %E %E %E %E %E %E\n" % (Disp[k,0], Disp[k,1], Disp[k,2], Vels[k,0], Vels[k,1], Vels[k,2], Accel[k,0], Accel[k,1], Accel[k,2]))
    DRMfile.close()

def SVbackground2Dfield(Values, t, X, X0, Xmin, nt):
    """
    """
    #Compute the 2D Field Components.
    U = np.zeros(nt)
    V = np.zeros(nt)

    return U,V

def RHbackground2Dfield(Values, t, X, X0, Xmin, nt):
    """
    """
    #Compute the 2D Field Components.
    U = np.zeros(nt)
    V = np.zeros(nt)
    
    return U,V

def SHbackground3Dfield(Values, t, X, X0, Xmin, nt):
    """
    """
    #Compute the 3D Field Components.
    U = np.zeros(nt)
    V = np.zeros(nt)
    W = np.zeros(nt)

    return U, V, W

def SVbackground3Dfield(Values, t, X, X0, Xmin, nt):
    """
    """
    #Compute the 3D Field Components.
    U = np.zeros(nt)
    V = np.zeros(nt)
    W = np.zeros(nt)

    return U, V, W

def RHbackground3Dfield(Values, t, X, X0, Xmin, nt):
    """
    """
    #Compute the 3D Field Components.
    U = np.zeros(nt)
    V = np.zeros(nt)
    W = np.zeros(nt)

    return U, V, W

def GenerateDRMFiles():
    """
    This function generates the domain reduction files for a plane-wave case.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2021

    Returns
    -------
    None
    """
    for lTag in Entities['Loads']:
        if 'fun' in Entities['Loads'][lTag]['attributes']:
            fTag = Entities['Loads'][lTag]['attributes']['fun']
            if Entities['Loads'][lTag]['attributes']['type'].upper() == 'PLANEWAVE':
                #Gets the Domain Reduction Method Information.
                nodes, conditions, time, Disp, Vels, Accel, dt, option = ParseDRMFile(Entities['Functions'][fTag])

                #Computes the input displacement, velocities and accelerations.
                GenerateTimeSeries(Disp, Vels, Accel, dt, option)

                #Creates the DRM Directory
                dirName = Options['path'] + '/' + 'DRM'
                if not os.path.exists(dirName):
                    os.mkdir(dirName)

                #Computes and Generates the Domain Reduction files for each node.
                x0 = Entities['Functions'][fTag]['attributes']['x0']
                xmin = Entities['Functions'][fTag]['attributes']['xmin']
                funName = Entities['Functions'][fTag]['name']
                funOption = Entities['Functions'][fTag]['attributes']['option']

                nt = len(time)
                nd = Options['dimension']
                U = np.zeros((nt, nd))
                V = np.zeros((nt, nd))
                A = np.zeros((nt, nd))
                '''
                if nd == 2:
                    if funOption.upper() == 'SV':
                        for k, n in enumerate(nodes):
                            x = Entities['Nodes'][n]['coords']
                            #Compute the SV-wave DRM fields
                            U[:,0] , U[:,1] = SVbackground2Dfield(Disp, time, x, x0, xmin, nt)
                            V[:,0] , V[:,1] = SVbackground2Dfield(Vels, time, x, x0, xmin, nt)
                            A[:,0] , A[:,1] = SVbackground2Dfield(Accel, time, x, x0, xmin, nt)
                            WriteDRMFile(dirName, funName, fTag, U, V, A, nt, 6, n, conditions[k])
                    elif funOption.upper() == 'Rayleigh':
                        for k, n in enumerate(nodes):
                            x = Entities['Nodes'][n]['coords']
                            #Compute the Rayleigh DRM fields
                            U[:,0] , U[:,1] = RHbackground2Dfield(Disp, time, x, x0, xmin, nt)
                            V[:,0] , V[:,1] = RHbackground2Dfield(Vels, time, x, x0, xmin, nt)
                            A[:,0] , A[:,1] = RHbackground2Dfield(Accel, time, x, x0, xmin, nt)
                            WriteDRMFile(dirName, funName, fTag, U, V, A, nt, 6, n, conditions[k])
                    else:
                        print('\x1B[31m ERROR \x1B[0m: The specified PLANEWAVE (2D) option (=%s) is not recognized' % funOption)
                elif nd == 3:
                    if funOption.upper() == 'SH':
                        for k, n in enumerate(nodes):
                            x = Entities['Nodes'][n]['coords']
                            #Compute the SH-wave DRM fields
                            U[:,0] , U[:,1] , U[:,2] = SHbackground3Dfield(Disp, time, x, x0, xmin, nt)
                            V[:,0] , V[:,1] , V[:,2] = SHbackground3Dfield(Vels, time, x, x0, xmin, nt)
                            A[:,0] , A[:,1] , A[:,2] = SHbackground3Dfield(Accel, time, x, x0, xmin, nt)
                            WriteDRMFile(dirName, funName, fTag, U, V, A, nt, 9, n, conditions[k])
                    elif funOption.upper() == 'SV':
                        for k, n in enumerate(nodes):
                            x = Entities['Nodes'][n]['coords']
                            #Compute the SV-wave DRM fields
                            U[:,0] , U[:,1] , U[:,2] = SVbackground3Dfield(Disp, time, x, x0, xmin, nt)
                            V[:,0] , V[:,1] , V[:,2] = SVbackground3Dfield(Vels, time, x, x0, xmin, nt)
                            A[:,0] , A[:,1] , A[:,2] = SVbackground3Dfield(Accel, time, x, x0, xmin, nt)
                            WriteDRMFile(dirName, funName, fTag, U, V, A, nt, 9, n, conditions[k])
                    elif funOption.upper() == 'Rayleigh':
                        for k, n in enumerate(nodes):
                            x = Entities['Nodes'][n]['coords']
                            #Compute the Rayleigh DRM fields
                            U[:,0] , U[:,1] , U[:,2] = RHbackground3Dfield(Disp, time, x, x0, xmin, nt)
                            V[:,0] , V[:,1] , V[:,2] = RHbackground3Dfield(Vels, time, x, x0, xmin, nt)
                            A[:,0] , A[:,1] , A[:,2] = RHbackground3Dfield(Accel, time, x, x0, xmin, nt)
                            WriteDRMFile(dirName, funName, fTag, U, V, A, nt, 9, n, conditions[k])
                    else:
                        print('\x1B[31m ERROR \x1B[0m: The specified PLANEWAVE (3D) option (=%s) is not recognized' % funOption)
                else:
                    print('\x1B[31m ERROR \x1B[0m: The specified dimension (=%d) is not possible for DRM' % Options['dimensions'])
                '''

                #Update the domain reduction time series path where files are located.
                Entities['Functions'][fTag]['attributes']['file'] = dirName + "/" + funName + "-" + str(fTag) + ".$.drm"