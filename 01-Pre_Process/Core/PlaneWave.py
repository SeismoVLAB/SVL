#!/usr/bin/python3
# -*- coding: Utf-8 -*-

import os
import time
import numpy as np
from scipy import signal
import concurrent.futures

from scipy import interpolate

import scipy.sparse as sps
from scipy.sparse import linalg as sla

from scipy import integrate
from Core.Utilities import *
from Core.Definitions import *

def GetDerivative(f, dx):
    """
    This function computes the numerical derivative of a time series\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2021

    Parameters
    ----------
    f  : array
        Vector (function) that contains the time series values
    dx : float 
        The time step for the given time series 

    Returns
    -------
    dfdx : array 
        The derivative of the given function
    """
    dfdx = np.gradient(f, dx)
    return dfdx

def GetIntegration(f, dx):
    """
    This function computes the numerical integration of a time series\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2021

    Parameters
    ----------
    f  : array
        Vector (function) that contains the time series values
    dx : float 
        The time step for the given time series 

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
    This function computes the missing time series for the displacement, 
    velocity and acceleration depending on the option provied by the user\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2021

    Parameters
    ----------
    disp   : array
        The provided displacement time series
    vels   : array
        The provided velocity time series
    accel  : array
        The provided acceleration time series
    dt     : float
        The time step for the given time series 
    option : str
        User's time series data, option=DISP, VEL, ACCEL, or ALL

    Returns
    -------
    Disp, Vels, Accel the time series
    """
    if  option.upper() == 'DISP':
        Disp  = disp
        Vels  = GetDerivative(disp, dt)
        Accel = GetDerivative(Vels, dt)
    elif option.upper() == 'VEL':
        Disp  = GetIntegration(vels, dt)
        Vels  = vels
        Accel = GetDerivative (vels, dt)
    elif option.upper() == 'ACCEL':
        Vels  = GetIntegration(accel, dt)
        Disp  = GetIntegration(Vels, dt)
        Accel = accel

    return Disp, Vels, Accel

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
        The displacement time series (generated or provided)
    vels   : array
        The velocity time series (generated or provided)
    accel  : array
        The acceleration time series (generated or provided)
    dt     : float
        The time step for the given time series 
    option : str
        User's time series data, option=DISP, VEL, ACCEL, or ALL
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

def Ricker(parameters, option=""):
    """
    This function generates a Ricker pulse provided with some parameters.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2021

    Parameters
    ----------
    parameters : dict
        'to'  : float
            Time where the velocity is maximum
        'f0'  : float
            Central frequency of the Ricker pulse
        'dt'  : float
            Time step to discretize the signal
        'Ts'  : float
            Duration of the time signal
        'Ap'  : float
            Amplitude of the Ricker pulse for velocity
    option: str
        string that specify the return vectors: 'DISP', 'VEL', 'ACCEL'

    Returns
    -------
    disp  : array
        Displacement time series
    Vels  : array
        Velocity time series
    accel  : array
        Acceleration time series
    t  : array
        Vector with time values 
    """
    #Unpack Ricker parameters
    to = parameters['to']
    f0 = parameters['f0']
    dt = parameters['dt']
    Ts = parameters['Ts']
    Ap = parameters['Ap']

    #Time domain discretization
    nt = int(Ts/dt + 1)
    t = np.linspace(0.0, Ts, nt)

    #Input Signal:
    beta = np.square(np.pi*f0*(t - to))
    vels = Ap*np.multiply(1 - 2.0*beta, np.exp(-beta))
    disp = integrate.cumtrapz(vels, t, initial=0) 
    accel = -Ap*np.multiply(np.multiply(2*np.pi**2*f0**2*(t - to), 2 + (1 - 2*beta)), np.exp(-beta))

    if option.upper() == 'DISP':
        return t, disp
    elif option.upper() == 'VEL':
        return t, vels
    elif option.upper() == 'ACCEL':
        return t, accel

    return t, disp, vels, accel

def WritePlaneWaveFile(DRM):
    """
    """
    #Unpack DRM information
    dt   = DRM['dt']
    t    = DRM['t']
    vals = DRM['signal']
    field = DRM['field']
    intDRM = DRM['Interior']
    extDRM = DRM['Exterior']
    filename = DRM['filename']

    #The file name
    filepath = Options['path'] + '/' + filename

    #Writes the DRM Plane-Wave file
    nt = len(t)
    nn = len(intDRM) + len(extDRM)
    with open(filepath, "w+") as DRMfile: 
        DRMfile.write("%d %E %d %s\n" % (nt, dt, nn, field.upper()))

        #Interior DRM Nodes
        for k in intDRM:
            DRMfile.write("%d 0\n" % k)

        #Exterior DRM Nodes
        for k in extDRM:
            DRMfile.write("%d 1\n" % k)

        if field.upper() == 'ALL':
            for k in range(nt):
                DRMfile.write("%E %E %E %E\n" % (t[k], vals[0,k], vals[1,k], vals[2,k]))
        else:
            for k in range(nt):
                DRMfile.write("%E %E\n" % (t[k], vals[k]))

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
        The displacement time series (generated or provided)
    Vels   : array
        The velocity time series (generated or provided)
    Accel  : array
        The acceleration time series (generated or provided)
    nt  : int
        The number of time steps in the time series
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

def GetKofLayer(k,p,s,h,mu,aSP):
    """
    This function calculates the K00 and K01 components of the stiffness matrix
    of a layer with finite thickness.\n
    
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Kien T. Nguyen 2021, ORCID: 0000-0001-5761-3156
    
    Parameters
    ----------
    k  : float
        The horizontal wavenumber (along x-axis)
    p, s  : complex
        Complex coefficients
    h  : float
        The thickness of soil layer
    mu  : float
        The shear modulus of soil
    aSP  : float
        The ratio between shear wave velocity and dilatational wave velocity

    Returns
    -------
    K00  : array
        2x2 matrix, block component of stiffness matrix of a layer 
    K01  : array
        2x2 matrix, block component of stiffness matrix of a layer 
    """
    K00 = np.zeros((2,2), dtype=complex)
    K01 = np.zeros((2,2), dtype=complex)
    
    a = np.real(k*p*h)
    b = np.imag(k*p*h)
    c = np.real(k*s*h)
    d = np.imag(k*s*h)
    
    cb = np.cos(b)
    sb = np.sin(b)
    
    plusA  = 1.0/2.0*(1.0+np.exp(-2.0*a))
    minusA = 1.0/2.0*(1.0-np.exp(-2.0*a))
    
    C1 = plusA*cb + 1j*minusA*sb
    S1 = minusA*cb + 1j*plusA*sb
    
    cd = np.cos(d)
    sd = np.sin(d)
    
    plusC  = 1.0/2.0*(1.0+np.exp(-2.0*c))
    minusC = 1.0/2.0*(1.0-np.exp(-2.0*c))
    
    C2 = plusC*cd + 1j*minusC*sd
    S2 = minusC*cd + 1j*plusC*sd
    
    D0 = 2.0*(np.exp(-a-c)-C1*C2)+(1.0/p/s+p*s)*S1*S2
    
    K00[0,0] = (1.0-s*s)/2.0/s*(C1*S2-p*s*C2*S1)/D0
    K00[0,1] = (1.0-s*s)/2.0*(np.exp(-a-c)-C1*C2+p*s*S1*S2)/D0 + (1.0+s*s)/2.0
    K00[1,0] = K00[0,1]
    K00[1,1] = (1.0-s*s)/2.0/p*(C2*S1-p*s*C1*S2)/D0
    K00 *= 2.0*k*mu
    
    K01[0,0] = 1.0/s*(p*s*S1*np.exp(-c) - S2*np.exp(-a))/D0
    K01[0,1] = (C1*np.exp(-c) - C2*np.exp(-a))/D0
    K01[1,0] = -K01[0,1]
    K01[1,1] = 1.0/p*(p*s*S2*np.exp(-a) - S1*np.exp(-c))/D0
    K01 *= 2.0*k*mu*(1.0-s*s)/2.0
    
    return K00, K01

def GetKofHalfSpace(k,p,s,mu):
    """
    This function calculates the stiffness matrix of half space.\n
    
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Kien T. Nguyen 2021, ORCID: 0000-0001-5761-3156
    
    Parameters
    ----------
    k  : float
        The horizontal wavenumber (along x-axis)
   p, s  : complex
        Complex coefficients
    mu  : float
        The shear modulus of soil

    Returns
    -------
    K  : array
        2x2 matrix, stiffness matrix of half space
    """
   
    K = np.zeros((2,2), dtype=complex)
    coef = (1.0-s*s)/2.0/(1.0-p*s)
    K[0,0] = coef*p
    K[0,1] = -coef + 1.0
    K[1,0] = K[0,1]
    K[1,1] = coef*s
    K *= 2.0*k*mu
    
    return K

def GetKofFullSpace(k,p,s,mu):
    """
    This function calculates the stiffness matrix of full space.\n
    
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Kien T. Nguyen 2021, ORCID: 0000-0001-5761-3156
    
    Parameters
    ----------
    k  : float
        The horizontal wavenumber (along x-axis)
    p, s  : complex
        Complex coefficients
    mu  : float
        The shear modulus of soil

    Returns
    -------
    K  : array
        2x2 matrix, stiffness matrix of full space
    """
    K = np.zeros((2,2), dtype=complex)
    coef = 2.0*k*mu*(1.0-s*s)/(1.0-p*s)
    K[0,0] = coef*p
    K[1,1] = coef*s
    
    return K

def GetDisplacementAtInteriorLayer(y,yTop,yBot,uTop,uBot,k,p,s,mu,aSP):
    """
    This function calculates the displacement at interior points of a layer.\n
    
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Kien T. Nguyen 2021, ORCID: 0000-0001-5761-3156
    
    Parameters
    ----------
    y  : float
        The y-coordinate of the query point
    yTop, yBot  : float
        The y-coordinate of the top and bottom of the parent layer containing query point
    uTop, uBot  : complex
        The displacement at the top and bottom of the parent layer containing query point 
    k  : float
        The horizontal wavenumber (along x-axis)
    p, s  : complex
        Complex coefficients
    h  : float
        The thickness of parent layer
    mu  : float
        The shear modulus of parent layer
    aSP  : float
        The ratio between shear wave velocity and dilatational wave velocity of parent layer
        
    Returns
    -------
    uz  : array
        displacement at the specific query point
    """
    xi = yTop - y
    eta = y - yBot
    [K00xi, K01xi]   = GetKofLayer(k,p,s,xi,mu,aSP) 
    [K00eta, K01eta] = GetKofLayer(k,p,s,eta,mu,aSP) 
    A = K00eta + K00xi*np.array([[1.0,-1.0],[-1.0,1.0]])
    b = -(np.dot(K01xi.T,uTop) + np.dot(K01eta,uBot))
    uz = np.linalg.solve(A, b)
    
    return uz

def DataPreprocessing(Disp, Vels, Accel, Layers, beta, rho, nu, angle, yDRMmin, nt, dt, fun):
    """
    This function performs the data pre-processing for the wave propagation
    problem, such as adding an imaginary layer at the bottom of DRM nodes 
    (if necessary), zero padding, transform wave signal into frequency domain.
    Particularly, it also calculates the motion at the position of half-space surface, 
    assuming that wave is propagating in the full space having same properties 
    as the half space. This full-space motion is subsequently used to 
    generate the force vector while calculating the response of soil interface 
    based on substructure technique.\n
    
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Kien T. Nguyen 2021, ORCID: 0000-0001-5761-3156
    
    Parameters
    ----------
    Disp, Vels, Accel  : array
        Displacement, Velocity, and Acceleration time series of incoming wave
    Layers  : array
        The y-coordinate of soil layer interfaces, from free ground surface to half-space surface 
    beta, rho, nu  : array
        Shear wave velocity, Mass density, and Poisson's ratio of soil layers, from top layer to half space 
    angle  : float
        Angle of incoming wave, with respect to vertical axis (y-axis), from 0 to 90 degrees
    yDRMmin  : float
        The minimum of y-coordinates among DRM nodes
    nt  : int
        The original number of time step of the wave signal
    dt  : float
        The time step increment 
    fun  : dict
        The dictionary that contains the information of DRM
        
    Returns
    -------
    ufull, vfull, afull  : array
        Displacement, Velocity, and Acceleration at the position of half-space surface
        Note that the vertical components are multiplied with -1j
    Layers, beta, rho, nu: array
        Similar as above explaination
    wVec  : array
        Angular frequency spectrum 
    p, s  : array
        Complex coefficients
    h  : array
        The thickness of soil layers
    mu  : array
        The shear modulus of soil layers
    aSP  : array
        The ratio between shear wave velocity and dilatational wave velocity of soil layers
    phaseVelIn  : float
        The phase velocity of incoming wave in the half space underneath
    sinTheta  : float
        Sine of the incoming angle
    N  : int
        Number of soil layers, including imaginary layer and half space
    Nt  : int
        Length of the Disp, Vels, Accel after zero padding
    """
    CutOffFrequency = fun['CutOffFrequency']
    waveType = fun['option']
    df = fun['df']

    #DEALING WITH ANGLE = 0, 90
    if np.isclose(angle, 0.0, rtol=1e-05):
        angle += 0.0001
    elif np.isclose(angle, 90.0, rtol=1e-05):
        angle -= 0.0001

    #ADDING IMAGINARY LAYER AT yDRMmin IF NECESSARY
    if yDRMmin < Layers[-1]:
        Layers = np.append(Layers, np.array([yDRMmin]))
        beta = np.append(beta, [beta[-1]])
        rho = np.append(rho, [rho[-1]])
        nu = np.append(nu, [nu[-1]])

    N = len(Layers)

    #PADDING ZERO TO HAVE DESIRED DISCRETIZED FREQUENCY STEP
    stepsNeeded = max(nt,1.0/dt/df)
    nextPowerOfTwo = int(np.ceil(np.log2(stepsNeeded)))

    Nt = 2**nextPowerOfTwo

    Disp = np.concatenate([Disp, np.zeros(Nt-nt)])
    Vels = np.concatenate([Vels, np.zeros(Nt-nt)])
    Accel = np.concatenate([Accel, np.zeros(Nt-nt)])

    n = np.arange(1,Nt/2+1)
    wVec = 2.0*np.pi/dt/Nt*n
    wVec = wVec[wVec<=2.0*np.pi*CutOffFrequency]

    FdispIn = np.fft.rfft(Disp)
    FdispIn = FdispIn[1:len(wVec)+1]

    FvelIn = np.fft.rfft(Vels)
    FvelIn = FvelIn[1:len(wVec)+1]

    FaccelIn = np.fft.rfft(Accel)
    FaccelIn = FaccelIn[1:len(wVec)+1]

    sinTheta = np.sin(angle/180.0*np.pi)
    cosTheta = np.cos(angle/180.0*np.pi)

    alpha = beta*np.sqrt(2.0*(1.0 - nu)/(1.0 - 2.0*nu))
    aSP = beta/alpha
    mu = rho*beta*beta

    #Polarization angle for P or SV wave
    if waveType=="SV":
        phaseVelIn = beta[-1]
        polarization = np.array([cosTheta, -sinTheta])
    elif waveType=="P":
        phaseVelIn = alpha[-1]
        polarization = np.array([sinTheta, cosTheta])

    p = np.lib.scimath.sqrt(1.0-np.square(phaseVelIn/alpha/sinTheta))
    s = np.lib.scimath.sqrt(1.0-np.square(phaseVelIn/beta/sinTheta))
    h = -np.diff(Layers)

    zrelHalfSpace = Layers[-1] - yDRMmin

    #FULL SPACE PROPAGATION OF SV WAVE IN FREQUENCY DOMAIN, RESULT AT HALF-SPACE INTERFACE y=yN
    ufullx = FdispIn*polarization[0]*np.exp(-1j*wVec*cosTheta/phaseVelIn*zrelHalfSpace)
    ufullz = FdispIn*polarization[1]*np.exp(-1j*wVec*cosTheta/phaseVelIn*zrelHalfSpace)
    ufull = np.vstack((ufullx, -1j*ufullz))

    vfullx = FvelIn*polarization[0]*np.exp(-1j*wVec*cosTheta/phaseVelIn*zrelHalfSpace)
    vfullz = FvelIn*polarization[1]*np.exp(-1j*wVec*cosTheta/phaseVelIn*zrelHalfSpace)
    vfull = np.vstack((vfullx, -1j*vfullz)) 

    afullx = FaccelIn*polarization[0]*np.exp(-1j*wVec*cosTheta/phaseVelIn*zrelHalfSpace)
    afullz = FaccelIn*polarization[1]*np.exp(-1j*wVec*cosTheta/phaseVelIn*zrelHalfSpace)
    afull = np.vstack((afullx, -1j*afullz)) 

    return ufull, vfull, afull, Layers, beta, rho, nu, wVec, p, s, h, mu, aSP, phaseVelIn, sinTheta, N, Nt

def SoilInterfaceResponse(ufull, wVec, p, s, h, mu, aSP, phaseVelIn, sinTheta, N):
    """
    This function calculates the displacement time series at soil interface positions
    in wave propagation problem, in which the SV or P wave incoming from the half space
    underneath under arbitrary incident angle from 0 to 90 degrees and propagating 
    through stratified soil domain. 
    Note: 
    [1] The coordinate system and displacement positive axes:

        y(V) ^
             |
             |
             o-----> x(U)
             
    [2] At each frequency, the horizontal and vertical displacements are calculated 
        based on Eduardo Kausel's Stiffness Matrix Method, in "Fundamental 
        Solutions in Elastodynamics, A Compendium", chap. 10, pp. 140--159 

    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Kien T. Nguyen 2021, ORCID: 0000-0001-5761-3156
    
    Parameters
    ----------
    ufull  : array
        Displacement of full-space problem at the position of half-space surface 
        Note that the vertical components are multiplied with -1j
    wVec  : array
        Angular frequency spectrum 
    p, s  : array
        Complex coefficients
    h  : array
        The thickness of soil layers
    mu  : array
        The shear modulus of soil layers
    aSP  : array
        The ratio between shear wave velocity and dilatational wave velocity of soil layers
    phaseVelIn  : float
        The phase velocity of incoming wave in the half space underneath
    sinTheta  : float
        Sine of the incoming angle
    N  : int
        Number of soil layers, including imaginary layer and half space
        
    Returns
    -------
    uInterface  : array
        Displacements at the soil layer interface positions, in frequency domain
    """
    nfi = len(wVec)
    uInterface = np.zeros((2*N,2*N,nfi), dtype=complex)

    for fi in range(nfi):
        w = wVec[fi]
        k = w*sinTheta/phaseVelIn
        Kglobal = np.zeros((2*N,2*N), dtype=complex)

        #Assemble each layer
        for i in range(N-1):
            [K00, K01] = GetKofLayer(k, p[i], s[i], h[i], mu[i], aSP[i])
            Kglobal[2*i:2*i+2,2*i:2*i+2]     += K00
            Kglobal[2*i:2*i+2,2*i+2:2*i+4]   += K01
            Kglobal[2*i+2:2*i+4,2*i:2*i+2]   += K01.T
            Kglobal[2*i+2:2*i+4,2*i+2:2*i+4] += K00*np.array([[1.0,-1.0],[-1.0,1.0]])

        #Assemble each layer 
        Khalfspace = GetKofHalfSpace(k, p[-1], s[-1], mu[-1])
        Kglobal[2*N-2:2*N,2*N-2:2*N]  += Khalfspace
    
        #Assembel force vector
        forceVec = np.zeros((2*N,1), dtype=complex)
        Kfull = GetKofFullSpace(k, p[-1], s[-1], mu[-1])
        forceVec[2*N-2:2*N,:] = (Kfull.dot(ufull[:,fi])).reshape(2,1) 
    
        #Displacement at interface
        uInterface[:,:,fi] = np.linalg.solve(Kglobal, forceVec)

    return uInterface

def GetRayleighFFTfields(Disp, Vels, Accel, endFrequency, dt, df, nt):
    '''
    '''
    stepsNeeded    = max(nt, 1.0/dt/df)
    nextPowerOfTwo = int(np.ceil(np.log2(stepsNeeded)))
    Nt = 2**nextPowerOfTwo

    Disp  = np.concatenate([Disp, np.zeros(Nt-nt)])
    Vels  = np.concatenate([Vels, np.zeros(Nt-nt)])
    Accel  = np.concatenate([Accel, np.zeros(Nt-nt)])

    n       = np.arange(1, Nt/2+1)
    wVec    = 2.0*np.pi/dt/Nt*n
    wVec    = wVec[wVec <= 2.0*np.pi*endFrequency]

    FFTdisp = np.fft.rfft(Disp)
    FFTdisp = FFTdisp[1:len(wVec)+1]

    FFTvels = np.fft.rfft(Vels)
    FFTvels = FFTvels[1:len(wVec)+1]

    FFTaccel = np.fft.rfft(Accel)
    FFTaccel = FFTaccel[1:len(wVec)+1]

    df = 1.0/dt/Nt
    startFrequency = df
    endFrequency  += 0.001*df #to include the endFre if there is a multiple of df, fDispersion = np.arange(startFre,endFre,df)

    return wVec, FFTdisp, FFTvels, FFTaccel, df, Nt, startFrequency, endFrequency

def Assemble(a,b,k):
    """
    This function appends b to the end of a with k elements overlapped.
    If k<0, abs(k) zeros will be added to the end of a before appending b.\n
    
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Kien T. Nguyen 2021, ORCID: 0000-0001-5761-3156
    
    Parameters
    ----------
    a, b  : ndarray
        ndarray with shape (1,*)
    k  : int
        Number of overlapped elements, k<0 means abs(k) zeros are added to the end of a

    Returns
    -------
    c  : ndarray
        The concatenated ndarray with shape (1,a.size+b.size)
    """
    if a.size==0:
        c = b
    else:
        la = a.size
        lb = b.size
        if k<0:
            c = np.concatenate((a,np.zeros((1,-k)),b),1)
        elif k==0:
            c = np.concatenate((a,b),1)
        else:
            c = np.concatenate((a[:,0:la-k],a[:,la-k:la] + b[:,0:k],b[:,k:lb]),1)
    
    return c
    
def GetRayleighDispersionAndModeShape(mode, intY, beta, rho, nu, dy1, yDRMmin, nepw, startFre, endFre, df, depthFactor):
    """
    This function calculates the dispersion curves and mode shapes of Rayleigh
    wave in stratified soil.\n
    
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Kien T. Nguyen 2021, ORCID: 0000-0001-5761-3156
    
    Parameters
    ----------
    mode  : int
        The chosen Rayleigh mode shape. Currently, only mode=0 for fundamental mode is available
    intY  : array
        The y-coordinate of soil layer interfaces, from free ground surface downwards 
    beta, rho, nu  : array
        Shear wave velocity, Mass density, and Poisson's ratio of soil layers, from top layer to half space 
    dy1  : float
        y-grid spacing of points used to interpolate the mode shape, from ground surface to yDRMmin
    yDRMmin  : float
        Minimum of y-coordinates among DRM nodes
    nepw  : int
        Number of element per wavelength in the extended region below yDRMmin
    startFre, endFre, df  : float
        Start, end, and increment of frequency
    depthFactor  :float
        Soil domain is extended up to depthFactor*wavelength to mimic the clamp condition
        
    Returns
    -------
    fDispersion, phaseVelDispersion  : array
        The frequency and phase velocity of the dispersion curve of chosen Rayleigh mode
    yGridModeShape  : array
        The y-coordinate grid points which is subsequently used to interpolate mode shape at DRM nodes
    uModeShape, vModeShape  : ndarray
        The mode shapes of horizontal and vertical displacements of chosen Rayleigh mode at yGridModeShape, 
        at each frequency fDispersion
    """
    Vr = np.zeros(len(intY))
    for jj in range(len(intY)):
        Vr[jj] = GetRayleighVelocity(beta[jj],nu[jj])
    
    Vlower = Vr.min()
    Vupper = beta[-1]
    
    fDispersion = np.arange(startFre, endFre, df)
    phaseVelDispersion = np.zeros(len(fDispersion))
    modeShape = []
    
    for indexFre, f in enumerate(fDispersion):
        A0 = np.array([[]])
        A2 = np.array([[]])
        B1 = np.array([[]])
        B3 = np.array([[]])
        G0 = np.array([[]])
        G2 = np.array([[]])
        M0 = np.array([[]])
        M2 = np.array([[]])
        
        #Calculate wavelength and element size in the extended region
        wavelengthMax = Vupper/f
        wavelengthMin = Vlower/f
        dy2 = wavelengthMin/nepw
        
        #Add imaginary layer at yDRMmin and at depth depthFactor*wavelength
        depth = min(intY[-1],yDRMmin,intY[0]-depthFactor*wavelengthMax)
        intAux = np.unique(np.concatenate([intY,[yDRMmin,depth]]))
        intYNew = intAux[::-1]

        indices = np.where(np.in1d(intYNew, intY))[0]
        repeats = np.hstack((np.diff(indices),[np.size(intYNew)-indices[-1]]))
        betaNew = np.repeat(beta,repeats)
        rhoNew = np.repeat(rho,repeats)
        nuNew = np.repeat(nu,repeats)
    
        muNew = rhoNew*betaNew*betaNew
        Lame1New = 2.0*muNew*nuNew/(1.0-2.0*nuNew)
        
        N = len(intYNew) #number of interfaces (including imaginary interface at yDRMmin and at depthFactor*wavelength)
        h = -np.diff(intYNew) #thickness of each soil layer
        
        #Mesh density in 2 region, from 0 to yDRMmin and yDRMmin to depthFactor*wavelength
        idx = np.where(intYNew == yDRMmin)[0][0]
        dy = np.concatenate([dy1*np.ones(idx),dy2*np.ones(N-idx-1)])
        ns = np.ceil(h/dy).astype(int) #number of sublayer in each soil layer
        dy = h/ns
        
        for ii in range(N-1):
            [A0e, A2e, B1e, B3e, G0e, G2e, M0e, M2e] = GetLayerStiffnessComponents(ns[ii],dy[ii],Lame1New[ii],muNew[ii],rhoNew[ii])     
            A0 = Assemble(A0,A0e,2)
            A2 = Assemble(A2,A2e,0)
            B1 = Assemble(B1,B1e,1)
            B3 = Assemble(B3,B3e,-1)
            G0 = Assemble(G0,G0e,2)
            G2 = Assemble(G2,G2e,0)
            M0 = Assemble(M0,M0e,2)
            M2 = Assemble(M2,M2e,0)
            
        #Add last 2x2 block for last interface
        A0 = Assemble(A0,A0e[0:1,0:2],2);
        B1 = Assemble(B1,B1e[0:1,0:1],1);
        G0 = Assemble(G0,G0e[0:1,0:2],2);
        M0 = Assemble(M0,M0e[0:1,0:2],2);
    
        Ns = np.sum(ns) #total number of sub layers
        noDoF = 2*Ns+2 #total number of degree of freedom
        A = sps.spdiags(np.vstack((np.concatenate((A2,[[0.0,0.0]]),1),A0,np.concatenate(([[0.0,0.0]],A2),1))),np.array([-2,0,2]),noDoF,noDoF)
        B = sps.spdiags(np.vstack((np.concatenate((B3,[[0.0,0.0,0.0]]),1),np.concatenate((B1,[[0.0]]),1),np.concatenate(([[0.0]],B1),1),np.concatenate(([[0.0,0.0,0.0]],B3),1))),np.array([-3,-1,1,3]),noDoF,noDoF)
        G = sps.spdiags(np.vstack((np.concatenate((G2,[[0.0,0.0]]),1),G0,np.concatenate(([[0.0,0.0]],G2),1))),np.array([-2,0,2]),noDoF,noDoF)
        M = sps.spdiags(np.vstack((np.concatenate((M2,[[0.0,0.0]]),1),M0,np.concatenate(([[0.0,0.0]],M2),1))),np.array([-2,0,2]),noDoF,noDoF)
        
        #Solve for eigenvalues and eigenvectors
        Ngrid = 2*(np.sum(ns[0:idx])+1) #number of degrees of freedom of the predefined grid used for DRM nodes interpolation
        w = 2.0*np.pi*f
        
        Ain = sps.vstack([sps.hstack([sps.csr_matrix((noDoF, noDoF), dtype = np.float), sps.identity(noDoF,dtype=np.float)]), sps.hstack([w*w*M-G, -B])])
        Min = sps.vstack([sps.hstack([sps.identity(noDoF,dtype=np.float),sps.csr_matrix((noDoF, noDoF), dtype = np.float)]), sps.hstack([sps.csr_matrix((noDoF, noDoF), dtype = np.float),A])])
        
        eigk, eigShape = sla.eigs(Ain, k=mode+1, M = Min, sigma=w/Vlower) #note that k is the number of modes desired
        
        wavenumber = np.real(eigk[mode])
        phaseVelDispersion[indexFre] = w/wavenumber
        modeShape.append(eigShape[0:Ngrid,mode])
    
    modeShapeReshape = np.transpose(np.vstack(modeShape)) #reshape the matrix, column is indexFre, row is interlacing of u and v of each point     
    uModeShape = modeShapeReshape[0::2,:]
    vModeShape = modeShapeReshape[1::2,:]
    
    yGridModeShape = [intYNew[0]]
    for nn in range(idx):
        yGridModeShape.append(np.linspace(intYNew[nn]-dy[nn],intYNew[nn+1],ns[nn]))
    yGridModeShape = np.hstack(yGridModeShape)
    
    return fDispersion, phaseVelDispersion, yGridModeShape, uModeShape, vModeShape

def GetLayerStiffnessComponents(ns,h,Lame1,mu,rho):
    """
    This function calculates the diagonal components of the stiffness matrices
    A, B, G, M for 1 soil layer with multiple soil sublayers in Thin Layer Method.\n
    
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Kien T. Nguyen 2021, ORCID: 0000-0001-5761-3156
    
    Parameters
    ----------
    ns  : int
        Number of soil sublayers in the currently-considered soil layer
    h  : float
        Uniform thickness of the sublayer
    Lame1, mu, rho  : float
        The first Lame parameter, shear modulus, and mass density of the soil layer

    Returns
    -------
    A0, A2, B1, B3, G0, G2, M0, M2  : ndarray
        The diagonal components of stiffness matrices A, B, G, M 
        0, 1, and 2 mean main, 1st, and 2nd diagonal 
    """
    
    coef1 = Lame1 + 2.0*mu
    coef2 = Lame1 + mu
    coef3 = Lame1 - mu
    
    A0 = h/6.0*np.concatenate(([[2.0*coef1,2.0*mu]],4.0*np.tile([[coef1,mu]],(1,ns-1)),[[2.0*coef1,2.0*mu]]),1)
    A2 = h/6.0*np.tile([[coef1,mu]],(1,ns))
    
    B1 = 1.0/2.0*np.concatenate(([[coef3,coef2]],np.tile([[0.0,coef2]],(1,ns-1)),[[-coef3]]),1)
    B3 = 1.0/2.0*np.concatenate((np.tile([[-coef2,0.0]],(1,ns-1)),[[-coef2]]),1)
    
    G0 = 1.0/h*np.concatenate(([[mu,coef1]],2.0*np.tile([[mu,coef1]],(1,ns-1)),[[mu,coef1]]),1)
    G2 = 1.0/h*np.tile([[-mu,-coef1]],(1,ns))

    M0 = rho*h/6.0*np.concatenate(([[2.0,2.0]],4.0*np.ones((1,2*ns-2)),[[2.0,2.0]]),1)
    M2 = rho*h/6.0*np.ones((1,2*ns))
    
    return A0, A2, B1, B3, G0, G2, M0, M2

def PSVbackground2Dfield(us, Layers, wVec, p, s, h, mu, aSP, phaseVelIn, sinTheta, N, Nt, x0, x, y):
    """
    This function calculates the displacement time series in 2D wave propagation
    problem, in which the SV or P wave incoming from the half space
    underneath under arbitrary incident angle from 0 to 90 degrees and propagating 
    through stratified soil domain. 
    Note: 
    [1] The coordinate system and displacement positive axes:

        y(V) ^
             |
             |
             o-----> x(U)
             
    [2] At each frequency, the horizontal and vertical displacements are calculated 
        based on Eduardo Kausel's Stiffness Matrix Method, in "Fundamental 
        Solutions in Elastodynamics, A Compendium", chap. 10, pp. 140--159 

    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Kien T. Nguyen 2021, ORCID: 0000-0001-5761-3156
    
    Parameters
    ----------
    us  : array
        Displacement at the layer interfaces in frequency domain
    Layers  : array
        The y-coordinate of soil layer interfaces, including imaginary and half-space interfaces
    wVec  : array
        Angular frequency spectrum 
    p, s  : array
        Complex coefficients
    h  : array
        The thickness of soil layers
    mu  : array
        The shear modulus of soil layers
    aSP  : array
        The ratio between shear wave velocity and dilatational wave velocity of soil layers
    phaseVelIn  : float
        The phase velocity of incoming wave in the half space underneath
    sinTheta  : float
        Sine of the incoming angle
    N  : int
        Number of soil layers, including imaginary layer and half space
    Nt  : int
        Length of the Disp, Vels, Accel after zero padding
    x0  :float
        x-coordinate of the reference point (where the incoming signal time series is prescribed)
    x, y :float
        x- and y-coordinate of the query point
        
    Returns
    -------
    Z  : array
        Displacement time series at the query point
    """
    nfi = len(wVec)
    U_fft = np.zeros(nfi, dtype=complex)
    V_fft = np.zeros(nfi, dtype=complex) 

    for fi in range(nfi):
        w = wVec[fi]
        k = w*sinTheta/phaseVelIn
   
        #Displacement at interface
        uInterface = us[:,:, fi]

        #find parent layer where yTopLayer>=y>yBotLayer
        parentLayer = N - 1 - np.searchsorted(Layers[::-1], y, side = "left") 
            
        if parentLayer == (N-1):
            if y==Layers[-1]:
                uz = uInterface[2*N-2:2*N,:]
        elif parentLayer < (N-1):
            yTop = Layers[parentLayer]
            yBot = Layers[parentLayer+1]
            uTop = uInterface[2*parentLayer:2*parentLayer+2,:]
            uBot = uInterface[2*parentLayer+2:2*parentLayer+4,:]
            if y < yTop:
                uz = GetDisplacementAtInteriorLayer(y,yTop,yBot,uTop,uBot,k,p[parentLayer],s[parentLayer],mu[parentLayer],aSP[parentLayer])
            elif np.isclose(y, yTop, rtol=1e-05):
                uz = uTop
                    
        U_fft[fi] = uz[0,0]*np.exp(-1j*k*(x-x0))
        V_fft[fi] = 1j*uz[1,0]*np.exp(-1j*k*(x-x0))

    #Add 0 for zero frequency and frequency larger than cutOffFrequency
    U_fft = np.concatenate((np.zeros(1, dtype=complex), U_fft, np.zeros(int(Nt/2)-nfi, dtype=complex)), axis=0)
    V_fft = np.concatenate((np.zeros(1), V_fft, np.zeros(int(Nt/2)-nfi)), axis=0)

    U = np.real(np.fft.irfft(U_fft, Nt, axis=0))
    V = np.real(np.fft.irfft(V_fft, Nt, axis=0))

    #Gathers the field components
    Z = np.stack([U,V], axis=1)

    return Z

def RHbackground2Dfield(FdispIn,wVec,interpDispersion,interpuMmodeShape,interpvMmodeShape,x,y,Nt,x0,y0):
    """
    This function calculates the displacements at a specific point in time domain for Rayleigh
    wave in stratified soil.\n
    
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Kien T. Nguyen 2021, ORCID: 0000-0001-5761-3156
    
    Parameters
    ----------
    wVec  : array
        The angular frequency
    interpDispersion  : object of class scipy.interpolate.interp1d
        The 1d interpolation function to interpolate the phase velocity at chosen angular frequency 
    interpuMmodeShape, interpvMmodeShape  : objects of class scipy.interpolate.RectBivariateSpline
        The 2d interpolation function to interpolate the horizontal and vertical displacement mode shape
        at the query point. Functions of y and angular frequency 
    x, y  : float
        x- and y-coordinate of the query point
    Nt  : int
        Length of the Disp, Vels, Accel after zero padding
    FdispIn  : array
        FFT of the input displacement
    x0, y0  : float
        x- and y-coordinate of the reference point (where the incoming signal time series is prescribed)
        
    Returns
    -------
    U, V  : array
        The horizontal and vertical displacements at the query point
    """
    U_fft = np.zeros(int(Nt/2), dtype=complex)
    V_fft = np.zeros(int(Nt/2), dtype=complex)
    
    for indexFre, w in enumerate(wVec):
        vr = interpDispersion(w)
        k = w/vr
        uRef = interpuMmodeShape(y0,w)
        Ufi = interpuMmodeShape(y,w)/uRef*FdispIn[indexFre]
        Vfi = interpvMmodeShape(y,w)/uRef*FdispIn[indexFre]
        
        U_fft[indexFre+1] = Ufi[0,0]*np.exp(-1j*k*(x-x0)) #since zero frequency is not accounted for, so it start from index [1]
        V_fft[indexFre+1] = -1j*Vfi[0,0]*np.exp(-1j*k*(x-x0))
        
    U = np.real(np.fft.irfft(U_fft, Nt, axis=0))
    V = np.real(np.fft.irfft(V_fft, Nt, axis=0))

    #Gathers the field components
    Z = np.stack([U,V], axis=1)
    
    return Z

def PSVbackground3Dfield(us, Layers, wVec, p, s, h, mu, aSP, phaseVelIn, di, sinTheta, N, Nt, x0, x1, x2, x3):
    """
    This function calculates the displacement time series in 3D wave propagation
    problem, in which the SV or P wave incoming from the half space
    underneath under arbitrary incident angle from 0 to 90 degrees and propagating 
    through stratified soil domain.

    Note: 
    [1] The coordinate system and displacement positive axes:

        z(W) ^  y(V)
             | /
             |/
             o-----> x(U)

    [2] At each frequency, the displacement components are calculated 
        based on Eduardo Kausel's Stiffness Matrix Method, in "Fundamental 
        Solutions in Elastodynamics, A Compendium", chap. 10, pp. 140--159 

    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Kien T. Nguyen 2021, ORCID: 0000-0001-5761-3156
    
    Parameters
    ----------
    us  : array
        Displacement at the layer interfaces in frequency domain
    Layers  : array
        The y-coordinate of soil layer interfaces, including imaginary and half-space interfaces
    wVec  : array
        Angular frequency spectrum 
    p, s  : array
        Complex coefficients
    h  : array
        The thickness of soil layers
    mu  : array
        The shear modulus of soil layers
    aSP  : array
        The ratio between shear wave velocity and dilatational wave velocity of soil layers
    phaseVelIn  : float
        The phase velocity of incoming wave in the half space underneath
    di  : array
        Polarization of the propagation direction with respect to horizontal axis (x-axis)
    sinTheta  : float
        Sine of the incoming angle
    N  : int
        Number of soil layers, including imaginary layer and half space
    Nt  : int
        Length of the Disp, Vels, Accel after zero padding
    x0  :float
        x-coordinate of the reference point (where the incoming signal time series is prescribed)
    x1, x2, x3 :float
        Cartesian coordinates of the query point
        
    Returns
    -------
    Z  : array
        Displacement time series at the query point
    """
    nfi = len(wVec)
    U_fft = np.zeros(nfi, dtype=complex)
    V_fft = np.zeros(nfi, dtype=complex) 

    x = x1*di[0] + x2*di[1]
    y = x3

    for fi in range(nfi):
        w = wVec[fi]
        k = w*sinTheta/phaseVelIn

        #Displacement at interface
        uInterface = us[:, :, fi]

        #find parent layer where yTopLayer>=y>yBotLayer
        parentLayer = N - 1 - np.searchsorted(Layers[::-1], y, side = "left") 
            
        if parentLayer == (N-1):
            if y==Layers[-1]:
                uz = uInterface[2*N-2:2*N,:]
        elif parentLayer < (N-1):
            yTop = Layers[parentLayer]
            yBot = Layers[parentLayer+1]
            uTop = uInterface[2*parentLayer:2*parentLayer+2,:]
            uBot = uInterface[2*parentLayer+2:2*parentLayer+4,:]
            if y < yTop:
                uz = GetDisplacementAtInteriorLayer(y,yTop,yBot,uTop,uBot,k,p[parentLayer],s[parentLayer],mu[parentLayer],aSP[parentLayer])
            elif np.isclose(y, yTop, rtol=1e-06):
                uz = uTop
                    
        U_fft[fi] = uz[0,0]*np.exp(-1j*k*(x-x0))
        V_fft[fi] = 1j*uz[1,0]*np.exp(-1j*k*(x-x0))

    #Add 0 for zero frequency and frequency larger than cutOffFrequency
    U_fft = np.concatenate((np.zeros(1, dtype=complex), U_fft, np.zeros(int(Nt/2)-nfi, dtype=complex)), axis=0)
    V_fft = np.concatenate((np.zeros(1), V_fft, np.zeros(int(Nt/2)-nfi)), axis=0)

    U = di[0]*np.real(np.fft.irfft(U_fft, Nt, axis=0))
    V = di[1]*np.real(np.fft.irfft(U_fft, Nt, axis=0))
    W = np.real(np.fft.irfft(V_fft, Nt, axis=0))

    #Gathers the field components
    Z = np.stack([U,V,W], axis=1)

    return Z

def SHbackground3Dfield(Values, t, X, X0, Xmin, di, nt, fTag):
    """
    """
    #TODO: This is not working and needs to be corrected
    #Compute the 3D Field Components.
    U = np.zeros(nt)
    V = np.zeros(nt)
    W = np.zeros(nt)

    #Gathers the field components
    Z = np.stack([U,V,W], axis=1)

    return Z

def RHbackground3Dfield(FdispIn, wVec, interpDispersion, interpuMmodeShape, interpvMmodeShape, di, x1, y1, z1, Nt, xmin, ymin, zmin):
    """
    This function calculates the displacements at a specific point in time domain for Rayleigh
    wave in stratified soil.\n
    
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Kien T. Nguyen 2021, ORCID: 0000-0001-5761-3156
    
    Parameters
    ----------
    wVec  : array
        The angular frequency
    interpDispersion  : object of class scipy.interpolate.interp1d
        The 1d interpolation function to interpolate the phase velocity at chosen angular frequency 
    interpuMmodeShape, interpvMmodeShape  : objects of class scipy.interpolate.RectBivariateSpline
        The 2d interpolation function to interpolate the horizontal and vertical displacement mode shape
        at the query point. Functions of y and angular frequency 
    x, y  : float
        x- and y-coordinate of the query point
    Nt  : int
        Length of the Disp, Vels, Accel after zero padding
    FdispIn  : array
        FFT of the input displacement
    x0, y0  : float
        x- and y-coordinate of the reference point (where the incoming signal time series is prescribed)
        
    Returns
    -------
    U, V  : array
        The horizontal and vertical displacements at the query point
    """
    x = x1*di[0] + y1*di[1]
    y = z1

    x0 = xmin*di[0] + ymin*di[1]
    y0 = zmin

    U_fft = np.zeros(int(Nt/2), dtype=complex)
    V_fft = np.zeros(int(Nt/2), dtype=complex)
    
    for indexFre, w in enumerate(wVec):
        vr = interpDispersion(w)
        k = w/vr
        uRef = interpuMmodeShape(y0,w)
        Ufi = interpuMmodeShape(y,w)/uRef*FdispIn[indexFre]
        Vfi = interpvMmodeShape(y,w)/uRef*FdispIn[indexFre]
        
        U_fft[indexFre+1] = Ufi[0,0]*np.exp(-1j*k*(x-x0)) #since zero frequency is not accounted for, so it start from index [1]
        V_fft[indexFre+1] = -1j*Vfi[0,0]*np.exp(-1j*k*(x-x0))
        
    U = di[0]*np.real(np.fft.irfft(U_fft, Nt, axis=0))
    V = di[1]*np.real(np.fft.irfft(U_fft, Nt, axis=0))
    W = np.real(np.fft.irfft(V_fft, Nt, axis=0))

    #Gathers the field components
    Z = np.stack([U,V,W], axis=1)

    return Z

def GetRayleighVelocity(Vs, nu):
    """
    This function calculates the Rayleigh wave phase velocity in homogeneous
    elastic domain.
    
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Kien T. Nguyen 2021, ORCID: 0000-0001-5761-3156
    
    Parameters
    ----------
    Vs  : float
        Shear wave velocity of the homogeneous soil domain 
    nu  : float
        The Poisson's ratio of soil
        
    Returns
    -------
    ratio*Vs  : float
        Rayleigh wave phase velocity
    """
    #Compute the Rayleigh wave velocity.
    lambda0 = 2.0*1.0*nu/(1.0 - 2.0*nu)
    hk = np.sqrt(1.0/(lambda0 + 2.0*1.0))
    
    p = [-16.0*(1.0 - hk**2), 24.0 - 16.0*hk**2, -8.0, 1.0]
    x = np.sort(np.roots(p))
    ratio = 1.0/np.sqrt(np.real(x[2]))
    
    return ratio*Vs

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
                print(" Generating DRM files. This may take a few minutes...")
                start_time = time.time()

                #Gets the Domain Reduction Method Information.
                nodes, conditions, t, Disp, Vels, Accel, dt, option = ParseDRMFile(Entities['Functions'][fTag])

                #Computes the input displacement, velocities and accelerations.
                Disp, Vels, Accel = GenerateTimeSeries(Disp, Vels, Accel, dt, option)

                #Creates the DRM Directory
                dirName = Options['path'] + '/' + 'DRM'
                if not os.path.exists(dirName):
                    os.mkdir(dirName)

                #Computes and Generates the Domain Reduction files for each node.
                x0 = Entities['Functions'][fTag]['attributes']['x0']
                xmin = Entities['Functions'][fTag]['attributes']['xmin']
                funName = Entities['Functions'][fTag]['name']
                funOption = Entities['Functions'][fTag]['attributes']['option']
                waveType = funOption.upper()

                #Layer material information
                nmat = len(Entities['Functions'][fTag]['attributes']['material'])
                beta = np.zeros((nmat,))
                rho = np.zeros((nmat,))
                nu = np.zeros((nmat,))
                for k, mTag in enumerate(Entities['Functions'][fTag]['attributes']['material']):
                    material = Entities['Materials'][mTag]['attributes']
                    beta[k] = np.sqrt(material['E']/2.0/material['rho']/(1.0 + material['nu']))
                    rho[k] = material['rho']
                    nu[k] = material['nu']

                nt = len(t)
                
                #Computes the DRM field depending on the option name (P,SV,SH,RH)
                if Options['dimension'] == 2:
                    if waveType == 'P' or waveType == 'SV':
                        fun = Entities['Functions'][fTag]['attributes']
                        if 'theta' not in fun:
                            fun['theta'] = 0.0
                        if 'CutOffFrequency' not in fun:
                            fun['CutOffFrequency'] = 30.0
                        if 'df' not in fun:
                            fun['df'] = 0.2

                        #Unpack Layer information
                        angle = fun['theta']
                        layers = fun['layer']
                        df = fun['df']
                        CutOffFrequency = fun['CutOffFrequency']

                        ufull, vfull, afull, layers, beta, rho, nu, wVec, p, s, h, mu, aSP, phaseVelIn, sinTheta, N, Nt = DataPreprocessing(Disp, Vels, Accel, layers, beta, rho, nu, angle, xmin[1], nt, dt, fun)

                        #Compute Interface responses
                        uInterface = SoilInterfaceResponse(ufull, wVec, p, s, h, mu, aSP, phaseVelIn, sinTheta, N)
                        vInterface = SoilInterfaceResponse(vfull, wVec, p, s, h, mu, aSP, phaseVelIn, sinTheta, N)
                        aInterface = SoilInterfaceResponse(afull, wVec, p, s, h, mu, aSP, phaseVelIn, sinTheta, N)

                        x0 = xmin[0]
                        with concurrent.futures.ProcessPoolExecutor() as executor:
                            for k, n in enumerate(nodes):
                                x = Entities['Nodes'][n]['coords']
                                U = PSVbackground2Dfield(uInterface, layers, wVec, p, s, h, mu, aSP, phaseVelIn, sinTheta, N, Nt, x0, x[0], x[1])
                                V = PSVbackground2Dfield(vInterface, layers, wVec, p, s, h, mu, aSP, phaseVelIn, sinTheta, N, Nt, x0, x[0], x[1])
                                A = PSVbackground2Dfield(aInterface, layers, wVec, p, s, h, mu, aSP, phaseVelIn, sinTheta, N, Nt, x0, x[0], x[1])
                                executor.submit(WriteDRMFile, dirName, funName, fTag, U, V, A, nt, 6, n, conditions[k])
                    elif waveType == 'RH':
                        #Unpack Layer information
                        fun    = Entities['Functions'][fTag]['attributes']
                        layers = fun['layer']
                        df     = fun['df']
                        endFrequency = fun['CutOffFrequency'] # energy of frequency larger than this is zero

                        #Consider until depth = depthFactor*wavelength to ensure fix end approximation
                        depthFactor = 3.0

                        #Y-grid spacing of points used to interpolate the mode shape, from ground surface to yDRMmin
                        dy1  = np.min(beta/endFrequency/16.0)

                        #Number of element per wavelength in the extended region, from yDRMmin to depth
                        nepw = 40

                        #Currently only work with mode = 0: fundamental mode of Rayleigh wave      
                        mode = 0

                        #Compute the FFT for displacement, velocity and acceleration fields
                        wVec, FFTdisp, FFTvels, FFTaccel, df, Nt, startFrequency, endFrequency = GetRayleighFFTfields(Disp, Vels, Accel, endFrequency, dt, df, nt)

                        #Computes Mode Shape and Phase velocity dispersion for generation of the interpolation functions: uModeShape, vModeShape, and yGridModeShape
                        fDispersion, phaseVelDispersion, yGridModeShape, uModeShape, vModeShape = GetRayleighDispersionAndModeShape(mode, layers, beta, rho, nu, dy1, xmin[1], nepw, startFrequency, endFrequency, df, depthFactor)

                        #
                        yGridModeShape = np.flipud(yGridModeShape)
                        uModeShape = np.flipud(np.real(uModeShape))
                        vModeShape = np.flipud(np.real(vModeShape))

                        #Interpolation function required to obtain displacement, velocity and acceleration fields in nodes inside the soil layers
                        interpDispersion = interpolate.interp1d(2.0*np.pi*fDispersion, phaseVelDispersion,kind='linear', fill_value='extrapolate')
                        interpuMmodeShape = interpolate.RectBivariateSpline(yGridModeShape,2.0*np.pi*fDispersion, uModeShape)
                        interpvMmodeShape = interpolate.RectBivariateSpline(yGridModeShape,2.0*np.pi*fDispersion, vModeShape)

                        with concurrent.futures.ProcessPoolExecutor() as executor:
                            for k, n in enumerate(nodes):
                                x = Entities['Nodes'][n]['coords']  
                                U = RHbackground2Dfield(FFTdisp, wVec, interpDispersion, interpuMmodeShape, interpvMmodeShape, x[0], x[1], Nt, xmin[0], x0[1])
                                V = RHbackground2Dfield(FFTvels, wVec, interpDispersion, interpuMmodeShape, interpvMmodeShape, x[0], x[1], Nt, xmin[0], x0[1])
                                A = RHbackground2Dfield(FFTaccel, wVec, interpDispersion, interpuMmodeShape, interpvMmodeShape, x[0], x[1], Nt, xmin[0], x0[1])
                                executor.submit(WriteDRMFile, dirName, funName, fTag, U, V, A, Nt, 6, n, conditions[k])
                    else:
                        print('\x1B[31m ERROR \x1B[0m: The specified PLANEWAVE (2D) option (=%s) is not recognized' % funOption)
                elif Options['dimension'] == 3:
                    if 'phi' not in Entities['Functions'][fTag]['attributes']:
                        Entities['Functions'][fTag]['attributes']['phi'] = 0.0

                    #The wave direction on horizontal plane
                    azimut  = Entities['Functions'][fTag]['attributes']['phi']
                    phi     = azimut*np.pi/180.0
                    di      = np.array([np.cos(phi), np.sin(phi)]) 

                    if waveType == 'P' or waveType == 'SV':
                        fun = Entities['Functions'][fTag]['attributes']
                        if 'theta' not in fun:
                            fun['theta'] = 0.0
                        if 'CutOffFrequency' not in fun:
                            fun['CutOffFrequency'] = 30.0
                        if 'df' not in fun:
                            fun['df'] = 0.2

                        #Unpack Layer information
                        angle = fun['theta']
                        layers = fun['layer']
                        df = fun['df']
                        CutOffFrequency = fun['CutOffFrequency']

                        ufull, vfull, afull, layers, beta, rho, nu, wVec, p, s, h, mu, aSP, phaseVelIn, sinTheta, N, Nt = DataPreprocessing(Disp, Vels, Accel, layers, beta, rho, nu, angle, xmin[2], nt, dt, fun)

                        #Compute Interface responses
                        uInterface = SoilInterfaceResponse(ufull, wVec, p, s, h, mu, aSP, phaseVelIn, sinTheta, N)
                        vInterface = SoilInterfaceResponse(vfull, wVec, p, s, h, mu, aSP, phaseVelIn, sinTheta, N)
                        aInterface = SoilInterfaceResponse(afull, wVec, p, s, h, mu, aSP, phaseVelIn, sinTheta, N)

                        x0 = xmin[0]*di[0] + xmin[0]*di[1]
                        with concurrent.futures.ProcessPoolExecutor() as executor:
                            for k, n in enumerate(nodes):
                                x = Entities['Nodes'][n]['coords']  
                                U = PSVbackground3Dfield(uInterface, layers, wVec, p, s, h, mu, aSP, phaseVelIn, di, sinTheta, N, Nt, x0, x[0], x[1], x[2])
                                V = PSVbackground3Dfield(vInterface, layers, wVec, p, s, h, mu, aSP, phaseVelIn, di, sinTheta, N, Nt, x0, x[0], x[1], x[2])
                                A = PSVbackground3Dfield(aInterface, layers, wVec, p, s, h, mu, aSP, phaseVelIn, di, sinTheta, N, Nt, x0, x[0], x[1], x[2])
                                executor.submit(WriteDRMFile, dirName, funName, fTag, U, V, A, nt, 9, n, conditions[k])
                    elif waveType == 'SH':
                        #TODO: Complete SH case in 3D
                        with concurrent.futures.ProcessPoolExecutor() as executor: 
                            for k, n in enumerate(nodes):
                                x = Entities['Nodes'][n]['coords']
                                U = SHbackground3Dfield(Disp, t, x, x0, xmin, di, nt, fTag)
                                V = SHbackground3Dfield(Vels, t, x, x0, xmin, di, nt, fTag)
                                A = SHbackground3Dfield(Accel, t, x, x0, xmin, di, nt, fTag)
                                executor.submit(WriteDRMFile, dirName, funName, fTag, U, V, A, nt, 9, n, conditions[k])
                    elif waveType == 'RH':
                        #Unpack Layer information
                        fun    = Entities['Functions'][fTag]['attributes']
                        layers = fun['layer']
                        df     = fun['df']
                        endFrequency = fun['CutOffFrequency'] # energy of frequency larger than this is zero

                        #Consider until depth = depthFactor*wavelength to ensure fix end approximation
                        depthFactor = 3.0

                        #Y-grid spacing of points used to interpolate the mode shape, from ground surface to yDRMmin
                        dy1  = np.min(beta/endFrequency/16.0) 

                        #Number of element per wavelength in the extended region, from yDRMmin to depth
                        nepw = 40

                        #Currently only work with mode = 0: fundamental mode of Rayleigh wave      
                        mode = 0

                        #Compute the FFT for displacement, velocity and acceleration fields
                        wVec, FFTdisp, FFTvels, FFTaccel, df, Nt, startFrequency, endFrequency = GetRayleighFFTfields(Disp, Vels, Accel, endFrequency, dt, df, nt)

                        #Computes Mode Shape and Phase velocity dispersion for generation of the interpolation functions: uModeShape, vModeShape, and yGridModeShape
                        fDispersion, phaseVelDispersion, yGridModeShape, uModeShape, vModeShape = GetRayleighDispersionAndModeShape(mode, layers, beta, rho, nu, dy1, xmin[2], nepw, startFrequency, endFrequency, df, depthFactor)

                        #
                        yGridModeShape = np.flipud(yGridModeShape)
                        uModeShape = np.flipud(np.real(uModeShape))
                        vModeShape = np.flipud(np.real(vModeShape))

                        #Interpolation function required to obtain displacement, velocity and acceleration fields in nodes inside the soil layers
                        interpDispersion = interpolate.interp1d(2.0*np.pi*fDispersion, phaseVelDispersion,kind='linear', fill_value='extrapolate')
                        interpuMmodeShape = interpolate.RectBivariateSpline(yGridModeShape,2.0*np.pi*fDispersion, uModeShape)
                        interpvMmodeShape = interpolate.RectBivariateSpline(yGridModeShape,2.0*np.pi*fDispersion, vModeShape)

                        with concurrent.futures.ProcessPoolExecutor() as executor:
                            for k, n in enumerate(nodes):
                                x = Entities['Nodes'][n]['coords']
                                U = RHbackground3Dfield(FFTdisp, wVec, interpDispersion, interpuMmodeShape, interpvMmodeShape, di, x[0], x[1], x[2], Nt, xmin[0], xmin[1], x0[2])
                                V = RHbackground3Dfield(FFTvels, wVec, interpDispersion, interpuMmodeShape, interpvMmodeShape, di, x[0], x[1], x[2], Nt, xmin[0], xmin[1], x0[2])
                                A = RHbackground3Dfield(FFTaccel, wVec, interpDispersion, interpuMmodeShape, interpvMmodeShape, di, x[0], x[1], x[2], Nt, xmin[0], xmin[1], x0[2])
                                executor.submit(WriteDRMFile, dirName, funName, fTag, U, V, A, nt, 9, n, conditions[k])
                    else:
                        print('\x1B[31m ERROR \x1B[0m: The specified PLANEWAVE (3D) option (=%s) is not recognized' % funOption)
                else:
                    print('\x1B[31m ERROR \x1B[0m: The specified dimension (=%d) is not possible for DRM' % Options['dimensions'])

                #Update the load type after files were successfully created
                Entities['Loads'][lTag]['attributes']['type'] = 'GENERALWAVE'

                #Update the domain reduction time series path where files are located.
                Entities['Functions'][fTag]['attributes']['file'] = dirName + "/" + funName + "-" + str(fTag) + ".$.drm"

                end_time = time.time()
                print(" DRM files (",len(nodes),") were created in ", end_time - start_time, "s\n")
