#!/usr/bin/python3
# -*- coding: Utf-8 -*-

import os
import numpy as np
from scipy import signal
import concurrent.futures
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
    None
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

def SVbackground2Dfield(Values, t, X, X0, Xmin, nt, fTag):
    """
    This function generates the 2D DRM field to be specified at a particular node. The 
    displacement,, velocity and acceleration fields are computed assuming homogenous
    half-space.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2021
    """
    #Soil properties
    mTag = Entities['Functions'][fTag]['attributes']['material']
    Es  = Entities['Materials'][mTag]['attributes']['E']
    nu  = Entities['Materials'][mTag]['attributes']['nu']
    rho = Entities['Materials'][mTag]['attributes']['rho']
    Vs  = np.sqrt(Es/2.0/rho/(1.0 + nu)) 
    Vp  = Vs*np.sqrt(2.0*(1.0 - nu)/(1.0 - 2.0*nu))

    #Signal properties
    angle   = Entities['Functions'][fTag]['attributes']['theta']
    theta_s = angle*np.pi/180.0
    theta_p = np.arcsin(Vp/Vs*np.sin(theta_s))

    #Reference time and coordinates for signal
    x_rela =  X[0] - X0[0]
    y_rela =  X[1] - X0[1]
    t_init = (X0[0] - Xmin[0])/Vs*np.sin(theta_s) + (X0[1] - Xmin[1])/Vs*np.cos(theta_s)

    #Time shift according to reference point X0.
    t1 = -x_rela/Vs*np.sin(theta_s) - y_rela/Vs*np.cos(theta_s) + t - t_init
    t2 = -x_rela/Vs*np.sin(theta_s) + y_rela/Vs*np.cos(theta_s) + t - t_init
    t3 = -x_rela/Vp*np.sin(theta_p) + y_rela/Vp*np.cos(theta_p) + t - t_init

    #Incident and Reflected amplitudes.
    k = Vp/Vs
    U_si = 1.000
    U_pr = U_si*(-2*k*np.sin(2*theta_s)*np.cos(2*theta_s))/(np.sin(2*theta_s)*np.sin(2*theta_p) + k**2*np.cos(2*theta_s)**2)
    U_sr = U_si*(np.sin(2*theta_s)*np.sin(2*theta_p) - k**2*np.cos(2*theta_s)**2)/(np.sin(2*theta_s)*np.sin(2*theta_p) + k**2*np.cos(2*theta_s)**2)

    #Compute the 2D Field Components.
    U = np.zeros(nt)
    V = np.zeros(nt)

    U += U_si*np.cos(theta_s)*np.interp(t1, t, Values, left=0.0, right=0.0)
    V -= U_si*np.sin(theta_s)*np.interp(t1, t, Values, left=0.0, right=0.0)

    U -= U_sr*np.cos(theta_s)*np.interp(t2, t, Values, left=0.0, right=0.0)
    V -= U_sr*np.sin(theta_s)*np.interp(t2, t, Values, left=0.0, right=0.0)

    U -= U_pr*np.sin(theta_p)*np.interp(t3, t, Values, left=0.0, right=0.0)
    V += U_pr*np.cos(theta_p)*np.interp(t3, t, Values, left=0.0, right=0.0)

    Z = np.stack([U,V], axis=1)

    return Z

def SVbackground2DfieldCritical(Values, SP, SS, Vp, Vs, cj, sj, p, w, xp, xmin):
    """
    The script is to calculate the displacement, velocity or acceleration components 
    of a 2D rayleigh wave by using real FFT and inverse real FFT, accounting for both 
    amplitude and phase or the incident angle larger than critical one.
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Kien T. Nguyen 2021
    """
    #TODO: This is not working and needs to be corrected
    x     = xp[0]
    y     = xp[1]
    x0    = xmin[0]
    y0    = xmin[1]
    xrel  = x - x0
    yrel  = y0 - y 
    FSVin = np.fft.rfft(Values)
    
    #Horizontal U1, U2, U3 due to upgoing SV, downgoing P, downgoing SV
    U1 = cj*FSVin*np.exp(-1j*(p*xrel-cj/Vs*yrel)*w)
    U2 = Vp*p*SP*np.exp(np.sqrt(p*p-1.0/Vp/Vp)*y*w)*FSVin*np.exp(-1j*(p*xrel-cj/Vs*y0)*w)
    U3 = cj*SS*FSVin*np.exp(-1j*(p*xrel-cj/Vs*(y+y0))*w)
    
    #Vertical V1, V2, V3 due to upgoing SV, downgoing P, downgoing SV
    V1 = -sj*FSVin*np.exp(-1j*(p*xrel-cj/Vs*yrel)*w)
    V2 = -1j*np.sqrt(Vp*Vp*p*p-1.0)*SP*np.exp(np.sqrt(p*p-1.0/Vp/Vp)*y*w)*FSVin*np.exp(-1j*(p*xrel-cj/Vs*y0)*w)
    V3 = sj*SS*FSVin*np.exp(-1j*(p*xrel-cj/Vs*(y+y0))*w)
        
    U_fft = U1 + U2 + U3
    V_fft = V1 + V2 + V3
    
    U = np.real(np.fft.irfftn(U_fft))
    V = np.real(np.fft.irfftn(V_fft))

    #Gathers the field components
    Z = np.stack([U,V], axis=1)

    return Z

def RHbackground2Dfield(Values, t, Vp, Vs, Vr, w, kR, qR, sR, x0, xmin, x1, x2):
    """
    The script is to calculate the displacement, velocity or acceleration components 
    of a 2D rayleigh wave by using real FFT and inverse real FFT, accounting for both 
    amplitude and phase
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Kien T. Nguyen 2021
    """
    x = x1
    y = x0[1] - x2
    Bn = np.fft.rfft(Values)
    An = 2.0*(Vs/Vr)**2*Bn # scale factor An = A*kR in I.9 Viktorov

    #Compute the 2D Field Components.
    U_fft = An*(np.exp(-qR*y)-(1.0-Vr**2.0/2.0/Vs**2)*np.exp(-sR*y))*np.exp(-1j*x/Vr*w)
    V_fft = np.exp(1j*(-np.pi/2.0))*An*np.sqrt(1.0-Vr**2.0/Vp**2.0)*(np.exp(-qR*y)-1.0/(1.0-Vr**2/2.0/Vs**2.0)*np.exp(-sR*y))*np.exp(-1j*x/Vr*w)
    
    U =  np.real(np.fft.irfftn(U_fft))
    V = -np.real(np.fft.irfftn(V_fft))

    #Gathers the field components
    Z = np.stack([U,V], axis=1)

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

def SVbackground3Dfield(Values, t, X, X0, Xmin, di, nt, fTag):
    """
    This function generates the 3D DRM field to be specified at a particular node. The 
    displacement,, velocity and acceleration fields are computed assuming homogenous
    half-space.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2021
    """
    #Soil properties
    mTag = Entities['Functions'][fTag]['attributes']['material']
    Es  = Entities['Materials'][mTag]['attributes']['E']
    nu  = Entities['Materials'][mTag]['attributes']['nu']
    rho = Entities['Materials'][mTag]['attributes']['rho']
    Vs  = np.sqrt(Es/2.0/rho/(1.0 + nu)) 
    Vp  = Vs*np.sqrt(2.0*(1.0 - nu)/(1.0 - 2.0*nu))

    #Signal properties
    angle   = Entities['Functions'][fTag]['attributes']['theta']
    theta_s = angle*np.pi/180.0
    theta_p = np.arcsin(Vp/Vs*np.sin(theta_s))

    #Reference time and coordinates for signal
    x_rela =  X[0] - X0[0]
    y_rela =  X[1] - X0[1]
    z_rela =  X[2] - X0[2]
    t_init = (X0[0] - Xmin[0])/Vs*np.sin(theta_s)*di[0] + (X0[1] - Xmin[1])/Vs*np.sin(theta_s)*di[1] + (X0[2] - Xmin[2])/Vs*np.cos(theta_s)

    #Time shift according to reference point X0.
    t1 = -x_rela/Vs*np.sin(theta_s)*di[0] - y_rela/Vs*np.sin(theta_s)*di[1] - z_rela/Vs*np.cos(theta_s) + t - t_init
    t2 = -x_rela/Vs*np.sin(theta_s)*di[0] - y_rela/Vs*np.sin(theta_s)*di[1] + z_rela/Vs*np.cos(theta_s) + t - t_init
    t3 = -x_rela/Vp*np.sin(theta_p)*di[0] - y_rela/Vp*np.sin(theta_p)*di[1] + z_rela/Vp*np.cos(theta_p) + t - t_init

    #Incident and Reflected amplitudes.
    k = Vp/Vs
    U_si = 1.000
    U_pr = U_si*(-2*k*np.sin(2*theta_s)*np.cos(2*theta_s))/(np.sin(2*theta_s)*np.sin(2*theta_p) + k**2*np.cos(2*theta_s)**2)
    U_sr = U_si*(np.sin(2*theta_s)*np.sin(2*theta_p) - k**2*np.cos(2*theta_s)**2)/(np.sin(2*theta_s)*np.sin(2*theta_p) + k**2*np.cos(2*theta_s)**2)

    #Compute the 3D Field Components.
    U = np.zeros(nt)
    V = np.zeros(nt)
    W = np.zeros(nt)

    A  = np.interp(t1, t, Values, left=0.0, right=0.0)
    U += U_si*di[0]*np.cos(theta_s)*A
    V += U_si*di[1]*np.cos(theta_s)*A
    W -= U_si*np.sin(theta_s)*A

    A  = np.interp(t2, t, Values, left=0.0, right=0.0)
    U -= U_sr*di[0]*np.cos(theta_s)*A
    V -= U_sr*di[1]*np.cos(theta_s)*A
    W -= U_sr*np.sin(theta_s)*A

    A  = np.interp(t3, t, Values, left=0.0, right=0.0)
    U -= U_pr*di[0]*np.sin(theta_p)*A
    V -= U_pr*di[1]*np.sin(theta_p)*A
    W += U_pr*np.cos(theta_p)*A

    #Gathers the field components
    Z = np.stack([U,V,W], axis=1)

    return Z

def SVbackground3DfieldCritical(Values, SP, SS, Vp, Vs, cj, sj, p, w, di, xp, xmin, nt):
    """
    The script is to calculate the displacement, velocity or acceleration components 
    of a 3D rayleigh wave by using real FFT and inverse real FFT, accounting for both 
    amplitude and phase or the incident angle larger than critical one.
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Kien T. Nguyen 2021
    """
    #TODO: This is not working and needs to be corrected
    #Compute the 3D Field Components.
    U = np.zeros(nt)
    V = np.zeros(nt)
    W = np.zeros(nt)

    #Gathers the field components
    Z = np.stack([U,V,W], axis=1)

    return Z

def RHbackground3Dfield(Values, Vp, Vs, Vr, w, kR, qR, sR, x0, xmin, di, x1, x2, x3):
    """
    The script is to calculate the displacement, velocity or acceleration components 
    of a 3D rayleigh wave by using real FFT and inverse real FFT, accounting for both 
    amplitude and phase
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Kien T. Nguyen 2021
    """
    x = x1*di[0] + x2*di[1]
    y = x0[2] - x3
    Bn = np.fft.rfft(Values)
    An = 2.0*(Vs/Vr)**2*Bn # scale factor An = A*kR in I.9 Viktorov

    #Compute the 2D Field Components.
    U_fft = An*(np.exp(-qR*y)-(1.0-Vr**2.0/2.0/Vs**2)*np.exp(-sR*y))*np.exp(-1j*x/Vr*w)
    V_fft = np.exp(1j*(-np.pi/2.0))*An*np.sqrt(1.0-Vr**2.0/Vp**2.0)*(np.exp(-qR*y)-1.0/(1.0-Vr**2/2.0/Vs**2.0)*np.exp(-sR*y))*np.exp(-1j*x/Vr*w)

    #Compute the 3D Field Components.
    U = di[0]*np.real(np.fft.irfftn(U_fft))
    V = di[1]*np.real(np.fft.irfftn(U_fft))
    W = -np.real(np.fft.irfftn(V_fft))

    #Gathers the field components
    Z = np.stack([U,V,W], axis=1)

    return Z

def ComputeSVCriticalAngle(fTag):
    """
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Kien T. Nguyen 2021
    """
    mTag = Entities['Functions'][fTag]['attributes']['material']
    Es  = Entities['Materials'][mTag]['attributes']['E']
    nu  = Entities['Materials'][mTag]['attributes']['nu']
    rho = Entities['Materials'][mTag]['attributes']['rho']
    Vs  = np.sqrt(Es/2.0/rho/(1.0 + nu)) 
    Vp  = Vs*np.sqrt(2.0*(1.0 - nu)/(1.0 - 2.0*nu))
    Angle = Entities['Functions'][fTag]['attributes']['theta']*np.pi/180.0
    Theta = Vp/Vs*np.sin(Angle)

    return Theta, Vs, Vp

def GetRayleighVelocity(Vs, nu):
    """
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Kien T. Nguyen 2021
    """
    #Compute the Rayleigh wave velocity.
    lambda0 = 2.0*1.0*nu/(1.0-2.0*nu)
    hk = np.sqrt(1.0/(lambda0+2.0*1.0))
    
    p =[-16.0*(1.0-hk**2), 24.0-16.0*hk**2, -8.0, 1.0]
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

                nt = len(t)
                nd = Options['dimension']

                #Tapper fucntion
                shift = 0 if (nt % 2) == 0 else 1
                window = signal.tukey(nt, alpha = 0.15)
                Disp *= window
                Vels *= window
                Accel *= window
                
                if nd == 2:
                    if funOption.upper() == 'SV':
                        if 'theta' not in Entities['Functions'][fTag]['attributes']:
                            Entities['Functions'][fTag]['attributes']['theta'] = 0.0
                        thetacr, Vs, Vp = ComputeSVCriticalAngle(fTag)

                        if thetacr < 1.0:
                            with concurrent.futures.ProcessPoolExecutor() as executor: 
                                for k, n in enumerate(nodes):
                                    x = Entities['Nodes'][n]['coords']
                                    U = SVbackground2Dfield(Disp, t, x, x0, xmin, nt, fTag)
                                    V = SVbackground2Dfield(Vels, t, x, x0, xmin, nt, fTag)
                                    A = SVbackground2Dfield(Accel, t, x, x0, xmin, nt, fTag)        
                                    executor.submit(WriteDRMFile, dirName, funName, fTag, U, V, A, nt, 6, n, conditions[k])
                        else:
                            angle = Entities['Functions'][fTag]['attributes']['theta']
                            sj = np.sin(angle/180.0*np.pi)
                            cj = np.cos(angle/180.0*np.pi)
                            p  = sj/Vs
                            ci = 1j*np.sqrt(Vp**2*p**2-1.0)

                            cons1 = 1.0/Vs/Vs-2.0*p*p
                            cons2 = 4.0*p*p*ci*cj/Vs/Vp

                            SP = 4.0*Vs/Vp*p*cj/Vs*cons1/(cons1*cons1 + cons2)
                            SS = (cons1*cons1-cons2)/(cons1*cons1 + cons2)
                            nn = np.arange(np.floor_divide(nt,2) + shift)
                            w  = 2*np.pi/t[-1]*nn

                            with concurrent.futures.ProcessPoolExecutor() as executor: 
                                for k, n in enumerate(nodes):
                                    x = Entities['Nodes'][n]['coords']
                                    U = SVbackground2DfieldCritical(Disp, SP, SS, Vp, Vs, cj, sj, p, w, x, xmin)
                                    V = SVbackground2DfieldCritical(Vels, SP, SS, Vp, Vs, cj, sj, p, w, x, xmin)
                                    A = SVbackground2DfieldCritical(Accel,SP, SS, Vp, Vs, cj, sj, p, w, x, xmin)        
                                    executor.submit(WriteDRMFile, dirName, funName, fTag, U, V, A, nt, 6, n, conditions[k])
                    elif funOption.upper() == 'RAYLEIGH':
                        mTag = Entities['Functions'][fTag]['attributes']['material']
                        Es  = Entities['Materials'][mTag]['attributes']['E']
                        nu  = Entities['Materials'][mTag]['attributes']['nu']
                        rho = Entities['Materials'][mTag]['attributes']['rho']
                        Vs  = np.sqrt(Es/2.0/rho/(1.0 + nu)) 

                        Vp = Vs*np.sqrt(2.0*(1.0 - nu)/(1.0 - 2.0*nu))
                        Vr = GetRayleighVelocity(Vs,nu)

                        nn = np.arange(np.floor_divide(nt,2) + shift)
                        w  = 2*np.pi/t[-1]*nn
                        kR = 1.0/Vr*w
                        qR = np.sqrt(1.0-(Vr/Vp)**2)*kR
                        sR = np.sqrt(1.0-(Vr/Vs)**2)*kR

                        with concurrent.futures.ProcessPoolExecutor() as executor:
                            for k, n in enumerate(nodes):
                                x = Entities['Nodes'][n]['coords']  
                                U = RHbackground2Dfield(Disp, t, Vp, Vs, Vr, w, kR, qR, sR, x0, xmin, x[0], x[1])
                                V = RHbackground2Dfield(Vels, t, Vp, Vs, Vr, w, kR, qR, sR, x0, xmin, x[0], x[1])
                                A = RHbackground2Dfield(Accel, t, Vp, Vs, Vr, w, kR, qR, sR, x0, xmin, x[0], x[1])
                                executor.submit(WriteDRMFile, dirName, funName, fTag, U, V, A, nt, 6, n, conditions[k])
                    else:
                        print('\x1B[31m ERROR \x1B[0m: The specified PLANEWAVE (2D) option (=%s) is not recognized' % funOption)
                elif nd == 3:
                    if 'phi' not in Entities['Functions'][fTag]['attributes']:
                        Entities['Functions'][fTag]['attributes']['phi'] = 0.0

                    #The wave direction on horizontal plane
                    azimut  = Entities['Functions'][fTag]['attributes']['phi']
                    phi     = azimut*np.pi/180.0
                    di      = np.array([np.cos(phi), np.sin(phi)]) 

                    if funOption.upper() == 'SH':
                        with concurrent.futures.ProcessPoolExecutor() as executor: 
                            for k, n in enumerate(nodes):
                                x = Entities['Nodes'][n]['coords']
                                U = SHbackground3Dfield(Disp, t, x, x0, xmin, di, nt, fTag)
                                V = SHbackground3Dfield(Vels, t, x, x0, xmin, di, nt, fTag)
                                A = SHbackground3Dfield(Accel, t, x, x0, xmin, di, nt, fTag)
                                executor.submit(WriteDRMFile, dirName, funName, fTag, U, V, A, nt, 9, n, conditions[k])
                    elif funOption.upper() == 'SV':
                        if 'theta' not in Entities['Functions'][fTag]['attributes']:
                            Entities['Functions'][fTag]['attributes']['theta'] = 0.0
                        if 'phi' not in Entities['Functions'][fTag]['attributes']:
                            Entities['Functions'][fTag]['attributes']['phi'] = 0.0
                        thetacr, Vs, Vp = ComputeSVCriticalAngle(fTag)

                        if thetacr < 1.0:
                            with concurrent.futures.ProcessPoolExecutor() as executor: 
                                for k, n in enumerate(nodes):
                                    x = Entities['Nodes'][n]['coords']
                                    U = SVbackground3Dfield(Disp, t, x, x0, xmin, di, nt, fTag)
                                    V = SVbackground3Dfield(Vels, t, x, x0, xmin, di, nt, fTag)
                                    A = SVbackground3Dfield(Accel, t, x, x0, xmin, di, nt, fTag)        
                                    executor.submit(WriteDRMFile, dirName, funName, fTag, U, V, A, nt, 9, n, conditions[k])
                        else:
                            angle = Entities['Functions'][fTag]['attributes']['theta']
                            sj = np.sin(angle/180.0*np.pi)
                            cj = np.cos(angle/180.0*np.pi)
                            p  = sj/Vs
                            ci = 1j*np.sqrt(Vp**2*p**2-1.0)

                            cons1 = 1.0/Vs/Vs-2.0*p*p
                            cons2 = 4.0*p*p*ci*cj/Vs/Vp

                            SP = 4.0*Vs/Vp*p*cj/Vs*cons1/(cons1*cons1 + cons2)
                            SS = (cons1*cons1-cons2)/(cons1*cons1 + cons2)
                            nn = np.arange(np.floor_divide(nt,2) + shift)
                            w  = 2*np.pi/t[-1]*nn

                            with concurrent.futures.ProcessPoolExecutor() as executor: 
                                for k, n in enumerate(nodes):
                                    x = Entities['Nodes'][n]['coords']
                                    U = SVbackground3DfieldCritical(Disp, SP, SS, Vp, Vs, cj, sj, p, w, di, x, xmin, nt)
                                    V = SVbackground3DfieldCritical(Vels, SP, SS, Vp, Vs, cj, sj, p, w, di, x, xmin, nt)
                                    A = SVbackground3DfieldCritical(Accel,SP, SS, Vp, Vs, cj, sj, p, w, di, x, xmin, nt)        
                                    executor.submit(WriteDRMFile, dirName, funName, fTag, U, V, A, nt, 9, n, conditions[k])
                    elif funOption.upper() == 'RAYLEIGH':
                        mTag = Entities['Functions'][fTag]['attributes']['material']
                        Es  = Entities['Materials'][mTag]['attributes']['E']
                        nu  = Entities['Materials'][mTag]['attributes']['nu']
                        rho = Entities['Materials'][mTag]['attributes']['rho']
                        Vs  = np.sqrt(Es/2.0/rho/(1.0 + nu)) 

                        Vp = Vs*np.sqrt(2.0*(1.0 - nu)/(1.0 - 2.0*nu))
                        Vr = GetRayleighVelocity(Vs,nu)

                        nn = np.arange(np.floor_divide(nt,2) + shift)
                        w  = 2*np.pi/t[-1]*nn
                        kR = 1.0/Vr*w
                        qR = np.sqrt(1.0-(Vr/Vp)**2)*kR
                        sR = np.sqrt(1.0-(Vr/Vs)**2)*kR

                        with concurrent.futures.ProcessPoolExecutor() as executor:
                            for k, n in enumerate(nodes):
                                x = Entities['Nodes'][n]['coords']
                                U = RHbackground3Dfield(Disp, Vp, Vs, Vr, w, kR, qR, sR, x0, xmin, di, x[0], x[1], x[2])
                                V = RHbackground3Dfield(Vels, Vp, Vs, Vr, w, kR, qR, sR, x0, xmin, di, x[0], x[1], x[2])
                                A = RHbackground3Dfield(Accel, Vp, Vs, Vr, w, kR, qR, sR, x0, xmin, di, x[0], x[1], x[2])
                                executor.submit(WriteDRMFile, dirName, funName, fTag, U, V, A, nt, 9, n, conditions[k])
                    else:
                        print('\x1B[31m ERROR \x1B[0m: The specified PLANEWAVE (3D) option (=%s) is not recognized' % funOption)
                else:
                    print('\x1B[31m ERROR \x1B[0m: The specified dimension (=%d) is not possible for DRM' % Options['dimensions'])

                #Update the load type after files were successfully created
                Entities['Loads'][lTag]['attributes']['type'] = 'GENERALWAVE'

                #Update the domain reduction time series path where files are located.
                Entities['Functions'][fTag]['attributes']['file'] = dirName + "/" + funName + "-" + str(fTag) + ".$.drm"
