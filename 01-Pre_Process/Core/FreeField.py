import numpy as np
import matplotlib.pyplot as plt
from Core.PlaneWave import *

def GetTimeSeries(vals, dt, option):
    '''
    This function computes the missing time series for the displacement, 
    velocity and acceleration depending on the option provied by the user\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2021

    Parameters
    ----------
    vals   : array
        The provided time series
    dt     : float
        The time step for the given time series 
    option : str
        User's time series data, option=DISP, VEL, or ACCEL

    Returns
    -------
    Disp, Vels, Accel the time series
    '''
    #Compute other fields according to option
    if  option.upper() == 'DISP':
        Disp  = vals
        Vels  = GetDerivative(vals, dt)
        Accel = GetDerivative(Vels, dt)
    elif option.upper() == 'VEL':
        Disp  = GetIntegration(vals, dt)
        Vels  = vals
        Accel = GetDerivative (vals, dt)
    elif option.upper() == 'ACCEL':
        Vels  = GetIntegration(vals, dt)
        Disp  = GetIntegration(Vels, dt)
        Accel = vals
    
    return Disp, Vels, Accel

def Compute2DFreeFieldBoundaries(x, y, signal, dt, xmin, layers, beta, rho, nu, option, df=0.2, cof=30.0):
    '''
    This function calculates the free-field forces required to be applied at the boundaries.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Feiruo (Flora) Xia, Eugene Loh, and Yaozhong Shi
    
    Parameters
    ----------
    x  : array
        Vector with the x-coordinates to compute the free-field forces
    y  : array
        Vector with the y-coordinates to compute the free-field forces
    signal  : array
        Vector with the time series values from which the forces will be computed
    dt  : float
        The time step of the time-series
    TODO: COMPLETE ACCORDINGLY
    '''
    #Unpack predefines SV-Wave 
    angle = 0.0001
    y0 = xmin[1]
    x0 = xmin[0]

    #Creates the fun dictionary required for DataPreprocessing
    fun = { 'option': 'SV', 'df': df, 'CutOffFrequency': cof}

    #Computes the input displacement, velocities and accelerations.
    Disp, Vels, Accel = GetTimeSeries(signal, dt, option)

    #Transform time-series in frequency domain and computes variables required to compute reponses layer interface
    nt = len(signal)
    ufull, vfull, afull, layers, beta, rho, nu, wVec, p, s, h, mu, aSP, phaseVelIn, sinTheta, N, Nt = DataPreprocessing(Disp, Vels, Accel, layers, beta, rho, nu, angle, y0, nt, dt, fun)

    #Compute Interface responses
    uInterface = SoilInterfaceResponse(ufull, wVec, p, s, h, mu, aSP, phaseVelIn, sinTheta, N)
    vInterface = SoilInterfaceResponse(vfull, wVec, p, s, h, mu, aSP, phaseVelIn, sinTheta, N)
    aInterface = SoilInterfaceResponse(afull, wVec, p, s, h, mu, aSP, phaseVelIn, sinTheta, N)

    for (xk, yk) in zip(x,y):
        #Compute Displacement, Velocity and acceleration in time domain at (x,y) coordinate
        U = PSVbackground2Dfield(uInterface, layers, wVec, p, s, h, mu, aSP, phaseVelIn, sinTheta, N, Nt, x0, xk, yk)
        V = PSVbackground2Dfield(vInterface, layers, wVec, p, s, h, mu, aSP, phaseVelIn, sinTheta, N, Nt, x0, xk, yk)
        A = PSVbackground2Dfield(aInterface, layers, wVec, p, s, h, mu, aSP, phaseVelIn, sinTheta, N, Nt, x0, xk, yk)

        #Computes the Reaction forces at (x,y) coordinate 
        #TODO: Use U,V,A for this

        #Saves the Time-Series to be loaded as PointLoad
        #TODO:

        #Computes the dashpot coefficients associated to this coordinate
        print('running!!!')

def Compute3DFreeFieldBoundaries(x, y, signal, dt, xmin, layers, beta, rho, nu, option, df=0.2, cof=30.0):
    '''
    '''
    pass

def Find2DBoundaries(xleft, xright, ybottom):
    '''
    '''
    Boundaries = {'left': [], 'right': [], 'bottom': []}
    return Boundaries