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

def Get2DFreeFieldMotion(x,):
    '''
    This function computes the displecement, velocity and acceleration at a
    given coordinate.
    '''
    pass

def Get2DFreeFieldReaction():
    '''
    This function computes the reaction forces at a given coordinate.
    '''
    pass 

def Compute2DFreeFieldBoundaries(x, y, vals, dt, xmin, layers, beta, rho, nu,fun):

    #plt.plot(t, vals)
    #plt.show()

    #Computes the input displacement, velocities and accelerations.
    Disp, Vels, Accel = GetTimeSeries(vals, dt, option)

    #plt.plot(t, Disp)
    #plt.show()

    #plt.plot(t, Accel)
    #plt.show()

    angle = 0.0001
    nt = len(vals)
    y0 = xmin[1]

    ufull, vfull, afull, layers, beta, rho, nu, wVec, p, s, h, mu, aSP, phaseVelIn, sinTheta, N, Nt = DataPreprocessing(Disp, Vels, Accel, layers, beta, rho, nu, angle, y0, nt, dt, fun)

    #Compute Interface responses
    uInterface = SoilInterfaceResponse(ufull, wVec, p, s, h, mu, aSP, phaseVelIn, sinTheta, N)
    vInterface = SoilInterfaceResponse(vfull, wVec, p, s, h, mu, aSP, phaseVelIn, sinTheta, N)
    aInterface = SoilInterfaceResponse(afull, wVec, p, s, h, mu, aSP, phaseVelIn, sinTheta, N)

    x0 = xmin[0]
    t = np.linspace(0.0, (Nt-1)*dt, Nt)

    for (xk, yk) in zip(x,y):
        U = PSVbackground2Dfield(uInterface, layers, wVec, p, s, h, mu, aSP, phaseVelIn, sinTheta, N, Nt, x0, xk, yk)
        V = PSVbackground2Dfield(vInterface, layers, wVec, p, s, h, mu, aSP, phaseVelIn, sinTheta, N, Nt, x0, xk, yk)
        A = PSVbackground2Dfield(aInterface, layers, wVec, p, s, h, mu, aSP, phaseVelIn, sinTheta, N, Nt, x0, xk, yk)
        #Use this info to get the velocities for example
        #compute reaction
        plt.plot(t, V[:,0])
        plt.show()

option = 'VEL'
dt = 0.001
parameters = {
    'to': 1.5,
    'f0': 1.0,
    'dt': dt,
    'Ts': 10.0,
    'Ap': 1.0
}

       

fun = {
    'option'  : 'SV',
    'df'      : 0.2,
    'CutOffFrequency': 20.0
}

layers = np.array([0.0, -20.0, -50.0]) 
beta = np.array([100.0, 200.0, 500.0]) 
rho = np.array([1500.0,1500.0,1500.0]) 
nu  = np.array([0.33,0.33,0.33]) 
xmin = [0.0,-100.0]
x = [0.0,  0.0,  0.0,  0.0]
y = [0.0,-20.0,-50.0,-80.0]
Compute2DFreeFieldBoundaries(x,y,vals,dt,xmin,layers,beta,rho,nu,fun)

'''

#
ufull, vfull, afull, Layers, beta, rho, nu, wVec, p, s, h, mu, aSP, phaseVelIn, sinTheta, N, Nt = DataPreprocessing(Disp, Vels, Accel, Layers, beta, rho, nu, angle, xmin[1], nt, dt, fun)

#Compute Interface responses
uInterface = SoilInterfaceResponse(ufull, wVec, p, s, h, mu, aSP, phaseVelIn, sinTheta, N)
vInterface = SoilInterfaceResponse(vfull, wVec, p, s, h, mu, aSP, phaseVelIn, sinTheta, N)
aInterface = SoilInterfaceResponse(afull, wVec, p, s, h, mu, aSP, phaseVelIn, sinTheta, N)

nmat = len(fun['material'])
beta = np.zeros((nmat,))
rho = np.zeros((nmat,)) 
nu = np.zeros((nmat,))


fun = {
    'option'  : 'SV',
    'layer'   : [0.0, -20.0, -50.0],
    'material': [1,2,3],
    'theta'   : 0.0,
    'phi'     : 0.0,
    'df'      : 0.2,
    'CutOffFrequency': 20.0
}

#SoilInterfaceResponse, DataPreprocessing,PSVbackground2Dfield, PSVbackground3Dfield
'''
