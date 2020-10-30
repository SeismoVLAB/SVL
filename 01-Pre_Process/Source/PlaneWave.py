#!/usr/bin/python3

import os
import numpy as np
from scipy import integrate

def GetDerivative(x, dt):
    """
    """
    dx = np.gradient(x, dt)

    return dx

def GetIntegration(x, dt):
    """
    """
    n  = len(x)
    t  = np.linspace(0.0, (n-1)*dt, n)
    dx = integrate.cumtrapz(x, t, initial=0)

    return dx

def ComputeField(disp, vel, accel, dt, option):
    """
    This function computes the absent time history signal for the displacement, 
    velocity and acceleration depending on the option=ALL,DISP,VEL,ACCEL provied 
    by the user.
    """
    if  option.upper() == 'DISP':
        U = disp
        V = GetDerivative(U, dt)
        A = GetDerivative(V, dt)
    elif option.upper() == 'VEL':
        V = vel
        U = GetIntegration(V, dt)
        A = GetDerivative (V, dt)
    elif option.upper() == 'ACCEL':
        A = accel
        V = GetIntegration(A, dt)
        U = GetIntegration(V, dt)
    else:
        U = disp
        V = vel
        A = accel

    return U, V, A

def Compute2DBackgroundField(Values, t, Vs, Vp, theta_s, theta_p, X, X0, Xmin, nt):
    """
    """
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

    return U, V

def Compute3DBackgroundField(Values, t, Vs, Vp, theta_s, theta_p, di, X, X0, Xmin, nt):
    """
    This function generates the 3D DRM field to be specified at a particular node.
    """
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

    return U, V, W

def WriteFile(filepath, filename, Disp, Vels, Accel, nt, nc, n, option):
    """
    """
    #The output file.
    path = filepath + "/" + filename + "." + str(n) + ".drm"

    #Domain Reduction file format.
    DRMfile = open(path, "w+")
    DRMfile.write("%d %d %d\n" % (nt, nc, option))

    if nc == 6:
        for k in range(nt):
            DRMfile.write("%1.15e %1.15e %1.15e %1.15e %1.15e %1.15e\n" % (Disp[k,0], Disp[k,1], Vels[k,0], Vels[k,1], Accel[k,0], Accel[k,1]))
    elif nc == 9:
        for k in range(nt):
            DRMfile.write("%1.15e %1.15e %1.15e %1.15e %1.15e %1.15e %1.15e %1.15e %1.15e\n" % (Disp[k,0], Disp[k,1], Disp[k,2], Vels[k,0], Vels[k,1], Vels[k,2], Accel[k,0], Accel[k,1], Accel[k,2]))
    DRMfile.close()

def ParseDRMFile(Function):
    """
    This function parses the DRM input signal information provided in the *FUNCTION. 
    It reads the displacement, velocity, or accelertion input signal depending on 
    the option=ALL,DISP,VEL,ACCEL provied by the user.
    """
    #Open the provided file
    path = Function['PATH']

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
        conditions = np.zeros(nn, dtype=int)

        #Parse DRM node information.
        m = 0
        for k in range(1, nn+1):
            line = list(filter(None, lines[k].strip().split()))
            nodes[m], conditions[m] = int(line[0]), int(line[1])
            m += 1

        #Parse DRM input signal information.
        m = 0
        if option.upper() == 'ALL':
            for k in range(nn+1, nn+nt+1):
                line = list(filter(None, lines[k].strip().split()))
                time[m], disp[m], vels[m], accel[m] = float(line[0]), float(line[1]), float(line[2]), float(line[3])
                m += 1
        elif option.upper() == 'DISP':
            for k in range(nn+1, nn+nt+1):
                line = list(filter(None, lines[k].strip().split()))
                time[m], disp[m] = float(line[0]), float(line[1])
                m += 1
        elif option.upper() == 'VEL':
            for k in range(nn+1, nn+nt+1):
                line = list(filter(None, lines[k].strip().split()))
                time[m], vels[m] = float(line[0]), float(line[1])
                m += 1
        elif option.upper() == 'ACCEL':
            for k in range(nn+1, nn+nt+1):
                line = list(filter(None, lines[k].strip().split()))
                time[m], accel[m] = float(line[0]), float(line[1])
                m += 1

    return nodes, conditions, time, disp, vels, accel, dt, option

def Driver(User, Point, Function, time, Disp, Vels, Accel, nodes, conditions):
    """
    This function writes (generates) the DRM input files to be used in the 
    Run-Analysis. The function writes the input information regarding 
    displacement, velocity and acceleration at each DRM node.
    """
    #Planar wave in homogeneous half-space parameters
    angle  = Function['THETA'] 
    azimut = Function['PHI'] 
    Vs     = Function['VS'] 
    nu     = Function['NU'] 
    x0     = Function['X0']
    xmin   = Function['XMIN']

    phi     = azimut*np.pi/180.0
    Vp      = Vs*np.sqrt(2.0*(1.0 - nu)/(1.0 - 2.0*nu))
    theta_s = angle*np.pi/180.0
    theta_p = np.arcsin(Vp/Vs*np.sin(theta_s))

    #Creates the DRM Directory
    dirName = User['FOLDER'] + 'DRM'
    if not os.path.exists(dirName):
        os.mkdir(dirName)

    #Generate the input files for Domain Reduction Method.
    if User['DIMENSION'] == 2:
        nt = len(time)
        U = np.zeros((nt,2))
        V = np.zeros((nt,2))
        A = np.zeros((nt,2))

        if 'UNDEFINED' in Function:
            for k in Function['UNDEFINED']:
                WriteFile(dirName, Function['NAME'], U, V, A, nt, 6, k, 0)

        for k, n in enumerate(nodes):
            x  = Point[n]['COORDINATES']
   
            U[:,0] , U[:,1] = Compute2DBackgroundField(Disp, time, Vs, Vp, theta_s, theta_p, x, x0, xmin, nt)
            V[:,0] , V[:,1] = Compute2DBackgroundField(Vels, time, Vs, Vp, theta_s, theta_p, x, x0, xmin, nt)
            A[:,0] , A[:,1] = Compute2DBackgroundField(Accel, time, Vs, Vp, theta_s, theta_p, x, x0, xmin, nt)

            WriteFile(dirName, Function['NAME'], U, V, A, nt, 6, n, conditions[k])
    elif User['DIMENSION'] == 3:
        nt = len(time)
        U = np.zeros((nt,3))
        V = np.zeros((nt,3))
        A = np.zeros((nt,3))
        d = np.array([np.cos(phi), np.sin(phi)]) 

        if 'UNDEFINED' in Function:
            for k in Function['UNDEFINED']:
                WriteFile(dirName, Function['NAME'], U, V, A, nt, 9, k, 0)

        for k, n in enumerate(nodes):
            x  = Point[n]['COORDINATES']

            U[:,0] , U[:,1] , U[:,2] = Compute3DBackgroundField(Disp, time, Vs, Vp, theta_s, theta_p, d, x, x0, xmin, nt)
            V[:,0] , V[:,1] , V[:,2] = Compute3DBackgroundField(Vels, time, Vs, Vp, theta_s, theta_p, d, x, x0, xmin, nt)
            A[:,0] , A[:,1] , A[:,2] = Compute3DBackgroundField(Accel, time, Vs, Vp, theta_s, theta_p, d, x, x0, xmin, nt)

            WriteFile(dirName, Function['NAME'], U, V, A, nt, 9, n, conditions[k])

    #Update the path where files are located.
    Function['PATH'] = dirName + "/" + Function['NAME'] + ".$.drm"
