import os
import json
import subprocess
import numpy as np
import collections.abc
from json import JSONEncoder

#https://stackoverflow.com/questions/50472095/how-can-i-solve-a-system-of-linear-equations-in-python-faster-than-with-numpy-li
import time
def lstsq(A, b):
    AA = A.T @ A
    bA = b @ A
    D, U = np.linalg.eigh(AA)
    Ap = (U * np.sqrt(D)).T
    bp = bA @ U / np.sqrt(D)
    return np.linalg.lstsq(Ap, bp, rcond=None)

class NumpyArrayEncoder(JSONEncoder):
    """
    This class is a Custom JSON Encoder to Serialize NumPy ndarray\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2021
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

def ParametersToScreen(un, iteration):
    """
    This function prints on screen the variable values during each 
    Ensemble Kalman Inversion Updating iteration\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2021

    Parameters
    ----------
    un : list
        The vector with identified values
    iteration : int
        The EnKI-FEM iteration step

    Returns
    -------
    None
    """
    ToScreen = ' [' + str(iteration) + '] Parameter values:'
    nParameters = len(un)
    for k in range(nParameters):
        ToScreen += ' u[' + str(k) + '] = ' + str(np.exp(un[k]))
        if k != nParameters-1:
            ToScreen += ','
    print(ToScreen)

def FunctionWrapper(function, uk, ensemble, EnsembleOutputFile, k):
    """
    This functions is a wrapper that takes a function defined by the user
    and change the ensemble values\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2021

    Parameters
    ----------
    uk : list
        The vector with identified values
    ensemble: dict
        The dictionary with the ensemble information to be updated
    EnsembleOutputFile: str
        The full file path where the ensemble JSON model will be stored
    k : int
        The EnKI-FEM iteration step

    Returns
    -------
    None
    """
    function(uk, ensemble, EnsembleOutputFile, k)

def update(d, u):
    """
    This function update the values of a dictionary using another dictionary\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2021

    Parameters
    ----------
    d : dict
        The dictionary to be updated
    u : dict
        The dictionary that holds the key, value to update

    Returns
    -------
    None
    """
    for k, v in u.items():
        if isinstance(v, collections.abc.Mapping):
            d[k] = update(d.get(k, {}), v)
        else:
            d[k] = v
    return d

def GetMeasurements(FilePath, OutputFile, dofs, nSkip, k=0):
    """
    This function reads the solution from a file and return the values in a vector\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2021

    Parameters
    ----------
    FilePath : str
        The full path of the file to be opened
    OutputFile : str
        The file name to be opened
    dofs : list
        The list of column indexes to be read
    nSkip : int
        The number of lines to skip when opening the file
    k : int
        The ensemble (particle) number 

    Returns
    -------
    Measurements: array
        The vectorized array that contains the requested information
    """
    MeasurementFile = FilePath + '/' + OutputFile.replace('$', str(k))
    Measurements = np.loadtxt(MeasurementFile, dtype='float', skiprows=nSkip)
    Measurements = Measurements[:,dofs]
    Measurements = Measurements.flatten('F')
    return Measurements

def RunForwardModel(SVLFolder, FilePath, FileName, k):
    """
    This function runs SeismoVLAB using a certain file\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2021

    Parameters
    ----------
    SVLFolder : str
        The full path where SeismoVLAB.exe is located
    FilePath : str
        The full path where the input JSON is located
    FileName : str
        The name of the input JSON file to be executed
    k : int
        The ensemble (particle) number 

    Returns
    -------
    None
    """
    ToRun = SVLFolder + '/SeismoVLAB.exe -dir ' + FilePath + ' -file ' + FileName.replace('$', str(k))
    subprocess.check_output(ToRun, shell=True)

def WriteEnsemble(ensemble, FileFolder, oldFileName, newFileName, k):
    """
    This function writes the SVL model in JSON format\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2021

    Parameters
    ----------
    ensemble : dict
        The dictionary that holds the key, value of the ensemble to update the model
    FileFolder : str
        The full path where the input JSON (background model) is located
    oldFileName : str
        The name of the input JSON file (original model) to be loaded
    newFileName : str
        The name of the output JSON file (new ensemble) to be executed
    k : int
        The ensemble (particle) number 

    Returns
    -------
    None
    """
    #Opens the SVL model
    FullFilePath = FileFolder + '/' + oldFileName
    with open(FullFilePath, 'r') as myfile:
        JSONdata = myfile.read()
    json_object = json.loads(JSONdata)
    myfile. close()

    #Modifies the requested field
    update(json_object, ensemble)

    #Serializing json and Writes the file in JSON format
    JSONdata = json.dumps(json_object, cls=NumpyArrayEncoder, indent=4)

    FullFilePath = FileFolder + '/' + newFileName
    FullFilePath = FullFilePath.replace('$', str(k))
    with open(FullFilePath, "w") as outfile: 
        outfile.write(JSONdata)

def EnsembleKalmanInversion(Options):
    """
    This function performs the Ensemble Kalman Inversion Finite Element Model 
    Updating framework for parameter estimation\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2021

    Parameters
    ----------
    Options : dict
        The dictionary that holds the user's defined configuration parameters. The 
        dictionary must have the following fields:
            "nLines"    : Number of lines to skip when reading the observation (true responses) file
            "nTimeStep" : Number of time steps that the observation (true responses) file has
            "nParticles": Number of ensembles or particules to be used
            "threshold" : The cut-off value when generating the covariance matrix
            "maxiter"   : Maximum number of iterations o be used in the EnKI-FEM framework 
            "dofs"      : Columns (indexes) to considered when reading the observation (true responses) file
            "u0"        : The vector of initial guess for the parameters
            "cov"       : The vector of coefficient of variation for each parameters
            "ensemble"  : A dictionary that contains the field to be changed
            "function"  : A user defined function that takes the ensemble an update its values
            "files":
                "BackgroundFile"  : The JSON file that contains the original model
                "ObservationFile" : The files where the true responses are going to be loaded
                "SeismoVLAB"      : The full path where the SeismoVLAB.exe is located
                "SolutionFolder"  : The full folder path where the solution will be stored
                "EnsembleFolder"  : The full folder path where the ensemble JSON files will be stored

    Returns
    -------
    None
    """
    #UNPACK THE ENSEMBLE KALMAN INVERSION PARAMETERS
    print(' Starting Ensemble Kalman Inversion Finite Element Updating Algorithm:')

    #Read file parameters
    dofs      = Options['dofs']
    nLines    = Options['nLines']
    threshold = Options['threshold']

    #Generation of Ensembles
    nTimeStep  = Options['nTimeStep']
    nParticles = Options['nParticles']
    ensemble   = Options['ensemble']
    function   = Options['function']

    #Vector of Parameter to be identified
    u0  = Options['u0']
    cov = Options['cov']
    maxiter = Options['maxiter']

    #Files and Folders where information will be stored
    EnsembleFile = 'Ensemble_$.0.json'
    BackgroundFile = Options['files']['BackgroundFile']
    ObservationFile = Options['files']['ObservationFile']
    EnsembleOutputFile = ObservationFile.replace('.', '-$.', 1) 

    SeismoVLABFolder = Options['files']['SeismoVLAB']
    EnsembleFolder = Options['files']['EnsembleFolder']
    SolutionFolder = Options['files']['SolutionFolder']
    ObservationFolder = Options['files']['SolutionFolder']

    #THE ENSEMBLE KALMAN INVERSION FINITE ELEMENT UPDATING FRAMEWORK
    nControl = len(dofs)
    nParameters = len(u0)
    N = nControl*(nTimeStep - 1)

    #Generation of Ensembles
    Particle = np.random.uniform(-1, 1, (nParameters,nParticles))
    for k in range(nParameters):
        Particle[k,:] = u0[k]*(1.0 + cov[k]*Particle[k,:])

    #Gets the observation data 
    y_measured = GetMeasurements(ObservationFolder, ObservationFile, dofs, nLines)

    #Covariance Matrix
    Rmatrix = np.square(y_measured)
    Rmatrix[Rmatrix < threshold] = threshold 
    Rmatrix = np.diag(Rmatrix)

    #Allocate memory for Ensemble
    uhat = np.log(Particle)
    what = np.zeros((N, nParticles))

    #Starts the Ensemble Kalman Inversion Algorithm
    s = 1.0
    iteration = 0
    while iteration < maxiter:
        #Observation Operator for each Ensamble.
        for k in range(nParticles):
            uk = np.exp(uhat[:,k])

            #FunctionWrapper(function, uk, ensemble, EnsembleOutputFile, k)
            #WriteEnsemble(ensemble, EnsembleFolder, BackgroundFile, EnsembleFile, k)
            #RunForwardModel(SeismoVLABFolder, EnsembleFolder, EnsembleFile, k)

            what[:,k] = GetMeasurements(SolutionFolder, EnsembleOutputFile, dofs, nLines, k)

        #Perturbation on the Measurements to Construct Artificial Data
        yperturb = np.zeros((N, nParticles)) #np.random.multivariate_normal(np.zeros(N), Rmatrix, nParticles).T

        #Computes the averages for Parameters and Observation Operator
        wnp1 = what.mean(axis=1)
        unp1 = uhat.mean(axis=1)

        #Computes the Empirical Covariance Matrix
        Cuw = np.zeros((nParameters, N))
        Cww = np.zeros((N, N))

        for k in range(nParticles):
            Cuw += 1.0/(nParticles - 1.0)*np.outer(uhat[:,k] - unp1, what[:,k] - wnp1)
            Cww += 1.0/(nParticles - 1.0)*np.outer(what[:,k] - wnp1, what[:,k] - wnp1)

        #Kalman Gain Matrix
        print('starting')
        start = time.time()
        K = np.linalg.solve(Cww + Rmatrix, Cuw.T).T
        #K, *info_acc = lstsq(Cww + Rmatrix, Cuw.T).T
        end = time.time()
        print('ending in' + str(end - start) + ' [s]')

        #Update the Ensembles for next iteration
        for k in range(nParticles):
            dy = y_measured + s*yperturb[:,k] - what[:,k]
            uhat[:,k] += np.dot(K,dy)

        #Prints the identified parameter values for this iteration
        ParametersToScreen(unp1, iteration)

        iteration += 1
