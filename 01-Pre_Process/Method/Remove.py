#!/usr/bin/python3
# -*- coding: Utf-8 -*-

import numpy as np
from Core.Utilities import debugInfo
from Core.Definitions import Entities, Options

def delNode(tag=None):
    """
    Deletes an specified Node\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of the Node to be deleted, i.e., tag > -1

    Returns
    -------
    bool
        Whether the deletion was successful (True) of failed (False)
    """
    #Deletes the requested node in Entities
    if tag in Entities['Nodes']:
        del Entities['Nodes'][tag]
        return True
    else:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Node[%s] cannot be deleted.' %(info.filename,info.lineno,tag))
        return False

def delMass(tag=None):
    """
    Deletes mass associated to a Node\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of the Node to add mass, i.e., tag > -1

    Returns
    -------
    bool
        Whether the mass deletion was successful (True) of failed (False)
    """
    #Deletes the requested node in Entities
    if tag in Entities['Masses']:
        del Entities['Masses'][tag]
        return True
    else:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Mass associated to Node[%s] cannot be deleted.' %(info.filename,info.lineno,tag))
        return False

def delRestrain(tag=None, dof=None):
    """
    Deletes an specified restrain\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of the node, i.e., tag > -1
    dof : int or list
        Number or list of degree of freedom to be freed.

    Returns
    -------
    bool
        Whether the constraint/s was/were removed (True) of failed (False)
    """
    if not tag or not dof:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d the variable delRestrain(tag=?,dof=?) must be specified.' %(info.filename,info.lineno))
        return False

    #Transform single value into a list
    if isinstance(dof, int):
        dof = list([dof])

    #Transform to zero base index
    dof = [k-1 for k in dof]

    #Assign restriction depending on one-value or list
    ndof = Entities['Nodes'][tag]['ndof']

    for n in dof:
        if n > ndof:
            print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d dof[%d] > %d (Node[%d]).' %(info.filename,info.lineno,n+1,ndof,tag))
            return False
        Entities['Nodes'][tag]['freedof'][n] = 0
    return True

def delConstraint(tag=None):
    """
    Deletes Constraint associated to a Node\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of the Constraint applied to Node, i.e., tag < -1

    Returns
    -------
    bool
        Whether the constraint deletion was successful (True) of failed (False)
    """
    #Deletes the requested constraint in Entities
    if tag in Entities['Constraints']:
        #Constrained degree of freedom in now free
        ntag = Entities['Constraints'][tag]['stag']
        dof  = Entities['Constraints'][tag]['sdof']
        Entities['Nodes'][ntag]['freedof'][dof] = 0

        del Entities['Constraints'][tag]
        return True
    else:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Constraint associated to Node[%s] cannot be deleted.' %(info.filename,info.lineno,tag))
        return False

def delRigidLink(tag=None):
    """
    Deletes an specified rigid link constraint\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2021

    Parameters
    ----------
    tag : int
        The identifier of the rigid link, i.e., tag > 1

    Returns
    -------
    bool
        Whether the deletion was successful (True) of failed (False)
    """
    #Check whether the Material exists
    if tag in Entities['RigidLinks']:
        del Entities['RigidLinks'][tag]
        return True
    else:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d RigidLinks[%s] cannot be deleted.' %(info.filename,info.lineno,tag))
        return False

def delDiaphragm(tag=None):
    """
    Deletes an specified diaphragm constraint\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of the diaphragm, i.e., tag > 1 (different from Nodes)

    Returns
    -------
    bool
        Whether the deletion was successful (True) of failed (False)
    """
    #Check whether the Material exists
    if tag in Entities['Diaphragms']:
        del Entities['Diaphragms'][tag]
        return True
    else:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Diaphragm[%s] cannot be deleted.' %(info.filename,info.lineno,tag))
        return False

def delRigidBody(tag=None):
    """
    Deletes an specified rigid body constraint\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of the rigid body, i.e., tag > 1 (different from Nodes)

    Returns
    -------
    bool
        Whether the deletion was successful (True) of failed (False)
    """
    #Check whether the Material exists
    if tag in Entities['RigidBodies']:
        del Entities['RigidBodies'][tag]
        return True
    else:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d RigidBodies[%s] cannot be deleted.' %(info.filename,info.lineno,tag))
        return False

def delSupportMotion(tag=None):
    """
    Deletes an specified diaphragm constraint\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of the support motion, i.e., tag > 1

    Returns
    -------
    bool
        Whether the deletion was successful (True) of failed (False)
    """
    #Check whether the Material exists
    if tag in Entities['Supports']:
        del Entities['Supports'][tag]
        return True
    else:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Supports[%d] cannot be deleted.' %(info.filename,info.lineno,tag))
        return False

def delMaterial(tag=None):
    """
    Deletes an specified material in the model\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of the Material to be deleted, i.e., tag > -1

    Returns
    -------
    int
        Whether the deletion was successful (True) of failed (False)
    """
    #Deletes the requested material in Entities
    if tag in Entities['Materials']:
        del Entities['Materials'][tag]
        return True
    else:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Element[%s] cannot be deleted.' %(info.filename,info.lineno,tag))
        return False

def delSection(tag=None):
    """
    Deletes an specified section in the model\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of the Section to be deleted, i.e., tag > -1

    Returns
    -------
    int
        Whether the deletion was successful (True) of failed (False)
    """
    #Deletes the requested material in Entities
    if tag in Entities['Sections']:
        del Entities['Sections'][tag]
        return True
    else:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Element[%d] cannot be deleted.' %(info.filename,info.lineno,tag))
        return False

def delElement(tag=None):
    """
    Deletes an specified Element in the model\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of the Element to be deleted, i.e., tag > -1

    Returns
    -------
    int
        Whether the deletion was successful (True) of failed (False)
    """
    #Deletes the requested element in Entities
    if tag in Entities['Elements']:
        del Entities['Elements'][tag]
        return True
    else:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Element[%s] cannot be deleted.' %(info.filename,info.lineno,tag))
        return False

def delSurface(tag=None):
    """
    Deletes an specified element's surface\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of the Surface to be deleted, i.e., tag > -1

    Returns
    -------
    int
        Whether the deletion was successful (True) of failed (False)
    """
    #Deletes the requested material in Entities
    if tag in Entities['Surfaces']:
        del Entities['Surfaces'][tag]
        return True
    else:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Surface[%s] cannot be deleted.' %(info.filename,info.lineno,tag))
        return False

def delDamping(tag=None):
    """
    Deletes an specified damping model\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of the Damping to be deleted, i.e., tag > -1

    Returns
    -------
    int
        Whether the deletion was successful (True) of failed (False)
    """
    #Deletes the requested Damping in Entities
    if tag in Entities['Dampings']:
        del Entities['Dampings'][tag]
        return True
    else:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Damping[%s] cannot be deleted.' %(info.filename,info.lineno,tag))
        return False

def delFunction(tag=None):
    """
    Deletes the specified Function\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of the Load to be deleted, i.e., tag > -1

    Returns
    -------
    int
        Whether the deletion was successful (True) of failed (False)
    """
    #Deletes the requested Loads in Entities
    if tag in Entities['Functions']:
        del Entities['Functions'][tag]
        return True
    else:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Function[%s] cannot be deleted.' %(info.filename,info.lineno,tag))
        return False

def delLoad(tag=None):
    """
    Deletes the specified Load\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of the Load to be deleted, i.e., tag > -1

    Returns
    -------
    int
        Whether the deletion was successful (True) of failed (False)
    """
    #Deletes the requested Loads in Entities
    if tag in Entities['Loads']:
        del Entities['Loads'][tag]
        return True
    else:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Load[%s] cannot be deleted.' %(info.filename,info.lineno,tag))
        return False

def delCombinationCase(tag=None):
    """
    Deletes the specified Combination\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of the Combination to be deleted, i.e., tag > -1

    Returns
    -------
    bool
        Whether the deletion was successful (True) of failed (False)
    """
    #Deletes the requested Combinations in Entities
    if tag in Entities['Combinations']:
        del Entities['Combinations'][tag]
        return True
    else:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Combination[%s] cannot be deleted.' %(info.filename,info.lineno,tag))
        return False

def delRecorder(tag=None):
    """
    Deletes the specified Recorder\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of this recorder, i.e., tag > -1

    Returns
    -------
    bool
        Whether the addition was successful (True) of failed (False)
    """
    #Deletes the Solver in Entities
    if tag in Entities['Recorders']:
        del Entities['Recorders'][tag]
        return True
    else:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Recorder[%s] cannot be deleted.' %(info.filename,info.lineno,tag))

def delSolver(tag=None):
    """
    Deletes the specified Solver\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of this solver, i.e., tag > -1

    Returns
    -------
    bool
        Whether the addition was successful (True) of failed (False)
    """
    #Deletes the Solver in Entities
    if tag in Entities['Solvers']:
        del Entities['Solvers'][tag]
        return True
    else:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Solver[%s] cannot be deleted.' %(info.filename,info.lineno,tag))
        return False

def delAlgorithm(tag=None):
    """
    Deletes the specified Algorithm\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of this algorithm, i.e., tag > -1

    Returns
    -------
    bool
        Whether the addition was successful (True) of failed (False)
    """
    #Deletes the Solver in Entities
    if tag in Entities['Algorithms']:
        del Entities['Algorithms'][tag]
        return True
    else:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Algorithm[%s] cannot be deleted.' %(info.filename,info.lineno,tag))
        return False

def delIntegrator(tag=None):
    """
    Deletes the specified Integrator\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of this integrator, i.e., tag > -1

    Returns
    -------
    bool
        Whether the addition was successful (True) of failed (False)
    """
    #Deletes the Solver in Entities
    if tag in Entities['Integrators']:
        del Entities['Integrators'][tag]
        return True
    else:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Integrator[%s] cannot be deleted.' %(info.filename,info.lineno,tag))
        return False

def delAnalysis(tag=None):
    """
    Deletes the specified Analysis\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of this analysis, i.e., tag > -1

    Returns
    -------
    bool
        Whether the addition was successful (True) of failed (False)
    """
    #Deletes the Solver in Entities
    if tag in Entities['Analyses']:
        del Entities['Analyses'][tag]
        return True
    else:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Analysis[%s] cannot be deleted.' %(info.filename,info.lineno,tag))
        return False

def delSimulation(tag=None):
    """
    Deletes the specified Simulation\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of the Simulation to be deleted, i.e., tag > -1

    Returns
    -------
    bool
        Whether the deletion was successful (True) of failed (False)
    """
    #Deletes the Simulation in Entities
    if tag in Entities['Simulations']:
        del Entities['Simulations'][tag]
        return True
    else:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Simulation[%s] cannot be deleted.' %(info.filename,info.lineno,tag))
        return False