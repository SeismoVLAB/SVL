#!/usr/bin/python3
# -*- coding: Utf-8 -*-

import numpy as np
from Method.Remove import delConstraint
from Core.Utilities import debugInfo
from Core.Definitions import Entities, Options, ConvergeTest, SolverOption

def addNode(tag=np.nan, ndof=0, coords=[], freedof=[], totaldof=[]):
    """
    Appends a new Node to the model\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of this node, i.e., tag > -1
    ndof : int
        Number of degree of freedom, i.e., ndof > 0
    coords : list
        Coordinates of the node, i.e., (x,y) or (x,y,z)
    freedof : list
        The free (with restrain/constaint) degree of freedom numbering of this node
    totaldof : list
        The total (without restrain/constaint) degree of freedom numbering of this node

    Returns
    -------
    bool
        Whether the addition was successful (True) of failed (False)
    """
    #Check whether the Node exists
    if tag not in Entities['Nodes']:
        #Create the total/free list
        free  = np.array(freedof, dtype=int)  if freedof  else np.zeros(ndof, dtype=int)
        total = np.array(totaldof, dtype=int) if totaldof else np.zeros(ndof, dtype=int)

        Entities['Nodes'][tag] = {'ndof': ndof, 'freedof': free, 'totaldof': total, 'coords': np.array(coords)}
        return True
    else:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Node[%d] has been already defined.' %(info.filename,info.lineno,tag))
        return False

def addMass(tag=np.nan, dof=[], vals=[]):
    """
    Add mass to a Node in specific degree-of-freedom\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of the Node to add mass, i.e., tag > -1
    dof : int or list
        Degree of freedom number where the mass will be added, i.e., ndof > 0
    vals : list
        Mass values to be added to specific degree of freedom

    Returns
    -------
    bool
        Whether the mass addition was successful (True) of failed (False)
    """
    #Transform single value into a list
    if isinstance(dof, int):
        dof = list([dof])
        vals = list([vals])

    #Check whether the Node exists
    if tag in Entities['Nodes']:
        ndof = Entities['Nodes'][tag]['ndof']

        if ndof > 0:
            mass = np.zeros(ndof, dtype=float)

            #Construct the mass vector
            for k, n in enumerate(dof):
                mass[k] = vals[k]

            #Check whether the Mass has been defined
            if tag in Entities['Masses']:
                Entities['Masses'][tag]['mass'] += mass
            else:
                Entities['Masses'][tag] = {'ndof': ndof, 'mass': mass}
            return True
        else:
            info = debugInfo(2) 
            print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Mass[%d] is applied to Node[%d] with ndof=%d.' %(info.filename,info.lineno,tag,tag,ndof))
            return False
    else:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Node[%d] has not been created.' %(info.filename,info.lineno,tag))
        return False

def addRestrain(tag=np.nan, dof=[]):
    """
    Specifies a restrain applied to a Node\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of the node, i.e., tag > -1
    dof : int or list
        Number or list of degree of freedom to be restrained.

    Returns
    -------
    bool
        Whether the Restrain/s was/were successful (True) of failed (False)
    """
    #Transform single value into a list
    if isinstance(dof, int):
        dof = list([dof])

    #Transform to zero base index
    dof = [k-1 for k in dof]

    #Assign restriction depending on one-value or list
    ndof = Entities['Nodes'][tag]['ndof']

    for n in dof:
        if n > ndof:
            info = debugInfo(2)
            print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d for Node[%d]: dof[%d] > %d (out-of-bound).' %(info.filename,info.lineno,n,ndof,tag))
            return False
        #The constraint is not applied if degree-of-freedom was restrained
        if n < -1:
            info = debugInfo(2)
            print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Node[%d] had Constraint[%d] applied, the restrain will be enforced.' %(info.filename,info.lineno,tag,n))
            delConstraint(tag=n)
        Entities['Nodes'][tag]['freedof'][n] = -1
    return True

def addConstraint(tag=np.nan, name='Unknown', attributes={}):
    """
    Specifies a constraint applied to a degree of freedom\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of the constraint, i.e., tag < -1
    name : str
        The constraint type that is being applied.
    attributes : dict
        Dictionary containing the master node(s) information. The fields are
        'stag'  : (int) the slave node(s) tags to be constrained
        'sdof'  : (int) the slave degree-of-freedom to be constrained
        'mtag'  : (int or list) the master node(s) tags to be employed
        'mdof'  : (int or list) the master degree-of-freedom to be combined
        'factor': (int or list) the combinational factors to be combined

    Returns
    -------
    bool
        Whether the constraint/s was/were successful (True) of failed (False)
    """
    #Transform single value into a list
    if isinstance(attributes['mtag'], int):
        attributes['mtag'] = [attributes['mtag']]
        attributes['mdof'] = [attributes['mdof']]
    if 'factor' not in attributes and name.upper() == 'EQUAL':
        attributes['factor'] = [1.00]

    #Transform to zero base index
    attributes['sdof'] -= 1
    attributes['mdof'] = [k-1 for k in attributes['mdof']]

    #Check whether the Constraint exists and is well-defined
    if tag not in Entities['Constraints']:
        if tag < -1:
            stag = attributes['stag']
            sdof = attributes['sdof']
            if stag in Entities['Nodes']:
                #The constraint is not applied if slave node was restrained
                if Entities['Nodes'][stag]['freedof'][sdof] == -1:
                    info = debugInfo(2)
                    print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Constraint[%d] applied to Node[%d] on dof=%d was restrained.' %(info.filename,info.lineno,tag,stag,sdof+1))
                    return False
                #Rise an ALERT if master node(s) is/are restrained
                for mnode, mdof in zip(attributes['mtag'], attributes['mdof']):
                    if Entities['Nodes'][mnode]['freedof'][mdof] < 0:
                        info = debugInfo(2)
                        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Constraint[%d] has master Node[%d] on dof=%d that is restrained/constrained.' %(info.filename,info.lineno,tag,mnode,mdof+1))
                        return False
                Entities['Constraints'][tag] = {'name': name.upper(), 'stag': stag, 'sdof': sdof, 'mtag': attributes['mtag'], 'mdof': attributes['mdof'], 'factor': attributes['factor']}
                Entities['Nodes'][stag]['freedof'][sdof] = tag
                return True
            else:
                info = debugInfo(2)
                print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Constraints[%d] could not be assign to Node[%d].' %(info.filename,info.lineno,tag,stag))
                return False
        else:
            info = debugInfo(2) 
            print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d the variable tag must be less than -1.' %(info.filename,info.lineno))
            return False
    else:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Constraints[%d] has been already defined.' %(info.filename,info.lineno,tag))
        return False

def addRigidLink(tag=np.nan, attributes={}):
    """
    Specifies a rigid link applied to a set of nodes\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2021

    Parameters
    ----------
    tag : int
        The identifier of the rigid body, i.e., tag > 1 (defferent from Nodes)
    attributes : dict
        Dictionary containing the diaphragm action and Nodes
        'tag'   : (int) the node tag that the rigid body node will take (must be different from 'Nodes')
        'ndof'  : (int) the number of degree of freedom of the rigid body i.e., 3 (2D) and 6 (3D)
        'list'  : (int or list) list of Nodes that belong to the rigid body
        'center' : (str) center of rotation of the rigid body

    Returns
    -------
    bool
        Whether the rigid link/s was/were successful (True) of failed (False)
    """
    #TODO: 
    return False

def addDiaphragm(tag=np.nan, attributes={}):
    """
    Specifies a rigid diaphragm applied to a set of nodes\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of the diaphragm, i.e., tag > 1
    attributes : dict
        Dictionary containing the diaphragm action and Nodes
        'tag'   : (int) the node tag that the diaphragm node will take (must be different from 'Nodes')
        'list'  : (int or list) list of Nodes that belong to the diaphragm
        'axis'  : (str) axis that is perpendicular to the diaphragm

    Returns
    -------
    bool
        Whether the diaphragm/s was/were successful (True) of failed (False)
    """
    if 'axis' in attributes:
        attributes['axis'] = attributes['axis'].upper()
        if attributes['axis'] not in ['X','Y','Z']:
            info = debugInfo(2)
            print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d in Diaphragm[%d] \'axis\'=%s is invalid.' %(info.filename,info.lineno,tag,attributes['axis']))
            return False
    else:
        info = debugInfo(2)
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d in Diaphragm[%d] the attributes[\'axis\'] must be defined.' %(info.filename,info.lineno,tag))
        return False

    #Check whether the Diaphragm exists
    if tag not in Entities['Diaphragms']:
        Entities['Diaphragms'][tag] = attributes
        return True
    else:
        info = debugInfo(2)
        if attributes['axis'] == Entities['Diaphragms'][tag]['axis']:
            indeces = list(set( Entities['Diaphragms'][tag]['list'] + attributes['list'] ))
            Entities['Diaphragms'][tag]['list'] = indeces
            print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Diaphragm[%d] has been already defined, but list was appended.' %(info.filename,info.lineno,tag))
        else: 
            print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Diaphragm[%d] has been already defined.' %(info.filename,info.lineno,tag))
            return False

def addRigidBody(tag=np.nan, attributes={}):
    """
    Specifies a rigid diaphragm applied to a set of nodes\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of the rigid body, i.e., tag > 1 (defferent from Nodes)
    attributes : dict
        Dictionary containing the diaphragm action and Nodes
        'tag'   : (int) the node tag that the rigid body node will take (must be different from 'Nodes')
        'ndof'  : (int) the number of degree of freedom of the rigid body i.e., 3 (2D) and 6 (3D)
        'list'  : (int or list) list of Nodes that belong to the rigid body
        'center' : (str) center of rotation of the rigid body

    Returns
    -------
    bool
        Whether the rigid body/s was/were successful (True) of failed (False)
    """
    if 'center' not in attributes:
        info = debugInfo(2)
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d in RigidBodies[%d] the attributes[\'center\'] must be defined.' %(info.filename,info.lineno,tag))
        return False
    if 'ndof' not in attributes:
        if Options['dimension'] == 2:
            attributes['ndof'] = 3
        elif Options['dimension'] == 3:
            attributes['ndof'] = 6

    #Check whether the Rigid Body exists
    if tag not in Entities['RigidBodies']:
        Entities['RigidBodies'][tag] = attributes
        return True
    else:
        info = debugInfo(2)
        indeces = list(set( Entities['RigidBodies'][tag]['list'] + attributes['list'] ))
        Entities['RigidBodies'][tag]['list'] = indeces
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d RigidBodies[%d] has been already defined, but list was appended.' %(info.filename,info.lineno,tag))

def addSupportMotion(tag=np.nan, attributes={}):
    """
    Specifies a support motion to be applied to a degree-of-freedom\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of the restrained node, i.e., tag > 1 to apply the support motion 
    attributes : dict
        Dictionary containing the support motion action on a Nodes in a direction
        'type'  : (str) Constant or TimeSerie
        'value' : (float) displacement value if type=static 
        'file'  : (str) displacement time serie file if type=dynamic
        'dof'   : (int) degree of freedom to apply the support motion

    Returns
    -------
    bool
        Whether the diaphragm/s was/were successful (True) of failed (False)
    """
    if 'type' in attributes:
        attributes['type'] = attributes['type'].upper()
    if isinstance(attributes['dof'], list):
        attributes['dof'] = [k-1 for k in attributes['dof']]
    if isinstance(attributes['dof'], int):
        attributes['dof'] = [ attributes['dof'] - 1 ]
    if 'value' in attributes:
        if isinstance(attributes['value'], float):
            attributes['value'] = [ attributes['value'] ]
    if 'file' in attributes:
        if isinstance(attributes['file'], str):
            attributes['file'] = [ attributes['file'] ]
        elif isinstance(attributes['file'], list):
            for k, fname in enumerate(attributes['file']):
                attributes['file'][k] = fname
    
    #Check whether the Support exists
    if tag not in Entities['Supports']:
        Entities['Supports'][tag] = attributes
        return True
    else:
        for m,dof in enumerate(attributes['dof']):
            if dof not in Entities['Supports'][tag]['dof']:
                Entities['Supports'][tag]['dof'].append(dof)
                if 'value' in attributes:
                    Entities['Supports'][tag]['value'].append(attributes['value'][m])
                elif 'file' in attributes:
                    Entities['Supports'][tag]['file'].append(attributes['file'][m])
            else:
                info = debugInfo(2) 
                print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Supports[%d] has been already defined.' %(info.filename,info.lineno,tag))
        return True

def addMaterial(tag=np.nan, name='Unknown', attributes={}):
    """
    Appends a new Material to the model\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of this material, i.e., tag > -1
    name : str
        Seismo-VLAB material class name
    attributes : dict
        Specific properties for the created material

    Returns
    -------
    bool
        Whether the addition was successful (True) of failed (False)
    """
    if name == 'Unknown':
        info = debugInfo(2)
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Section[%d] \'name\' must be specified.' %(info.filename,info.lineno,tag))
        
    #Check whether the Material exists
    if tag not in Entities['Materials']:
        Entities['Materials'][tag] = {'name': name.upper(), 'attributes': attributes}
        return True
    else:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Material[%d] has been already defined.' %(info.filename,info.lineno,tag))
        return False

def addSection(tag=np.nan, name='Unknown', model='Unknown', attributes={}):
    """
    Appends a new Section to the model\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of this section, i.e., tag > -1
    name : str
        Seismo-VLAB section class name
    model : str
        model of the section, i.e., model=plain,fiber
    attributes : dict
        Specific properties for the created section

    Returns
    -------
    bool
        Whether the addition was successful (True) of failed (False)
    """
    #Check if insetion point and section rotation are provided 
    if 'theta' not in attributes:
        attributes['theta'] = 0.0
    if 'ip' not in attributes:
        attributes['ip'] = 10
    if name == 'Unknown':
        info = debugInfo(2)
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Section[%d] \'name\' must be specified.' %(info.filename,info.lineno,tag))
    if model == 'Unknown':
        info = debugInfo(2)
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Section[%d] \'model\' must be specified.' %(info.filename,info.lineno,tag))

    #Check whether the Section exists
    if tag not in Entities['Sections']:
        Entities['Sections'][tag] = {'name': name.upper(), 'model': model.upper(), 'attributes': attributes}
        return True
    else:
        info = debugInfo(2)
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Section[%d] has been already defined.' %(info.filename,info.lineno,tag))
        return False

def addElement(tag=np.nan, name='Unknown', conn=[], attributes={}):
    """
    Appends a new Element to the model\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of this element, i.e., tag > -1
    conn : list
        Connectivity array of this element
    name : str
        Seismo-VLAB element class name
    attributes : dict
        Specific properties for the created element

    Returns
    -------
    bool
        Whether the addition was successful (True) of failed (False)
    """
    if 'model' in attributes:
        attributes['model'] = attributes['model'].upper()
    if 'dir' in attributes:
        attributes['dir'] -= 1
    if 'rule' in attributes:
        attributes['rule'] = attributes['rule'].upper()
        
    #Check whether the Element exists
    if tag not in Entities['Elements']:
        Entities['Elements'][tag] = {'name': name.upper(), 'conn': conn, 'attributes': attributes}
        return True
    else:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Element[%s] has been already defined.' %(info.filename,info.lineno,tag))
        return False

def addSurface(tag=np.nan, etag=-1, conn=[]):
    """
    Appends a surface for an specific element\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of this surface, i.e., tag > 0
    etag : int
        The identifier of this element, i.e., tag > 0
    conn : list
        Connectivity array of this element surface

    Returns
    -------
    bool
        Whether the addition was successful (True) of failed (False)
    """
    #Check whether the Surface exists
    if tag not in Entities['Surfaces']:
        Entities['Surfaces'][tag] = {'etag': etag, 'conn': conn}
        return True
    else:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Surface[%d] has been already defined.' %(info.filename,info.lineno,tag))
        return False

def addDamping(tag=np.nan, name='Unknown', attributes={}):
    """
    Appends a damping formulation for a group of elements\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of this damping model, i.e., tag > -1
    name : str
        The damping class name
    attributes : dict
        Specific properties for the created damping, for example
        'list' : The element identifier that share this damping
        'ak'   : The stiffness proportinal factor used in Rayleigh Damping
        'am'   : The mass proportinal factor used in Rayleigh Damping

    Returns
    -------
    bool
        Whether the addition was successful (True) of failed (False)
    """
    #Check list format
    if 'list' in attributes:
        if isinstance(attributes['list'], int):
            attributes['list'] = [ attributes['list'] ]
    
    #Check whether the Damping exists
    if tag not in Entities['Dampings']:
        Entities['Dampings'][tag] = {'name': name.upper(), 'attributes': attributes}
        return True
    else:
        info = debugInfo(2)
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Damping[%d] has been already defined.' %(info.filename,info.lineno,tag))
        return False

def addFunction(tag=np.nan, name='Unknown', attributes={}):
    """
    Appends a function to describe a load behavior\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of this function, i.e., tag > -1
    name : str
        The Seismo-VLAB function class name: name=TimeSeries,Constant,PlaneWave,GeneralWave
    attributes : dict
        Specific properties for the created function
        'dir'  : (list) The direction of the acting function
        'mag'  : (float) The magnitude of the load
        'file' : (str) The file where the load is going to be loaded. The file location is with respect to the main python script being executed

    Returns
    -------
    bool
        Whether the addition was successful (True) of failed (False)
    """
    #Generates the global path for file
    if 'file' in attributes:
        attributes['file'] = attributes['file']
    
    #Check whether the Function exists
    if tag not in Entities['Functions']:
        Entities['Functions'][tag] = {'name': name.upper(), 'attributes': attributes}
        return True
    else:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Function[%d] has been already defined.' %(info.filename,info.lineno,tag))
        return False

def addLoad(tag=np.nan, name='Unknown', attributes={}):
    """
    Appends a damping formulation for a group of elements\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of this load, i.e., tag > -1
    name : str
        The Seismo-VLAB load class name: name=PointLoad, ElementLoad, SupportMotion
    attributes : dict
        Specific properties for the created load
        'fun' (int) The function identifier
        'list' (list) The identifier where this load acts upon
        'type' (str) The type of element load:
            type=Constant,TimeSeries (PointLoad)
            type=Surcace,Body,PlaneWave,GeneralWave (ElementLoad)
 
    Returns
    -------
    bool
        Whether the addition was successful (True) of failed (False)
    """
    #Check the load type:
    if 'list' in attributes:
        if isinstance(attributes['list'], int):
            attributes['list'] = [ attributes['list'] ]
    if 'type' in attributes:
        attributes['type'] = attributes['type'].upper()
    else:
        if name.upper() != 'SUPPORTMOTION':
            info = debugInfo(2) 
            print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d for Load[%d] the attributes[\'type\'] must be defined.' %(info.filename,info.lineno,tag))
            return False

    #Check whether the Load exists
    if tag not in Entities['Loads']:
        #Check whether the associated function exists
        if 'fun' in attributes:
            ftag = attributes['fun']
            if ftag not in Entities['Functions']:
                info = debugInfo(2) 
                print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d for Load[%d] the Function[%d] has not been defined.' %(info.filename,info.lineno,tag,ftag))
                return False

        Entities['Loads'][tag] = {'name': name.upper(), 'attributes': attributes}
        return True
    else:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Load[%d] has been already defined.' %(info.filename,info.lineno,tag))
        return False

def addCombinationCase(tag=np.nan, name='Unknown', attributes={}):
    """
    Appends a Combination case for an analysis\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of this combination, i.e., tag > -1
    name : str
        The combination case name
    attributes : dict
        Specific properties for the created combination
        'folder' : (str) Name of the folder to store solution
        'load'   : (list) Load identifiers to be combined
        'factor' : (list) Combinational factors

    Returns
    -------
    bool
        Whether the addition was successful (True) of failed (False)
    """
    if 'folder' not in attributes:
        attributes['folder'] = name
    if isinstance( attributes['load'], int):
        attributes['load'] = [ attributes['load'] ]
        attributes['factor'] = [ attributes['factor'] ]

    #Check whether the Combination exists
    if tag not in Entities['Combinations']:
        Entities['Combinations'][tag] = {'name': name, 'attributes': attributes}
        return True
    else:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Combination[%d] has been already defined.' %(info.filename,info.lineno,tag))
        return False

def addRecorder(tag=np.nan, attributes={}):
    """
    Appends a Recorder to output some quantity\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of this recorder, i.e., tag > -1
    attributes : dict
        Specific properties for the created recorder
        'name'  : (str) The recorder type: 'Node', 'Element', 'Material', 'Section', 'ParaView'
        'file'  : (str) The output's file name
        'ndps'  : (int) Number of digit of precision for solution
        'nsamp' : (int) Number of point to skip when sampling the response
        'resp'  : (str) The response associated to 'name' to be redorded
            Node     = disp, vel, accel, reaction
            Section  = strain, stress
            Element  = strain, stress, strainrate, internalforce
        'pos'   : (list) coordinate where to measure stress/strain in section
        'list'  : (list) The identifier list to be recorded

    Returns
    -------
    bool
        Whether the addition was successful (True) of failed (False)
    """
    #Check the recorder type is correct
    if 'name' in attributes:
        attributes['name'] = attributes['name'].upper()
    if 'nsamp' not in attributes:
        attributes['nsamp'] = 1
    if attributes['name'] not in ['NODE', 'SECTION', 'ELEMENT', 'PARAVIEW']:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d addRecorder(attributes=?) attributes[%s] is not recognized.' %(info.filename,info.lineno,attributes['name']))
        return False
    if 'list' in attributes:
        if isinstance(attributes['list'], int):
            attributes['list'] = [ attributes['list'] ]

    #Check whether the Recorder exists
    if tag not in Entities['Recorders']:
        Entities['Recorders'][tag] = attributes
        return True
    else:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Recorder[%d] has been already defined.' %(info.filename,info.lineno,tag))
        return False

def addSolver(tag=np.nan, attributes={}):
    """
    Appends a Solver to be performed\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of this solver, i.e., tag > -1
    attributes : dict
        Specific properties for the created solver
        'name'   : (str) The solver's name
        'update' : (str) ON for nonlinear, OFF for linear
        'option' : (str) Solver's matrix structure, or iterative algorithm name
        'tol'    : (float) Iteration tolerance for iterative solver

    Returns
    -------
    bool
        Whether the addition was successful (True) of failed (False)
    """
    if 'name' in attributes:
        attributes['name'] = attributes['name'].upper()
        if attributes['name'] == 'PETSC':
            Options['allocation'] = 'YES'
    else:
        info = debugInfo(2)
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d attributes[%s] must be specified.' %(info.filename,info.lineno,'name')) 
        return False
    if 'update' in attributes:
        UPDATE = attributes['update'].upper()
        attributes['update'] = 1  if UPDATE == 'OFF'  else 0
    if 'option' in attributes:
        OPTION = attributes['option'].upper()
        attributes['option'] = SolverOption[OPTION]

    #Check whether the Solver exists
    if tag not in Entities['Solvers']:
        Entities['Solvers'][tag] = attributes
        return True
    else:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Solver[%d] has been already defined.' %(info.filename,info.lineno,tag))
        return False

def addAlgorithm(tag=np.nan, attributes={}):
    """
    Appends an Algorithm to be performed\n
    @visit  https://github.com/SeismoVLAB/SVL\n\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of this algorithm, i.e., tag > -1
    attributes : dict
        Specific properties for the created algorithm
        'name'     : (str) The algorithm's name
        'nstep'    : (int) The algorithm's number of increment steps
        'cnvgtol'  : (float) Tolerance to accept the solution has converged
        'cnvgtest' : (str) The converge test to apply, these can be
            'UnbalanceForce'
            'IncrementalDisplacement'
            'RelativeUnbalanceForce'
            'RelativeIncrementalDisplacement'
            
    Returns
    -------
    bool
        Whether the addition was successful (True) of failed (False)
    """
    #Sets default values
    if 'nstep' not in attributes:
        attributes['nstep'] = 1
    if 'cnvgtol' not in attributes:
        attributes['cnvgtol'] = 1E-3
    if 'cnvgtest' not in attributes:
        attributes['cnvgtest'] = ConvergeTest['UNBALANCEFORCE']
    else:
        CNVGTEST = attributes['cnvgtest'].upper()
        attributes['cnvgtest'] = ConvergeTest[CNVGTEST] 

    #Upper case Algorithm name
    if 'name' in attributes:
        attributes['name'] = attributes['name'].upper()

    #Check whether the Combination exists
    if tag not in Entities['Algorithms']:
        Entities['Algorithms'][tag] = attributes
        return True
    else:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Algorithm[%d] has been already defined.' %(info.filename,info.lineno,tag))
        return False

def addIntegrator(tag=np.nan, attributes={}):
    """
    Appends an Integrator to be performed\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of this integrator, i.e., tag > -1
    attributes : dict
        Specific properties for the created integrator
        'name' : (str) The integrator's name
        'dt'   : (float) The integrator's time step
        'ktol' : (float) Tolerance to set a stiffness matrix value kij = 0.0
        'mtol' : (float) Tolerance to set a mass matrix value  mij = 0.0
        'ftol' : (float) Tolerance to set a force vector value fj = 0.0

    Returns
    -------
    bool
        Whether the addition was successful (True) of failed (False)
    """
    #Sets default values
    if 'ktol' not in attributes:
        attributes['ktol'] = 1E-12
    if 'mtol' not in attributes:
        attributes['mtol'] = 1E-12
    if 'ftol' not in attributes:
        attributes['ftol'] = 1E-12
    if 'dt' not in attributes:
        attributes['dt'] = 0.0
    
    #Upper case Integrator name
    if 'name' in attributes:
        attributes['name'] = attributes['name'].upper()

    #Check whether the Combination exists
    if tag not in Entities['Integrators']:
        Entities['Integrators'][tag] = attributes
        return True
    else:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Integrator[%d] has been already defined.' %(info.filename,info.lineno,tag))
        return False

def addAnalysis(tag=np.nan, attributes={}):
    """
    Appends an Analysis to be performed\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of this analysis, i.e., tag > -1
    attributes : dict
        Specific properties for the created analysis
        'name' : (str) The analysis's name
        'nt'   : (int) The number of time step to be evolved

    Returns
    -------
    bool
        Whether the addition was successful (True) of failed (False)
    """
    #Upper case Analysis name
    if 'name' in attributes:
        attributes['name'] = attributes['name'].upper()

    #Check whether the Combination exists
    if tag not in Entities['Analyses']:
        Entities['Analyses'][tag] = attributes
        return True
    else:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Analysis[%d] has been already defined.' %(info.filename,info.lineno,tag))
        return False

def addSimulation(tag=np.nan, combo=None, attributes={}):
    """
    Appends a Simulation to be performed\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    tag : int
        The identifier of this simulation, i.e., tag > -1
    combo : int
        The combination identifier for this simulation, i.e., ctag > -1
    attributes : dict
        Specific properties for the created simulation
        'analysis'   : (int) The analysis identifier to be used
        'algorithm'  : (int) The algorithm identifier to be used
        'integrator' : (int) The integrator identifier to be used
        'solver'     : (int) The solver identifier to be used

    Returns
    -------
    bool
        Whether the addition was successful (True) of failed (False)
    """
    #Check whether the Combination exists
    if tag not in Entities['Simulations']:
        Entities['Simulations'][tag] = {'combo': combo, 'attributes': attributes}
        return True
    else:
        info = debugInfo(2) 
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d Simulation[%d] has been already defined.' %(info.filename,info.lineno,tag))
        return False
