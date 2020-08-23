#!/usr/bin/python3

import numpy as np
import scipy.sparse as sps
import matplotlib.pylab as plt
from scipy.sparse import find

def PetscAllocation(User, Point, Element, Constraint):
    """
    """
    N = User['NPART']
    M  = User['NFREEDOF']

    nr = int(M/N);
    n_nz = np.zeros(N, dtype='int')
    o_nz = np.zeros(N, dtype='int')

    #Assembles the matrix.
    A = FormSparseMatrix(User, Point, Element, Constraint)

    #Re-define matrix.
    (I,J,V) = find(A)
    V = 0.0*V + 1.0;
    A = sps.coo_matrix((V,(I, J)), shape=(M,M))
    A = A.tocsr()

    #Loop over the partitions.
    for k in range(0,N):
        #Matrix partition.
        inf = k*nr
        sup = (k + 1)*nr

        #Sub-Matrix stripe. 
        S = A[inf:sup,:]

        #Sub-Matrix block-diagonal.
        D = S[:,inf:sup]

        n_nz[k] = np.max(D.sum(axis=1))
        o_nz[k] = np.max(S.sum(axis=1) - D.sum(axis=1))

    #Sets the closest multipple of five.
    n_nz = np.floor( np.divide(n_nz,5) + 1) * 5
    o_nz = np.floor( np.divide(o_nz,5) + 1) * 5

    User['D_NZ'] = n_nz
    User['O_NZ'] = o_nz

def FormSparseMatrix(User, Point, Element, Constraint):
    """
    This function compute/emulates the Sparse Matrix Pattern to be employed 
    in the user's defined ordering scheme.

    Parameters
    ----------
    User : dict 
           The dictionary containing all the user's relevant information.
    Point : dict
           The dictionary containing all point information.
    Element : dict
           The dictionary containing all element information.
    Constraint : dict
           The dictionary containing all constraint information.

    Output
    -------
    Modifies User['NFREEDOF'], User['NTOTALDOF'], Point['FREE'] and 
    Point['TOTAL'] fields
    """
    #Emulates the Element Assembly Pattern 
    N = User['NTOTALDOF'  ]
    M = User['NCONSISTENT']

    I = np.zeros(M, dtype='int')
    J = np.zeros(M, dtype='int')
    V = np.zeros(M, dtype='float')

    m = 0
    for eTag in Element:
        #Element's node connectivity array
        connection = Element[eTag]['CONNECTIVITY']

        #Element's degree-of-freedom list
        total = []
        for nTag in connection:
            total = np.concatenate((total, Point[nTag]['TOTAL']))

        for idof in total:
            for jdof in total:
                I[m] = idof
                J[m] = jdof
                V[m] = 1.00
                m += 1

    A = sps.coo_matrix((V,(I, J)), shape=(N,N))

    #Emulates the Transformation matrix 
    M = User['NFREEDOF' ]
    N = User['NTOTALDOF']
    S = User['NFREEDOF' ] + User['NCONSTRAINTDOF']

    I = np.zeros(S, dtype='int')
    J = np.zeros(S, dtype='int')
    V = np.zeros(S, dtype='float')

    m = 0
    for nTag in Point:
        nDOF  = Point[nTag]['NDOF' ]
        free  = Point[nTag]['FREE' ]
        total = Point[nTag]['TOTAL']

        for k in range(nDOF):
            #The Free Degree-Of-Freedom is unconstrained
            if free[k] > -1:
                I[m] = total[k]
                J[m] = free[k]
                V[m] = 1.00
                m += 1

            #The Free Degree-Of-Freedom is constrained
            if free[k] < -1:
                SlaveNode = Constraint[free[k]]['SLAVENODE']
                SlaveDOF  = Constraint[free[k]]['SLAVEDOF' ]
                Slave     = Point[SlaveNode]['TOTAL'][SlaveDOF]

                MasterNode = Constraint[free[k]]['MASTERNODE']
                MasterDOF  = Constraint[free[k]]['MASTERDOF' ]
                CombCoeffs = Constraint[free[k]]['FACTORS']

                for i in range(len(MasterNode)):
                    Master = Point[MasterNode[i]]['FREE'][MasterDOF[i]]
                    I[m] = Slave
                    J[m] = Master
                    V[m] = 1.00 
                    m += 1

    T = sps.coo_matrix((V,(I, J)), shape=(N,M))

    #The Final Element Assembly pattern 
    A = T.transpose()*A*T
    A = A.tocsr()

    #plt.spy(A, markersize=2)
    #plt.savefig("OrderMatrix.png")
    #plt.grid(True)
    #plt.close()

    return A

def PlainScheme(User, Point, Element, Constraint):
    """
    This function assign a Plain Scheme to the degree of freedom of each Point, 
    i.e., the node numbering is assigned consecutively from the lowest to the 
    highest Point identifier. 

    Parameters
    ----------
    User : dict 
           The dictionary containing all the user's relevant information.
    Point : dict
           The dictionary containing all point information.
    Element : dict
           The dictionary containing all element information.

    Output
    -------
    Modifies User['NFREEDOF'], User['NTOTALDOF'], Point['FREE'] and 
    Point['TOTAL'] fields
    """
    #Counter Variables
    count0 = 0
    count1 = 0
    nConstraintDofs = 0

    #Assign consecutively degree of freedom numbering
    for node in Point:
        nDofs = Point[node]['NDOF' ]
        free  = Point[node]['FREE' ]
        total = Point[node]['TOTAL']

        #Total degree-of-freedom numbering.
        for k in range(nDofs):
            total[k] = count0
            count0 += 1

        #Free degree-of-freedom numbering.
        for k in range(nDofs):
            if free[k] > -1:
                free[k] = count1
                count1 += 1
            elif free[k] < -1:
                MasterNodes = Constraint[free[k]]['MASTERNODE'] 
                nConstraintDofs += len(MasterNodes)

        #Assign ne Free/Total degree-of-freedom numbering.
        Point[node]['FREE' ] = free
        Point[node]['TOTAL'] = total

    #Write number of Free/Total/Constrained degree-of-freedom.
    User['NFREEDOF'      ] = count1
    User['NTOTALDOF'     ] = count0
    User['NCONSTRAINTDOF'] = nConstraintDofs

def MinimumDegreeScheme(User, Point, Element, Constraint):
    """
    This function performs the Minimum-Degree ordering Scheme to the free
    degree of freedom. First compute the Sparse Matrix Pattern, and then
    computes the permutation vector. Such vercor is finally employed to
    re-label the free-degree-of-freedom

    Parameters
    ----------
    User : dict 
           The dictionary containing all the user's relevant information.
    Point : dict
           The dictionary containing all point information.
    Element : dict
           The dictionary containing all element information.
    Constraint : dict
           The dictionary containing all constraint information.

    Output
    -------
    Modifies Point['FREE'] field
    """
    print('\x1B[33m ALERT \x1B[0m: The degree-of-freedom ordering MAXIMUM-DEGREE is not implemented yet.')
    print(' The degree-of-freedom numbering will assume a Plain Scheme.')

def CutHillMcKeeScheme(User, Point, Element, Constraint):
    """
    This function performs the CutHill-McKee ordering Scheme to the free
    degree of freedom. First compute the Sparse Matrix Pattern, and then
    computes the permutation vector. Such vercor is finally employed to
    re-label the free-degree-of-freedom

    Parameters
    ----------
    User : dict 
           The dictionary containing all the user's relevant information.
    Point : dict
           The dictionary containing all point information.
    Element : dict
           The dictionary containing all element information.
    Constraint : dict
           The dictionary containing all constraint information.

    Output
    -------
    Modifies Point['FREE'] field
    """
    #Gets the Sparse Matrix to Perform Permutation
    A = FormSparseMatrix(User, Point, Element, Constraint)

    #Computes the Permutation Vector for the Matrix
    perm = sps.csgraph.reverse_cuthill_mckee(A, True)

    #The New Free-Degree-Of-Freedom Numbering
    N       = len(perm)
    I       = np.zeros(N, dtype='int')
    I[perm] = np.arange(0, N, 1)

    #Transform the degree of freedom numbering form Plain to CutHill-McKee
    for nTag in Point:
        nDofs = Point[nTag]['NDOF']
        free  = Point[nTag]['FREE']

        #New-Free degree-of-freedom numbering.
        for k in range(nDofs):
            if free[k] > -1:
                free[k] = I[free[k]]

        #Assign the Free degree-of-freedom numbering.
        Point[nTag]['FREE'] = free

def FindDefectiveNodes(User, Point, Constraint, Element):
    """
    This function loops over the Element and identify Point that does not
    belong to them. Then, it fix such Point that are defective. Constraints 
    of type Diaphragm or General are not considered defective. 

    Parameters
    ----------
    User : dict 
           The dictionary containing all the user's relevant information.
    Point : dict
           The dictionary containing all point information.
    Constraint : dict
           The dictionary containing all constraints information.
    Element : dict
           The dictionary containing all element information.

    Output
    -------
    Modifies Point['FREE'] field
    """
    #Memory Storage Variables
    nLumpedStorage     = 0
    nConsistentStorage = 0

    #Identifies Nodes which don't belong to Element.
    for eTag in Element:
        nDOFelem   = 0
        connection = Element[eTag]['CONNECTIVITY']
        for conn in connection:
            Point[conn]['DEFECTIVE'] = False
            nDOFelem += Point[conn]['NDOF']

        nLumpedStorage     += nDOFelem
        nConsistentStorage += nDOFelem*nDOFelem

    #Diaphragm, Equal, and General Constraints Nodes are considered non-defective.
    for cTag in Constraint:
        mNodes = Constraint[cTag]['MASTERNODE']
        for mTag in mNodes:
            Point[mTag]['DEFECTIVE'] = False

    #Fix all Detected Defective Point
    for nTag in Point:
        if Point[nTag]['DEFECTIVE']:
            nDOFs = Point[nTag]['NDOF']
            Point[nTag]['FREE'] = np.full(nDOFs,  -1, dtype='int')

            print('\x1B[33m ALERT \x1B[0m: The Point[' + str(nTag) + '] does not belong to an element. The Point will be fixed.')

    #Saves in User's file the Matrix Memory Storage 
    User['NLUMPED'    ] = nLumpedStorage
    User['NCONSISTENT'] = nConsistentStorage

def SetDegreeOfFreedom(User, Point, Element, Constraint):
    """
    This function assigns the degree of freedom numbering for each Point 
    according to the User's numbering pattern.

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
    Constraint : dict
           The dictionary containing all constraints information.

    Output
    -------
    Modifies Point['FREE'] and Point['TOTAL'] fields
    """
    #Detects Point that does not belong to Element
    FindDefectiveNodes(User, Point, Constraint, Element)

    #Assign Plain Numbering Scheme
    PlainScheme(User, Point, Element, Constraint)

    #Impose User's Define Numbering Scheme
    if User['ORDERING'] == 'PLAIN':
       pass
    elif User['ORDERING'] == 'CUTHILL-MCKEE':
        CutHillMcKeeScheme(User, Point, Element, Constraint)
    elif User['ORDERING'] == 'MAXIMUM-DEGREE':
        MinimumDegreeScheme(User, Point, Element, Constraint)
    else:
        print('\x1B[33m ALERT \x1B[0m: The degree-of-freedom ordering ' + User['ORDERING'] + ' is not identified.')
        print(' The degree-of-freedom numbering will assume a Plain Scheme.')

    #Compute memory storage for Petsc matrix.
    if User['PETSC'] == 'YES':
        PetscAllocation(User, Point, Element, Constraint)
