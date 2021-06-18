#!/usr/bin/env python3
# -*- coding: Utf-8 -*-

Entities = {
    #SVL dictionary with associated geometry and analysis entities
    # https://github.com/SeismoVLAB/SVL
    # @author Danilo S. Kusanovic 2020
    'Nodes'       : {}, 
    'Masses'      : {},
    'RigidLinks'  : {},
    'Diaphragms'  : {}, 
    'RigidBodies' : {},
    'Constraints' : {},
    'Supports'    : {}, 
    'Surfaces'    : {}, 
    'Materials'   : {}, 
    'Sections'    : {}, 
    'Elements'    : {},
    'Dampings'    : {}, 
    'Functions'   : {}, 
    'Loads'       : {}, 
    'Combinations': {},
    'Recorders'   : {},
    'Algorithms'  : {},
    'Solvers'     : {},
    'Integrators' : {},
    'Analyses'    : {},
    'Simulations' : {}
}

Options = {
    #SVL dictionary with user's defined variables
    # https://github.com/SeismoVLAB/SVL
    # @author Danilo S. Kusanovic 2020
    'run'         : '',
    'path'        : '',
    'file'        : 'SeismoVLAB',
    'description' : '\n',
    'execfile'    : '',
    'execpath'    : '',
    'preanalysis' : '',
    'runanalysis' : '',
    'allocation'  : 'NO',
    'numbering'   : 'Plain',
    'metispath'   : '',
    'solution'    : 'Sequential',
    'massform'    : 'Consistent',
    'wasChecked'  : False,
    'nparts'      :  1,
    'dimension'   :  0,
    'nfree'       :  0,
    'ntotal'      :  0,
    'nconstraint' :  0,
    'nlumped'     :  0,
    'nconsistent' :  0,
    'nparaview'   :  0,
    'nfeatures'   :  0,
    'd_nz'        : [],
    'o_nz'        : [],
    'partition'   : []
}

ConvergeTest = {
    'UNBALANCEFORCE'                 : 1,
    'INCREMENTALDISPLACEMENT'        : 2,
    'RELATIVEUNBALANCEFORCE'         : 3,
    'RELATIVEINCREMENTALDISPLACEMENT': 4
}

SolverOption = {
    'SPD'     : 0,
    'SYM'     : 1,
    'USYM'    : 2,
    'KSPCG'   : 0,
    'KSPBCGS' : 1,
    'KSPCGS'  : 2,
    'KSPBICG' : 3
}

SVLclasses = {
    #VTK dictionary with the SVL associated material, section, and element classes
    # https://github.com/SeismoVLAB/SVL
    # @author Danilo S. Kusanovic 2021
    #  'Materials' : Stores the material classes defined in SVL
    #  'Sections'  : Stores the section classes defined in SVL
    #  'Elements'  : Stores the elements classes defined in SVL
    'Materials': {
        'ELASTIC1DLINEAR': {
            'dim'  : [1,2,3]
        },
        'HERTZIAN1DLINEAR': {
            'dim': [1,2,3]
        },
        'VISCOUS1DLINEAR': {
            'dim': [1,2,3]
        },
        'PLASTIC1DJ2': {
            'dim': [1,2,3]
        },
        'ELASTIC2DPLANESTRAIN': {
            'dim': [2,3]
        },
        'ELASTIC2DPLANESTRESS': {
            'dim': [2,3]
        },
        'PLASTICPLANESTRAINJ2': {
            'dim': [2]
        },
        'PLASTICPLANESTRAINBA': {
            'dim': [2]
        },
        'ELASTIC3DLINEAR': {
            'dim': [3]
        },
        'PLASTIC3DJ2': {
            'color': [237,237,237], 
            'dim': [3]
        },
        'PLASTIC3DBA': {
            'dim': [3]
        },
        'ELASTIC1DFIBER': {
            'dim': [1,2,3]
        },
        'STEEL1DFIBER': {
            'dim': [1,2,3]
        },
        'CONCRETE1DFIBER': {
            'dim': [1,2,3]
        },
        'ELASTIC1DGAP': {
            'dim': [1,2,3]
        },
        'PLASTIC1DGAP': {
            'dim': [1,2,3]
        }
    },
    'Sections': {
        'LIN2DANGLE': {
            'dim': [2]
        },
        'LIN2DCHANNEL': {
            'dim': [2]
        },
        'LIN2DCIRCULAR': {
            'dim': [2]
        },
        'LIN2DCIRCULARTUBE': {
            'dim': [2]
        },
        'LIN2DRECTANGULAR': {
            'dim': [2]
        },
        'LIN2DRECTANGULARTUBE' : {
            'dim': [2]
        },
        'LIN2DTEE': {
            'dim': [2]
        },
        'LIN2DWIDEFLANGE': {
            'dim': [2]
        },
        'LIN3DANGLE': {
            'dim': [3]
        },
        'LIN3DCHANNEL': {
            'dim': [3]
        },
        'LIN3DCIRCULAR': {
            'color': [255,0,255], 
            'dim': [3]
        },
        'LIN3DCIRCULARTUBE': {
            'dim': [3]
        },
        'LIN3DRECTANGULAR': {
            'dim': [3]
        },
        'LIN3DRECTANGULARTUBE': {
            'dim': [3]
        },
        'LIN3DTEE': {
            'dim': [3]
        },
        'LIN3DTHINAREA': {
            'dim': [3]
        },
        'LIN3DWIDEFLANGE': {
            'dim': [3]
        },
        'FIB3DLINESECTION': {
            'dim': [3]
        },
        'FIB3DAREASECTION': {
            'dim': [3]
        },
    },
    'Elements': {
        'ZEROLENGTH1D': {
            'type': 'LINE',
            'group': 31,
            'paraview': 3, 
            'dim': [1,2,3]
        },
        'LIN2DTRUSS2': {
            'type': 'LINE', 
            'group': 10,
            'paraview': 3, 
            'dim': [2]
        },
        'KIN2DTRUSS2': {
            'type': 'LINE',
            'group': 10,
            'paraview': 3, 
            'dim': [2]
        },
        'LIN3DTRUSS2': {
            'type': 'LINE' , 
            'group': 10,
            'paraview': 3, 
            'dim': [3]
        },
        'KIN3DTRUSS2': {
            'type': 'LINE', 
            'group': 10,
            'paraview': 3, 
            'dim': [3]
        },
        'LIN2DTRUSS3': {
            'type': 'LINE', 
            'group': 10,
            'paraview': 21, 
            'dim': [2]
        },
        'LIN3DTRUSS3': {
            'type': 'LINE', 
            'group': 10,
            'paraview': 21, 
            'dim': [3]
        },
        'LIN2DTRIA3': {
            'type': 'TRIA' , 
            'group': 11,
            'paraview': 5, 
            'dim': [2]
        },
        'LIN2DTRIA6': {
            'type': 'TRIA' , 
            'group': 11,
            'paraview': 22, 
            'dim': [2]
        },
        'LIN2DQUAD4': {
            'type': 'QUAD', 
            'group': 12,
            'paraview': 9, 
            'dim': [2]
        },
        'KIN2DQUAD4': {
            'type': 'QUAD' , 
            'group': 12,
            'paraview': 9, 
            'dim': [2]
        },
        'LIN2DQUAD8': {
            'type': 'QUAD', 
            'group': 12,
            'paraview': 23, 
            'dim': [2]
        },
        'PML2DQUAD4': {
            'type': 'QUAD', 
            'group': 12,
            'paraview': 9, 
            'dim': [2]
        },
        'PML2DQUAD8': {
            'type': 'QUAD', 
            'group': 32,
            'paraview': 23, 
            'dim': [2]
        },
        'LIN2DFRAME2': {
            'type': 'LINE', 
            'group': 20,
            'paraview': 3, 
            'dim': [2]
        },
        'KIN2DFRAME2': {
            'type': 'LINE', 
            'group': 20,
            'paraview': 3, 
            'dim': [2]
        },
        'LIN3DFRAME2': {
            'type': 'LINE', 
            'group': 20,
            'paraview': 3, 
            'dim': [3]
        },
        'LIN3DSHELL3': {
            'type': 'TRIA', 
            'group': 21,
            'paraview': 5, 
            'dim': [3]
        },
        'LIN3DSHELL4': {
            'type': 'QUAD', 
            'group': 21,
            'paraview': 9, 
            'dim': [3]
        },
        'LIN3DTETRA4': {
            'type': 'TETRA', 
            'group': 13,
            'paraview': 10, 
            'dim': [3]
        },
        'LIN3DTETRA10': {
            'type': 'TETRA', 
            'group': 13,
            'paraview': 24, 
            'dim': [3]
        },
        'LIN3DHEXA8': {
            'type': 'HEXA', 
            'group': 14,
            'paraview': 12, 
            'dim': [3]
        },
        'KIN3DHEXA8': {
            'type': 'HEXA', 
            'group': 14,
            'paraview': 12, 
            'dim': [3]
        },
        'LIN3DHEXA20': {
            'type': 'HEXA', 
            'group': 14,
            'paraview': 25, 
            'dim': [3]
        },
        'PML3DHEXA8': {
            'type': 'HEXA', 
            'group': 14,
            'paraview': 12, 
            'dim': [3]
        },
        'PML3DHEXA20': {
            'type': 'HEXA', 
            'group': 14,
            'paraview': 25, 
            'dim': [3]
        },
        'NULL2DFRAME2': {
            'type': 'LINE', 
            'group': 50,
            'paraview': 3, 
            'dim': [2]
        },
        'NULL3DFRAME2': {
            'type': 'LINE', 
            'group': 50,
            'paraview': 3, 
            'dim': [3]
        },
        'EQLIN2DQUAD4': {
            'type': 'QUAD', 
            'group': 12,
            'paraview': 9, 
            'dim': [2]
        },
        'TIEQLIN2DQUAD4': {
            'type': 'QUAD', 
            'group': 12,
            'paraview': 9, 
            'dim': [2]
        },
        'UNXBOUCWEN2DLINK': {
            'type': 'LINE', 
            'group': 41,
            'paraview': 3, 
            'dim': [2]
        },
        'UNXBOUCWEN3DLINK': {
            'type': 'LINE', 
            'group': 41,
            'paraview': 3, 
            'dim': [3]
        },
        'HDRBYAMAMOTO2DLINK': {
            'type': 'LINE', 
            'group': 42,
            'paraview': 3, 
            'dim': [2]
        },
        'HDRBYAMAMOTO3DLINK': {
            'type': 'LINE', 
            'group': 42, 
            'paraview': 3, 
            'dim': [3]
        }
    }
}
