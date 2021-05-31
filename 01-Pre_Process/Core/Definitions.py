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
            'color': [0,255,255], 
            'dim'  : [1,2,3]
        },
        'HERTZIAN1DLINEAR': {
            'color': [29,249,20], 
            'dim': [1,2,3]
        },
        'VISCOUS1DLINEAR': {
            'color': [255,0,255], 
            'dim': [1,2,3]
        },
        'PLASTIC1DJ2': {
            'color': [146,110,174], 
            'dim': [1,2,3]
        },
        'ELASTIC2DPLANESTRAIN': {
            'color': [0,255,0], 
            'dim': [2,3]
        },
        'ELASTIC2DPLANESTRESS': {
            'color': [255,0,0], 
            'dim': [2,3]
        },
        'PLASTICPLANESTRAINJ2': {
            'color': [59,176,143], 
            'dim': [2]
        },
        'PLASTICPLANESTRAINBA': {
            'color': [165,105,79], 
            'dim': [2]
        },
        'ELASTIC3DLINEAR': {
            'color': [0,0,255], 
            'dim': [3]
        },
        'PLASTIC3DJ2': {
            'color': [237,237,237], 
            'dim': [3]
        },
        'PLASTIC3DBA': {
            'color': [255,182, 83], 
            'dim': [3]
        },
        'ELASTIC1DFIBER': {
            'color': [153,76,0], 
            'dim': [1,2,3]
        },
        'STEEL1DFIBER': {
            'color': [0,153,0], 
            'dim': [1,2,3]
        },
        'CONCRETE1DFIBER': {
            'color': [128,128, 128], 
            'dim': [1,2,3]
        },
        'ELASTIC1DGAP': {
            'color': [183, 149, 11], 
            'dim': [1,2,3]
        },
        'PLASTIC1DGAP': {
            'color': [44, 62, 80], 
            'dim': [1,2,3]
        }
    },
    'Sections': {
        'LIN2DANGLE': {
            'color': [0,255,255], 
            'dim': [2]
        },
        'LIN2DCHANNEL': {
            'color': [29,249,20], 
            'dim': [2]
        },
        'LIN2DCIRCULAR': {
            'color': [255,0,255], 
            'dim': [2]
        },
        'LIN2DCIRCULARTUBE': {
            'color': [146,110,174], 
            'dim': [2]
        },
        'LIN2DRECTANGULAR': {
            'color': [255,0,0], 
            'dim': [2]
        },
        'LIN2DRECTANGULARTUBE' : {
            'color': [0,255,0], 
            'dim': [2]
        },
        'LIN2DTEE': {
            'color': [165,105,79], 
            'dim': [2]
        },
        'LIN2DWIDEFLANGE': {
            'color': [255,182,83], 
            'dim': [2]
        },
        'LIN3DANGLE': {
            'color': [0,255,255], 
            'dim': [3]
        },
        'LIN3DCHANNEL': {
            'color': [29,249,20], 
            'dim': [3]
        },
        'LIN3DCIRCULAR': {
            'color': [255,0,255], 
            'dim': [3]
        },
        'LIN3DCIRCULARTUBE': {
            'color': [146,110,174], 
            'dim': [3]
        },
        'LIN3DRECTANGULAR': {
            'color': [255,0,0], 
            'dim': [3]
        },
        'LIN3DRECTANGULARTUBE': {
            'color': [0,255,0], 
            'dim': [3]
        },
        'LIN3DTEE': {
            'color': [165,105,79], 
            'dim': [3]
        },
        'LIN3DTHINAREA': {
            'color': [0,0,255], 
            'dim': [3]
        },
        'LIN3DWIDEFLANGE': {
            'color': [255,182,83], 
            'dim': [3]
        },
        'FIB3DLINESECTION': {
            'color': [255,0,255], 
            'dim': [3]
        },
        'FIB3DAREASECTION': {
            'color': [221,160,221], 
            'dim': [3]
        },
    },
    'Elements': {
        'ZEROLENGTH1D': {
            'type': 'LINE' , 
            'color': [29,249, 20], 
            'paraview': 3, 
            'dim': [1,2,3]
        },
        'LIN2DTRUSS2': {
            'type': 'LINE', 
            'color': [25,25,25], 
            'paraview': 3, 
            'dim': [2]
        },
        'KIN2DTRUSS2': {
            'type': 'LINE', 
            'color': [25,25,25], 
            'paraview': 3, 
            'dim': [2]
        },
        'LIN3DTRUSS2': {
            'type': 'LINE' , 
            'color': [25,25,25], 
            'paraview': 3, 
            'dim': [3]
        },
        'KIN3DTRUSS2': {
            'type': 'LINE', 
            'color': [25,25,25], 
            'paraview': 3, 
            'dim': [3]
        },
        'LIN2DTRUSS3': {
            'type': 'LINE', 
            'color': [ 25, 25, 25], 
            'paraview': 4, 
            'dim': [2]
        },
        'LIN3DTRUSS3': {
            'type': 'LINE', 
            'color': [25,25,25], 
            'paraview': 4, 
            'dim': [3]
        },
        'LIN2DTRIA3': {
            'type': 'TRIA' , 
            'color': [146,110,174], 
            'paraview': 4, 
            'dim': [2]
        },
        'LIN2DTRIA6': {
            'type': 'TRIA' , 
            'color': [146,110,174], 
            'paraview': 7, 
            'dim': [2]
        },
        'LIN2DQUAD4': {
            'type': 'QUAD', 
            'color': [ 59,176,143], 
            'paraview': 5, 
            'dim': [2]
        },
        'KIN2DQUAD4': {
            'type': 'QUAD' , 
            'color': [ 59,176,143], 
            'paraview': 5, 
            'dim': [2]
        },
        'LIN2DQUAD8': {
            'type': 'QUAD', 
            'color': [ 59,176,143], 
            'paraview': 9, 
            'dim': [2]
        },
        'PML2DQUAD4': {
            'type': 'QUAD', 
            'color': [255,182, 83], 
            'paraview': 5, 
            'dim': [2]
        },
        'PML2DQUAD8': {
            'type': 'QUAD', 
            'color': [255,182, 83], 
            'paraview': 9, 
            'dim': [2]
        },
        'LIN2DFRAME2': {
            'type': 'LINE', 
            'color': [165,105, 79], 
            'paraview': 3, 
            'dim': [2]
        },
        'KIN2DFRAME2': {
            'type': 'LINE', 
            'color': [165,105, 79], 
            'paraview': 3, 'dim': [2]
        },
        'LIN3DFRAME2': {
            'type': 'LINE', 
            'color': [165,105, 79], 
            'paraview': 3, 
            'dim': [3]
        },
        'LIN3DSHELL3': {
            'type': 'TRIA', 
            'color': [128,218,235], 
            'paraview': 4, 
            'dim': [3]
        },
        'LIN3DSHELL4': {
            'type': 'QUAD', 
            'color': [128,218,235], 
            'paraview': 5, 
            'dim': [3]
        },
        'LIN3DTETRA4': {
            'type': 'TETRA', 
            'color': [255,155,170], 
            'paraview': 5, 
            'dim': [3]
        },
        'LIN3DTETRA10': {
            'type': 'TETRA', 
            'color': [255,155,170], 
            'paraview': 11, 
            'dim': [3]
        },
        'LIN3DHEXA8': {
            'type': 'HEXA', 
            'color': [255,255,159], 
            'paraview': 9, 
            'dim': [3]
        },
        'KIN3DHEXA8': {
            'type': 'HEXA', 
            'color': [255,255,159], 
            'paraview': 9, 
            'dim': [3]
        },
        'LIN3DHEXA20': {
            'type': 'HEXA', 
            'color': [255,255,159], 
            'paraview': 21, 
            'dim': [3]
        },
        'PML3DHEXA8': {
            'type': 'HEXA', 
            'color': [176,183,198], 
            'paraview': 9, 
            'dim': [3]
        },
        'PML3DHEXA20': {
            'type': 'HEXA', 
            'color': [176,183,198], 
            'paraview': 21, 
            'dim': [3]
        },
        'NULL2DFRAME2': {
            'type': 'LINE', 
            'color': [29,249, 20], 
            'paraview': 3, 
            'dim': [2]
        },
        'NULL3DFRAME2': {
            'type': 'LINE', 
            'color': [29,249, 20], 
            'paraview': 3, 
            'dim': [3]
        },
        'EQLIN2DQUAD4': {
            'type': 'QUAD', 
            'color': [59,176,143], 
            'paraview': 5, 
            'dim': [2]
        },
        'TIEQLIN2DQUAD4': {
            'type': 'QUAD', 
            'color': [59,176,143], 
            'paraview': 5, 
            'dim': [2]
        },
        'UNXBOUCWEN2DLINK': {
            'type': 'LINE', 
            'color': [0,0,0], 
            'paraview': 3, 
            'dim': [2]
        },
        'UNXBOUCWEN3DLINK': {
            'type': 'LINE', 
            'color': [0,0,0], 
            'paraview': 3, 
            'dim': [3]
        },
        'HDRBYAMAMOTO2DLINK': {
            'type': 'LINE', 
            'color': [0,0,0], 
            'paraview': 3, 
            'dim': [2]
        },
        'HDRBYAMAMOTO3DLINK': {
            'type': 'LINE', 
            'color': [0,0,0], 
            'paraview': 3, 
            'dim': [3]
        }
    }
}