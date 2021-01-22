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
    'format'      : 'SVL',
    'description' : '\n',
    'execfile'    : '',
    'execpath'    : '',
    'allocation'  : 'NO',
    'numbering'   : 'Plain',
    'metispath'   : '',
    'massform'    : 'consistent',
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
    # @author Danilo S. Kusanovic 2020
    #  'Materials' : Stores the material classes defined in SVL
    #  'Sections'  : Stores the section classes defined in SVL
    #  'Elements'  : Stores the elements classes defined in SVL
    'Materials': [
        'ELASTIC1DLINEAR',
        'HERTZIAN1DLINEAR',
        'VISCOUS1DLINEAR',
        'PLASTIC1DJ2',
        'ELASTIC2DPLANESTRAIN',
        'ELASTIC2DPLANESTRESS',
        'PLASTICPLANESTRAINJ2',
        'PLASTICPLANESTRAINBA',
        'ELASTIC3DLINEAR',
        'PLASTIC3DJ2',
        'PLASTIC3DBA'
        ],
    'Sections': [
        'LIN2DANGLE',
        'LIN2DCHANNEL',
        'LIN2DCIRCULAR',
        'LIN2DCIRCULARTUBE',
        'LIN2DRECTANGULAR',
        'LIN2DRECTANGULARTUBE',
        'LIN2DTEE',
        'LIN2DWIDEFLANGE',
        'LIN3DANGLE',
        'LIN3DCHANNEL',
        'LIN3DCIRCULAR',
        'LIN3DCIRCULARTUBE',
        'LIN3DRECTANGULAR',
        'LIN3DRECTANGULARTUBE',
        'LIN3DTEE',
        'LIN3DTHINAREA',
        'LIN3DWIDEFLANGE',
        ],
    'Elements': [
        'ZEROLENGTH1D',
        'LIN2DTRUSS2',
        'KIN2DTRUSS2',
        'LIN3DTRUSS2',
        'KIN3DTRUSS2',
        'LIN2DTRUSS3',
        'LIN3DTRUSS3',
        'LIN2DTRIA3',
        'LIN2DTRIA6',
        'LIN2DQUAD4',
        'KIN2DQUAD4',
        'LIN2DQUAD8',
        'PML2DQUAD4',
        'PML2DQUAD8',
        'LIN2DFRAME2',
        'KIN2DFRAME2',
        'LIN3DFRAME2',
        'LIN3DSHELL3',
        'LIN3DSHELL4',
        'LIN3DTETRA4',
        'LIN3DTETRA10',
        'LIN3DHEXA8',
        'KIN3DHEXA8',
        'LIN3DHEXA20',
        'PML3DHEXA8',
        'PML3DHEXA20',
        'UNXBOUCWEN2DLINK',
        'UNXBOUCWEN3DLINK',
        'HDRBYAMAMOTO2DLINK',
        'HDRBYAMAMOTO3DLINK',
        'EQLIN2DQUAD4',
        'TIEQLIN2DQUAD4',
        'NULL2DFRAME2',
        'NULL3DFRAME2'
        ]
}

VTKcolors = {
    #Default VTK colors associated to an option
    # https://github.com/SeismoVLAB/SVL
    # @author Danilo S. Kusanovic 2020
    #  'part' : can handle up to 23 partitions
    #  'elem' : defines the color for element type
    'part': {
        0      : [  0,  0,255],
        1      : [  0,255,  0],
        2      : [255,  0,  0],
        3      : [  0,255,255],
        4      : [255,  0,255],
        5      : [255,255,  0],
        6      : [ 29,249, 20], 
        7      : [237,237,237],
        8      : [165,105, 79],
        9      : [146,110,174],
        10     : [ 59,176,143],
        11     : [255,182, 83],
        12     : [128,218,235],
        13     : [255,155,170],
        14     : [255,255,159],
        15     : [176,183,198],
        16     : [ 25, 25, 25],
        17     : [  0,  0,128],
        18     : [  0,128,  0],
        19     : [128,  0,  0],
        20     : [  0,128,128],
        21     : [128,  0,128],
        22     : [128,128,  0]
    },
    'elem': {
        'NONE' : [ 29,249, 20], 
        'TRUSS': [ 25, 25, 25],
        'FRAME': [165,105, 79],
        'TRIA' : [146,110,174],
        'QUAD' : [ 59,176,143],
        'QPML' : [255,182, 83],
        'SHELL': [128,218,235],
        'TETRA': [255,155,170],
        'HEXA' : [255,255,159],
        'HPML' : [176,183,198],
        'LINK' : [ 25, 25, 25]
    }
}

VTKelems = {
    #VTK dictionary with associated element properties
    # https://github.com/SeismoVLAB/SVL
    # @author Danilo S. Kusanovic 2020
    #  'VTKname'  : stores the basic element type (LINE, QUAD, HEXA) for the SVL (element class) 
    #  'VTKtype'  : stores the VTK cell element type for rendering
    #  'VTKcolor' : stores the color to be applied to each element in domain
    #  'VTKsvl'   : defines a map between the SVL (element class) and VTKcolors (dictionary)
    'VTKname'  : {},
    'VTKtype'  : {},
    'VTKcolor' : {},
    'VTKsvl'   : {
        'ZEROLENGTH1D': 'LINE',
        'LIN2DTRUSS2' : 'LINE',
        'LIN2DTRUSS3' : 'LINE',
        'KIN2DTRUSS2' : 'LINE',
        'LIN3DTRUSS2' : 'LINE',
        'KIN3DTRUSS2' : 'LINE',
        'LIN3DTRUSS3' : 'LINE',
        'LIN2DTRIA3'  : 'TRIA',
        'LIN2DTRIA6'  : 'TRIA',
        'LIN2DQUAD4'  : 'QUAD',
        'KIN2DQUAD4'  : 'QUAD',
        'LIN2DQUAD8'  : 'QUAD',
        'PML2DQUAD4'  : 'QUAD',
        'PML2DQUAD8'  : 'QUAD',
        'LIN2DFRAME2' : 'LINE',
        'KIN2DFRAME2' : 'LINE',
        'LIN3DFRAME2' : 'LINE',
        'LIN3DSHELL3' : 'TRIA',
        'LIN3DSHELL4' : 'QUAD',
        'LIN3DTETRA4' : 'TETRA',
        'LIN3DTETRA10': 'TETRA',
        'LIN3DHEXA8'  : 'HEXA',
        'KIN3DHEXA8'  : 'HEXA',
        'LIN3DHEXA20' : 'HEXA',
        'PML3DHEXA8'  : 'HEXA',
        'PML3DHEXA20' : 'HEXA',
        'EQLIN2DQUAD4'  : 'QUAD',
        'TIEQLIN2DQUAD4': 'QUAD',
        'NULL2DFRAME2'  : 'LINE',
        'NULL3DFRAME2'  : 'LINE',
        'UNXBOUCWEN2DLINK'  : 'LINE',
        'UNXBOUCWEN3DLINK'  : 'LINE',
        'HDRBYAMAMOTO2DLINK': 'LINE',
        'HDRBYAMAMOTO3DLINK': 'LINE'
    }
}
