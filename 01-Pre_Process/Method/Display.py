#!/usr/bin/env python3
# -*- coding: Utf-8 -*-

import vtk
from Core.Utilities import debugInfo
from Core.Definitions import Entities, Options, VTKelems, VTKcolors

def setVTKtype(option=None):
    """
    This function classifies each element depending on its name and
    associate to it color, type, and name for VTK.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    option : str
        The color to be used in display: 
        option = None, Elements, Partition

    Returns
    -------
    None
    """
    count = 0
    for eTag in Entities['Elements']:
        name = Entities['Elements'][eTag]['name']
        if name == 'ZEROLENGTH1D':
            VTKelems['VTKtype'][eTag] = 3
            VTKelems['VTKname'][eTag] = VTKelems['VTKsvl'][name]
            if option.upper() == 'ELEMENT':
                VTKelems['VTKcolor'][eTag] = VTKcolors['elem']['TRUSS']
        elif name == 'LIN2DTRUSS2':
            VTKelems['VTKtype'][eTag] = 3
            VTKelems['VTKname'][eTag] = VTKelems['VTKsvl'][name]
            if option.upper() == 'ELEMENT':
                VTKelems['VTKcolor'][eTag] = VTKcolors['elem']['TRUSS']
        elif name == 'KIN2DTRUSS2':
            VTKelems['VTKtype'][eTag] = 3
            VTKelems['VTKname'][eTag] = VTKelems['VTKsvl'][name]
            if option.upper() == 'ELEMENT':
                VTKelems['VTKcolor'][eTag] = VTKcolors['elem']['TRUSS']
        elif name == 'LIN3DTRUSS2':
            VTKelems['VTKtype'][eTag] = 3
            VTKelems['VTKname'][eTag] = VTKelems['VTKsvl'][name]
            if option.upper() == 'ELEMENT':
                VTKelems['VTKcolor'][eTag] = VTKcolors['elem']['TRUSS']
        elif name == 'KIN3DTRUSS2':
            VTKelems['VTKtype'][eTag] = 3
            VTKelems['VTKname'][eTag] = VTKelems['VTKsvl'][name]
            if option.upper() == 'ELEMENT':
                VTKelems['VTKcolor'][eTag] = VTKcolors['elem']['TRUSS']
        elif name == 'LIN3DTRUSS3':
            VTKelems['VTKtype'][eTag] = 21
            VTKelems['VTKname'][eTag] = VTKelems['VTKsvl'][name]
            if option.upper() == 'ELEMENT':
                VTKelems['VTKcolor'][eTag] = VTKcolors['elem']['TRUSS']
        elif name == 'LIN2DTRIA3':
            VTKelems['VTKtype'][eTag] = 5
            VTKelems['VTKname'][eTag] = VTKelems['VTKsvl'][name]
            if option.upper() == 'ELEMENT':
                VTKelems['VTKcolor'][eTag] = VTKcolors['elem']['TRIA']
        elif name == 'LIN2DTRIA6':
            VTKelems['VTKtype'][eTag] = 22
            VTKelems['VTKname'][eTag] = VTKelems['VTKsvl'][name]
            if option.upper() == 'ELEMENT':
                VTKelems['VTKcolor'][eTag] = VTKcolors['elem']['TRIA']
        elif name == 'LIN2DQUAD4':
            VTKelems['VTKtype'][eTag] = 9
            VTKelems['VTKname'][eTag] = VTKelems['VTKsvl'][name]
            if option.upper() == 'ELEMENT':
                VTKelems['VTKcolor'][eTag] = VTKcolors['elem']['QUAD']
        elif name == 'LIN2DQUAD8':
            VTKelems['VTKtype'][eTag] = 23
            VTKelems['VTKname'][eTag] = VTKelems['VTKsvl'][name]
            if option.upper() == 'ELEMENT':
                VTKelems['VTKcolor'][eTag] = VTKcolors['elem']['QUAD']
        elif name == 'PML2DQUAD4':
            VTKelems['VTKtype'][eTag] = 9
            VTKelems['VTKname'][eTag] = VTKelems['VTKsvl'][name]
            if option.upper() == 'ELEMENT':
                VTKelems['VTKcolor'][eTag] = VTKcolors['elem']['QPML']
        elif name == 'PML2DQUAD8':
            VTKelems['VTKtype'][eTag] = 23
            VTKelems['VTKname'][eTag] = VTKelems['VTKsvl'][name]
            if option.upper() == 'ELEMENT':
                VTKelems['VTKcolor'][eTag] = VTKcolors['elem']['QPML']
        elif name == 'LIN2DFRAME2':
            VTKelems['VTKtype'][eTag] = 3
            VTKelems['VTKname'][eTag] = VTKelems['VTKsvl'][name]
            if option.upper() == 'ELEMENT':
                VTKelems['VTKcolor'][eTag] = VTKcolors['elem']['FRAME']
        elif name == 'KIN2DFRAME2':
            VTKelems['VTKtype'][eTag] = 3
            VTKelems['VTKname'][eTag] = VTKelems['VTKsvl'][name]
            if option.upper() == 'ELEMENT':
                VTKelems['VTKcolor'][eTag] = VTKcolors['elem']['FRAME']
        elif name == 'LIN3DFRAME2':
            VTKelems['VTKtype'][eTag] = 3
            VTKelems['VTKname'][eTag] = VTKelems['VTKsvl'][name]
            if option.upper() == 'ELEMENT':
                VTKelems['VTKcolor'][eTag] = VTKcolors['elem']['FRAME']
        elif name == 'LIN3DSHELL3':
            VTKelems['VTKtype'][eTag] = 5
            VTKelems['VTKname'][eTag] = VTKelems['VTKsvl'][name]
            if option.upper() == 'ELEMENT':
                VTKelems['VTKcolor'][eTag] = VTKcolors['elem']['SHELL']
        elif name == 'LIN3DSHELL4':
            VTKelems['VTKtype'][eTag] = 9
            VTKelems['VTKname'][eTag] = VTKelems['VTKsvl'][name]
            if option.upper() == 'ELEMENT':
                VTKelems['VTKcolor'][eTag] = VTKcolors['elem']['SHELL']
        elif name == 'LIN3DTETRA4':
            VTKelems['VTKtype'][eTag] = 10
            VTKelems['VTKname'][eTag] = VTKelems['VTKsvl'][name]
            if option.upper() == 'ELEMENT':
                VTKelems['VTKcolor'][eTag] = VTKcolors['elem']['TETRA']
        elif name == 'LIN3DTETRA10':
            VTKelems['VTKtype'][eTag] = 24
            VTKelems['VTKname'][eTag] = VTKelems['VTKsvl'][name]
            if option.upper() == 'ELEMENT':
                VTKelems['VTKcolor'][eTag] = VTKcolors['elem']['TETRA']
        elif name == 'LIN3DHEXA8':
            VTKelems['VTKtype'][eTag] = 12
            VTKelems['VTKname'][eTag] = VTKelems['VTKsvl'][name]
            if option.upper() == 'ELEMENT':
                VTKelems['VTKcolor'][eTag] = VTKcolors['elem']['HEXA']
        elif name == 'LIN3DHEXA20':
            VTKelems['VTKtype'][eTag] = 25
            VTKelems['VTKname'][eTag] = VTKelems['VTKsvl'][name]
            if option.upper() == 'ELEMENT':
                VTKelems['VTKcolor'][eTag] = VTKcolors['elem']['HEXA']
        elif name == 'PML3DHEXA8':
            VTKelems['VTKtype'][eTag] = 12
            VTKelems['VTKname'][eTag] = VTKelems['VTKsvl'][name]
            if option.upper() == 'ELEMENT':
                VTKelems['VTKcolor'][eTag] = VTKcolors['elem']['HPML']
        elif name == 'PML3DHEXA20':
            VTKelems['VTKtype'][eTag] = 25
            VTKelems['VTKname'][eTag] = VTKelems['VTKsvl'][name]
            if option.upper() == 'ELEMENT':
                VTKelems['VTKcolor'][eTag] = VTKcolors['elem']['HPML']
        elif name == 'UNXBOUCWEN2DLINK':
            VTKelems['VTKtype'][eTag] = 3
            VTKelems['VTKname'][eTag] = VTKelems['VTKsvl'][name]
            if option.upper() == 'ELEMENT':
                VTKelems['VTKcolor'][eTag] = VTKcolors['elem']['LINK']
        elif name == 'UNXBOUCWEN3DLINK':
            VTKelems['VTKtype'][eTag] = 3
            VTKelems['VTKname'][eTag] = VTKelems['VTKsvl'][name]
            if option.upper() == 'ELEMENT':
                VTKelems['VTKcolor'][eTag] = VTKcolors['elem']['LINK']
        elif name == 'HDRBYAMAMOTO2DLINK':
            VTKelems['VTKtype'][eTag] = 3
            VTKelems['VTKname'][eTag] = VTKelems['VTKsvl'][name]
            if option.upper() == 'ELEMENT':
                VTKelems['VTKcolor'][eTag] = VTKcolors['elem']['LINK']
        elif name == 'HDRBYAMAMOTO3DLINK':
            VTKelems['VTKtype'][eTag] = 3
            VTKelems['VTKname'][eTag] = VTKelems['VTKsvl'][name]
            if option.upper() == 'ELEMENT':
                VTKelems['VTKcolor'][eTag] = VTKcolors['elem']['LINK']
        elif name == 'EQLIN2DQUAD4':
            VTKelems['VTKtype'][eTag] = 9
            VTKelems['VTKname'][eTag] = VTKelems['VTKsvl'][name]
            if option.upper() == 'ELEMENT':
                VTKelems['VTKcolor'][eTag] = VTKcolors['elem']['QUAD']
        elif name == 'TIEQLIN2DQUAD4':
            VTKelems['VTKtype'][eTag] = 9
            VTKelems['VTKname'][eTag] = VTKelems['VTKsvl'][name]
            if option.upper() == 'ELEMENT':
                VTKelems['VTKcolor'][eTag] = VTKcolors['elem']['QUAD']
        elif name == 'NULL2DFRAME2':
            VTKelems['VTKtype'][eTag] = 3
            VTKelems['VTKname'][eTag] = VTKelems['VTKsvl'][name]
            if option.upper() == 'ELEMENT':
                VTKelems['VTKcolor'][eTag] = VTKcolors['elem']['NONE']
        elif name == 'NULL3DFRAME2':
            VTKelems['VTKtype'][eTag] = 3
            VTKelems['VTKname'][eTag] = VTKelems['VTKsvl'][name]
            if option.upper() == 'ELEMENT':
                VTKelems['VTKcolor'][eTag] = VTKcolors['elem']['NONE']
        elif name == 'LINK':
            VTKelems['VTKtype'][eTag] = 3
            VTKelems['VTKname'][eTag] = VTKelems['VTKsvl'][name]
            if option.upper() == 'ELEMENT':
                VTKelems['VTKcolor'][eTag] = VTKcolors['elem']['NONE']

        if option.upper() == 'PARTITION':
            pTag = Options['partition'][count]
            VTKelems['VTKcolor'][eTag] = VTKcolors['part'][pTag]
        count += 1

def renderData(option=None):
    """
    This function render a 3D view of the model created so far.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    option : str
        The color to be used in display: 
        option = None, Elements, Partition

    Returns
    -------
    None
    """
    #Checks minimum parameters are provided
    info = debugInfo(2) 
    if not option:
        option = 'element'
    if option.upper() not in ['ELEMENT', 'PARTITION']:
        print('\x1B[33m ALERT \x1B[0m: In file=\'%s\' at line=%d renderData(option=?) is not recognize.' %(info.filename,info.lineno))
        option = 'element'
    
    #Sets color scheme and VTK types
    setVTKtype(option)

    #Import VTK color schemes
    colors = vtk.vtkNamedColors()

    renderer = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.SetSize(800, 800)
    renderWindow.SetWindowName("Seismo-VLAB: Model Render Windows")
    renderWindow.AddRenderer(renderer)

    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)
    axes = vtk.vtkAxesActor()

    #Creates the Coordinate System
    rgba = [0, 0, 0, 0]
    colors.GetColor("Carrot", rgba)

    widget = vtk.vtkOrientationMarkerWidget()
    widget.SetOutlineColor(rgba[0], rgba[1], rgba[2])
    widget.SetOrientationMarker(axes)
    widget.SetInteractor(renderWindowInteractor)
    widget.SetViewport(0.0, 0.0, 0.35, 0.35)
    widget.SetEnabled(1)
    widget.InteractiveOn()

    #Color list for elements
    Colors = vtk.vtkUnsignedCharArray()
    Colors.SetNumberOfComponents(3)

    #Generates the Node list for VTK
    points = vtk.vtkPoints()

    count = 0
    VTKmap = dict()
    for nTag in Entities['Nodes']:
        VTKmap[nTag] = count
        coordinates = Entities['Nodes'][nTag]['coords']
        if Options['dimension'] == 1:
            coordinates = [coordinates[0], 0.0, 0.0]
        if Options['dimension'] == 2:
            coordinates = [coordinates[0], coordinates[1], 0.0]
        elif Options['dimension'] == 3:
            coordinates = [coordinates[0], coordinates[1], coordinates[2]]
        points.InsertPoint(count, coordinates)
        count += 1

    #Generate the VTK cell list
    mesh = vtk.vtkUnstructuredGrid()

    nElems = len(Entities['Elements'])
    mesh.Allocate(nElems)

    for eTag in Entities['Elements']:
        name = VTKelems['VTKname'][eTag]
        conn = Entities['Elements'][eTag]['conn']
        if name == 'LINE':
            connection = [VTKmap[conn[0]], VTKmap[conn[1]]]
            mesh.InsertNextCell(vtk.VTK_LINE, 2, connection)
        elif name == 'TRIA':
            connection = [VTKmap[conn[0]], VTKmap[conn[1]], VTKmap[conn[2]]]
            mesh.InsertNextCell(vtk.VTK_TRIANGLE, 3, connection)
        elif name == 'QUAD':
            connection = [VTKmap[conn[0]], VTKmap[conn[1]], VTKmap[conn[2]], VTKmap[conn[3]]]
            mesh.InsertNextCell(vtk.VTK_QUAD, 4, connection)
        elif name == 'TETRA':
            connection = [VTKmap[conn[0]], VTKmap[conn[1]], VTKmap[conn[2]], VTKmap[conn[3]]]
            mesh.InsertNextCell(vtk.VTK_TETRA, 4, connection)
        elif name == 'HEXA':
            connection = [VTKmap[conn[0]], VTKmap[conn[1]], VTKmap[conn[2]], VTKmap[conn[3]], VTKmap[conn[4]], VTKmap[conn[5]], VTKmap[conn[6]], VTKmap[conn[7]]]
            mesh.InsertNextCell(vtk.VTK_HEXAHEDRON, 8, connection)
        Colors.InsertNextTuple3(VTKelems['VTKcolor'][eTag][0], VTKelems['VTKcolor'][eTag][1], VTKelems['VTKcolor'][eTag][2])
    mesh.SetPoints(points)

    #Assign the colors scheme
    mesh.GetCellData().SetScalars(Colors)
    mesh.Modified()
    
    meshMapper = vtk.vtkDataSetMapper()
    meshMapper.SetInputData(mesh)

    meshActor = vtk.vtkActor()
    meshActor.SetMapper(meshMapper)
    meshActor.GetProperty().EdgeVisibilityOn()
    meshActor.GetProperty().SetPointSize(5)
    meshActor.GetProperty().SetLineWidth(0.5)
    meshActor.GetProperty().SetSpecularColor(1.0, 1.0, 1.0)
    meshActor.GetProperty().SetSpecularPower(100)

    #Creates the Renderer
    renderer.AddActor(meshActor)
    renderer.SetBackground2(colors.GetColor3d('LightBlue'))
    renderer.SetBackground(colors.GetColor3d('LavenderBlush'))
    renderer.GradientBackgroundOn()

    #Loacte the Camera
    renderer.ResetCamera()
    renderer.GetActiveCamera().Elevation(95.0)
    renderer.GetActiveCamera().Azimuth(30.0)
    renderer.GetActiveCamera().Dolly(1.2)

    #Starts and Display the Finite Element Mesh
    renderWindowInteractor.Initialize()
    renderWindow.Render()
    renderWindowInteractor.Start()