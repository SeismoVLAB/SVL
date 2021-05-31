#!/usr/bin/env python3
# -*- coding: Utf-8 -*-

from Core.Utilities import debugInfo, setFileExtension
from Core.Definitions import Entities, Options, SVLclasses

def GetNumberOfFeatures():
    #Gets the elements for this partition
    nparaview = 0
    for eTag in Entities['Elements']:
        nparaview += (len(Entities['Elements'][eTag]['conn']) + 1)
    return nparaview

def renderData(filename=''):
    """
    This function render a 3D view of the model created so far.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2020

    Parameters
    ----------
    filename : str
        The file name for the ParaView VTK file

    Returns
    -------
    The ParaView VTK file
    """
    #Pre-define variables
    ThereAreNodes = False
    ThereAreElems = False
    ThereIsPartition = False

    #Pre-defined VTK ParaView file name
    filename = setFileExtension(filename, '.vtk')
    filepath = Options['path'] + '/' + filename

    print(" Generating the VTK file : " + filename)

    #Opens the ParaView File 
    Paraviewfile = open(filepath, "w+")
    Paraviewfile.write("# vtk DataFile Version 4.0\n")
    Paraviewfile.write("SeismoVLAB: ParaView Visualization\n")
    Paraviewfile.write("ASCII\n")
    Paraviewfile.write("DATASET UNSTRUCTURED_GRID\n\n")

    #Write Node Information in ParaView
    count = 0
    options = 0
    ParaviewMap = dict()
    if Entities['Nodes']:
        ThereAreNodes = True
        Paraviewfile.write("POINTS %d float\n" % len(Entities['Nodes']))
        for nTag in Entities['Nodes']:
            coords = Entities['Nodes'][nTag]['coords']
            if Options['dimension'] == 2:
                coords = [coords[0], coords[1], 0.0]
            elif Options['dimension'] == 1:
                coords = [coords[0], 0.0, 0.0]
            #Write the Node coordinates
            Paraviewfile.write("%E %E %E\n" % (coords[0], coords[1], coords[2]))
            #Creates the Node map
            ParaviewMap[nTag] = count 
            count += 1
        Paraviewfile.write("\n")

    #Write Element Information in ParaView
    if ThereAreNodes and Entities['Elements']:
        options += 3 
        ThereAreElems = True
        nfeatures = Options['nfeatures']
        if Options['nfeatures'] == 0:
            nfeatures = GetNumberOfFeatures()

        Paraviewfile.write("CELLS %d %d\n" % (len(Entities['Elements']), nfeatures))
        for eTag in Entities['Elements']:
            connection = Entities['Elements'][eTag]['conn']
            Paraviewfile.write("%d" % len(connection))
            for nTag in connection:
                Paraviewfile.write(" %d" % ParaviewMap[nTag])
            Paraviewfile.write("\n")
        Paraviewfile.write("\n")

        #Element Rendering Type.
        Paraviewfile.write("CELL_TYPES %d\n" % len(Entities['Elements']))
        for eTag in Entities['Elements']:
            name = Entities['Elements'][eTag]['name']
            cell = SVLclasses['Elements'][name]['paraview']
            Paraviewfile.write("%d\n" % cell)
        Paraviewfile.write("\n")
    elif ThereAreNodes and not ThereAreElems:
        nNodes = len(Entities['Nodes'])
        Paraviewfile.write("CELLS %d %d\n" % (nNodes, 2*nNodes))
        for nTag in Entities['Nodes']:
            Paraviewfile.write("1 %d\n" % ParaviewMap[nTag])
        Paraviewfile.write("\n")

        Paraviewfile.write("CELL_TYPES %d\n" % nNodes)
        for nTag in Entities['Nodes']:
            Paraviewfile.write("1\n")
    elif not ThereAreNodes and not ThereAreElems:
        print(" Entities['Nodes'] and Entities['Elements'] are empty")
        Paraviewfile.close()
        return

    #Scalar Attributes (Node Tags) applied to Nodes.
    Paraviewfile.write("POINT_DATA %d\n" % len(Entities['Nodes']))
    Paraviewfile.write("SCALARS NodeIds int 1\n")
    Paraviewfile.write("LOOKUP_TABLE pointIds\n")
    for nTag in Entities['Nodes']:
        Paraviewfile.write("%d\n" % nTag)
    Paraviewfile.write("\n")

    #Scalar Attributes (Partition) applied to Element.
    if len(Options['partition']) > 0:
        options += 1 
        ThereIsPartition = True

    if ThereAreElems:
        Paraviewfile.write("CELL_DATA %d\n" % len(Entities['Elements']))
        Paraviewfile.write("FIELD attributes %d\n" % options)

        #Mesh Partition
        if ThereIsPartition:
            Paraviewfile.write("Partition 1 %d int\n" % len(Entities['Elements']))
            for pid in Options['partition']:
                Paraviewfile.write("%d\n" % pid)
            Paraviewfile.write("\n")

        #Element's Material
        Paraviewfile.write("Materials 1 %d int\n" % len(Entities['Elements']))
        for eTag in Entities['Elements']:
            if 'material' in Entities['Elements'][eTag]['attributes']:
                Paraviewfile.write("%d\n" % Entities['Elements'][eTag]['attributes']['material'])
            else:
                Paraviewfile.write("-1\n")
        Paraviewfile.write("\n")

        #Element's Section
        Paraviewfile.write("Sections 1 %d int\n" % len(Entities['Elements']))
        for eTag in Entities['Elements']:
            if 'section' in Entities['Elements'][eTag]['attributes']:
                Paraviewfile.write("%d\n" % Entities['Elements'][eTag]['attributes']['section'])
            else:
                Paraviewfile.write("-1\n")
        Paraviewfile.write("\n")

        #Mesh Partition
        Paraviewfile.write("elemIds 1 %d int\n" % len(Entities['Elements']))
        for eTag in Entities['Elements']:
            Paraviewfile.write("%d\n" % eTag)
    Paraviewfile.close()
    print(" The VTK file can now be open in ParaView!\n")
