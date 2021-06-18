#!/usr/bin/env python3
# -*- coding: Utf-8 -*-

import numpy as np
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
    @author Danilo S. Kusanovic 2021

    Parameters
    ----------
    filename : str
        The file name for the ParaView VTK file

    Returns
    -------
    The ParaView VTK file
    """
    #Pre-defined VTK ParaView file name
    filename = setFileExtension(filename, '.vtu')
    filepath = Options['path'] + '/' + filename
    nNodes   = len(Entities['Nodes'])
    nElems   = len(Entities['Elements'])

    print(" Generating the VTK file : " + filename)

    #Opens the ParaView File 
    Paraviewfile = open(filepath, "w+")

    #Header line skips for XML file format
    head1 = " "*2
    head2 = " "*4
    head3 = " "*6
    head4 = " "*8
    
    Paraviewfile.write("<?xml version=\"1.0\"?>\n")
    Paraviewfile.write("<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"BigEndian\">\n%s" % head1)
    Paraviewfile.write("<UnstructuredGrid>\n%s" % head2)
    
    #Write Mesh date in ParaView VTU format
    if Entities['Nodes'] and Entities['Elements']:
        #BOTH NODES AND ELEMENTS HAVE BEEN PROVIDED
        Paraviewfile.write("<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n" % (nNodes, nElems))
        Paraviewfile.write("%s<Points>\n%s" % (head3, head4))
        Paraviewfile.write("<DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n")

        count = 0
        ParaviewMap = dict()
        for nTag in Entities['Nodes']:
            coords = Entities['Nodes'][nTag]['coords']
            if Options['dimension'] == 2:
                coords = [coords[0], coords[1], 0.0]
            elif Options['dimension'] == 1:
                coords = [coords[0], 0.0, 0.0]
            Paraviewfile.write("%s  %E %E %E\n" % (head3, coords[0], coords[1], coords[2]))
            ParaviewMap[nTag] = count 
            count += 1

        Paraviewfile.write("%s</DataArray>\n%s" % (head4, head3))
        Paraviewfile.write("</Points>\n%s" % head3)
        Paraviewfile.write("<Cells>\n%s" % head4)
        Paraviewfile.write("<DataArray type=\"Int64\" Name=\"connectivity\" Format=\"ascii\">")

        for eTag in Entities['Elements']:
            connection = Entities['Elements'][eTag]['conn']
            Paraviewfile.write("\n%s " % head3)
            for conn in connection:
                Paraviewfile.write(" %d" % ParaviewMap[conn])

        Paraviewfile.write("\n%s" % head4)
        Paraviewfile.write("</DataArray>\n%s" % head4)
        Paraviewfile.write("<DataArray type=\"Int64\" Name=\"offsets\" Format=\"ascii\">\n%s " % head3)

        nOffsets = 0
        for eTag in Entities['Elements']:
            nOffsets += len(Entities['Elements'][eTag]['conn'])
            Paraviewfile.write(" %d" % nOffsets)
        Paraviewfile.write("\n%s" % head4)
        Paraviewfile.write("</DataArray>\n%s" % head4)
        Paraviewfile.write("<DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n%s " % head3)

        for eTag in Entities['Elements']:
            name = Entities['Elements'][eTag]['name']
            cell = SVLclasses['Elements'][name]['paraview']
            Paraviewfile.write(" %d" % cell)

        Paraviewfile.write("\n%s" % head4)
        Paraviewfile.write("</DataArray>\n%s" % head3)
        Paraviewfile.write("</Cells>\n%s" % head3)
        Paraviewfile.write("<PointData>\n%s" % head4)
        Paraviewfile.write("<DataArray type=\"Int64\" Name=\"GlobalNodeId\" format=\"ascii\">\n%s" % head4)

        for nTag in Entities['Nodes']:
            Paraviewfile.write("%d " % nTag)

        Paraviewfile.write("\n%s" % head4)
        Paraviewfile.write("</DataArray>\n%s" % head4)
        Paraviewfile.write("<DataArray type=\"Int32\" Name=\"DegreeOfFreedom\" format=\"ascii\">\n%s" % head4)

        for nTag in Entities['Nodes']:
            Paraviewfile.write("%d " % Entities['Nodes'][nTag]['ndof'])

        Paraviewfile.write("\n%s" % head4)
        Paraviewfile.write("</DataArray>\n%s" % head4)

        bc = np.zeros(nNodes, dtype=int)
        for nTag in Entities['Nodes']:
            freedof = Entities['Nodes'][nTag]['freedof']
            for dof in freedof:
                if dof == -1:
                    bc[ParaviewMap[nTag]] = 1
                    break
                elif dof < -1:
                    bc[ParaviewMap[nTag]] = -1
                    break

        Paraviewfile.write("<DataArray type=\"Int32\" Name=\"NodeConditions\" format=\"ascii\">\n%s" % head4)

        for n in bc:
            Paraviewfile.write("%d " % n)

        Paraviewfile.write("\n%s" % head4)
        Paraviewfile.write("</DataArray>\n%s" % head3)
        Paraviewfile.write("</PointData>\n%s" % head3)
        Paraviewfile.write("<CellData>\n%s" % head4)
        Paraviewfile.write("<DataArray type=\"Int64\" Name=\"GlobalElementId\" format=\"ascii\">\n%s" % head4)

        for eTag in Entities['Elements']:
            Paraviewfile.write("%d " % eTag)

        Paraviewfile.write("\n%s" % head4)
        Paraviewfile.write("</DataArray>\n%s" % head4)

        Paraviewfile.write("<DataArray type=\"Int32\" Name=\"ElementGroup\" format=\"ascii\">\n%s" % head4)

        for eTag in Entities['Elements']:
            name = Entities['Elements'][eTag]['name']
            Paraviewfile.write("%d " % SVLclasses['Elements'][name]['group'])

        Paraviewfile.write("\n%s" % head4)
        Paraviewfile.write("</DataArray>\n%s" % head4)

        if len(Options['partition']) != 0:
            Paraviewfile.write("<DataArray type=\"Int32\" Name=\"DomainPartition\" format=\"ascii\">\n%s" % head4)

            for k in Options['partition']:
                Paraviewfile.write("%d " % k)
            Paraviewfile.write("\n%s" % head4)
            Paraviewfile.write("</DataArray>\n%s" % head3)

        Paraviewfile.write("</CellData>\n%s" % head2)

    elif Entities['Nodes'] and not Entities['Elements']:
        #ONLY NODES HAVE BEEN PROVIDED
        Paraviewfile.write("<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n" % (nNodes, nNodes))
        Paraviewfile.write("%s<Points>\n%s" % (head3, head4))
        Paraviewfile.write("<DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n")

        count = 0
        ParaviewMap = dict()
        for nTag in Entities['Nodes']:
            coords = Entities['Nodes'][nTag]['coords']
            if Options['dimension'] == 2:
                coords = [coords[0], coords[1], 0.0]
            elif Options['dimension'] == 1:
                coords = [coords[0], 0.0, 0.0]
            Paraviewfile.write("%s  %E %E %E\n" % (head3, coords[0], coords[1], coords[2]))
            ParaviewMap[nTag] = count 
            count += 1

        Paraviewfile.write("%s</DataArray>\n%s" % (head4, head3))
        Paraviewfile.write("</Points>\n%s" % head3)
        Paraviewfile.write("<Cells>\n%s" % head4)
        Paraviewfile.write("<DataArray type=\"Int64\" Name=\"connectivity\" Format=\"ascii\">\n")

        for nTag in Entities['Nodes']:
            Paraviewfile.write("%s  %d\n" % (head3, ParaviewMap[nTag]))
        Paraviewfile.write("%s</DataArray>\n%s" % (head4, head4))
        Paraviewfile.write("<DataArray type=\"Int64\" Name=\"offsets\" Format=\"ascii\">\n%s " % head3)

        nOffsets = 0
        for nTag in Entities['Nodes']:
            nOffsets += 1
            Paraviewfile.write(" %d" % nOffsets)
        Paraviewfile.write("\n%s" % head4)
        Paraviewfile.write("</DataArray>\n%s" % head4)

        Paraviewfile.write("<DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n%s " % head3)

        for nTag in Entities['Nodes']:
            Paraviewfile.write(" 1")
        Paraviewfile.write("\n%s" % head4)
        Paraviewfile.write("</DataArray>\n%s" % head3)
        Paraviewfile.write("</Cells>\n%s" % head2)

    elif not Entities['Nodes'] and not Entities['Elements']:
        #NEITHER NODES NOR ELEMENTS HAVE BEEN PROVIDED
        print(" Entities['Nodes'] and Entities['Elements'] are empty")
        Paraviewfile.write("<Piece NumberOfPoints=\"1\" NumberOfCells=\"1\">\n")
        Paraviewfile.write("%s<Points>\n%s" % (head3, head4))
        Paraviewfile.write("<DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n")
        Paraviewfile.write("%s  0.0 0.0 0.0\n" % head3)
        Paraviewfile.write("%s</DataArray>\n%s" % (head4, head3))
        Paraviewfile.write("</Points>\n%s" % head3)
        Paraviewfile.write("<Cells>\n%s" % head4)
        Paraviewfile.write("<DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n")
        Paraviewfile.write("%s  0\n" % head3)
        Paraviewfile.write("%s</DataArray>\n%s" % (head4, head4))
        Paraviewfile.write("<DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n%s " % head3)
        Paraviewfile.write(" 1\n%s" % head4)
        Paraviewfile.write("</DataArray>\n%s" % head4)
        Paraviewfile.write("<DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n%s " % head3)
        Paraviewfile.write(" 1\n%s" % head4)
        Paraviewfile.write("</DataArray>\n%s" % head3)
        Paraviewfile.write("</Cells>\n")
        Paraviewfile.write("%s</Piece>\n" % head2)
        Paraviewfile.write("%s</UnstructuredGrid>\n" % head1)
        Paraviewfile.write("</VTKFile>")

        Paraviewfile.close()
        return

    #Close the ParaView VTU file
    Paraviewfile.write("</Piece>\n%s" % head1)
    Paraviewfile.write("</UnstructuredGrid>\n")
    Paraviewfile.write("</VTKFile>")

    Paraviewfile.close()
    print(" The VTK file can now be open in ParaView!")