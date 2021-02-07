#!/usr/bin/python3
# -*- coding: Utf-8 -*-

import os
import sys
import subprocess

def main():
    """
    This function generates a PDF report of all debugging cases.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2021

    Parameters
    ----------
    files : list
        The list of files (names) to be run 

    Returns
    -------
    None
    """
    #The current working path.
    cwd = os.path.abspath(os.path.dirname(sys.argv[0]))

    #List of Debugging Seismo-VLAB File to be Run.
    files = []
    files.append("A01-DY_Lin_2D_Elastic_ZeroLength")
    files.append("A02-DY_Lin_2D_Viscous_ZeroLength")
    files.append("A03-DY_Lin_2D_Plastic_ZeroLength")
    files.append("A04-DY_Lin_1DPointMass_Elastic_ZeroLength")
    files.append("A05-DY_Lin_Hertzian_Contact_ZeroLength")
    files.append("A06-DY_Lin_1DPointMass_Elastic_ZeroLength_SupportMotion")
    files.append("A08-DY_2D_UniAxial_BoucWen_Link_Param1")
    files.append("A10-DY_3D_UniAxial_BoucWen_Link_Param1")
    files.append("A12-DY_2D_UniAxial_YamamotoHDRB_Link")
    files.append("A13-DY_3D_UniAxial_YamamotoHDRB_Link")
    files.append("C01-ST_Lin_2DAxial_Elastic_Truss2")
    files.append("C02-ST_Lin_2DRoof_Elastic_Truss2")
    files.append("C03-ST_Lin_3DAxial_Elastic_Truss2")
    files.append("C04-ST_Lin_3DAxial_Plastic_Truss2")
    files.append("C05-ST_Lin_3DPiramid_Elastic_Truss2")
    ####files.append("C06-DY_Lin_3DConstrained_Elastic_Truss2") #To Be Done!!
    ####files.append("C07-ST_Lin_3DCantilever_Elastic_Truss2") #To Be Done!!
    files.append("C08-ST_kin_2DCantilever_Elastic_Truss2")
    files.append("C09-ST_kin_3DCantilever_Elastic_Truss2")
    files.append("C11-ST_Lin_2DSurface_Elastic_Truss2")
    files.append("C12-ST_Lin_3DSurface_Elastic_Truss2")
    files.append("C15-ST_Lin_2DSurface_Elastic_Truss3")
    files.append("C16-ST_Lin_3DSurface_Elastic_Truss3")
    files.append("D01-ST_Lin_2DBernoulli_Elastic_Frame2")
    files.append("D02-ST_Lin_2DTimoshenko_Elastic_Frame2")
    files.append("D03-ST_Lin_3DBernoulli_Elastic_Frame2")
    files.append("D04-ST_Lin_3DTimoshenko_Elastic_Frame2")
    files.append("D05-ST_Lin_2DBernoulliArc_Elastic_Frame2")
    files.append("D06-ST_Lin_2DTimoshenkoArc_Elastic_Frame2")
    files.append("D07-ST_Lin_2DConstrainedBuilding_Elastic_Frame2")
    files.append("D08-ST_Lin_2DVolForce_Elastic_Frame2")
    files.append("D09-ST_Lin_3DVolForce_Elastic_Frame2")
    files.append("D10-ST_kin_2DPointLoad_Bernoulli_Elastic_Frame2")
    files.append("D11-ST_kin_2DMomentBernoulli_Elastic_Frame2")
    files.append("D12-ST_Lin_2DSurfaceHorizontal_Elastic_Frame2")
    files.append("D13-ST_Lin_3DSurfaceHorizontal_Elastic_Frame2")
    files.append("D16-DY_Free_Rectangular_3DPointLoad_Elastic_Frame2")
    files.append("D18-DY_Free_Circular_2DPointLoad_Elastic_Frame2")
    files.append("D19-DY_Damped_Angle_2DPointLoad_Elastic_Frame2")
    files.append("D20-DY_Free_Rectangular_BodyLoad_Elastic_Frame2")
    files.append("D21-DY_Damped_Rectangular_BodyLoad_Elastic_Frame2")
    files.append("D23-DY_Rectangular_SupportMotion_Elastic_Frame2")
    files.append("D17-DY_Damped_WideFlange_3DPointLoad_Elastic_Frame2")
    files.append("D22-ST_Rectangular_SupportMotion_Elastic_Frame2")
    files.append("E01-ST_Lin_2DPointLoad_Elastic_Quad4")
    files.append("E02-ST_Lin_2DPointLoad_Elastic_Quad8")
    files.append("E03-ST_Lin_2DSurfaceLoad_Elastic_Quad4")
    files.append("E04-ST_Lin_2DSurfaceLoad_Elastic_Quad8")
    files.append("E05-ST_Lin_2DRigidLink_Elastic_Frame2")
    files.append("F01-ST_Lin_2DPointLoad_ElasticPStrain_Quad4")
    files.append("F02-DY_Lin_2DPointLoad_ElasticPStrain_Quad4")
    files.append("F03-DY_Lin_2DPointLoad_J2PStrain_Quad4")
    files.append("F04-DY_Lin_2DPointLoad_BAPStrain_Quad4")
    files.append("F06-DY_Lin_2DSoilColumn_ElasticPStrain_Quad4")
    files.append("F07-DY_Lin_2DSoilColumn_J2PStrain_Quad4")
    files.append("F08-DY_Lin_2DSoilColumn_BAPStrain_Quad4")
    ####files.append("F10-ST_Lin_2DPointLoad_ElasticPStrain_Quad8") #To Be Done!!
    files.append("F11-DY_Lin_2DPMLSoilColumn_ElasticPStrain_Quad4")
    files.append("F14-ST_Lin_2DConstrainedSSI_Elastic_Quad4_Frame2")
    ####files.append("F16-DY_Lin_2DSoilColumn_ElasticPStrain_Quad8") #To Be Done!!
    files.append("H01-ST_Lin_3DinPlanePointLoad_ElasticPStress_Shell4")
    files.append("H02-ST_Lin_3DoutPlanePointLoad_ElasticPStress_Shell4")
    files.append("H03-ST_Lin_3DSlabPointLoad_ElasticPStress_Shell4")
    files.append("H04-ST_Lin_3DSlabBodyLoad_ElasticPStress_Shell4")
    files.append("H05-ST_Lin_3DBuildingDiaphragm_ElasticPStress_Frame2_Shell4")
    files.append("H06-ST_Lin_3DSlabSurfaceLoad_ElasticPStress_Shell4")
    files.append("H07-ST_Lin_3DCantileverSurfaceLoad_Shell4")
    files.append("H08-DY_Damped_3DPointLoad_Plate_Elastic_Shell4")
    files.append("I01-ST_Lin_3DPointLoad_Elastic_Hexa8")
    files.append("I02-ST_Lin_3DSurfaceLoad_Elastic_Hexa8")
    files.append("I03-ST_Lin_3DPointLoad_Elastic_Hexa20")
    files.append("I04-ST_Lin_3DSurfaceLoad_Elastic_Hexa20")
    files.append("I05-ST_Lin_3DBodyLoad_Elastic_Hexa8")
    files.append("I06-ST_Lin_3DBodyLoad_Elastic_Hexa20")
    files.append("J02-DY_Lin_3DPointLoad_Elastic_Hexa8")
    files.append("J05-DY_Lin_3DSoilColumn_Elastic_Hexa8")
    ####files.append("J06-DY_Lin_3DSoilColumn_Elastic_Hexa20") #To Be Done!!
    files.append("J12-DY_Axial_Load_Long_Rod_PML3D")

    #The Global LaTeX files to be Included.
    LaTeXFiles = []

    n = len(files)
    for k in range(n):
        LaTeXFiles.append(cwd + "/../01-Debugging/" + files[k] + "/LaTeX/LaTeXFile.tex")

    #Run all the validation cases.
    print('Running all the validation cases')
    RunValidationCases(files)

    #Generate the Validation LaTeX Files.
    print('Generating the Validation LaTeX Files')
    GenValidationResults(files)

    #Generates the LaTeX Main Document.
    print('Generating the LaTeX Main Document')
    MainTeX(LaTeXFiles, files)

    #Generates the Final Report in PDF.
    print('Generating the Final Report in PDF')
    CompileTeX()

    print('Process Completed Successfully!')

def RunValidationCases(files):
    """
    This function runs all provided debugging cases in files.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2021

    Parameters
    ----------
    files : list
        The list of files (names) to be run 

    Returns
    -------
    None
    """
    #The current working path:
    cwd = os.path.abspath(os.path.dirname(sys.argv[0]))

	#Number of Debugging Cases.
    n = len(files)
    for k in range(n):
        #Excecutes the Pre-Analysis and Generate Files.
        cmdline = "python3 " + cwd + "/../01-Debugging/" + files[k] + "/" + files[k] + ".py"
        subprocess.check_output(cmdline, shell=True)

        #Excecutes the Run-Analysis and Generate Files.
        cmdline = cwd + "/../../02-Run_Process/SeismoVLAB.exe -dir " + cwd + "/../01-Debugging/" + files[k] + "/Partition -file Debugging_" + files[k][0:3] + ".$.json"
        subprocess.check_output(cmdline, shell=True)

        #Removes the unnecessary Files.
        cmdline = "rm -r " + cwd + "/../01-Debugging/" + files[k] + "/Partition" 
        os.system(cmdline)

        #Removes the unnecessary Files.
        cmdline = "rm -r " + cwd + "/../01-Debugging/" + files[k] + "/Paraview" 
        os.system(cmdline)

def GenValidationResults(files):
    """
    This function generates comparison between SVL and other results.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2021

    Parameters
    ----------
    files : list
        The list of files (names) to be run 

    Returns
    -------
    None
    """
    #The current working path.
    cwd = os.path.abspath(os.path.dirname(sys.argv[0]))

	#Number of Debugging Cases.
    n = len(files)
    for k in range(n):
        #Obtains the debugging folder path.
        ChangePath = cwd + "/../01-Debugging/" + files[k]      

        #Enters the Debug's Case Folder.
        os.chdir(ChangePath)

        #Excecutes the Debugging Case:
        exec(open('./LaTeX/cmpResults.py').read())

        #Returns to the Main Folder.
        os.chdir(cwd)

        #Removes the unnecessary Files.
        cmdline = "rm -r " + cwd + "/../01-Debugging/" + files[k] + "/Solution"
        os.system(cmdline)

def MainTeX(LaTeXFiles, files):
    """
    This function writes the LaTeX main document.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2021

    Parameters
    ----------
    files : list
        The list of the python (.py) files to be run 
    LaTeXFiles : list
        The list of the latex (.tex) files to be run 

    Returns
    -------
    Verification.tex file located at: 03-Validations/03-Report/LaTeX
    """
	#The document's abstract:
    abstract = """\\textbf{Seismo-VLAB} (\\textbf{SVL}) is a simple, fast, and extendable C++ multi-platform finite 
    element software designed to run large-scale simulations of dynamic, nonlinear soil-structure interaction problems. 
    In this report several debugging cases are presented in order to verify the accuracy and well-behaviour of the 
    implemented features. The DEBUG CASE's names are as follows: 
    \\textbf{L01-Analysis\_Formulation\_Comment\_Material\_Element}, where \\textbf{L} 
    is a letter that denotes complexity, \\textbf{Analysis} can be ST=Static, or DY=Dynamic, \\textbf{Formulation} can 
    be Lin=Linearized or Kin=kinematics, \\textbf{Comment} is a description, \\textbf{Material} and \\textbf{Element} are 
    the \\textbf{SVL}'s class.\n"""

    #Create the LateX File:
    LaTeXfile = open("./LaTeX/Verification.tex", "w+")

    LaTeXfile.write("\\documentclass[fleqn]{article}\n")
    LaTeXfile.write("\\usepackage{graphicx}\n")
    LaTeXfile.write("\\usepackage{subfigure}\n")
    LaTeXfile.write("\\usepackage{hyperref}\n")
    LaTeXfile.write("\\usepackage{amsthm}\n")
    LaTeXfile.write("\\usepackage{color}\n")
    LaTeXfile.write("\\usepackage{float}\n")
    LaTeXfile.write("\\usepackage{caption}\n")
    LaTeXfile.write("\\usepackage{array}\n")
    LaTeXfile.write("\\usepackage{multirow}\n")
    LaTeXfile.write("\\usepackage[left=2.00cm,right=2.00cm, bottom=2.0cm, top=1.5cm]{geometry} \n")
    LaTeXfile.write("\\usepackage{amsmath,amsfonts,amssymb, bbm}\n")
    LaTeXfile.write("\n")
    LaTeXfile.write("\\usepackage{tikz}\n")
    LaTeXfile.write("\\newcommand{\\redline}{\\raisebox{2pt}{\\tikz{\draw[-,red, line width = 1pt](0,0) -- (7mm,0);}}}\n")
    LaTeXfile.write("\n")
    LaTeXfile.write("\\title{\\bf Debugging Report for Verification Test Cases of SeismoVLAB} \n")
    LaTeXfile.write("\\author{Danilo S. Kusanovic, Elnaz E. Seylabi, and Domniki Asimaki.} \n")
    LaTeXfile.write("\n")
    LaTeXfile.write("\\begin{document}\n")
    LaTeXfile.write("\n")
    LaTeXfile.write("\\maketitle \n")
    LaTeXfile.write("\n")
    LaTeXfile.write("\\begin{abstract}\n")
    LaTeXfile.write(abstract)
    LaTeXfile.write("\\end{abstract}\n")
    LaTeXfile.write("\n")
    LaTeXfile.write("\n")

    #Writes the Debugging cases results.
    n = len(LaTeXFiles)
    for k in range(n):
        LaTeXfile.write("\\subsection*{DEBUG CASE : " + files[k].replace("_", "\_") + "} \n")
        LaTeXfile.write("\\input{" + LaTeXFiles[k] + "} \n") 

    LaTeXfile.write("\n")
    LaTeXfile.write("\\end{document}")

    LaTeXfile.close()

def CompileTeX():
    """
    This fucntion compiles the LaTeX main document.\n
    @visit  https://github.com/SeismoVLAB/SVL\n
    @author Danilo S. Kusanovic 2021

    Parameters
    ----------
    None 

    Returns
    -------
    PDF file located at: 03-Validations/03-Report/LaTeX
    """
    #Get the current path.
    cwd = os.getcwd() + "/"

    #Obtains the debugging folder path.
    ChangePath = cwd + "LaTeX"
    
    #Enters the Debug's Case Folder.
    os.chdir(ChangePath)

    #Enters the Debug's Case Folder.
    s = subprocess.check_output("pdflatex Verification.tex", shell=True)

    #Returns to the Main Folder.
    os.chdir(cwd)

if __name__ == "__main__":
    main()
