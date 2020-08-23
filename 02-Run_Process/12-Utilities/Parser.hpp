//==============================================================================
//
//                    (Seismo)s (V)irtual (Lab)oratory
//             Module for Serial and Parallel Analysis of seismic 
//         wave propagation and soil-structure interaction simulation
//         Copyright (C) 2018, The California Institute of Technology
//                           All Rights Reserved.
//
// Commercial use of this program without express permission of the California
// Institute of Technology, is strictly  prohibited. See  file "COPYRIGHT"  in
// main  directory  for  information on  usage  and  redistribution, and for a
// DISCLAIMER OF ALL WARRANTIES.
//
//==============================================================================
//
// Written by:
//   Danilo S. Kusanovic (dkusanov@caltech.edu)
//   Elnaz E. Seylabi    (elnaze@unr.edu)
//
// Supervised by:
//   Domniki M. Asimaki  (domniki@caltech.edu)
//
// References : 
//   [1]
//
// Description:
///This is a parallel parser object that takes the program's format and 
///generates the Mesh and Analysis objects.
//------------------------------------------------------------------------------

#ifndef _PARSER_HPP_
#define _PARSER_HPP_

#include <map>
#include <vector>
#include <string>
#include <memory>
#include <iostream>
#include <Eigen/Dense>

#include "Node.hpp"
#include "Load.hpp"
#include "Material.hpp"
#include "Section.hpp"
#include "Element.hpp"
#include "Damping.hpp"
#include "Constraint.hpp"
#include "LoadCombo.hpp"
#include "Recorder.hpp"

#include "Mesh.hpp"
#include "Analysis.hpp"

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      December 18, 2018
/// @version   1.0
/// @file      Parser.hpp
/// @class     Parser
/// @see       Mesh.hpp Analysis.hpp
/// @brief     Class that translate an input file into Seismo-VLAB classes
class Parser {

    public:
        ///Creates a Parser object to load information from file.
        ///@param np The number of processors.
        ///@param folder Folder address where the files will be loaded.
        ///@param file The file name where Mesh components are loaded.
        ///@param execute The file name where Analysis components are loaded.
        ///@note More details can be found at @ref linkParser.
        ///@see Parser::Folder, Parser::File, Parser::Execute.
        Parser(unsigned int np, std::string folder, std::string file);

        ///Destroys this Parser object.
        ~Parser();

        ///Parse the User's Input Model.
        ///@param theMesh Pointer to the Mesh where objects are loaded.
        ///@param theAnalysis Pointer to the analysis case to be simulated.
        ///@return Whether or not the information was successfully loaded.
        ///@see Mesh, Analysis.
        bool GetFromFile(std::shared_ptr<Mesh> &theMesh, std::unique_ptr<Analysis> &theAnalysis);

    protected:
        ///Parse the Input Mesh File.
        ///@param theMesh Pointer to the Mesh where objects are loaded.
        ///@param theAnalysis Pointer to the analysis case to be simulated.
        ///@param theFile The file name where Mesh will be populated.
        ///@return Whether or not the Mesh information was successfully loaded.
        ///@note More details can be found at @ref linkParser.
        ///@see Mesh, Element, Node, Material, Section, Load, Constraint.
        bool ParseSVLFile(std::shared_ptr<Mesh> &theMesh, std::unique_ptr<Analysis> &theAnalysis, std::string theFile);

        ///Fix blank spaces provided by user in path.
        ///@param theFile String with full path.
        ///@param toReplace The string pattern to be replaced with. 
        ///@return String with the replaced pattern.
        std::string GetSpacedName(std::string theFile, std::string toReplace);

        ///Creates the Subdomain File Name.
        ///@param file String with the file name to be modified.
        ///@param k The process number to be replaced in file. 
        ///@param cond Whether or not to consider the full path. 
        ///@return String with the modified pattern.
        std::string GetPartitionName(std::string file, int k, bool cond);

    private:
        ///Mass Formulation for elements.
        bool MassForm;

        ///Include file flag.
        bool Flag;

        ///Partition Variables.
        int nProcessors;

        ///Model's Dimension.
        unsigned int nDimensions;

        ///Folder where the file is loaded.
        std::string Folder; 

        ///File name to be loaded.
        std::string File;

        ///Vector of Recorders to store solution. 
        std::vector<std::shared_ptr<Recorder> > Recorders;

        ///Load combination cases.
        std::map<unsigned int, std::shared_ptr<LoadCombo> > LoadCombos;

        ///Obtains a Point from Input File.
        ///@param InputFile The input file where Node objects are loaded.
        ///@param theNode Pointer to the Node objects where information is stored.
        ///@see Node, Mesh.
        unsigned int CreatePoint(std::ifstream& InputFile, std::shared_ptr<Node> &theNode);

        ///Obtains a Material from Input File.
        ///@param InputFile The input file where objects are loaded.
        ///@param theMaterial Pointer to the Material objects where information is stored.
        ///@see Material, Mesh.
        unsigned int CreateMaterial(std::ifstream& InputFile, std::unique_ptr<Material> &theMaterial);

        ///Obtains a Section from Input File.
        ///@param InputFile The input file where the section objects is loaded.
        ///@param theMesh Pointer to the Mesh where section objects are loaded.
        ///@param theSection Pointer to the section object where information is stored.
        unsigned int CreateSection(std::ifstream& InputFile, std::shared_ptr<Mesh> &theMesh, std::unique_ptr<Section> &theSection);

        ///Obtains a Constraint from Input File.
        ///@param InputFile The input file where Constraint objects is loaded.
        ///@param theConstraint Pointer to the constraint where information is stored.
        ///@see Constraint.
        int CreateConstraint(std::ifstream& InputFile, std::unique_ptr<Constraint> &theConstraint);

        ///Obtains an Element from Input File.
        ///@param InputFile The input file where the Element objects is loaded.
        ///@param theMesh Pointer to the Mesh where objects are loaded.
        ///@param theElement Pointer to the Element where objects is generated.
        ///@see Mesh, Element.
        unsigned int CreateElement(std::ifstream& InputFile, std::shared_ptr<Mesh> &theMesh, std::shared_ptr<Element> &theElement);

        ///Obtains the Support Motion from Input File.
        ///@param InputFile The input file where the support motion is loaded.
        ///@param Xo The vector with the support motion values in time.
        ///@param dof The degree of freedom where the support motion is applied.
        ///@see Node::SupportMotion, Node::SetSupportMotion().
        unsigned int CreateSupportMotion(std::ifstream& InputFile, std::vector<double>& Xo, unsigned int &dof);

        ///Obtains a Point Load from Input File.
        ///@param InputFile The input file where the Load objects are loaded.
        ///@param theLoad Pointer to the Load object where information is stored. 
        ///@see Node, Element, Load.
        unsigned int CreatePointLoad(std::ifstream& InputFile, std::shared_ptr<Load> &theLoad);

        ///Obtains an Element Load from Input File.
        ///@param InputFile The input file where element load objects is loaded.
        ///@param theMesh Pointer to the Mesh where objects are loaded.
        ///@param theLoad Pointer to the Load object where information is stored.
        ///@see Element, Load, Mesh.
        unsigned int CreateElementLoad(std::ifstream& InputFile, std::shared_ptr<Mesh> &theMesh, std::shared_ptr<Load> &theLoad);

        ///Obtains a Damping objec from Input File.
        ///@param InputFile The input file where Damping objects are loaded.
        ///@param eList The list of \@ Element that share this damping model.
        ///@param theDamping Pointer to the damping object where information is loaded.
        ///@see Element, Damping, Mesh.
        unsigned int CreateDamping(std::ifstream& InputFile, std::vector<unsigned int>& eList, std::shared_ptr<Damping> &theDamping);

        ///Obtains a Combination object from Input File.
        ///@param InputFile The input file where objects are loaded.
        ///@param theCombo The load combination map to be performed by each analysis.
        ///@see LoadCombo.
        unsigned int CreateCombination(std::ifstream& InputFile, std::shared_ptr<LoadCombo> &theCombo);

        ///Obtains a Recorder object from Input File.
        ///@param InputFile The input file where objects are loaded.
        ///@param theRecorder Pointer to the recorder where solution will be stored.
        ///@see Analysis, Recorder
        void CreateRecorder(std::ifstream& InputFile, std::shared_ptr<Recorder> &theRecorder);

        ///Obtains a Analysis object from Input File.
        ///@param InputFile The input file where objects are loaded.
        ///@param theMesh Pointer to the Mesh where objects are loaded.
        ///@param theAnalysis Pointer to the analysis case to be simulated.
        ///@param theCombo The load combination map to be performed by each analysis.
        ///@see Mesh, Analysis, LoadCombo
        void CreateAnalysis(std::ifstream& InputFile, std::shared_ptr<Mesh> &theMesh, std::unique_ptr<Analysis> &theAnalysis, std::map<unsigned int, std::shared_ptr<LoadCombo> > &theCombo);
};

#endif
