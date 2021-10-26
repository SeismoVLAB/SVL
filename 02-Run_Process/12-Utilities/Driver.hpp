//==============================================================================
//
//                       Seismo Virtual Laboratory
//             Module for Serial and Parallel Analysis of seismic 
//         wave propagation and soil-structure interaction simulation
//         Copyright (C) 2018-2021, The California Institute of Technology
//                         All Rights Reserved.
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
//  [1] 
//
// Description:
//This file contains the "parsing" functions to create, modify and 
//analyze a finite element model provided in JSON format.
//------------------------------------------------------------------------------

#ifndef _DRIVER_HPP_
#define _DRIVER_HPP_

#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <fstream>
#include <iostream>
#include <algorithm>

#include "Viscous1DLinear.hpp"
#include "Elastic1DLinear.hpp"
#include "Hertzian1DLinear.hpp"
#include "Elastic2DPlaneStress.hpp"
#include "Elastic2DPlaneStrain.hpp"
#include "PlasticPlaneStrainJ2.hpp"
#include "PlasticPlaneStrainBA.hpp"
#include "Elastic3DLinear.hpp"
#include "Plastic1DJ2.hpp"
#include "Plastic3DJ2.hpp"
#include "Plastic3DBA.hpp"

#include "Elastic1DGap.hpp"
#include "Plastic1DGap.hpp"
#include "Elastic1DFiber.hpp"
#include "Steel1DFiber.hpp"
#include "Concrete1DFiber.hpp"

#include "Lin2DRectangular.hpp"
#include "Lin3DRectangular.hpp"
#include "Lin2DAngle.hpp"
#include "Lin3DAngle.hpp"
#include "Lin2DChannel.hpp"
#include "Lin3DChannel.hpp"
#include "Lin2DTee.hpp"
#include "Lin3DTee.hpp"
#include "Lin2DWideFlange.hpp"
#include "Lin3DWideFlange.hpp"
#include "Lin2DCircular.hpp"
#include "Lin3DCircular.hpp"
#include "Lin2DRectangularTube.hpp"
#include "Lin3DRectangularTube.hpp"
#include "Lin2DCircularTube.hpp"
#include "Lin3DCircularTube.hpp"
#include "Lin3DThinArea.hpp"
#include "Lin2DUserDefined.hpp"
#include "Lin3DUserDefined.hpp"

#include "Fib3DLineSection.hpp"

#include "ZeroLength1D.hpp"
#include "lin2DTruss2.hpp"
#include "kin2DTruss2.hpp"
#include "lin3DTruss2.hpp"
#include "kin3DTruss2.hpp"
#include "lin2DTruss3.hpp"
#include "lin3DTruss3.hpp"
#include "lin2DFrame2.hpp"
#include "kin2DFrame2.hpp"
#include "lin3DFrame2.hpp"
#include "lin2DTria3.hpp"
#include "lin2DTria6.hpp"
#include "lin2DQuad4.hpp"
#include "lin2DQuad8.hpp"
#include "PML2DQuad4.hpp"
#include "kin2DQuad4.hpp"
#include "PML2DQuad8.hpp"
#include "lin3DShell4.hpp"
#include "lin3DTetra4.hpp"
#include "lin3DTetra10.hpp"
#include "lin3DHexa8.hpp"
#include "kin3DHexa8.hpp"
#include "PML3DHexa8.hpp"
#include "lin3DHexa20.hpp"
#include "EQlin2DQuad4.hpp"
#include "TIEQlin2DQuad4.hpp"
#include "UnxBoucWen2DLink.hpp"
#include "UnxBoucWen3DLink.hpp"
#include "HDRBYamamoto2DLink.hpp"
#include "HDRBYamamoto3DLink.hpp"
#include "null2DFrame2.hpp"
#include "null3DFrame2.hpp"

#include "StaticAnalysis.hpp"
#include "DynamicAnalysis.hpp"

#include "Linear.hpp"
#include "NewtonRaphson.hpp"
//#include "GeneralDisplacementPath.hpp"

#include "QuasiStatic.hpp"
#include "NewmarkBeta.hpp"
#include "CompositeBathe.hpp"
#include "CentralDifference.hpp"
#include "ExtendedNewmarkBeta.hpp"

#include "EigenSolver.hpp"
#include "MumpsSolver.hpp"
#include "PetscSolver.hpp"

#include "RSJparser.hpp"
#include "Definitions.hpp"
#include "Profiler.hpp"

///This file contains the "parsing" functions to create, modify and analyze a finite element model provided in JSON format.
/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      June 18, 2021
/// @version   1.0
/// @file      Driver.hpp
/// @see       main.cpp

///Sorts a vector of tag number form smaller to larger
///@param v vector with the data to be sorted
template<typename T> 
void 
sortVector(std::vector<T>& v){
    std::sort(v.begin(), v.end());
}

///Performs union, difference or intersection operation between two vectors.
///@param v1 vector with the identifiers in json file
///@param v1 vector with the identifiers in Mesh object
///@param op operation to be performed between these two sets
///@return A vector with the Difference, Intersection, or Union of the two sets.
template<typename T> 
std::vector<T>
setOperation(std::vector<T>& v1, std::vector<T>& v2, std::string op){
    //The resulting vector
    std::vector<T> v;

    //The possible options for set operations   
    if( strcasecmp(op.c_str(),"DIFFERENCE") == 0) {
        std::set_difference(v1.begin(), v1.end(), v2.begin(), v2.end(), std::inserter(v, v.begin())); 
    }
    else if(strcasecmp(op.c_str(),"INTERSECTION") == 0){
        std::set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(v));
    }
    else if(strcasecmp(op.c_str(),"UNION") == 0){
        std::set_union(v1.begin(), v1.end(), v2.begin(), v2.end(), std::back_inserter(v));
    } 

    //Removes possible repeated values in result vector
    if(v.size() != 0){
        auto last = std::unique(v.begin(), v.end());
        v.erase(last, v.end());
    }

    return v;
}

///Query the JSON file and obtains the Identifier list
///@param jsonFile json file where mesh entities will be readden.
///@param Name name of the entity to be readden from json file.
///@return A vector with the identifiers in the json file.
template<typename T> 
std::vector<T>
GetIDsFromJSON(RSJresource& jsonFile, std::string Name){
    //Number of Entities in Property:
    unsigned int num = jsonFile[Name].as_object().size();
    std::vector<T> Tags(num);

    //Goes Over the identifiers of the JSON Property 
    unsigned int k = 0;
    for(auto it = jsonFile[Name].as_object().begin(); it != jsonFile[Name].as_object().end(); ++it){
        //Property identifier
        unsigned int Tag = std::stoi(it->first);
        Tags[k] = Tag;
        k++;
    }

    //Sorts the identifier list
    sortVector<T>(Tags);

    return Tags;
}

///Query the MESH object and obtains the Identifier list
///@param mesh Pointer to the Mesh container.
///@param Name name of the entity to be readden from json file.
///@return A vector with the identifiers in Mesh.
template<typename T> 
std::vector<T>
GetIDsFromMESH(std::shared_ptr<Mesh>& mesh, std::string Name){
    //The resulting ID vector
    std::vector<T> Tags = mesh->GetVectorIDs<T>(Name);

    //Sorts the identifier list
    sortVector<T>(Tags);

    return Tags;
}

///The Entities indexes to be updated
///@param mesh Pointer to the Mesh container.
///@param jsonFile json file where mesh entities will be readden.
///@param Name name of the entity to be readden from json file.
///@return A map with the added ("add"), deleted ("del") or modified ("mod") entities with respect to previous mesh.
template<typename T> 
std::map<std::string, std::vector<T> >
Entities2Update(std::shared_ptr<Mesh>& mesh, RSJresource& jsonFile, std::string Name){
    //Entity Objects To Update
    std::vector<T> v1 = GetIDsFromJSON<T>(jsonFile, Name);
    std::vector<T> v2 = GetIDsFromMESH<T>(mesh, Name);

    //Map with indexes for which Node Objects needs to be added, deleted, or modified
    std::map<std::string, std::vector<T> > Tags;
    Tags["add"] = setOperation<T>(v1,v2,"Difference");
    Tags["del"] = setOperation<T>(v2,v1,"Difference");
    Tags["mod"] = setOperation<T>(v1,v2,"Intersection");

    return Tags;
}

///Sets the partition subdomain tag number.
///@param theFile file that contains a character to be replaced.
///@param k the number to be replaced with.
///@param cond condition to return full path of the file or just the file name.
///@return the modified string where the preplacement is performed.
std::string 
GetPartitionName(std::string theFile, int k, bool cond){
    //Auxiliary variables.
    std::string auxName = theFile;
    std::string subDomainMesh;

    //Conver process number into string.
    std::stringstream process;
    process << k;

    //Replace strings.
    size_t pos;
    pos = auxName.find("$");
    if(pos !=std::string::npos)
        auxName.replace(pos, std::string("$").length(), process.str());

    //Modified File Name.
    if(cond){
        subDomainMesh = filePath + "/" + auxName;
    }
    else{
        subDomainMesh = auxName;
    }

    return subDomainMesh;
}

///Fix blank spaces provided by user in path (if any).
///@param theFile file that contains a character to be replaced.
///@param toReplace the character to be replaced with.
///@return the modified string where the preplacement is performed.
std::string 
GetSpacedName(std::string theFile, std::string toReplace){
    //Auxiliary variables.
    std::string auxName = theFile;
    std::string subDomainMesh;

    //Replace strings.
    size_t pos;
    pos = auxName.find("~");
    while(pos != std::string::npos){
        auxName.replace(pos, std::string("~").length(), toReplace);
        pos = auxName.find("~");
    }

    //Modified File Name.
    subDomainMesh = auxName;

    return subDomainMesh;
}

///Updates the Node Entities in Mesh Object 
///@param theMesh Pointer to the Mesh container.
///@param jsonFile json file where mesh entities will be readden.
void 
UpdateNodes(std::shared_ptr<Mesh> &theMesh, RSJresource& jsonFile){
    //Node Identifiers To Be Updated
    std::map<std::string, std::vector<unsigned int> > Tags = Entities2Update<unsigned int>(theMesh, jsonFile, "Nodes");

    //Nodes to be Removed from Mesh Object
    for(unsigned int k = 0; k < Tags["del"].size(); k++){
        unsigned int Tag = Tags["del"][k];
        theMesh->DelNode(Tag);
    }

    //Nodes to be Added To Mesh Object
    for(unsigned int k = 0; k < Tags["add"].size(); k++){
        unsigned int Tag = Tags["add"][k];

        //Node identifier and number of degree-of-freedom:
        std::string nTag = std::to_string(Tag);
        unsigned int nDofs = jsonFile["Nodes"][nTag]["ndof"].as<int>();

        //Allocate memory for Free/Total degree of freedom lists.
        std::vector<int> freeDofs(nDofs,0);
        std::vector<int> totalDofs(nDofs,0);

        //Allocate memory for the coordinate vector.
        Eigen::VectorXd coordinates(nDimensions);

        for(int k = 0; k < jsonFile["Nodes"][nTag]["freedof"].size(); ++k)
            freeDofs[k] = jsonFile["Nodes"][nTag]["freedof"][k].as<int>();
            
        for(int k = 0; k < jsonFile["Nodes"][nTag]["totaldof"].size(); ++k)
            totalDofs[k] = jsonFile["Nodes"][nTag]["totaldof"][k].as<int>();

        for(int k = 0; k < jsonFile["Nodes"][nTag]["coords"].size(); ++k)
            coordinates(k) = jsonFile["Nodes"][nTag]["coords"][k].as<double>();

        //Checks if the node is free/fixed.
        bool IsFixed = false;
        if(std::find(freeDofs.begin(), freeDofs.end(), -1) != freeDofs.end())
            IsFixed = true;

        //Creates a node object.
        std::shared_ptr<Node> theNode = std::make_shared<Node>(nDofs, coordinates, IsFixed);

        //Sets preliminary degree-of-freedom numbering.
        theNode->SetFreeDegreeOfFreedom(freeDofs);
        theNode->SetTotalDegreeOfFreedom(totalDofs);  

        //Stores the information in node container.
        theMesh->AddNode(Tag, theNode);
    }

    //Nodes to be Modified in Mesh Object
    if( Tags["mod"].size() > 0){
        std::map<unsigned int, std::shared_ptr<Node> > theNodes = theMesh->GetNodes();

        for(unsigned int k = 0; k < Tags["mod"].size(); k++){
            unsigned int Tag = Tags["mod"][k];

            //Node identifier and number of degree-of-freedom:
            std::string nTag = std::to_string(Tag);
            unsigned int nDofs = jsonFile["Nodes"][nTag]["ndof"].as<int>();

            //Allocate memory and get information from JSON file
            std::vector<int> freeDofs(nDofs,0);
            for(int k = 0; k < jsonFile["Nodes"][nTag]["freedof"].size(); ++k)
                freeDofs[k] = jsonFile["Nodes"][nTag]["freedof"][k].as<int>();
                
            std::vector<int> totalDofs(nDofs,0);
            for(int k = 0; k < jsonFile["Nodes"][nTag]["totaldof"].size(); ++k)
                totalDofs[k] = jsonFile["Nodes"][nTag]["totaldof"][k].as<int>();

            Eigen::VectorXd coordinates(nDimensions);
            for(int k = 0; k < jsonFile["Nodes"][nTag]["coords"].size(); ++k)
                coordinates(k) = jsonFile["Nodes"][nTag]["coords"][k].as<double>();

            //Checks if the node is free/fixed.
            bool IsFixed = false;
            if(std::find(freeDofs.begin(), freeDofs.end(), -1) != freeDofs.end())
                IsFixed = true;

            //Modify the Node's properties according to JSON file.
            theNodes[Tag]->SetAsFixed(IsFixed);
            theNodes[Tag]->SetCoordinates(coordinates);
            theNodes[Tag]->SetFreeDegreeOfFreedom(freeDofs);
            theNodes[Tag]->SetTotalDegreeOfFreedom(totalDofs);
        }
    }
}

///Updates the Masses in Mesh Object 
///@param theMesh Pointer to the Mesh container.
///@param jsonFile json file where mesh entities will be readden.
void 
UpdateMasses(std::shared_ptr<Mesh> &theMesh, RSJresource& jsonFile){
    //Mass Identifiers To Be Updated
    std::map<std::string, std::vector<unsigned int> > Tags = Entities2Update<unsigned int>(theMesh, jsonFile, "Masses");

    //Masses to be Removed from Mesh Object
    for(unsigned int k = 0; k < Tags["del"].size(); k++){
        unsigned int Tag = Tags["del"][k];
        theMesh->DelMass(Tag);
    }

    //Masses to be Added To Mesh Object
    for(unsigned int k = 0; k < Tags["add"].size(); k++){
        //Mass identifier and number of degree-of-freedom:
        unsigned int Tag = Tags["add"][k];
        std::string nTag = std::to_string(Tag);
        unsigned int nDofs = jsonFile["Masses"][nTag]["ndof"].as<int>();

        Eigen::VectorXd Mass(nDofs);
        for(unsigned int k = 0; k < nDofs; ++k)
            Mass(k) = jsonFile["Masses"][nTag]["mass"][k].as<double>();

        //Adds the nodal mass to the node.
        theMesh->AddMass(Tag, Mass);
    }

    //Masses to be Modified in Mesh Object 
    for(unsigned int k = 0; k < Tags["mod"].size(); k++){
        unsigned int Tag = Tags["mod"][k];
        std::string nTag = std::to_string(Tag);
        unsigned int nDofs = jsonFile["Masses"][nTag]["ndof"].as<int>();

        Eigen::VectorXd Mass(nDofs);
        for(unsigned int k = 0; k < nDofs; ++k)
            Mass(k) = jsonFile["Masses"][nTag]["mass"][k].as<double>();

        //Sets (replace) the new mass vector.
        theMesh->AddMass(Tag, Mass);
    }
}

///Updates the Constraint Entities in Mesh Object 
///@param theMesh Pointer to the Mesh container.
///@param jsonFile json file where mesh entities will be readden.
void 
UpdateConstraints(std::shared_ptr<Mesh> &theMesh, RSJresource& jsonFile){
    //Constraint Identifiers To Be Updated
    std::map<std::string, std::vector<int> > Tags = Entities2Update<int>(theMesh, jsonFile, "Constraints");

    //Constraint to be Removed from Mesh Object
    for(unsigned int k = 0; k < Tags["del"].size(); k++){
        int Tag = Tags["del"][k];
        theMesh->DelConstraint(Tag);
    }

    //Constraint to be Added To Mesh Object
    for(unsigned int k = 0; k < Tags["add"].size(); k++){
        //Node identifier and number of degree-of-freedom:
        int Tag = Tags["add"][k];
        std::string cTag = std::to_string(Tag);

        //Constraint and slave identifiers 
        unsigned int stag = jsonFile["Constraints"][cTag]["stag"].as<int>();
        unsigned int nCombs = jsonFile["Constraints"][cTag]["mtag"].size();

        //The master identifiers and combinational factors
        std::vector<double> factors(nCombs);
        std::vector<unsigned int> mtag(nCombs);

        for(unsigned int k = 0; k < nCombs; ++k){
            mtag[k] = jsonFile["Constraints"][cTag]["mtag"][k].as<int>();
            factors[k] = jsonFile["Constraints"][cTag]["factor"][k].as<double>();
        }

        //Creates a constaint object.
        std::shared_ptr<Constraint> theConstraint = std::make_shared<Constraint>(stag, mtag, factors);

        //Stores the information in constraint container.
        theMesh->AddConstraint(Tag, theConstraint);
    }

    //Constraint to be Modified in Mesh Object
    if( Tags["mod"].size() > 0){
        std::map<int, std::shared_ptr<Constraint> > theConstraints = theMesh->GetConstraints();

        for(unsigned int k = 0; k < Tags["mod"].size(); k++){
            int Tag = Tags["mod"][k];
            std::string cTag = std::to_string(Tag);

            //Constraint and slave identifiers 
            unsigned int stag = jsonFile["Constraints"][cTag]["stag"].as<int>();
            unsigned int nCombs = jsonFile["Constraints"][cTag]["mtag"].size();

            //The master identifiers and combinational factors
            std::vector<double> factors(nCombs);
            std::vector<unsigned int> mtag(nCombs);

            for(unsigned int k = 0; k < nCombs; ++k){
                mtag[k] = jsonFile["Constraints"][cTag]["mtag"][k].as<int>();
                factors[k] = jsonFile["Constraints"][cTag]["factor"][k].as<double>();
            }

            //Modify the Constraint's properties according to JSON file.
            theConstraints[Tag]->SetSlaveInformation(stag);
            theConstraints[Tag]->SetMasterInformation(mtag);
            theConstraints[Tag]->SetCombinationFactors(factors);
        }
    }
}

///Updates the Support Motions in Mesh Object 
///@param theMesh Pointer to the Mesh container.
///@param jsonFile json file where mesh entities will be readden.
void 
UpdateSupportMotion(std::shared_ptr<Mesh> &theMesh, RSJresource& jsonFile){
    //Support Motions Identifiers To Be Updated
    std::map<std::string, std::vector<int> > Tags = Entities2Update<int>(theMesh, jsonFile, "Supports");

    //Support Motions to be Removed from Mesh Object
    for(unsigned int k = 0; k < Tags["del"].size(); k++){
        unsigned int Tag = Tags["del"][k];
        theMesh->DelSupportMotion(Tag);
    }

    //Support Motions to be Added To Mesh Object
    for(unsigned int k = 0; k < Tags["add"].size(); k++){
        //Support motion name and identifier.
        unsigned int Tag = Tags["add"][k];
        std::string nTag = std::to_string(Tag);
        std::string Name = jsonFile["Supports"][nTag]["type"].as<std::string>();

        for(int n = 0; n < jsonFile["Supports"][nTag]["dof"].size(); ++n){
            unsigned int dof;
            std::vector<double> Xo;

            if(strcasecmp(Name.c_str(),"CONSTANT") == 0){
                Xo.resize(1);
                Xo[0] = jsonFile["Supports"][nTag]["value"][n].as<double>();
            }
            else if(strcasecmp(Name.c_str(),"TIMESERIES") == 0){
                std::string  File = jsonFile["Supports"][nTag]["file"][n].as<std::string>();
                //Loads the time history values into memory.
                std::ifstream motion(File.c_str());

                //The File is Opened and Ready to be Loaded.
                if(motion.is_open()){
                    //Number of time steps.
                    unsigned int nt;
                    motion >> nt;

                    //Time-history load vector.
                    Xo.resize(nt);
                    for(unsigned int j = 0; j < nt; j++)
                        motion >> Xo[j];
                }

                motion.close();
            }
            dof = jsonFile["Supports"][nTag]["dof"][n].as<int>();

            //Assign support motion to the node.
            theMesh->SetSupportMotion(Tag, dof, Xo);
        }
    }

    //No Support Motion Modification. Support should be erased from previous analysis
}

///Updates the Material Entities in Mesh Object 
///@param theMesh Pointer to the Mesh container.
///@param jsonFile json file where mesh entities will be readden.
void 
UpdateMaterials(std::shared_ptr<Mesh> &theMesh, RSJresource& jsonFile){
    //List of Materials already defined
    std::vector<unsigned int> Tags = theMesh->GetVectorIDs<unsigned int>("Materials");

    //Materials to be Added To Mesh Object
    for(auto it = jsonFile["Materials"].as_object().begin(); it != jsonFile["Materials"].as_object().end(); ++it){
        //Material name and identifier.
        unsigned int Tag = std::stoi(it->first);

        //Whether the material should be created 
        if (std::find(Tags.begin(), Tags.end(), Tag) == Tags.end()) {
            //Material name and identifier.
            unsigned int Tag = std::stoi(it->first);
            std::string Name = it->second["name"].as<std::string>();

            //Creates a material object.
            std::unique_ptr<Material> theMaterial;

            if(strcasecmp(Name.c_str(),"Elastic1DLinear") == 0){
                double E = it->second["attributes"]["E"].as<double>(0.0);
                double nu = it->second["attributes"]["nu"].as<double>(0.0);
                double rho = it->second["attributes"]["rho"].as<double>(0.0);

                //Instantiate the material object.
                theMaterial = std::make_unique<Elastic1DLinear>(E, nu, rho);
            }
            else if(strcasecmp(Name.c_str(),"Hertzian1DLinear") == 0){
                double k1 = it->second["attributes"]["k1"].as<double>(0.0);
                double k2 = it->second["attributes"]["k2"].as<double>(0.0);
                double k3 = it->second["attributes"]["k3"].as<double>(0.0);
                double rho = it->second["attributes"]["rho"].as<double>(0.0);

                //Instantiate the material object.
                theMaterial = std::make_unique<Hertzian1DLinear>(k1, k2, k3, rho);
            }
            else if(strcasecmp(Name.c_str(),"Viscous1DLinear") == 0){
                double eta = it->second["attributes"]["eta"].as<double>(0.0);

                //Instantiate the material object.
                theMaterial =std::make_unique<Viscous1DLinear>(eta);
            }
            else if(strcasecmp(Name.c_str(),"Elastic2DPlaneStrain") == 0){
                double E = it->second["attributes"]["E"].as<double>(0.0);
                double nu = it->second["attributes"]["nu"].as<double>(0.0);
                double rho = it->second["attributes"]["rho"].as<double>(0.0);

                //Instantiate the material object.
                theMaterial = std::make_unique<Elastic2DPlaneStrain>(E, nu, rho);
            }
            else if(strcasecmp(Name.c_str(),"Elastic2DPlaneStress") == 0){
                double E = it->second["attributes"]["E"].as<double>(0.0);
                double nu = it->second["attributes"]["nu"].as<double>(0.0);
                double rho = it->second["attributes"]["rho"].as<double>(0.0);

                //Instantiate the material object.
                theMaterial = std::make_unique<Elastic2DPlaneStress>(E, nu, rho);
            }
            else if(strcasecmp(Name.c_str(),"Elastic3DLinear") == 0){
                double E = it->second["attributes"]["E"].as<double>(0.0);
                double nu = it->second["attributes"]["nu"].as<double>(0.0);
                double rho = it->second["attributes"]["rho"].as<double>(0.0);

                //Instantiate the material object.
                theMaterial = std::make_unique<Elastic3DLinear>(E, nu, rho);
            }
            else if(strcasecmp(Name.c_str(),"Plastic1DJ2") == 0){
                double E = it->second["attributes"]["E"].as<double>(0.0);
                double nu = it->second["attributes"]["nu"].as<double>(0.0);
                double rho = it->second["attributes"]["rho"].as<double>(0.0);
                double K = it->second["attributes"]["k"].as<double>(0.0);
                double H = it->second["attributes"]["h"].as<double>(0.0);
                double Sy = it->second["attributes"]["Sy"].as<double>(0.0);
                        
                //Instantiate the material object.
                theMaterial = std::make_unique<Plastic1DJ2>(E, nu, rho, K, H, Sy);
            }
            else if(strcasecmp(Name.c_str(),"PlasticPlaneStrainJ2") == 0){
                double K = it->second["attributes"]["K"].as<double>(0.0);
                double G = it->second["attributes"]["G"].as<double>(0.0);
                double rho = it->second["attributes"]["rho"].as<double>(0.0);
                double H = it->second["attributes"]["h"].as<double>(0.0);
                double Sy = it->second["attributes"]["Sy"].as<double>(0.0);
                double beta = it->second["attributes"]["beta"].as<double>(0.0);

                //Instantiate the material object.
                theMaterial = std::make_unique<PlasticPlaneStrainJ2>(K, G, rho, H, beta, Sy);
            }
            else if(strcasecmp(Name.c_str(),"PlasticPlaneStrainBA") == 0){
                double K = it->second["attributes"]["K"].as<double>(0.0);
                double G = it->second["attributes"]["G"].as<double>(0.0);
                double rho = it->second["attributes"]["rho"].as<double>(0.0);
                double h0 = it->second["attributes"]["h0"].as<double>(0.0);
                double h = it->second["attributes"]["h"].as<double>(0.0);
                double m = it->second["attributes"]["m"].as<double>(0.0);
                double Su = it->second["attributes"]["Su"].as<double>(0.0);
                double beta = it->second["attributes"]["beta"].as<double>(0.0);

                //Instantiate the material object.
                theMaterial = std::make_unique<PlasticPlaneStrainBA>(K, G, rho, h0, h, m, Su, beta);
            }
            else if(strcasecmp(Name.c_str(),"Plastic3DJ2") == 0){
                double K = it->second["attributes"]["K"].as<double>(0.0);
                double G = it->second["attributes"]["G"].as<double>(0.0);
                double rho = it->second["attributes"]["rho"].as<double>(0.0);
                double h = it->second["attributes"]["h"].as<double>(0.0);
                double Sy = it->second["attributes"]["Sy"].as<double>(0.0);
                double beta = it->second["attributes"]["beta"].as<double>(0.0);

                //Instantiate the material object.
                theMaterial = std::make_unique<Plastic3DJ2>(K, G, rho, h, beta, Sy);
            }
            else if(strcasecmp(Name.c_str(),"Plastic3DBA") == 0){
                double K = it->second["attributes"]["K"].as<double>(0.0);
                double G = it->second["attributes"]["G"].as<double>(0.0);
                double rho = it->second["attributes"]["rho"].as<double>(0.0);
                double h0 = it->second["attributes"]["h0"].as<double>(0.0);
                double h = it->second["attributes"]["h"].as<double>(0.0);
                double m = it->second["attributes"]["m"].as<double>(0.0);
                double Su = it->second["attributes"]["Su"].as<double>(0.0);
                double beta = it->second["attributes"]["beta"].as<double>(0.0);

                //Instantiate the material object.
                theMaterial = std::make_unique<Plastic3DBA>(K, G, rho, h0, h, m, Su, beta);
            }
            else if(strcasecmp(Name.c_str(),"Elastic1DFiber") == 0){
                double E = it->second["attributes"]["E"].as<double>(0.0);
                double nu = it->second["attributes"]["nu"].as<double>(0.0);
                double rho = it->second["attributes"]["rho"].as<double>(0.0);

                //Instantiate the material (fiber) object.
                theMaterial = std::make_unique<Elastic1DFiber>(E, nu, rho);
            }
            else if(strcasecmp(Name.c_str(),"Elastic1DGap") == 0){
                double E = it->second["attributes"]["E"].as<double>(0.0);
                double gap = it->second["attributes"]["gap"].as<double>(0.0);
                double behavior = it->second["attributes"]["behavior"].as<bool>(false);

                //Instantiate the material (fiber) object.
                theMaterial = std::make_unique<Elastic1DGap>(E, gap, behavior);
            }
            else if(strcasecmp(Name.c_str(),"Plastic1DGap") == 0){
                double E = it->second["attributes"]["E"].as<double>(0.0);
                double fy = it->second["attributes"]["fy"].as<double>(0.0);
                double eta = it->second["attributes"]["ratio"].as<double>(0.0);
                double gap = it->second["attributes"]["gap"].as<double>(0.0);
                double behavior = it->second["attributes"]["behavior"].as<bool>(false);

                //Instantiate the material (fiber) object.
                theMaterial = std::make_unique<Plastic1DGap>(E, fy, gap, eta, behavior);
            }
            else if(strcasecmp(Name.c_str(),"Steel1DFiber") == 0){
                double E = it->second["attributes"]["E"].as<double>(0.0);
                double fy = it->second["attributes"]["fy"].as<double>(0.0);
                double b = it->second["attributes"]["b"].as<double>(0.0);
                double R0 = it->second["attributes"]["R0"].as<double>(15.00);
                double cR1 = it->second["attributes"]["cR1"].as<double>(0.925);
                double cR2 = it->second["attributes"]["cR2"].as<double>(0.150);
                double a1 = it->second["attributes"]["a1"].as<double>(0.0);
                double a2 = it->second["attributes"]["a2"].as<double>(1.0);
                double a3 = it->second["attributes"]["a3"].as<double>(0.0);
                double a4 = it->second["attributes"]["a4"].as<double>(1.0);
                double nu = it->second["attributes"]["nu"].as<double>(0.33);
                double rho = it->second["attributes"]["rho"].as<double>(0.0);

                //Instantiate the material (fiber) object.
                theMaterial = std::make_unique<Steel1DFiber>(fy, E, b, R0, cR1, cR2, a1, a2, a3, a4, nu, rho);
            }
            else if(strcasecmp(Name.c_str(),"Concrete1DFiber") == 0){
                double fc = it->second["attributes"]["fc"].as<double>();
                double ecc = it->second["attributes"]["ecc"].as<double>(-0.002);
                double fcu = it->second["attributes"]["fcu"].as<double>();
                double ecu = it->second["attributes"]["ecu"].as<double>(-0.012);
                double ratio = it->second["attributes"]["ratio"].as<double>(0.1);
                double ft = it->second["attributes"]["ft"].as<double>(0.0);
                double Et = it->second["attributes"]["Et"].as<double>(0.0);
                double nu = it->second["attributes"]["nu"].as<double>(0.25);
                double rho = it->second["attributes"]["rho"].as<double>(0.0);

                //Instantiate the material (fiber) object.
                theMaterial = std::make_unique<Concrete1DFiber>(fc, ecc, fcu, ecu, ratio, ft, Et, nu, rho);
            }

            //TODO: Add more material models here.

            //Stores the information in material container.
            theMesh->AddMaterial(Tag, theMaterial);
        }
    }
}

///Updates the Section Entities in Mesh Object 
///@param theMesh Pointer to the Mesh container.
///@param jsonFile json file where mesh entities will be readden.
void 
UpdateSections(std::shared_ptr<Mesh> &theMesh, RSJresource& jsonFile){
    //List of Sections already defined
    std::vector<unsigned int> Tags = theMesh->GetVectorIDs<unsigned int>("Sections");

    //Materials to be Added To Mesh Object
    for(auto it = jsonFile["Sections"].as_object().begin(); it != jsonFile["Sections"].as_object().end(); ++it){
        //Material name and identifier.
        unsigned int Tag = std::stoi(it->first);

        //Whether the material should be created 
        if (std::find(Tags.begin(), Tags.end(), Tag) == Tags.end()) {
            //Section name and identifier.
            unsigned int Tag = std::stoi(it->first);
            std::string Name = it->second["name"].as<std::string>();
            std::string Model = it->second["model"].as<std::string>();

            //Creates a section object.
            std::unique_ptr<Section> theSection;

            if(strcasecmp(Model.c_str(),"Plain") == 0){
                //Insertion point, section rotation and material identifier
                double theta = it->second["attributes"]["theta"].as<double>(0.0);
                unsigned int ip = it->second["attributes"]["ip"].as<int>(10);
                unsigned int matTag = it->second["attributes"]["material"].as<int>();

                if(strcasecmp(Name.c_str(),"Lin2DRectangular") == 0){
                    double h = it->second["attributes"]["h"].as<double>();
                    double b = it->second["attributes"]["b"].as<double>();

                    //Instantiate the section object.
                    theSection = std::make_unique<Lin2DRectangular>(h, b, theMesh->GetMaterial(matTag), theta, ip);
                }
                else if(strcasecmp(Name.c_str(),"Lin3DRectangular") == 0){
                    double h = it->second["attributes"]["h"].as<double>();
                    double b = it->second["attributes"]["b"].as<double>();

                    //Instantiate the section object.
                    theSection = std::make_unique<Lin3DRectangular>(h, b, theMesh->GetMaterial(matTag), theta, ip);
                }
                else if(strcasecmp(Name.c_str(),"Lin2DRectangularTube") == 0){
                    double h = it->second["attributes"]["h"].as<double>();
                    double b = it->second["attributes"]["b"].as<double>();
                    double tw = it->second["attributes"]["tw"].as<double>();
                    double tf = it->second["attributes"]["tf"].as<double>();

                    //Instantiate the section object.
                    theSection = std::make_unique<Lin2DRectangularTube>(h, b, tw, tf, theMesh->GetMaterial(matTag), theta, ip);
                }
                else if(strcasecmp(Name.c_str(),"Lin3DRectangularTube") == 0){
                    double h = it->second["attributes"]["h"].as<double>();
                    double b = it->second["attributes"]["b"].as<double>();
                    double tw = it->second["attributes"]["tw"].as<double>();
                    double tf = it->second["attributes"]["tf"].as<double>();

                    //Instantiate the section object.
                    theSection = std::make_unique<Lin3DRectangularTube>(h, b, tw, tf, theMesh->GetMaterial(matTag), theta, ip);
                }
                else if(strcasecmp(Name.c_str(),"Lin2DCircular") == 0){
                    double r = it->second["attributes"]["r"].as<double>();

                    //Instantiate the section object.
                    theSection = std::make_unique<Lin2DCircular>(r, theMesh->GetMaterial(matTag), theta, ip);
                }
                else if(strcasecmp(Name.c_str(),"Lin3DCircular") == 0){
                    double r = it->second["attributes"]["r"].as<double>();

                    //Instantiate the section object.
                    theSection = std::make_unique<Lin3DCircular>(r, theMesh->GetMaterial(matTag), theta, ip);
                }
                else if(strcasecmp(Name.c_str(),"Lin2DCircularTube") == 0){
                    double re = it->second["attributes"]["re"].as<double>();
                    double ri = it->second["attributes"]["ri"].as<double>();

                    //Instantiate the section object.
                    theSection = std::make_unique<Lin2DCircularTube>(re, ri, theMesh->GetMaterial(matTag), theta, ip);
                }
                else if(strcasecmp(Name.c_str(),"Lin3DCircularTube") == 0){
                    double re = it->second["attributes"]["re"].as<double>();
                    double ri = it->second["attributes"]["ri"].as<double>();

                    //Instantiate the section object.
                    theSection = std::make_unique<Lin3DCircularTube>(re, ri, theMesh->GetMaterial(matTag), theta, ip);
                }
                else if(strcasecmp(Name.c_str(),"Lin2DAngle") == 0){
                    double h = it->second["attributes"]["h"].as<double>();
                    double b = it->second["attributes"]["b"].as<double>();
                    double tw = it->second["attributes"]["tw"].as<double>();
                    double tf = it->second["attributes"]["tf"].as<double>();

                    //Instantiate the section object.
                    theSection = std::make_unique<Lin2DAngle>(h, b, tw, tf, theMesh->GetMaterial(matTag), theta, ip);
                }
                else if(strcasecmp(Name.c_str(),"Lin3DAngle") == 0){
                    double h = it->second["attributes"]["h"].as<double>();
                    double b = it->second["attributes"]["b"].as<double>();
                    double tw = it->second["attributes"]["tw"].as<double>();
                    double tf = it->second["attributes"]["tf"].as<double>();

                    //Instantiate the section object.
                    theSection = std::make_unique<Lin3DAngle>(h, b, tw, tf, theMesh->GetMaterial(matTag), theta, ip);
                }
                else if(strcasecmp(Name.c_str(),"Lin2DChannel") == 0){
                    double h = it->second["attributes"]["h"].as<double>();
                    double b = it->second["attributes"]["b"].as<double>();
                    double tw = it->second["attributes"]["tw"].as<double>();
                    double tf = it->second["attributes"]["tf"].as<double>();

                    //Instantiate the section object.
                    theSection = std::make_unique<Lin2DChannel>(h, b, tw, tf, theMesh->GetMaterial(matTag), theta, ip);
                }
                else if(strcasecmp(Name.c_str(),"Lin3DChannel") == 0){
                    double h = it->second["attributes"]["h"].as<double>();
                    double b = it->second["attributes"]["b"].as<double>();
                    double tw = it->second["attributes"]["tw"].as<double>();
                    double tf = it->second["attributes"]["tf"].as<double>();

                    //Instantiate the section object.
                    theSection = std::make_unique<Lin3DChannel>(h, b, tw, tf, theMesh->GetMaterial(matTag), theta, ip);
                }
                if(strcasecmp(Name.c_str(),"Lin2DTee") == 0){
                    double h = it->second["attributes"]["h"].as<double>();
                    double b = it->second["attributes"]["b"].as<double>();
                    double tw = it->second["attributes"]["tw"].as<double>();
                    double tf = it->second["attributes"]["tf"].as<double>();

                    //Instantiate the section object.
                    theSection = std::make_unique<Lin2DTee>(h, b, tw, tf, theMesh->GetMaterial(matTag), theta, ip);
                }
                if(strcasecmp(Name.c_str(),"Lin3DTee") == 0){
                    double h = it->second["attributes"]["h"].as<double>();
                    double b = it->second["attributes"]["b"].as<double>();
                    double tw = it->second["attributes"]["tw"].as<double>();
                    double tf = it->second["attributes"]["tf"].as<double>();

                    //Instantiate the section object.
                    theSection = std::make_unique<Lin3DTee>(h, b, tw, tf, theMesh->GetMaterial(matTag), theta, ip);
                }
                else if(strcasecmp(Name.c_str(),"Lin2DWideFlange") == 0){
                    double h = it->second["attributes"]["h"].as<double>();
                    double b = it->second["attributes"]["b"].as<double>();
                    double tw = it->second["attributes"]["tw"].as<double>();
                    double tf = it->second["attributes"]["tf"].as<double>();

                    //Instantiate the section object.
                    theSection = std::make_unique<Lin2DWideFlange>(h, b, tw, tf, theMesh->GetMaterial(matTag), theta, ip);
                }
                else if(strcasecmp(Name.c_str(),"Lin3DWideFlange") == 0){
                    double h = it->second["attributes"]["h"].as<double>();
                    double b = it->second["attributes"]["b"].as<double>();
                    double tw = it->second["attributes"]["tw"].as<double>();
                    double tf = it->second["attributes"]["tf"].as<double>();

                    //Instantiate the section object.
                    theSection = std::make_unique<Lin3DWideFlange>(h, b, tw, tf, theMesh->GetMaterial(matTag), theta, ip);
                }
                else if(strcasecmp(Name.c_str(),"Lin2DUserDefined") == 0){
                    std::vector<double> prop(3, 0.0);
                    prop[0] = it->second["attributes"]["A"].as<double>();
                    prop[1] = it->second["attributes"]["As2"].as<double>(0.0);
                    prop[2] = it->second["attributes"]["I33"].as<double>();

                    //Instantiate the section object.
                    theSection = std::make_unique<Lin2DUserDefined>(prop, theMesh->GetMaterial(matTag), theta);
                }
                else if(strcasecmp(Name.c_str(),"Lin3DUserDefined") == 0){
                    std::vector<double> prop(7, 0.0);
                    prop[0] = it->second["attributes"]["A"].as<double>();
                    prop[3] = it->second["attributes"]["J"].as<double>();
                    prop[1] = it->second["attributes"]["As2"].as<double>(0.0);
                    prop[2] = it->second["attributes"]["As3"].as<double>(0.0);
                    prop[4] = it->second["attributes"]["I22"].as<double>();
                    prop[5] = it->second["attributes"]["I33"].as<double>();
                    prop[6] = it->second["attributes"]["I23"].as<double>(0.0);

                    //Instantiate the section object.
                    theSection = std::make_unique<Lin3DUserDefined>(prop, theMesh->GetMaterial(matTag), theta);
                }
                else if(strcasecmp(Name.c_str(),"Lin3DThinArea") == 0){
                    double th = it->second["attributes"]["th"].as<double>();

                    //Instantiate the section object.
                    theSection = std::make_unique<Lin3DThinArea>(th, theMesh->GetMaterial(matTag));
                }
            }
            else if(strcasecmp(Model.c_str(),"Fiber") == 0){
                if(strcasecmp(Name.c_str(),"Fib3DLineSection") == 0){
                    //Initialize the vector with fiber coordinates, and area
                    unsigned int numOfFibers = it->second["attributes"]["fiber"].size();
                    std::vector<double> zi(numOfFibers);
                    std::vector<double> yi(numOfFibers);
                    std::vector<double> Ai(numOfFibers);
                    std::vector<std::unique_ptr<Material> > fibers(numOfFibers);

                    //Generates the Fiber section object
                    for(unsigned int k =0; k < numOfFibers; k++){
                        unsigned int fibTag = it->second["attributes"]["fiber"][k].as<int>();
                        zi[k] = it->second["attributes"]["zi"][k].as<double>();
                        yi[k] = it->second["attributes"]["yi"][k].as<double>();
                        Ai[k] = it->second["attributes"]["Ai"][k].as<double>();

                        theMesh->AddFiber(fibTag, fibers[k]);
                    }
                    //Cross section Shear factor (Ash/A): example 5/6 rectangle, 0.9 circular 
                    double kappa2 = it->second["attributes"]["kappa2"].as<double>(1.00);
                    double kappa3 = it->second["attributes"]["kappa3"].as<double>(1.00);

                    //Section Insertion Point bounding Box
                    double h = it->second["attributes"]["h"].as<double>();
                    double b = it->second["attributes"]["b"].as<double>();
                    unsigned int ip = it->second["attributes"]["ip"].as<int>(10);

                    //Section angle of rotation
                    double theta = it->second["attributes"]["theta"].as<double>(0.0);
                    
                    //Instantiate the section object.
                    theSection = std::make_unique<Fib3DLineSection>(h, b, fibers, zi, yi, Ai, kappa2, kappa3, theta, ip);
                }
                else if(strcasecmp(Name.c_str(),"Fib3DAreaSection") == 0){
                    //Instantiate the section object.
                    //theSection = std::make_unique<Fib3DAreaSection>(...);
                }
            }
            else if(strcasecmp(Model.c_str(),"General") == 0){
                //TODO: Parse general sections components.
            }

            //TODO: Add more section models here.

            //Stores the section in the mesh.
            theMesh->AddSection(Tag, theSection);
        }
    }
}

///Updates the Element Entities in Mesh Object 
///@param theMesh Pointer to the Mesh container.
///@param jsonFile json file where mesh entities will be readden.
void 
UpdateElements(std::shared_ptr<Mesh> &theMesh, RSJresource& jsonFile){
    //Element Identifiers To Be Updated
    std::map<std::string, std::vector<unsigned int> > Tags = Entities2Update<unsigned int>(theMesh, jsonFile, "Elements");

    //Element to be Removed from Mesh Object
    for(unsigned int k = 0; k < Tags["del"].size(); k++){
        unsigned int Tag = Tags["del"][k];
        theMesh->DelElement(Tag);
    }

    //Element to be Added To Mesh Object
    for(unsigned int k = 0; k < Tags["add"].size(); k++){
        unsigned int Tag = Tags["add"][k];

        //Node identifier and number of degree-of-freedom:
        std::string eTag = std::to_string(Tag);
        std::string Name = jsonFile["Elements"][eTag]["name"].as<std::string>();

        //Node connectivity
        unsigned int n = jsonFile["Elements"][eTag]["conn"].size();
        std::vector<unsigned int> nodes(n);
        for(unsigned int k = 0; k < n; ++k)
            nodes[k] = jsonFile["Elements"][eTag]["conn"][k].as<int>();

        //Creates an element object.
        std::shared_ptr<Element> theElement;

        if(strcasecmp(Name.c_str(),"lin2DTruss2") == 0){
            unsigned int matID = jsonFile["Elements"][eTag]["attributes"]["material"].as<int>();       
            double parameters = jsonFile["Elements"][eTag]["attributes"]["area"].as<double>();

            //Instantiate the lin2DTruss2 element.
            theElement = std::make_shared<lin2DTruss2>(nodes, theMesh->GetMaterial(matID), parameters);
        }
        else if(strcasecmp(Name.c_str(),"kin2DTruss2") == 0){  
            unsigned int matID = jsonFile["Elements"][eTag]["attributes"]["material"].as<int>();  
            double parameters = jsonFile["Elements"][eTag]["attributes"]["area"].as<double>();

            //Instantiate the kin2DTruss2 element.
            theElement = std::make_shared<kin2DTruss2>(nodes, theMesh->GetMaterial(matID), parameters);
        }
        else if(strcasecmp(Name.c_str(),"lin3DTruss2") == 0){   
            unsigned int matID = jsonFile["Elements"][eTag]["attributes"]["material"].as<int>();     
            double parameters = jsonFile["Elements"][eTag]["attributes"]["area"].as<double>();

            //Instantiate the lin3DTruss2 element.
            theElement = std::make_shared<lin3DTruss2>(nodes, theMesh->GetMaterial(matID), parameters);
        }
        else if(strcasecmp(Name.c_str(),"kin3DTruss2") == 0){
            unsigned int matID = jsonFile["Elements"][eTag]["attributes"]["material"].as<int>();
            double parameters = jsonFile["Elements"][eTag]["attributes"]["area"].as<double>();

            //Instantiate the kin3DTruss2 element.
            theElement = std::make_shared<kin3DTruss2>(nodes, theMesh->GetMaterial(matID), parameters);
        }
        else if(strcasecmp(Name.c_str(),"lin2DTruss3") == 0){      
            double parameters = jsonFile["Elements"][eTag]["attributes"]["area"].as<double>();
            unsigned int nGauss = jsonFile["Elements"][eTag]["attributes"]["np"].as<int>(3);
            unsigned int matID = jsonFile["Elements"][eTag]["attributes"]["material"].as<int>();
            std::string Quadrature = jsonFile["Elements"][eTag]["attributes"]["rule"].as<std::string>("GAUSS");

            //Instantiate the lin2DTruss3 element.
            theElement = std::make_shared<lin2DTruss3>(nodes, theMesh->GetMaterial(matID), parameters, Quadrature, nGauss);
        }
        else if(strcasecmp(Name.c_str(),"lin3DTruss3") == 0){   
            double parameters = jsonFile["Elements"][eTag]["attributes"]["area"].as<double>();
            unsigned int nGauss = jsonFile["Elements"][eTag]["attributes"]["np"].as<int>(3);
            unsigned int matID = jsonFile["Elements"][eTag]["attributes"]["material"].as<int>();
            std::string Quadrature = jsonFile["Elements"][eTag]["attributes"]["rule"].as<std::string>("GAUSS");

            //Instantiate the lin3DTruss3 element.
            theElement = std::make_shared<lin3DTruss3>(nodes, theMesh->GetMaterial(matID), parameters, Quadrature, nGauss);
        }
        else if(strcasecmp(Name.c_str(),"ZeroLength1D") == 0){
            unsigned int dir = jsonFile["Elements"][eTag]["attributes"]["dir"].as<int>();
            unsigned int matID = jsonFile["Elements"][eTag]["attributes"]["material"].as<int>();

            //Instantiate the zerolength1D element.
            theElement = std::make_shared<ZeroLength1D>(nodes, theMesh->GetMaterial(matID), dir);
        }
        else if(strcasecmp(Name.c_str(),"UnxBoucWen2DLink") == 0){
            std::vector<double> variables(2, 0.0);
            std::vector<double> parameters(4, 0.0);
            double tol = jsonFile["Elements"][eTag]["attributes"]["tol"].as<double>(1E-06);
            unsigned int nmax = jsonFile["Elements"][eTag]["attributes"]["nmax"].as<int>(50);
            unsigned int dir = jsonFile["Elements"][eTag]["attributes"]["dir"].as<int>();
            unsigned int dim = jsonFile["Elements"][eTag]["attributes"]["dim"].as<int>();

            variables[0] = jsonFile["Elements"][eTag]["attributes"]["Fy"].as<double>();
            variables[1] = jsonFile["Elements"][eTag]["attributes"]["k0"].as<double>();
            parameters[0] = jsonFile["Elements"][eTag]["attributes"]["alpha"].as<double>(1.0);
            parameters[1] = jsonFile["Elements"][eTag]["attributes"]["eta"].as<double>(1.0);
            parameters[2] = jsonFile["Elements"][eTag]["attributes"]["beta"].as<double>(0.5);
            parameters[3] = jsonFile["Elements"][eTag]["attributes"]["gamma"].as<double>(0.5);

            //Instantiate the UnxBoucWen2DLink element.
            theElement = std::make_shared<UnxBoucWen2DLink>(nodes, parameters, variables, dim, dir, tol, nmax);
        }
        else if(strcasecmp(Name.c_str(),"UnxBoucWen3DLink") == 0){
            std::vector<double> variables(2, 0.0);
            std::vector<double> parameters(4, 0.0);
            double tol = jsonFile["Elements"][eTag]["attributes"]["tol"].as<double>(1E-06);
            unsigned int nmax = jsonFile["Elements"][eTag]["attributes"]["nmax"].as<int>(50);
            unsigned int dir = jsonFile["Elements"][eTag]["attributes"]["dir"].as<int>();
            unsigned int dim = jsonFile["Elements"][eTag]["attributes"]["dim"].as<int>();

            variables[0] = jsonFile["Elements"][eTag]["attributes"]["Fy"].as<double>();
            variables[1] = jsonFile["Elements"][eTag]["attributes"]["k0"].as<double>();
            parameters[0] = jsonFile["Elements"][eTag]["attributes"]["alpha"].as<double>(1.0);
            parameters[1] = jsonFile["Elements"][eTag]["attributes"]["eta"].as<double>(1.0);
            parameters[2] = jsonFile["Elements"][eTag]["attributes"]["beta"].as<double>(0.5);
            parameters[3] = jsonFile["Elements"][eTag]["attributes"]["gamma"].as<double>(0.5);

            //Instantiate the UnxBoucWen3DLink element.
            theElement = std::make_shared<UnxBoucWen3DLink>(nodes, parameters, variables, dim, dir, tol, nmax);
        }
        else if(strcasecmp(Name.c_str(),"HDRBYamamoto2DLink") == 0){
            double De = jsonFile["Elements"][eTag]["attributes"]["De"].as<double>(1.3);
            double Di = jsonFile["Elements"][eTag]["attributes"]["Di"].as<double>(0.3);
            double Hr = jsonFile["Elements"][eTag]["attributes"]["Hr"].as<double>(0.261);
            unsigned int dim = jsonFile["Elements"][eTag]["attributes"]["dim"].as<int>();

            //Instantiate the HDRBYamamoto2DLink element.
            theElement = std::make_shared<HDRBYamamoto2DLink>(nodes, De, Di, Hr, dim);
        }
        else if(strcasecmp(Name.c_str(),"HDRBYamamoto3DLink") == 0){
            double De = jsonFile["Elements"][eTag]["attributes"]["De"].as<double>(1.3);
            double Di = jsonFile["Elements"][eTag]["attributes"]["Di"].as<double>(0.3);
            double Hr = jsonFile["Elements"][eTag]["attributes"]["Hr"].as<double>(0.261);
            unsigned int dim = jsonFile["Elements"][eTag]["attributes"]["dim"].as<int>();

            //Instantiate the HDRBYamamoto3DLink element.
            theElement = std::make_shared<HDRBYamamoto3DLink>(nodes, De, Di, Hr, dim);
        }
        else if(strcasecmp(Name.c_str(),"lin2DFrame2") == 0){
            unsigned int nGauss = jsonFile["Elements"][eTag]["attributes"]["np"].as<int>(3);
            unsigned int secID = jsonFile["Elements"][eTag]["attributes"]["section"].as<int>();
            std::string Quadrature = jsonFile["Elements"][eTag]["attributes"]["rule"].as<std::string>("GAUSS");
            std::string Formulation = jsonFile["Elements"][eTag]["attributes"]["formulation"].as<std::string>("Bernoulli");

            //Frame element formulation.
            bool Condition = false;
            if(strcasecmp(Formulation.c_str(),"Timoshenko") == 0)
                Condition = true;
                        
            //Instantiate the lin2DFrame2 element.
            theElement = std::make_shared<lin2DFrame2>(nodes, theMesh->GetSection(secID), Condition, Quadrature, nGauss);
        }
        else if(strcasecmp(Name.c_str(),"kin2DFrame2") == 0){
            unsigned int secID = jsonFile["Elements"][eTag]["attributes"]["section"].as<int>();
                        
            //Instantiate the kin2DFrame2 element.
            theElement = std::make_shared<kin2DFrame2>(nodes, theMesh->GetSection(secID));
        }
        else if(strcasecmp(Name.c_str(),"lin3DFrame2") == 0){
            unsigned int nGauss = jsonFile["Elements"][eTag]["attributes"]["np"].as<int>(3);
            unsigned int secID = jsonFile["Elements"][eTag]["attributes"]["section"].as<int>();
            std::string Quadrature = jsonFile["Elements"][eTag]["attributes"]["rule"].as<std::string>("GAUSS");
            std::string Formulation = jsonFile["Elements"][eTag]["attributes"]["formulation"].as<std::string>("Bernoulli");

            //Frame element formulation.
            bool Condition = false;
            if(strcasecmp(Formulation.c_str(),"Timoshenko") == 0)
                Condition = true;

            //Instantiate the lin3DFrame2 element.
            theElement = std::make_shared<lin3DFrame2>(nodes, theMesh->GetSection(secID), Condition, Quadrature, nGauss);
        }
        else if(strcasecmp(Name.c_str(),"lin2DTria3") == 0){
            double thickness = jsonFile["Elements"][eTag]["attributes"]["th"].as<double>();
            unsigned int nGauss = jsonFile["Elements"][eTag]["attributes"]["np"].as<int>(4);
            unsigned int matID = jsonFile["Elements"][eTag]["attributes"]["material"].as<int>();
            std::string Quadrature = jsonFile["Elements"][eTag]["attributes"]["rule"].as<std::string>("GAUSS");

            //Instantiate the lin2DTria3 element.
            theElement = std::make_shared<lin2DTria3>(nodes, theMesh->GetMaterial(matID), thickness, Quadrature, nGauss);
        }
        else if(strcasecmp(Name.c_str(),"lin2DTria6") == 0){
            double thickness = jsonFile["Elements"][eTag]["attributes"]["th"].as<double>();
            unsigned int nGauss = jsonFile["Elements"][eTag]["attributes"]["np"].as<int>(9);
            unsigned int matID = jsonFile["Elements"][eTag]["attributes"]["material"].as<int>();
            std::string Quadrature = jsonFile["Elements"][eTag]["attributes"]["rule"].as<std::string>("GAUSS");

            //Instantiate the lin2DTria6 element.
            theElement = std::make_shared<lin2DTria6>(nodes, theMesh->GetMaterial(matID), thickness, Quadrature, nGauss);
        }
        else if(strcasecmp(Name.c_str(),"lin2DQuad4") == 0){
            double thickness = jsonFile["Elements"][eTag]["attributes"]["th"].as<double>();
            unsigned int nGauss = jsonFile["Elements"][eTag]["attributes"]["np"].as<int>(4);
            unsigned int matID = jsonFile["Elements"][eTag]["attributes"]["material"].as<int>();
            std::string Quadrature = jsonFile["Elements"][eTag]["attributes"]["rule"].as<std::string>("GAUSS");

            //Instantiate the lin2DQuad4 element.
            theElement = std::make_shared<lin2DQuad4>(nodes, theMesh->GetMaterial(matID), thickness, Quadrature, nGauss);
        }
        else if(strcasecmp(Name.c_str(),"lin2DQuad8") == 0){
            double thickness = jsonFile["Elements"][eTag]["attributes"]["th"].as<double>();
            unsigned int nGauss = jsonFile["Elements"][eTag]["attributes"]["np"].as<int>(9);
            unsigned int matID = jsonFile["Elements"][eTag]["attributes"]["material"].as<int>();
            std::string Quadrature = jsonFile["Elements"][eTag]["attributes"]["rule"].as<std::string>("GAUSS");

            //Instantiate the lin2DQuad8 element.
            theElement = std::make_shared<lin2DQuad8>(nodes, theMesh->GetMaterial(matID), thickness, Quadrature, nGauss);
        }
        else if(strcasecmp(Name.c_str(),"PML2DQuad4") == 0){
            std::vector<double> parameters(8);
            parameters[0] = jsonFile["Elements"][eTag]["attributes"]["th"].as<double>();
            parameters[1] = jsonFile["Elements"][eTag]["attributes"]["n"].as<double>();
            parameters[2] = jsonFile["Elements"][eTag]["attributes"]["L"].as<double>();
            parameters[3] = jsonFile["Elements"][eTag]["attributes"]["R"].as<double>();
            parameters[4] = jsonFile["Elements"][eTag]["attributes"]["x0"][0].as<double>();
            parameters[5] = jsonFile["Elements"][eTag]["attributes"]["x0"][1].as<double>();
            parameters[6] = jsonFile["Elements"][eTag]["attributes"]["npml"][0].as<double>();
            parameters[7] = jsonFile["Elements"][eTag]["attributes"]["npml"][1].as<double>();

            unsigned int nGauss = jsonFile["Elements"][eTag]["attributes"]["np"].as<int>(4);
            unsigned int matID = jsonFile["Elements"][eTag]["attributes"]["material"].as<int>();
            std::string Quadrature = jsonFile["Elements"][eTag]["attributes"]["rule"].as<std::string>("GAUSS");

            //Instantiate the PML2DQuad4 element.
            theElement = std::make_shared<PML2DQuad4>(nodes, theMesh->GetMaterial(matID), parameters, Quadrature, nGauss);
        }
        else if(strcasecmp(Name.c_str(),"PML2DQuad8") == 0){
            std::vector<double> parameters(8);
            parameters[0] = jsonFile["Elements"][eTag]["attributes"]["th"].as<double>();
            parameters[1] = jsonFile["Elements"][eTag]["attributes"]["n"].as<double>();
            parameters[2] = jsonFile["Elements"][eTag]["attributes"]["L"].as<double>();
            parameters[3] = jsonFile["Elements"][eTag]["attributes"]["R"].as<double>();
            parameters[4] = jsonFile["Elements"][eTag]["attributes"]["x0"][0].as<double>();
            parameters[5] = jsonFile["Elements"][eTag]["attributes"]["x0"][1].as<double>();
            parameters[6] = jsonFile["Elements"][eTag]["attributes"]["npml"][0].as<double>();
            parameters[7] = jsonFile["Elements"][eTag]["attributes"]["npml"][1].as<double>();

            unsigned int nGauss = jsonFile["Elements"][eTag]["attributes"]["np"].as<int>(9);
            unsigned int matID = jsonFile["Elements"][eTag]["attributes"]["material"].as<int>();
            std::string Quadrature = jsonFile["Elements"][eTag]["attributes"]["rule"].as<std::string>("GAUSS");

            //Instantiate the PML2DQuad4 element.
            theElement = std::make_shared<PML2DQuad8>(nodes, theMesh->GetMaterial(matID), parameters, Quadrature, nGauss);
        }
        else if(strcasecmp(Name.c_str(),"kin2DQuad4") == 0){
            double thickness = jsonFile["Elements"][eTag]["attributes"]["th"].as<double>();
            unsigned int nGauss = jsonFile["Elements"][eTag]["attributes"]["np"].as<int>(4);
            unsigned int matID = jsonFile["Elements"][eTag]["attributes"]["material"].as<int>();
            std::string Quadrature = jsonFile["Elements"][eTag]["attributes"]["rule"].as<std::string>("GAUSS");

            //Instantiate the kin2DQuad4 element.
            theElement = std::make_shared<kin2DQuad4>(nodes, theMesh->GetMaterial(matID), thickness, Quadrature, nGauss);
        }
        else if(strcasecmp(Name.c_str(),"lin3DShell4") == 0){
            unsigned int nGauss = jsonFile["Elements"][eTag]["attributes"]["np"].as<int>(9);
            unsigned int secID = jsonFile["Elements"][eTag]["attributes"]["section"].as<int>();
            std::string Quadrature = jsonFile["Elements"][eTag]["attributes"]["rule"].as<std::string>("GAUSS");
                    
            //Instantiate the lin3DShell4 element.
            theElement = std::make_shared<lin3DShell4>(nodes, theMesh->GetSection(secID), Quadrature, nGauss);
        }
        else if(strcasecmp(Name.c_str(),"lin3DTetra4") == 0){
            unsigned int nGauss = jsonFile["Elements"][eTag]["attributes"]["np"].as<int>(4);
            unsigned int matID = jsonFile["Elements"][eTag]["attributes"]["material"].as<int>();
            std::string Quadrature = jsonFile["Elements"][eTag]["attributes"]["rule"].as<std::string>("GAUSS");

            //Instantiate the lin3DTetra4 element.
            theElement = std::make_shared<lin3DTetra4>(nodes, theMesh->GetMaterial(matID), Quadrature, nGauss);
        }
        else if(strcasecmp(Name.c_str(),"lin3DTetra10") == 0){
            unsigned int nGauss = jsonFile["Elements"][eTag]["attributes"]["np"].as<int>(11);
            unsigned int matID = jsonFile["Elements"][eTag]["attributes"]["material"].as<int>();
            std::string Quadrature = jsonFile["Elements"][eTag]["attributes"]["rule"].as<std::string>("GAUSS");

            //Instantiate the lin3DTetra10 element.
            theElement = std::make_shared<lin3DTetra10>(nodes, theMesh->GetMaterial(matID), Quadrature, nGauss);
        }
        else if(strcasecmp(Name.c_str(),"lin3DHexa8") == 0){
            unsigned int nGauss = jsonFile["Elements"][eTag]["attributes"]["np"].as<int>(8);
            unsigned int matID = jsonFile["Elements"][eTag]["attributes"]["material"].as<int>();
            std::string Quadrature = jsonFile["Elements"][eTag]["attributes"]["rule"].as<std::string>("GAUSS");

            //Instantiate the lin3DHexa8 element.
            theElement = std::make_shared<lin3DHexa8>(nodes, theMesh->GetMaterial(matID), Quadrature, nGauss);
        }
        else if(strcasecmp(Name.c_str(),"kin3DHexa8") == 0){
            unsigned int nGauss = jsonFile["Elements"][eTag]["attributes"]["np"].as<int>(8);
            unsigned int matID = jsonFile["Elements"][eTag]["attributes"]["material"].as<int>();
            std::string Quadrature = jsonFile["Elements"][eTag]["attributes"]["rule"].as<std::string>("GAUSS");

            //Instantiate the kin3DHexa8 element.
            theElement = std::make_shared<kin3DHexa8>(nodes, theMesh->GetMaterial(matID), Quadrature, nGauss);
        }
        else if(strcasecmp(Name.c_str(),"PML3DHexa8") == 0){
            std::vector<double> parameters(9);
            parameters[0] = jsonFile["Elements"][eTag]["attributes"]["n"].as<double>();
            parameters[1] = jsonFile["Elements"][eTag]["attributes"]["L"].as<double>();
            parameters[2] = jsonFile["Elements"][eTag]["attributes"]["R"].as<double>();
            parameters[3] = jsonFile["Elements"][eTag]["attributes"]["x0"][0].as<double>();
            parameters[4] = jsonFile["Elements"][eTag]["attributes"]["x0"][1].as<double>();
            parameters[5] = jsonFile["Elements"][eTag]["attributes"]["x0"][2].as<double>();
            parameters[6] = jsonFile["Elements"][eTag]["attributes"]["npml"][0].as<double>();
            parameters[7] = jsonFile["Elements"][eTag]["attributes"]["npml"][1].as<double>();
            parameters[8] = jsonFile["Elements"][eTag]["attributes"]["npml"][2].as<double>();

            unsigned int nGauss = jsonFile["Elements"][eTag]["attributes"]["np"].as<int>(8);
            unsigned int matID = jsonFile["Elements"][eTag]["attributes"]["material"].as<int>();
            std::string Quadrature = jsonFile["Elements"][eTag]["attributes"]["rule"].as<std::string>("GAUSS");

            //Instantiate the PML3DHexa8 element.
            theElement = std::make_shared<PML3DHexa8>(nodes, theMesh->GetMaterial(matID), parameters, Quadrature, nGauss);
        }
        else if(strcasecmp(Name.c_str(),"lin3DHexa20") == 0){
            unsigned int nGauss = jsonFile["Elements"][eTag]["attributes"]["np"].as<int>(27);
            unsigned int matID = jsonFile["Elements"][eTag]["attributes"]["material"].as<int>();
            std::string Quadrature = jsonFile["Elements"][eTag]["attributes"]["rule"].as<std::string>("GAUSS");

            //Instantiate the lin3DHexa20 element.
            theElement = std::make_shared<lin3DHexa20>(nodes, theMesh->GetMaterial(matID), Quadrature, nGauss);
        }
        else if(strcasecmp(Name.c_str(),"TIEQlin2DQuad4") == 0){
            double cf1 = jsonFile["Elements"][eTag]["attributes"]["cf1"].as<double>();
            double cf2 = jsonFile["Elements"][eTag]["attributes"]["cf2"].as<double>();
            double zref = jsonFile["Elements"][eTag]["attributes"]["zref"].as<double>();
            double eref = jsonFile["Elements"][eTag]["attributes"]["eref"].as<double>();
            double thickness = jsonFile["Elements"][eTag]["attributes"]["th"].as<double>();
            unsigned int matID = jsonFile["Elements"][eTag]["attributes"]["material"].as<int>();
            std::string eType = jsonFile["Elements"][eTag]["attributes"]["type"].as<std::string>("Darandelli");
            unsigned int nGauss = jsonFile["Elements"][eTag]["attributes"]["np"].as<int>(4);
            std::string Quadrature = jsonFile["Elements"][eTag]["attributes"]["rule"].as<std::string>("GAUSS");

            //Instantiate the TIEQlin2DQuad4 element.
            theElement = std::make_shared<TIEQlin2DQuad4>(nodes, theMesh->GetMaterial(matID), thickness, Quadrature, nGauss, eType, zref, cf1, cf2, eref);
        }
        else if(strcasecmp(Name.c_str(),"EQlin2DQuad4") == 0){
            double cf1 = jsonFile["Elements"][eTag]["attributes"]["cf1"].as<double>();
            double cf2 = jsonFile["Elements"][eTag]["attributes"]["cf2"].as<double>();
            double zref = jsonFile["Elements"][eTag]["attributes"]["zref"].as<double>();
            double thickness = jsonFile["Elements"][eTag]["attributes"]["th"].as<double>();
            unsigned int matID = jsonFile["Elements"][eTag]["attributes"]["material"].as<int>();
            std::string eType = jsonFile["Elements"][eTag]["attributes"]["type"].as<std::string>("Darandelli");
            unsigned int nGauss = jsonFile["Elements"][eTag]["attributes"]["np"].as<int>(4);
            std::string Quadrature = jsonFile["Elements"][eTag]["attributes"]["rule"].as<std::string>("GAUSS");

            //Instantiate the EQlin2DQuad4 element.
            theElement = std::make_shared<EQlin2DQuad4>(nodes, theMesh->GetMaterial(matID), thickness, Quadrature, nGauss, eType, zref, cf1, cf2);
        }
        else if(strcasecmp(Name.c_str(),"null2DFrame2") == 0){
            //Instantiate the null2DFrame2 element.
            theElement = std::make_shared<null2DFrame2>(nodes);
        }
        else if(strcasecmp(Name.c_str(),"null3DFrame2") == 0){
            //Instantiate the null3DFrame2 element.
            theElement = std::make_shared<null3DFrame2>(nodes);
        }

        //TODO: Add more elements models here.

        //Stores the information in element container.
        theMesh->AddElement(Tag, theElement);
    }
}

///Updates the Damping Entities in Mesh Object 
///@param theMesh Pointer to the Mesh container.
///@param jsonFile json file where mesh entities will be readden.
void 
UpdateDampings(std::shared_ptr<Mesh> &theMesh, RSJresource& jsonFile){
    //Damping Identifiers To Be Updated
    std::map<std::string, std::vector<unsigned int> > Tags = Entities2Update<unsigned int>(theMesh, jsonFile, "Dampings");

    //Damping to be Removed from Mesh Object
    for(unsigned int k = 0; k < Tags["add"].size(); k++){
        unsigned int Tag = Tags["add"][k];

        //Damping name and identifier.
        std::string dTag = std::to_string(Tag);
        std::string Name = jsonFile["Dampings"][dTag]["name"].as<std::string>();

        //Damping parameters
        std::vector<double> parameters;
        if(strcasecmp(Name.c_str(),"Free") == 0){
            parameters.resize(2);
            parameters[0] = 0.0;
            parameters[1] = 0.0;
        }
        else if(strcasecmp(Name.c_str(),"Rayleigh") == 0){
            parameters.resize(2);
            parameters[0] = jsonFile["Dampings"][dTag]["attributes"]["am"].as<double>(0.0);
            parameters[1] = jsonFile["Dampings"][dTag]["attributes"]["ak"].as<double>(0.0);
        }
        else if(strcasecmp(Name.c_str(),"Caughey") == 0){
            //TODO: Caughey-viscous damping is not implemented yet.
            //parameters.resize(4);
        }
        else if(strcasecmp(Name.c_str(),"Capped") == 0){
            //TODO: Capped-viscous damping is not implemented yet.
        }

        //Damping associated to elements
        unsigned int n = jsonFile["Dampings"][dTag]["attributes"]["list"].size();
        std::vector<unsigned int> eList(n);

        for(unsigned int k = 0; k < n; k++)
            eList[k] = jsonFile["Dampings"][dTag]["attributes"]["list"][k].as<int>();
            
        //Creates a damping object.
        std::shared_ptr<Damping> theDamping = std::make_shared<Damping>(Name, parameters);

        //Stores the damping in the mesh.
        theMesh->AddDamping(Tag, theDamping);

        //Assign damping to group of elements.
        theMesh->SetDamping(Tag, eList);
    }

    //Damping to be Modified in Mesh Object
    if( Tags["mod"].size() > 0){
        std::map<unsigned int, std::shared_ptr<Damping> > theDampings = theMesh->GetDampings();

        for(unsigned int k = 0; k < Tags["mod"].size(); k++){
            unsigned int Tag = Tags["mod"][k];

            //Damping name and identifier.
            std::string dTag = std::to_string(Tag);
            std::string Name = jsonFile["Dampings"][dTag]["name"].as<std::string>();

            //Damping parameters
            std::vector<double> parameters;
            if(strcasecmp(Name.c_str(),"Free") == 0){
                parameters.resize(2);
                parameters[0] = 0.0;
                parameters[1] = 0.0;
            }
            else if(strcasecmp(Name.c_str(),"Rayleigh") == 0){
                parameters.resize(2);
                parameters[0] = jsonFile["Dampings"][dTag]["attributes"]["am"].as<double>(0.0);
                parameters[1] = jsonFile["Dampings"][dTag]["attributes"]["ak"].as<double>(0.0);
            }
            else if(strcasecmp(Name.c_str(),"Caughey") == 0){
                //TODO: Caughey-viscous damping is not implemented yet.
                //parameters.resize(4);
            }
            else if(strcasecmp(Name.c_str(),"Capped") == 0){
                //TODO: Capped-viscous damping is not implemented yet.
            }

            //Update damping associated to elements
            unsigned int n = jsonFile["Dampings"][dTag]["attributes"]["list"].size();
            std::vector<unsigned int> eList(n);

            for(unsigned int k = 0; k < n; k++)
                eList[k] = jsonFile["Dampings"][dTag]["attributes"]["list"][k].as<int>();

            //Update at Damping and Element level.
            theDampings[Tag]->SetName(Name);
            theDampings[Tag]->SetParameters(parameters);
            theMesh->SetDamping(Tag, eList);
        }
    }
}

///Updates the Loads Entities in Mesh Object 
///@param theMesh Pointer to the Mesh container.
///@param jsonFile json file where mesh entities will be readden.
void 
UpdateLoads(std::shared_ptr<Mesh> &theMesh, RSJresource& jsonFile){
    //Load Identifiers To Be Updated
    std::map<std::string, std::vector<unsigned int> > Tags = Entities2Update<unsigned int>(theMesh, jsonFile, "Loads");

    //Load to be Removed from Mesh Object
    for(unsigned int k = 0; k < Tags["del"].size(); k++){
        int Tag = Tags["del"][k];
        theMesh->DelLoad(Tag);
    }

    //Load to be Added To Mesh Object
    for(unsigned int k = 0; k < Tags["add"].size(); k++){
        unsigned int Tag = Tags["add"][k];

        //Load name and identifier.
        std::string lTag = std::to_string(Tag);
        std::string Name = jsonFile["Loads"][lTag]["name"].as<std::string>();

        //Creates a load object.
        std::shared_ptr<Load> theLoad;

        if(strcasecmp(Name.c_str(),"PointLoad") == 0){
            //Node identifiers that share this load
            unsigned int n = jsonFile["Loads"][lTag]["attributes"]["list"].size();
            std::vector<unsigned int> nodes(n);
            for(unsigned int k = 0; k < n; k++)
                nodes[k] = jsonFile["Loads"][lTag]["attributes"]["list"][k].as<int>();

            //Direction of the applied load
            unsigned int ndirs = jsonFile["Loads"][lTag]["attributes"]["dir"].size();
            Eigen::VectorXd direction(ndirs);
            for(unsigned int k = 0; k < ndirs; k++)
                direction(k) = jsonFile["Loads"][lTag]["attributes"]["dir"][k].as<double>();

            //The assigned load pattern function
            std::vector<double> value;
            std::string pattern = jsonFile["Loads"][lTag]["attributes"]["type"].as<std::string>();
            std::string function = jsonFile["Loads"][lTag]["attributes"]["name"].as<std::string>();

            if(strcasecmp(pattern.c_str(),"CONCENTRATED") == 0){
                if (strcasecmp(function.c_str(),"CONSTANT") == 0){
                    //The constant force magnitude
                    value.resize(1);
                    value[0] = jsonFile["Loads"][lTag]["attributes"]["mag"].as<double>();

                    //Creates the constant point load.
                    theLoad = std::make_shared<Load>(direction, value, POINTLOAD_CONCENTRATED_CONSTANT);

                    //Add the nodes that shares this load.
                    theLoad->AddNodes(nodes);
                }
                else if(strcasecmp(function.c_str(),"TIMESERIES") == 0){
                    //Loads the time history values into memory
                    std::string pathfile = jsonFile["Loads"][lTag]["attributes"]["file"].as<std::string>();
                    std::ifstream load(pathfile.c_str());

                    //The File is Opened and Ready to be Loaded.
                    if(load.is_open()){
                        //Number of time steps.
                        unsigned int nt;
                        load >> nt;

                        //Time-history load vector.
                        value.resize(nt);
                        for(unsigned int j = 0; j < nt; j++)
                            load >> value[j];
                    }
                    load.close();

                    //Creates the timeseries point load.
                    theLoad = std::make_shared<Load>(direction, value, POINTLOAD_CONCENTRATED_DYNAMIC);

                    //Add the nodes that shares this load.
                    theLoad->AddNodes(nodes);
                }
            }
            else if(strcasecmp(pattern.c_str(),"BODY") == 0){
                if (strcasecmp(function.c_str(),"CONSTANT") == 0){
                    //The constant force magnitude
                    value.resize(1);
                    value[0] = jsonFile["Loads"][lTag]["attributes"]["mag"].as<double>();

                    //Creates the constant point load.
                    theLoad = std::make_shared<Load>(direction, value, POINTLOAD_BODY_CONSTANT);

                    //Add the nodes that shares this load.
                    theLoad->AddNodes(nodes);
                }
                else if(strcasecmp(function.c_str(),"TIMESERIES") == 0){
                    //Loads the time history values into memory
                    std::string pathfile = jsonFile["Loads"][lTag]["attributes"]["file"].as<std::string>();
                    std::ifstream load(pathfile.c_str());

                    //The File is Opened and Ready to be Loaded.
                    if(load.is_open()){
                        //Number of time steps.
                        unsigned int nt;
                        load >> nt;

                        //Time-history load vector.
                        value.resize(nt);
                        for(unsigned int j = 0; j < nt; j++)
                            load >> value[j];
                    }
                    load.close();

                    //Creates the timeseries point load.
                    theLoad = std::make_shared<Load>(direction, value, POINTLOAD_BODY_DYNAMIC);

                    //Add the nodes that shares this load.
                    theLoad->AddNodes(nodes);
                }
            }
        }
        else if(strcasecmp(Name.c_str(),"ElementLoad") == 0){
            //The function and pattern applied to these elements
            std::string pattern = jsonFile["Loads"][lTag]["attributes"]["type"].as<std::string>();
            std::string function = jsonFile["Loads"][lTag]["attributes"]["name"].as<std::string>();
                
            //Node identifiers that share this load
            unsigned int n = jsonFile["Loads"][lTag]["attributes"]["list"].size();
            std::vector<unsigned int> elements(n);
            for(unsigned int k = 0; k < n; k++)
                elements[k] = jsonFile["Loads"][lTag]["attributes"]["list"][k].as<int>();

            if(strcasecmp(pattern.c_str(),"SURFACE") == 0){
                //Direction of the applied surface load
                unsigned int ndirs = jsonFile["Loads"][lTag]["attributes"]["dir"].size();
                Eigen::VectorXd direction(ndirs);
                for(unsigned int k = 0; k < ndirs; k++)
                    direction(k) = jsonFile["Loads"][lTag]["attributes"]["dir"][k].as<double>();

                //Face (surfaces) identifiers that share this load
                unsigned int m = jsonFile["Loads"][lTag]["attributes"]["list"].size();
                std::vector<unsigned int> faces(m);
                    
                for(unsigned int k = 0; k < m; k++){
                    std::string Tag = std::to_string(elements[k]);
                    unsigned int sTag = jsonFile["Surfaces"][Tag]["face"].as<int>();
                    unsigned int eTag = jsonFile["Surfaces"][Tag]["element"].as<int>();

                    faces[k] = sTag;
                    elements[k] = eTag;
                }

                //The assigned load pattern function
                std::vector<double> value;
                if (strcasecmp(function.c_str(),"CONSTANT") == 0){
                    //The constant force magnitude
                    value.resize(1);
                    value[0] = jsonFile["Loads"][lTag]["attributes"]["mag"].as<double>();

                    //Creates the constant surface load.
                    theLoad = std::make_shared<Load>(direction, value, ELEMENTLOAD_SURFACE_CONSTANT);

                    //Add the element and face list that shares this load.
                    theLoad->AddFaces(faces);
                    theLoad->AddElements(elements);
                }
                else if(strcasecmp(function.c_str(),"TIMESERIES") == 0){
                    //TODO: Implement Dynamic Element Surface Load.
                }
            }
            else if(strcasecmp(pattern.c_str(),"BODY") == 0){
                //Direction of the applied body load
                unsigned int ndirs = jsonFile["Loads"][lTag]["attributes"]["dir"].size();
                Eigen::VectorXd direction(ndirs);
                for(unsigned int k = 0; k < ndirs; k++)
                    direction(k) = jsonFile["Loads"][lTag]["attributes"]["dir"][k].as<double>();

                //The assigned load pattern function
                std::vector<double> value;
                if (strcasecmp(function.c_str(),"CONSTANT") == 0){
                    //The constant force magnitude
                    value.resize(1);
                    value[0] = jsonFile["Loads"][lTag]["attributes"]["mag"].as<double>();

                    //Creates the constant body load.
                    theLoad = std::make_shared<Load>(direction, value, ELEMENTLOAD_BODY_CONSTANT);

                    //Add the element list that shares this load.
                    theLoad->AddElements(elements);
                }
                else if(strcasecmp(function.c_str(),"TIMESERIES") == 0){
                    //Loads the time history values into memory
                    std::string pathfile = jsonFile["Loads"][lTag]["attributes"]["file"].as<std::string>();
                    std::ifstream load(pathfile.c_str());

                    //The File is Opened and Ready to be Loaded.
                    if(load.is_open()){
                        //Number of time steps.
                        unsigned int nt;
                        load >> nt;

                        //Time-history load vector.
                        value.resize(nt);
                        for(unsigned int j = 0; j < nt; j++)
                            load >> value[j];
                    }
                    load.close();

                    //Creates the variable body load.
                    theLoad = std::make_shared<Load>(direction, value, ELEMENTLOAD_BODY_DYNAMIC);

                    //Add the element list that shares this load.
                    theLoad->AddElements(elements);
                }
            }
            else if(strcasecmp(pattern.c_str(),"GENERALWAVE") == 0){
                //Creates the domain reduction load.
                theLoad = std::make_shared<Load>(ELEMENTLOAD_DOMAIN_REDUCTION);
                theLoad->AddElements(elements);

                //Gets all nodes involved in elements.
                std::map<unsigned int, bool> nodes;
                std::vector<unsigned int> elemIDs = theLoad->GetElements();

                //Mesh node and element information.
                std::map<unsigned int, std::shared_ptr<Node> > theNodes = theMesh->GetNodes();
                std::map<unsigned int, std::shared_ptr<Element> > theElements = theMesh->GetElements();
                        
                for(unsigned int k = 0; k < elemIDs.size(); k++){
                    std::vector<unsigned int> nodeIDs = theElements[elemIDs[k]]->GetNodes();

                    for(unsigned int j = 0; j < nodeIDs.size(); j++)
                        nodes[nodeIDs[j]] = true;
                }

                //Reads/Copy the node information in specified folder.
                std::string pathfile = jsonFile["Loads"][lTag]["attributes"]["file"].as<std::string>();

                for(auto it : nodes){
                    auto &ind = it.first;
                    std::string LoadFile = GetPartitionName(pathfile, ind, false);
                                                
                    //Loads the time history values into memory.
                    std::ifstream load(LoadFile.c_str());

                    //The File is Opened and Ready to be Loaded.
                    if (load.is_open()){
                        //Number of time-steps and field-components.
                        bool cond;
                        unsigned int nt, nFields;
                        load >> nt >> nFields >> cond;

                        //Time-history displacement field.
                        Eigen::MatrixXd Signal(nt,nFields);

                        for(unsigned int i = 0; i < nt; i++){
                            for(unsigned int j = 0; j < nFields; j++)
                                load >> Signal(i,j);
                        }

                        //The node is exterior (change the sign).
                        if(cond)
                            Signal = -1.00*Signal;

                        theLoad->AddDRMCondition(ind,cond);
                        theNodes[ind]->SetDomainReductionMotion(Signal);
                    }
                    load.close();
                }
            }
        }
        else if(strcasecmp(Name.c_str(),"SupportMotion") == 0){
            //Creates the support motion as point load.
            theLoad = std::make_shared<Load>(POINTLOAD_SUPPORT_MOTION);

            unsigned int n = jsonFile["Loads"][lTag]["attributes"]["list"].size();
            std::vector<unsigned int> nodes(n);
            for(unsigned int k = 0; k < n; k++)
                nodes[k] = jsonFile["Loads"][lTag]["attributes"]["list"][k].as<int>();

            theLoad->AddNodes(nodes);
        }      

        //Stores the load in the mesh.
        theMesh->AddLoad(Tag, theLoad); 
    }
}

///Creates a new Analysis object from the json file.
///@param theAnalysis Pointer to the analysis to be  created from json file.
///@param Recorders Pointer to the recorder vector to be parsed.
///@param LoadCombos Pointer to the load combination to be parsed.
///@param InputFile json file where mesh entities will be readden.
bool 
UpdateAnalysis(std::shared_ptr<Mesh> &theMesh, std::unique_ptr<Analysis> &theAnalysis, std::vector<std::shared_ptr<Recorder> > &Recorders, std::map<unsigned int, std::shared_ptr<LoadCombo> > &LoadCombos, std::string InputFile){
    //Gets the partitioned mesh file name.
    std::string file2open = GetPartitionName(InputFile, rank, true);
    file2open = GetSpacedName(file2open, " ");

    //Opens the JSON file
    std::ifstream file2stream(file2open.c_str());

    if(file2stream.is_open()){
        RSJresource jsonFile(file2stream);

        //The pointes to define the simulation
        std::unique_ptr<LinearSystem> theSolver;
        std::shared_ptr<Algorithm> theAlgorithm;
        std::shared_ptr<Integrator> theIntegrator;

        //Creates the solver
        std::string solver = jsonFile["Simulations"]["attributes"]["solver"]["name"].as<std::string>();
        unsigned int Tag = jsonFile["Simulations"]["combo"].as<int>();

        if(strcasecmp(solver.c_str(),"EIGEN") == 0){
            bool update = jsonFile["Simulations"]["attributes"]["solver"]["update"].as<bool>(false);
            theSolver = std::make_unique<EigenSolver>(update);
        }
        else if(strcasecmp(solver.c_str(),"MUMPS") == 0){
            bool update = jsonFile["Simulations"]["attributes"]["solver"]["update"].as<bool>(false);
            unsigned int option = jsonFile["Simulations"]["attributes"]["solver"]["option"].as<int>(2);
            theSolver = std::make_unique<MumpsSolver>(option, update);
        }
        else if(strcasecmp(solver.c_str(),"PETSC") == 0){
            unsigned int dnz = jsonFile["Simulations"]["attributes"]["solver"]["d_nz"].as<int>();
            unsigned int onz = jsonFile["Simulations"]["attributes"]["solver"]["o_nz"].as<int>();
            unsigned int ksp = jsonFile["Simulations"]["attributes"]["solver"]["option"].as<int>(3);
            double tol = jsonFile["Simulations"]["attributes"]["solver"]["tol"].as<double>(1e-08);
            theSolver = std::make_unique<PetscSolver>(dnz, onz, tol, ksp);
        }

        //Creates the algorithm
        std::string algorithm = jsonFile["Simulations"]["attributes"]["algorithm"]["name"].as<std::string>();

        if(strcasecmp(algorithm.c_str(),"LINEAR") == 0){
            theAlgorithm = std::make_shared<Linear>(theSolver, theMesh);
        }
        else if(strcasecmp(algorithm.c_str(),"NEWTON") == 0){
            double ConvergenceTol = jsonFile["Simulations"]["attributes"]["algorithm"]["cnvgtol"].as<double>(1e-06);
            unsigned int nMaxIteration = jsonFile["Simulations"]["attributes"]["algorithm"]["nstep"].as<int>(50);
            unsigned int ConvergenceTest = jsonFile["Simulations"]["attributes"]["algorithm"]["cnvgtest"].as<int>(4);
            theAlgorithm = std::make_shared<NewtonRaphson>(theSolver, theMesh, ConvergenceTol, nMaxIteration, ConvergenceTest);
        }
        //else if(strcasecmp(algorithm.c_str(),"CONTROLPATH") == 0){
        //    theAlgorithm = std::make_shared<GeneralDisplacementPath>(theSolver, theIntegrator, theMesh); 
        //}

        //Creates the integrator
        std::string integrator = jsonFile["Simulations"]["attributes"]["integrator"]["name"].as<std::string>();
        double dt = jsonFile["Simulations"]["attributes"]["integrator"]["dt"].as<double>(0.0);
        double mtol = jsonFile["Simulations"]["attributes"]["integrator"]["mtol"].as<double>(1e-12);
        double ktol = jsonFile["Simulations"]["attributes"]["integrator"]["ktol"].as<double>(1e-12);
        double ftol = jsonFile["Simulations"]["attributes"]["integrator"]["ftol"].as<double>(1e-12);

        if(strcasecmp(integrator.c_str(),"STATIC") == 0){
            theIntegrator = std::make_shared<QuasiStatic>(theMesh, mtol, ktol, ftol);
        }
        else if(strcasecmp(integrator.c_str(),"NEWMARK") == 0){
            theIntegrator = std::make_shared<NewmarkBeta>(theMesh, dt, mtol, ktol, ftol);
        }
        else if(strcasecmp(integrator.c_str(),"BATHE") == 0){
            theIntegrator = std::make_shared<CompositeBathe>(theMesh, dt, mtol, ktol, ftol);
        }
        else if(strcasecmp(integrator.c_str(),"CENTRALDIFFERENCE") == 0){
            theIntegrator = std::make_shared<CentralDifference>(theMesh, dt, mtol, ktol, ftol);
        }
        else if(strcasecmp(integrator.c_str(),"EXTENDEDNEWMARK") == 0){
            theIntegrator = std::make_shared<ExtendedNewmarkBeta>(theMesh, dt, mtol, ktol, ftol);
        }

        //Creates Circular dependency between algorithm and integrator.
        theIntegrator->SetAlgorithm(theAlgorithm);
        theAlgorithm->SetIntegrator(theIntegrator);

        //Creates the analysis case.
        std::string analysis = jsonFile["Simulations"]["attributes"]["analysis"]["name"].as<std::string>();
        unsigned int nt = jsonFile["Simulations"]["attributes"]["analysis"]["nt"].as<int>(1);

        if(strcasecmp(analysis.c_str(),"STATIC") == 0){
            double Factor = 1.00 / (double)nt;
            theAlgorithm->SetLoadFactor(Factor);
            theAnalysis = std::make_unique<StaticAnalysis>(theMesh, theAlgorithm, theIntegrator, LoadCombos[Tag], nt);
        }
        else if(strcasecmp(analysis.c_str(),"DYNAMIC") == 0){
            theAlgorithm->SetLoadFactor(1.00);
            theAnalysis = std::make_unique<DynamicAnalysis>(theMesh, theAlgorithm, theIntegrator, LoadCombos[Tag], nt);
        }

        //Sets the Recorders for the analysis.
        for(unsigned int k = 0; k < Recorders.size(); k++)
            theAnalysis->SetRecorder(Recorders[k]);

        file2stream.close();
    }
    else{
        return true;
    }

    return false;
}

///Updates the Recorder Entities required for the Analysis 
///@param Recorders Pointer to the recorder vector to be parsed.
///@param InputFile json file where mesh entities will be readden.
void 
UpdateRecorders(std::vector<std::shared_ptr<Recorder> > &Recorders, std::string InputFile){
    //Gets the partitioned mesh file name.
    std::string file2open = GetPartitionName(InputFile, rank, true);
    file2open = GetSpacedName(file2open, " ");

    //Opens the JSON file
    std::ifstream file2stream(file2open.c_str());

    if(file2stream.is_open()){
        RSJresource jsonFile(file2stream);

        //Recorder Objects
        for(auto it = jsonFile["Recorders"].as_object().begin(); it != jsonFile["Recorders"].as_object().end(); ++it){
            //Recorder name and identifier.
            std::string name = it->second["name"].as<std::string>();
            std::string file = it->second["file"].as<std::string>();
            unsigned int precision = it->second["ndps"].as<int>();
            unsigned int nsample = it->second["nsamp"].as<int>();

            //Creates a Recorder Object.
            std::shared_ptr<Recorder> theRecorder;

            if(strcasecmp(name.c_str(),"PARAVIEW") == 0){
                unsigned int features = it->second["features"].as<int>();

                //Creates the Paraview Recorder object.
                theRecorder = std::make_shared<Recorder>(file, name, features, nsample, precision);
            }
            else if(strcasecmp(name.c_str(),"SECTION") == 0){
                unsigned int nlist = it->second["list"].size();
                unsigned int ndims = it->second["coords"].size();
                std::string response = it->second["resp"].as<std::string>();

                std::vector<double> coord(ndims);
                for(unsigned int k = 0; k < ndims; k++)
                    coord[k] = it->second["coords"][k].as<double>();

                std::vector<unsigned int> IDs(nlist);
                for(unsigned int k = 0; k < nlist; k++)
                    IDs[k] =  it->second["list"][k].as<int>();

                //Creates the Paraview Recorder object.
                theRecorder = std::make_shared<Recorder>(file, name, response, coord, IDs, nsample, precision);
            }
            else{
                std::string response = it->second["resp"].as<std::string>();
                unsigned int nlist = it->second["list"].size();
                std::vector<unsigned int> IDs(nlist);

                for(unsigned int k = 0; k < nlist; k++)
                    IDs[k] =  it->second["list"][k].as<int>();

                //Creates the Point/Element Recorder object.
                theRecorder = std::make_shared<Recorder>(file, name, response, IDs, nsample, precision);
            }

            //Adds the recorder to the analysis.
            Recorders.push_back(theRecorder);
        }

        file2stream.close();
    }
    else{
        std::cout << "\x1B[31m ERROR: \x1B[0mThe JSON file in \'Driver::UpdateRecorders()\' in Processor [" << rank << "] couldn't be opened. \n";
    }
}

///Updates the Combination Entities required for the Analysis 
///@param LoadCombos Pointer to the load combination to be parsed.
///@param InputFile json file where mesh entities will be readden.
void 
UpdateCombinations(std::map<unsigned int, std::shared_ptr<LoadCombo> > &LoadCombos, std::string InputFile){
    //Gets the partitioned mesh file name.
    std::string file2open = GetPartitionName(InputFile, rank, true);
    file2open = GetSpacedName(file2open, " ");

    //Opens the JSON file
    std::ifstream file2stream(file2open.c_str());

    if(file2stream.is_open()){
        RSJresource jsonFile(file2stream);

        //Combination Objects
        for(auto it = jsonFile["Combinations"].as_object().begin(); it != jsonFile["Combinations"].as_object().end(); ++it){
            //Combination name and identifier.
            unsigned int Tag = std::stoi(it->first);
            std::string Name = it->second["name"].as<std::string>();

            //Load identifiers and combinational factors
            std::vector<double> factors;
            std::vector<unsigned int> loads;

            if( it->second["attributes"]["load"].exists() ){
                unsigned int nCombs = it->second["attributes"]["load"].size();
                loads.resize(nCombs);
                factors.resize(nCombs);

                for(unsigned int k = 0; k < loads.size(); k++){
                    loads[k] = it->second["attributes"]["load"][k].as<int>();
                    factors[k] = it->second["attributes"]["factor"][k].as<double>();
                }
            }

            //Creates the load combination.
            std::shared_ptr<LoadCombo> theCombo = std::make_unique<LoadCombo>(Name, loads, factors);

            //Adds the combination to the analysis.
            LoadCombos[Tag] = theCombo; 
        }

        file2stream.close();
    }
    else{
        std::cout << "\x1B[31m ERROR: \x1B[0mThe JSON file in \'Driver::UpdateCombinations()\' in Processor [" << rank << "] couldn't be opened. \n";
    }
}

///Populate the Mesh object with the json provided entities
///@param theMesh Pointer to the Mesh container.
///@param InputFile json file where mesh entities will be readden.
///@return whether the mesh update was successful or not.
bool
UpdateMesh(std::shared_ptr<Mesh> &theMesh, std::string InputFile){
    //Gets the partitioned mesh file name.
    std::string file2open = GetPartitionName(InputFile, rank, true);
    file2open = GetSpacedName(file2open, " ");

    //Opens the JSON file
    std::ifstream file2stream(file2open.c_str());

    if(file2stream.is_open()){
        RSJresource jsonFile(file2stream);

        //Global Variables
        if( jsonFile["Global"].exists() ){
            nDimensions = jsonFile["Global"]["ndim"].as<int>();
            numberOfTotalDofs = jsonFile["Global"]["ntotal"].as<int>();
            numberOfFreeDofs = jsonFile["Global"]["nfree"].as<int>();
            UpdateOption = jsonFile["Global"]["update"].as<std::string>("RESTART");
            std::string MassForm = jsonFile["Global"]["massform"].as<std::string>("CONSISTENT");

            //Mass formulation for elements mass matrix
            if(strcasecmp(MassForm.c_str(),"LUMPED") == 0){
                MassFormulation = true;
            }
            else if(strcasecmp(MassForm.c_str(),"CONSISTENT") == 0){
                MassFormulation = false;
            }
        }

        //Node Objects
        UpdateNodes(theMesh, jsonFile);
        
        //Mass Objects
        UpdateMasses(theMesh, jsonFile);

        //Constraint Objects
        UpdateConstraints(theMesh, jsonFile);

        //Support Motion Objects
        UpdateSupportMotion(theMesh, jsonFile);

        //Material Objects
        UpdateMaterials(theMesh, jsonFile);

        //Section Objects
        UpdateSections(theMesh, jsonFile);

        //Element Objects
        UpdateElements(theMesh, jsonFile);

        //Damping Objects
        UpdateDampings(theMesh, jsonFile);

        //Load Objects
        UpdateLoads(theMesh, jsonFile);

        //TODO: Include File inside JSON

        file2stream.close();
    }
    else{
        return true;  
    }

    return false;
}

///Runs the User's JSON Input file.
void 
RunDriverFile(){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    if(driverFile){
        //Creates the Mesh subdomain Pointer.
        std::shared_ptr<Mesh> theMesh = std::make_shared<Mesh>();

        //Performs the one simulations after the other.
        for(unsigned int k = 0; k < fileName.size(); k++){
            //Vector of Recorders to store solution. 
            std::vector<std::shared_ptr<Recorder> > Recorders;

            //Load combination cases.
            std::map<unsigned int, std::shared_ptr<LoadCombo> > LoadCombos;

            //Update the Mesh object from json file
            bool FailMesh = UpdateMesh(theMesh, fileName[k]);

            if(!FailMesh){
                //Sets up memory allocation
                theMesh->Initialize();

                //Update the Combination object from json file
                UpdateCombinations(LoadCombos, fileName[k]);

                //Update the Combination object from json file
                UpdateRecorders(Recorders, fileName[k]);

                //Update the Analysis object from json file
                std::unique_ptr<Analysis> theAnalysis;
                bool FailAnalysis = UpdateAnalysis(theMesh, theAnalysis, Recorders, LoadCombos, fileName[k]);

                if(!FailAnalysis){
                    //Runs the current Analysis.
                    bool StopSimulation = theAnalysis->Analyze();

                    if(StopSimulation){
                        if(rank == 0)
                            std::cout << " **A PROBLEM WAS ENCOUNTERED. THE ANALYSIS WILL BE TERMINATED**\n\n";
                        return;
                    }
                }
                else{
                    std::cout << "\x1B[31m ERROR: \x1B[0mThe JSON file \'" << fileName[k]  << "\' in Driver::UpdateAnalysis() in Processor [" << rank << "] couldn't be opened. \n";
                    return;
                }
            }
            else{
                std::cout << "\x1B[31m ERROR: \x1B[0mThe JSON file \'" << fileName[k]  << "\' in Driver::UpdateMesh() in Processor [" << rank << "] couldn't be opened. \n";
                return;
            }
        }
    }
}

#endif
