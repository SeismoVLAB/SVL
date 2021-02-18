#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <iostream>
#include <fstream>

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
#include "lin2DQuad4.hpp"
#include "lin2DQuad8.hpp"
#include "PML2DQuad4.hpp"
#include "kin2DQuad4.hpp"
#include "PML2DQuad8.hpp"
#include "lin3DShell4.hpp"
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

//Compare the transformed string (json tag number) into integer value
bool 
compareStrings(std::string a, std::string b){
    return std::stoi(a) < std::stoi(b);
}

//Sets the partition subdomain tag number.
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

//Fix blank spaces provided by user in path (if any).
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

//Populate the Mesh object with the json provided entities
void
ParseMesh(std::shared_ptr<Mesh> &theMesh, std::vector<std::shared_ptr<Recorder> > &Recorders, std::map<unsigned int, std::shared_ptr<LoadCombo> > &LoadCombos, std::vector<std::string> &Simulations, std::string InputFile){
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
            std::string MassForm = jsonFile["Global"]["mass"].as<std::string>();

            //Mass formulation for elements mass matrix
            if(strcasecmp(MassForm.c_str(),"LUMPED") == 0){
                MassFormulation = true;
            }
            if(strcasecmp(MassForm.c_str(),"CONSISTENT") == 0){
                MassFormulation = false;
            }
        }

        //Node Objects
        for(auto it = jsonFile["Nodes"].as_object().begin(); it != jsonFile["Nodes"].as_object().end(); ++it){
            //Node identifier and number of degree-of-freedom:
            unsigned int Tag = std::stoi(it->first);
            unsigned int nDofs = it->second["ndof"].as<int>();

            //Allocate memory for Free/Total degree of freedom lists.
            std::vector<int> freeDofs(nDofs,0);
            std::vector<int> totalDofs(nDofs,0);

            //Allocate memory for the coordinate vector.
            Eigen::VectorXd coordinates(nDimensions);

            for(int k = 0; k < it->second["freedof"].size(); ++k)
                freeDofs[k] = it->second["freedof"][k].as<int>();
            
            for(int k = 0; k < it->second["totaldof"].size(); ++k)
                totalDofs[k] = it->second["totaldof"][k].as<int>();

            for(int k = 0; k < it->second["coords"].size(); ++k)
                coordinates(k) = it->second["coords"][k].as<double>();

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

        //Mass Objects
        for(auto it = jsonFile["Masses"].as_object().begin(); it != jsonFile["Masses"].as_object().end(); ++it){
            //Node identifier and number of degree-of-freedom:
            unsigned int Tag = std::stoi(it->first);
            unsigned int nDofs = it->second["ndof"].as<int>();

            Eigen::VectorXd Mass(nDofs);
            for(unsigned int k = 0; k < nDofs; ++k)
                Mass(k) = it->second["mass"][k].as<double>();

            //Adds the nodal mass to the node.
            theMesh->AddMass(Tag, Mass);
        }

        //Constraint Objects
        for(auto it = jsonFile["Constraints"].as_object().begin(); it != jsonFile["Constraints"].as_object().end(); ++it){
            //Constraint and slave identifiers 
            int Tag = std::stoi(it->first);
            unsigned int stag = it->second["stag"].as<int>();
            unsigned int nCombs = it->second["mtag"].size();

            //The master identifiers and combinational factors
            std::vector<double> factors(nCombs);
            std::vector<unsigned int> mtag(nCombs);

            for(unsigned int k = 0; k < nCombs; ++k){
                mtag[k] = it->second["mtag"][k].as<int>();
                factors[k] = it->second["factor"][k].as<double>();
            }

            //Creates a constaint object.
            std::unique_ptr<Constraint> theConstraint = std::make_unique<Constraint>(stag, mtag, factors);

            //Stores the information in constraint container.
            theMesh->AddConstraint(Tag, theConstraint);
        }

        //Support Motion Objects
        for(auto it = jsonFile["Supports"].as_object().begin(); it != jsonFile["Supports"].as_object().end(); ++it){
            //Support motion name and identifier.
            unsigned int Tag = std::stoi(it->first);
            std::string Name = it->second["type"].as<std::string>();

            for(int k = 0; k < it->second["dof"].size(); ++k){
                unsigned int dof;
                std::vector<double> Xo;

                if(strcasecmp(Name.c_str(),"CONSTANT") == 0){
                    Xo.resize(1);
                    Xo[0] = it->second["value"][k].as<double>();
                }
                else if(strcasecmp(Name.c_str(),"TIMESERIE") == 0){
                    std::string  File = it->second["file"][k].as<std::string>();
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
                dof = it->second["dof"][k].as<int>();

                //Assign support motion to the node.
                theMesh->SetSupportMotion(Tag, dof, Xo);
            }
        }

        //Material Objects
        for(auto it = jsonFile["Materials"].as_object().begin(); it != jsonFile["Materials"].as_object().end(); ++it){
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

            //TODO: Add more material models here.

            //Stores the information in material container.
            theMesh->AddMaterial(Tag, theMaterial);
        }

        //Section Objects
        for(auto it = jsonFile["Sections"].as_object().begin(); it != jsonFile["Sections"].as_object().end(); ++it){
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
                    double re = it->second["attributes"]["h"].as<double>();
                    double ri = it->second["attributes"]["b"].as<double>();

                    //Instantiate the section object.
                    theSection = std::make_unique<Lin2DCircularTube>(re, ri, theMesh->GetMaterial(matTag), theta, ip);
                }
                else if(strcasecmp(Name.c_str(),"Lin3DCircularTube") == 0){
                    double re = it->second["attributes"]["h"].as<double>();
                    double ri = it->second["attributes"]["b"].as<double>();

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
                else if(strcasecmp(Name.c_str(),"Lin3DThinArea") == 0){
                    double th = it->second["attributes"]["th"].as<double>();

                    //Instantiate the section object.
                    theSection = std::make_unique<Lin3DThinArea>(th, theMesh->GetMaterial(matTag));
                }
            }
            else if(strcasecmp(Model.c_str(),"Fiber") == 0){
                //TODO: Parse fiber sections components.
            }
            else if(strcasecmp(Model.c_str(),"General") == 0){
                //TODO: Parse general sections components.
            }

            //TODO: Add more section models here.

            //Stores the section in the mesh.
            theMesh->AddSection(Tag, theSection);
        }

        //Element Objects
        for(auto it = jsonFile["Elements"].as_object().begin(); it != jsonFile["Elements"].as_object().end(); ++it){
            //Element name and identifier.
            unsigned int Tag = std::stoi(it->first);
            std::string Name = it->second["name"].as<std::string>();

            //Node connectivity
            unsigned int n = it->second["conn"].size();
            std::vector<unsigned int> nodes(n);
            for(unsigned int k = 0; k < n; ++k)
                nodes[k] = it->second["conn"][k].as<int>();

            //Creates an element object.
            std::shared_ptr<Element> theElement;

            if(strcasecmp(Name.c_str(),"lin2DTruss2") == 0){
                unsigned int matID = it->second["attributes"]["material"].as<int>();       
                double parameters = it->second["attributes"]["area"].as<double>();

                //Instantiate the lin2DTruss2 element.
                theElement = std::make_shared<lin2DTruss2>(nodes, theMesh->GetMaterial(matID), parameters);
            }
            else if(strcasecmp(Name.c_str(),"kin2DTruss2") == 0){  
                unsigned int matID = it->second["attributes"]["material"].as<int>();  
                double parameters = it->second["attributes"]["area"].as<double>();

                //Instantiate the kin2DTruss2 element.
                theElement = std::make_shared<kin2DTruss2>(nodes, theMesh->GetMaterial(matID), parameters);
            }
            else if(strcasecmp(Name.c_str(),"lin3DTruss2") == 0){   
                unsigned int matID = it->second["attributes"]["material"].as<int>();     
                double parameters = it->second["attributes"]["area"].as<double>();

                //Instantiate the lin3DTruss2 element.
                theElement = std::make_shared<lin3DTruss2>(nodes, theMesh->GetMaterial(matID), parameters);
            }
            else if(strcasecmp(Name.c_str(),"kin3DTruss2") == 0){
                unsigned int matID = it->second["attributes"]["material"].as<int>();
                double parameters = it->second["attributes"]["area"].as<double>();

                //Instantiate the kin3DTruss2 element.
                theElement = std::make_shared<kin3DTruss2>(nodes, theMesh->GetMaterial(matID), parameters);
            }
            else if(strcasecmp(Name.c_str(),"lin2DTruss3") == 0){      
                double parameters = it->second["attributes"]["area"].as<double>();
                unsigned int nGauss = it->second["attributes"]["np"].as<int>(3);
                unsigned int matID = it->second["attributes"]["material"].as<int>();
                std::string Quadrature = it->second["attributes"]["rule"].as<std::string>("GAUSS");

                //Instantiate the lin2DTruss3 element.
                theElement = std::make_shared<lin2DTruss3>(nodes, theMesh->GetMaterial(matID), parameters, Quadrature, nGauss);
            }
            else if(strcasecmp(Name.c_str(),"lin3DTruss3") == 0){   
                double parameters = it->second["attributes"]["area"].as<double>();
                unsigned int nGauss = it->second["attributes"]["np"].as<int>(3);
                unsigned int matID = it->second["attributes"]["material"].as<int>();
                std::string Quadrature = it->second["attributes"]["rule"].as<std::string>("GAUSS");

                //Instantiate the lin3DTruss3 element.
                theElement = std::make_shared<lin3DTruss3>(nodes, theMesh->GetMaterial(matID), parameters, Quadrature, nGauss);
            }
            else if(strcasecmp(Name.c_str(),"ZeroLength1D") == 0){
                unsigned int dir = it->second["attributes"]["dir"].as<int>();
                unsigned int matID = it->second["attributes"]["material"].as<int>();

                //Instantiate the zerolength1D element.
                theElement = std::make_shared<ZeroLength1D>(nodes, theMesh->GetMaterial(matID), dir);
            }
            else if(strcasecmp(Name.c_str(),"UnxBoucWen2DLink") == 0){
                std::vector<double> variables(2, 0.0);
                std::vector<double> parameters(4, 0.0);
                double tol = it->second["attributes"]["tol"].as<double>(1E-06);
                unsigned int nmax = it->second["attributes"]["nmax"].as<int>(50);
                unsigned int dir = it->second["attributes"]["dir"].as<int>();
                unsigned int dim = it->second["attributes"]["dim"].as<int>();

                variables[0] = it->second["attributes"]["Fy"].as<double>();
                variables[1] = it->second["attributes"]["k0"].as<double>();
                parameters[0] = it->second["attributes"]["alpha"].as<double>(1.0);
                parameters[1] = it->second["attributes"]["eta"].as<double>(1.0);
                parameters[2] = it->second["attributes"]["beta"].as<double>(0.5);
                parameters[3] = it->second["attributes"]["gamma"].as<double>(0.5);

                //Instantiate the UnxBoucWen2DLink element.
                theElement = std::make_shared<UnxBoucWen2DLink>(nodes, parameters, variables, dim, dir, tol, nmax);
            }
            else if(strcasecmp(Name.c_str(),"UnxBoucWen3DLink") == 0){
                std::vector<double> variables(2, 0.0);
                std::vector<double> parameters(4, 0.0);
                double tol = it->second["attributes"]["tol"].as<double>(1E-06);
                unsigned int nmax = it->second["attributes"]["nmax"].as<int>(50);
                unsigned int dir = it->second["attributes"]["dir"].as<int>();
                unsigned int dim = it->second["attributes"]["dim"].as<int>();

                variables[0] = it->second["attributes"]["Fy"].as<double>();
                variables[1] = it->second["attributes"]["k0"].as<double>();
                parameters[0] = it->second["attributes"]["alpha"].as<double>(1.0);
                parameters[1] = it->second["attributes"]["eta"].as<double>(1.0);
                parameters[2] = it->second["attributes"]["beta"].as<double>(0.5);
                parameters[3] = it->second["attributes"]["gamma"].as<double>(0.5);

                //Instantiate the UnxBoucWen3DLink element.
                theElement = std::make_shared<UnxBoucWen3DLink>(nodes, parameters, variables, dim, dir, tol, nmax);
            }
            else if(strcasecmp(Name.c_str(),"HDRBYamamoto2DLink") == 0){
                double De = it->second["attributes"]["De"].as<double>(1.3);
                double Di = it->second["attributes"]["Di"].as<double>(0.3);
                double Hr = it->second["attributes"]["Hr"].as<double>(0.261);
                unsigned int dim = it->second["attributes"]["dim"].as<int>();

                //Instantiate the HDRBYamamoto2DLink element.
                theElement = std::make_shared<HDRBYamamoto2DLink>(nodes, De, Di, Hr, dim);
            }
            else if(strcasecmp(Name.c_str(),"HDRBYamamoto3DLink") == 0){
                double De = it->second["attributes"]["De"].as<double>(1.3);
                double Di = it->second["attributes"]["Di"].as<double>(0.3);
                double Hr = it->second["attributes"]["Hr"].as<double>(0.261);
                unsigned int dim = it->second["attributes"]["dim"].as<int>();

                //Instantiate the HDRBYamamoto3DLink element.
                theElement = std::make_shared<HDRBYamamoto3DLink>(nodes, De, Di, Hr, dim);
            }
            else if(strcasecmp(Name.c_str(),"lin2DFrame2") == 0){
                unsigned int nGauss = it->second["attributes"]["np"].as<int>(3);
                unsigned int secID = it->second["attributes"]["section"].as<int>();
                std::string Quadrature = it->second["attributes"]["rule"].as<std::string>("GAUSS");
                std::string Formulation = it->second["attributes"]["formulation"].as<std::string>("Bernoulli");

                //Frame element formulation.
                bool Condition = false;
                if(strcasecmp(Formulation.c_str(),"Timoshenko") == 0)
                    Condition = true;
                        
                //Instantiate the lin2DFrame2 element.
                theElement = std::make_shared<lin2DFrame2>(nodes, theMesh->GetSection(secID), Condition, Quadrature, nGauss);
            }
            else if(strcasecmp(Name.c_str(),"kin2DFrame2") == 0){
                unsigned int secID = it->second["attributes"]["section"].as<int>();
                        
                //Instantiate the kin2DFrame2 element.
                theElement = std::make_shared<kin2DFrame2>(nodes, theMesh->GetSection(secID));
            }
            else if(strcasecmp(Name.c_str(),"lin3DFrame2") == 0){
                unsigned int nGauss = it->second["attributes"]["np"].as<int>(3);
                unsigned int secID = it->second["attributes"]["section"].as<int>();
                std::string Quadrature = it->second["attributes"]["rule"].as<std::string>("GAUSS");
                std::string Formulation = it->second["attributes"]["formulation"].as<std::string>("Bernoulli");

                //Frame element formulation.
                bool Condition = false;
                if(strcasecmp(Formulation.c_str(),"Timoshenko") == 0)
                    Condition = true;

                //Instantiate the lin3DFrame2 element.
                theElement = std::make_shared<lin3DFrame2>(nodes, theMesh->GetSection(secID), Condition, Quadrature, nGauss);
            }
            else if(strcasecmp(Name.c_str(),"lin2DQuad4") == 0){
                double thickness = it->second["attributes"]["th"].as<double>();
                unsigned int nGauss = it->second["attributes"]["np"].as<int>(4);
                unsigned int matID = it->second["attributes"]["material"].as<int>();
                std::string Quadrature = it->second["attributes"]["rule"].as<std::string>("GAUSS");

                //Instantiate the lin2DQuad4 element.
                theElement = std::make_shared<lin2DQuad4>(nodes, theMesh->GetMaterial(matID), thickness, Quadrature, nGauss);
            }
            else if(strcasecmp(Name.c_str(),"lin2DQuad8") == 0){
                double thickness = it->second["attributes"]["th"].as<double>();
                unsigned int nGauss = it->second["attributes"]["np"].as<int>(9);
                unsigned int matID = it->second["attributes"]["material"].as<int>();
                std::string Quadrature = it->second["attributes"]["rule"].as<std::string>("GAUSS");

                //Instantiate the lin2DQuad8 element.
                theElement = std::make_shared<lin2DQuad8>(nodes, theMesh->GetMaterial(matID), thickness, Quadrature, nGauss);
            }
            else if(strcasecmp(Name.c_str(),"PML2DQuad4") == 0){
                std::vector<double> parameters(8);
                parameters[0] = it->second["attributes"]["th"].as<double>();
                parameters[1] = it->second["attributes"]["n"].as<double>();
                parameters[2] = it->second["attributes"]["L"].as<double>();
                parameters[3] = it->second["attributes"]["R"].as<double>();
                parameters[4] = it->second["attributes"]["x0"][0].as<double>();
                parameters[5] = it->second["attributes"]["x0"][1].as<double>();
                parameters[6] = it->second["attributes"]["npml"][0].as<double>();
                parameters[7] = it->second["attributes"]["npml"][1].as<double>();

                unsigned int nGauss = it->second["attributes"]["np"].as<int>(4);
                unsigned int matID = it->second["attributes"]["material"].as<int>();
                std::string Quadrature = it->second["attributes"]["rule"].as<std::string>("GAUSS");

                //Instantiate the PML2DQuad4 element.
                theElement = std::make_shared<PML2DQuad4>(nodes, theMesh->GetMaterial(matID), parameters, Quadrature, nGauss);
            }
            else if(strcasecmp(Name.c_str(),"PML2DQuad8") == 0){
                std::vector<double> parameters(8);
                parameters[0] = it->second["attributes"]["th"].as<double>();
                parameters[1] = it->second["attributes"]["n"].as<double>();
                parameters[2] = it->second["attributes"]["L"].as<double>();
                parameters[3] = it->second["attributes"]["R"].as<double>();
                parameters[4] = it->second["attributes"]["x0"][0].as<double>();
                parameters[5] = it->second["attributes"]["x0"][1].as<double>();
                parameters[6] = it->second["attributes"]["npml"][0].as<double>();
                parameters[7] = it->second["attributes"]["npml"][1].as<double>();

                unsigned int nGauss = it->second["attributes"]["np"].as<int>(9);
                unsigned int matID = it->second["attributes"]["material"].as<int>();
                std::string Quadrature = it->second["attributes"]["rule"].as<std::string>("GAUSS");

                //Instantiate the PML2DQuad4 element.
                theElement = std::make_shared<PML2DQuad8>(nodes, theMesh->GetMaterial(matID), parameters, Quadrature, nGauss);
            }
            else if(strcasecmp(Name.c_str(),"kin2DQuad4") == 0){
                double thickness = it->second["attributes"]["th"].as<double>();
                unsigned int nGauss = it->second["attributes"]["np"].as<int>(4);
                unsigned int matID = it->second["attributes"]["material"].as<int>();
                std::string Quadrature = it->second["attributes"]["rule"].as<std::string>("GAUSS");

                //Instantiate the kin2DQuad4 element.
                theElement = std::make_shared<kin2DQuad4>(nodes, theMesh->GetMaterial(matID), thickness, Quadrature, nGauss);
            }
            else if(strcasecmp(Name.c_str(),"lin3DShell4") == 0){
                unsigned int nGauss = it->second["attributes"]["np"].as<int>(9);
                unsigned int secID = it->second["attributes"]["section"].as<int>();
                std::string Quadrature = it->second["attributes"]["rule"].as<std::string>("GAUSS");
                    
                //Instantiate the lin3DShell4 element.
                theElement = std::make_shared<lin3DShell4>(nodes, theMesh->GetSection(secID), Quadrature, nGauss);
            }
            else if(strcasecmp(Name.c_str(),"lin3DHexa8") == 0){
                unsigned int nGauss = it->second["attributes"]["np"].as<int>(8);
                unsigned int matID = it->second["attributes"]["material"].as<int>();
                std::string Quadrature = it->second["attributes"]["rule"].as<std::string>("GAUSS");

                //Instantiate the lin3DHexa8 element.
                theElement = std::make_shared<lin3DHexa8>(nodes, theMesh->GetMaterial(matID), Quadrature, nGauss);
            }
            else if(strcasecmp(Name.c_str(),"kin3DHexa8") == 0){
                unsigned int nGauss = it->second["attributes"]["np"].as<int>(8);
                unsigned int matID = it->second["attributes"]["material"].as<int>();
                std::string Quadrature = it->second["attributes"]["rule"].as<std::string>("GAUSS");

                //Instantiate the kin3DHexa8 element.
                theElement = std::make_shared<kin3DHexa8>(nodes, theMesh->GetMaterial(matID), Quadrature, nGauss);
            }
            else if(strcasecmp(Name.c_str(),"PML3DHexa8") == 0){
                std::vector<double> parameters(9);
                parameters[0] = it->second["attributes"]["n"].as<double>();
                parameters[1] = it->second["attributes"]["L"].as<double>();
                parameters[2] = it->second["attributes"]["R"].as<double>();
                parameters[3] = it->second["attributes"]["x0"][0].as<double>();
                parameters[4] = it->second["attributes"]["x0"][1].as<double>();
                parameters[5] = it->second["attributes"]["x0"][2].as<double>();
                parameters[6] = it->second["attributes"]["npml"][0].as<double>();
                parameters[7] = it->second["attributes"]["npml"][1].as<double>();
                parameters[8] = it->second["attributes"]["npml"][2].as<double>();

                unsigned int nGauss = it->second["attributes"]["np"].as<int>(8);
                unsigned int matID = it->second["attributes"]["material"].as<int>();
                std::string Quadrature = it->second["attributes"]["rule"].as<std::string>("GAUSS");

                //Instantiate the PML3DHexa8 element.
                theElement = std::make_shared<PML3DHexa8>(nodes, theMesh->GetMaterial(matID), parameters, Quadrature, nGauss);
            }
            else if(strcasecmp(Name.c_str(),"lin3DHexa20") == 0){
                unsigned int nGauss = it->second["attributes"]["np"].as<int>(27);
                unsigned int matID = it->second["attributes"]["material"].as<int>();
                std::string Quadrature = it->second["attributes"]["rule"].as<std::string>("GAUSS");

                //Instantiate the lin3DHexa20 element.
                theElement = std::make_shared<lin3DHexa20>(nodes, theMesh->GetMaterial(matID), Quadrature, nGauss);
            }
            else if(strcasecmp(Name.c_str(),"TIEQlin2DQuad4") == 0){
                double cf1 = it->second["attributes"]["cf1"].as<double>();
                double cf2 = it->second["attributes"]["cf2"].as<double>();
                double zref = it->second["attributes"]["zref"].as<double>();
                double eref = it->second["attributes"]["eref"].as<double>();
                double thickness = it->second["attributes"]["th"].as<double>();
                unsigned int matID = it->second["attributes"]["material"].as<int>();
                std::string eType = it->second["attributes"]["type"].as<std::string>("Darandelli");
                unsigned int nGauss = it->second["attributes"]["np"].as<int>(4);
                std::string Quadrature = it->second["attributes"]["rule"].as<std::string>("GAUSS");

                //Instantiate the TIEQlin2DQuad4 element.
                theElement = std::make_shared<TIEQlin2DQuad4>(nodes, theMesh->GetMaterial(matID), thickness, Quadrature, nGauss, eType, zref, cf1, cf2, eref);
            }
            else if(strcasecmp(Name.c_str(),"EQlin2DQuad4") == 0){
                double cf1 = it->second["attributes"]["cf1"].as<double>();
                double cf2 = it->second["attributes"]["cf2"].as<double>();
                double zref = it->second["attributes"]["zref"].as<double>();
                double thickness = it->second["attributes"]["th"].as<double>();
                unsigned int matID = it->second["attributes"]["material"].as<int>();
                std::string eType = it->second["attributes"]["type"].as<std::string>("Darandelli");
                unsigned int nGauss = it->second["attributes"]["np"].as<int>(4);
                std::string Quadrature = it->second["attributes"]["rule"].as<std::string>("GAUSS");

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

        //Damping Objects
        for(auto it = jsonFile["Dampings"].as_object().begin(); it != jsonFile["Dampings"].as_object().end(); ++it){
            //Damping name and identifier.
            unsigned int Tag = std::stoi(it->first);
            std::string Name = it->second["name"].as<std::string>();

            //Damping parameters
            std::vector<double> parameters;
            if(strcasecmp(Name.c_str(),"Free") == 0){
                parameters.resize(2);
                parameters[0] = 0.0;
                parameters[1] = 0.0;
            }
            else if(strcasecmp(Name.c_str(),"Rayleigh") == 0){
                parameters.resize(2);
                parameters[0] = it->second["attributes"]["am"].as<double>(0.0);
                parameters[1] = it->second["attributes"]["ak"].as<double>(0.0);
            }
            else if(strcasecmp(Name.c_str(),"Caughey") == 0){
                //TODO: Caughey-viscous damping is not implemented yet.
                //parameters.resize(4);
            }
            else if(strcasecmp(Name.c_str(),"Capped") == 0){
                //TODO: Capped-viscous damping is not implemented yet.
            }

            //Damping associated to elements
            unsigned int n = it->second["attributes"]["list"].size();
            std::vector<unsigned int> eList(n);

            for(unsigned int k = 0; k < n; k++)
                eList[k] = it->second["attributes"]["list"][k].as<int>();
            
            //Creates a damping object.
            std::shared_ptr<Damping> theDamping = std::make_shared<Damping>(Name, parameters);

            //Stores the damping in the mesh.
            theMesh->AddDamping(Tag, theDamping);

            //Assign damping to group of elements.
            theMesh->SetDamping(Tag, eList);
        }

        //Load Objects
        for(auto it = jsonFile["Loads"].as_object().begin(); it != jsonFile["Loads"].as_object().end(); ++it){
            //Load name and identifier.
            unsigned int Tag = std::stoi(it->first);
            std::string Name = it->second["name"].as<std::string>();

            //Creates a load object.
            std::shared_ptr<Load> theLoad;

            if(strcasecmp(Name.c_str(),"PointLoad") == 0){
                //Node identifiers that share this load
                unsigned int n = it->second["attributes"]["list"].size();
                std::vector<unsigned int> nodes(n);
                for(unsigned int k = 0; k < n; k++)
                    nodes[k] = it->second["attributes"]["list"][k].as<int>();

                //Direction of the applied load
                unsigned int ndirs = it->second["attributes"]["dir"].size();
                Eigen::VectorXd direction(ndirs);
                for(unsigned int k = 0; k < ndirs; k++)
                    direction(k) = it->second["attributes"]["dir"][k].as<double>();

                //The assigned load pattern function
                std::vector<double> value;
                std::string function = it->second["attributes"]["name"].as<std::string>();
                if (strcasecmp(function.c_str(),"CONSTANT") == 0){
                    //The constant force magnitude
                    value.resize(1);
                    value[0] = it->second["attributes"]["mag"].as<double>();

                    //Creates the point load.
                    theLoad = std::make_shared<Load>(direction, value, 1);

                    //Add the nodes that shares this load.
                    theLoad->AddNodes(nodes);
                }
                else if(strcasecmp(function.c_str(),"TIMESERIE") == 0){
                    //Loads the time history values into memory
                    std::string pathfile = it->second["attributes"]["file"].as<std::string>();
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

                    //Creates the point load.
                    theLoad = std::make_shared<Load>(direction, value, 2);

                    //Add the nodes that shares this load.
                    theLoad->AddNodes(nodes);
                }
            }
            else if(strcasecmp(Name.c_str(),"ElementLoad") == 0){
                //The function and pattern applied to these elements
                std::string pattern = it->second["attributes"]["type"].as<std::string>();
                std::string function = it->second["attributes"]["name"].as<std::string>();
                
                //Node identifiers that share this load
                unsigned int n = it->second["attributes"]["list"].size();
                std::vector<unsigned int> elements(n);
                for(unsigned int k = 0; k < n; k++)
                    elements[k] = it->second["attributes"]["list"][k].as<int>();

                if(strcasecmp(pattern.c_str(),"SURFACE") == 0){
                    //Direction of the applied surface load
                    unsigned int ndirs = it->second["attributes"]["dir"].size();
                    Eigen::VectorXd direction(ndirs);
                    for(unsigned int k = 0; k < ndirs; k++)
                        direction(k) = it->second["attributes"]["dir"][k].as<double>();

                    //Face (surfaces) identifiers that share this load
                    unsigned int m = it->second["attributes"]["list"].size();
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
                        value[0] = it->second["attributes"]["mag"].as<double>();

                        //Creates the constant surface load.
                        theLoad = std::make_shared<Load>(direction, value, 3);

                        //Add the element and face list that shares this load.
                        theLoad->AddFaces(faces);
                        theLoad->AddElements(elements);
                    }
                    else if(strcasecmp(function.c_str(),"TIMESERIE") == 0){
                        //TODO: Implement Dynamic Element Surface Load.
                    }
                }
                else if(strcasecmp(pattern.c_str(),"BODY") == 0){
                    //Direction of the applied body load
                    unsigned int ndirs = it->second["attributes"]["dir"].size();
                    Eigen::VectorXd direction(ndirs);
                    for(unsigned int k = 0; k < ndirs; k++)
                        direction(k) = it->second["attributes"]["dir"][k].as<double>();

                    //The assigned load pattern function
                    std::vector<double> value;
                    if (strcasecmp(function.c_str(),"CONSTANT") == 0){
                        //The constant force magnitude
                        value.resize(1);
                        value[0] = it->second["attributes"]["mag"].as<double>();

                        //Creates the constant body load.
                        theLoad = std::make_shared<Load>(direction, value, 4);

                        //Add the element list that shares this load.
                        theLoad->AddElements(elements);
                    }
                    else if(strcasecmp(function.c_str(),"TIMESERIE") == 0){
                        //Loads the time history values into memory
                        std::string pathfile = it->second["attributes"]["file"].as<std::string>();
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
                        theLoad = std::make_shared<Load>(direction, value, 5);

                        //Add the element list that shares this load.
                        theLoad->AddElements(elements);
                    }
                }
                else if(strcasecmp(pattern.c_str(),"GENERALWAVE") == 0){
                    //Creates the domain reduction load.
                    theLoad = std::make_shared<Load>(7);
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
                    std::string pathfile = it->second["attributes"]["file"].as<std::string>();

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
                theLoad = std::make_shared<Load>(8);

                unsigned int n = it->second["attributes"]["list"].size();
                std::vector<unsigned int> nodes(n);
                for(unsigned int k = 0; k < n; k++)
                    nodes[k] = it->second["attributes"]["list"][k].as<int>();

                theLoad->AddNodes(nodes);
            }      

            //Stores the load in the mesh.
            theMesh->AddLoad(Tag, theLoad); 
        }

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

        //Recorder Objects
        for(auto it = jsonFile["Recorders"].as_object().begin(); it != jsonFile["Recorders"].as_object().end(); ++it){
            //Recorder name and identifier.
            std::string name = it->second["name"].as<std::string>();

            std::string file = it->second["file"].as<std::string>();
            unsigned int presicion = it->second["ndps"].as<int>();
            unsigned int nsample = it->second["nsamp"].as<int>();

            //Creates a Recorder Object.
            std::shared_ptr<Recorder> theRecorder;

            if(strcasecmp(name.c_str(),"PARAVIEW") == 0){
                unsigned int features = it->second["features"].as<int>();

                //Creates the Paraview Recorder object.
                theRecorder = std::make_shared<Recorder>(file, name, features, nsample, presicion);
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
                theRecorder = std::make_shared<Recorder>(file, name, response, coord, IDs, nsample, presicion);
            }
            else{
                std::string response = it->second["resp"].as<std::string>();
                unsigned int nlist = it->second["list"].size();
                std::vector<unsigned int> IDs(nlist);

                for(unsigned int k = 0; k < nlist; k++)
                    IDs[k] =  it->second["list"][k].as<int>();

                //Creates the Point/Element Recorder object.
                theRecorder = std::make_shared<Recorder>(file, name, response, IDs, nsample, presicion);
            }

            //Adds the recorder to the analysis.
            Recorders.push_back(theRecorder);
        }

        //TODO: Include File
        //if( jsonFile["Include"].exists() ){
            //std::string JsonFileName = jsonFile["Include"]["file"]
            //std::string fileInclude = filePath + "/" + JsonFileName;
            //ParseMesh(theMesh, Recorders, LoadCombos, Simulations, fileInclude);
        //}

        //Simulation Objects List
        for(auto it = jsonFile["Simulations"].as_object().begin(); it != jsonFile["Simulations"].as_object().end(); ++it)
            Simulations.push_back(it->first);
        std::sort(Simulations.begin(), Simulations.end(), compareStrings);

        file2stream.close();
    }
    else{
        std::cout << "\x1B[31mERROR: \x1B[0mThe JSON file in Processor [" << rank << "] couldn't be opened. \n";  
    }
}

//Populate the Analysis object with the json provided entities.
void
ParseAnalysis(std::shared_ptr<Mesh> &theMesh, std::unique_ptr<Analysis> &theAnalysis, std::vector<std::shared_ptr<Recorder> > &Recorders, std::map<unsigned int, std::shared_ptr<LoadCombo> > &LoadCombos, std::string sim, std::string InputFile){
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
        std::string solver = jsonFile["Simulations"][sim]["attributes"]["solver"]["name"].as<std::string>();
        unsigned int Tag = jsonFile["Simulations"][sim]["combo"].as<int>();

        if(strcasecmp(solver.c_str(),"EIGEN") == 0){
            bool update = jsonFile["Simulations"][sim]["attributes"]["solver"]["update"].as<bool>(false);
            theSolver = std::make_unique<EigenSolver>(update);
        }
        else if(strcasecmp(solver.c_str(),"MUMPS") == 0){
            bool update = jsonFile["Simulations"][sim]["attributes"]["solver"]["update"].as<bool>(false);
            unsigned int option = jsonFile["Simulations"][sim]["attributes"]["solver"]["option"].as<int>(2);
            theSolver = std::make_unique<MumpsSolver>(option, update);
        }
        else if(strcasecmp(solver.c_str(),"PETSC") == 0){
            unsigned int dnz = jsonFile["Simulations"][sim]["attributes"]["solver"]["d_nz"].as<int>();
            unsigned int onz = jsonFile["Simulations"][sim]["attributes"]["solver"]["o_nz"].as<int>();
            unsigned int ksp = jsonFile["Simulations"][sim]["attributes"]["solver"]["option"].as<int>(3);
            double tol = jsonFile["Simulations"][sim]["attributes"]["solver"]["tol"].as<double>(1e-08);
            theSolver = std::make_unique<PetscSolver>(dnz, onz, tol, ksp);
        }

        //Creates the algorithm
        std::string algorithm = jsonFile["Simulations"][sim]["attributes"]["algorithm"]["name"].as<std::string>();

        if(strcasecmp(algorithm.c_str(),"LINEAR") == 0){
            theAlgorithm = std::make_shared<Linear>(theSolver, theMesh);
        }
        else if(strcasecmp(algorithm.c_str(),"NEWTON") == 0){
            double ConvergenceTol = jsonFile["Simulations"][sim]["attributes"]["algorithm"]["cnvgtol"].as<double>(1e-06);
            unsigned int nMaxIteration = jsonFile["Simulations"][sim]["attributes"]["algorithm"]["nstep"].as<int>(50);
            unsigned int ConvergenceTest = jsonFile["Simulations"][sim]["attributes"]["algorithm"]["cnvgtest"].as<int>(4);
            theAlgorithm = std::make_shared<NewtonRaphson>(theSolver, theMesh, ConvergenceTol, nMaxIteration, ConvergenceTest);
        }
        //else if(strcasecmp(algorithm.c_str(),"CONTROLPATH") == 0){
        //    theAlgorithm = std::make_shared<GeneralDisplacementPath>(theSolver, theIntegrator, theMesh); 
        //}

        //Creates the integrator
        std::string integrator = jsonFile["Simulations"][sim]["attributes"]["integrator"]["name"].as<std::string>();
        double dt = jsonFile["Simulations"][sim]["attributes"]["integrator"]["dt"].as<double>(0.0);
        double mtol = jsonFile["Simulations"][sim]["attributes"]["integrator"]["mtol"].as<double>(1e-12);
        double ktol = jsonFile["Simulations"][sim]["attributes"]["integrator"]["ktol"].as<double>(1e-12);
        double ftol = jsonFile["Simulations"][sim]["attributes"]["integrator"]["ftol"].as<double>(1e-12);

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
        std::string analysis = jsonFile["Simulations"][sim]["attributes"]["analysis"]["name"].as<std::string>();
        unsigned int nt = jsonFile["Simulations"][sim]["attributes"]["analysis"]["nt"].as<int>(1);

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
    }
}

//Perform modification on the Mesh object as specified in json file
void 
PerformChanges(std::shared_ptr<Mesh>& UNUSED(theMesh), std::string UNUSED(InputFile)){
    //TODO: Parse instruction in JSON to modify mesh
}

//Runs the User's JSON Input file.
void 
RunFromJSON(){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Number of Simulations to be performed
    std::vector<std::string> Simulations;

    //Vector of Recorders to store solution. 
    std::vector<std::shared_ptr<Recorder> > Recorders;

    //Load combination cases.
    std::map<unsigned int, std::shared_ptr<LoadCombo> > LoadCombos;

    //Creates the Mesh subdomain Pointer.
    std::shared_ptr<Mesh> theMesh = std::make_shared<Mesh>();
    ParseMesh(theMesh, Recorders, LoadCombos, Simulations, fileName);

    //Performs the simulations sequentially.
    for(unsigned int sim = 0; sim < Simulations.size(); sim++){
        //Sets up memory allocation
        theMesh->Initialize();

        //Creates Analysis Pointer.
        std::unique_ptr<Analysis> theAnalysis;
        ParseAnalysis(theMesh, theAnalysis, Recorders, LoadCombos, Simulations[sim], fileName);

        //Runs the current Analysis.
        bool Stop = theAnalysis->Analyze();
        if(Stop){
            if(rank == 0)
                std::cout << " **A PROBLEM WAS ENCOUNTERED. THE ANALYSIS WILL BE TERMINATED** \n";
        }

        //Change the Mesh subdomain for next analysis
        PerformChanges(theMesh, fileName);
    }
}
