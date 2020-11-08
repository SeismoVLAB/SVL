#include <cstdlib>

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

#include "Parser.hpp"
#include "Definitions.hpp"
#include "Profiler.hpp"

//Overload Constructor.
Parser::Parser(unsigned int np, std::string folder, std::string file) : 
MassForm(false), Flag(true), nProcessors(np), Folder(folder), File(file){
    //Does nothing.
}

//Destructor.
Parser::~Parser(){
    //does nothing.
}

//Sets the partition subdomain tag.
std::string 
Parser::GetPartitionName(std::string theFile, int k, bool cond){
    //Auxiliar variables.
    std::string dumm0 = theFile;
    std::string subDomainMesh;

    //Conver process number into string.
    std::stringstream process;
    process << k;

    //Replace strings.
    size_t pos;
    pos = dumm0.find("$");
    if(pos !=std::string::npos)
        dumm0.replace(pos, std::string("$").length(), process.str());

    //Modified File Name.
    if(cond){
        subDomainMesh = Folder + "/" + dumm0;
    }
    else{
        subDomainMesh = dumm0;
    }

    return subDomainMesh;
}

//Fix blank spaces provided by user in path.
std::string 
Parser::GetSpacedName(std::string theFile, std::string toReplace){
    //Auxiliar variables.
    std::string dumm0 = theFile;
    std::string subDomainMesh;

    //Replace strings.
    size_t pos;
    pos = dumm0.find("~");
    while(pos != std::string::npos){
        dumm0.replace(pos, std::string("~").length(), toReplace);
        pos = dumm0.find("~");
    }

    //Modified File Name.
    subDomainMesh = dumm0;

    return subDomainMesh;
}

//Obtains a Point from Input File.
unsigned int 
Parser::CreatePoint(std::ifstream& InputFile, std::shared_ptr<Node> &theNode){
    //Auxiliar variable.
    unsigned int Tag, nDofs;

    //Parser the node object.
    InputFile >> Tag >> nDofs;

    //Allocate memory for the coordinate vector.
    Eigen::VectorXd coordinates;            
    coordinates.resize(nDimensions);

    //Allocate memory for Free/Total degree of freedom lists.
    std::vector<int> freeDofs(nDofs,0);
    std::vector<int> totalDofs(nDofs,0);

    for(unsigned int k = 0; k < nDofs; k++)
        InputFile >> totalDofs[k];

    for(unsigned int k = 0; k < nDofs; k++)
        InputFile >> freeDofs[k];

    for(unsigned int k = 0; k < nDimensions; k++)
        InputFile >> coordinates(k);

    //Checks if the node is free/fixed.
    bool IsFixed = false;
    if(std::find(freeDofs.begin(), freeDofs.end(), -1) != freeDofs.end())
        IsFixed = true;

    //Creates a node object.
    theNode = std::make_shared<Node>(nDofs, coordinates, IsFixed);

    //Sets preliminar degree-of-freedom numbering.
    theNode->SetFreeDegreeOfFreedom(freeDofs);
    theNode->SetTotalDegreeOfFreedom(totalDofs);        

    return Tag;
}

//Obtains a Material from Input File.
unsigned int 
Parser::CreateMaterial(std::ifstream& InputFile, std::unique_ptr<Material> &theMaterial){
    //Auxiliar variables.
    unsigned int Tag;
    std::string MatName;

    //Parser the material model.
    InputFile >> Tag >> MatName;

    if(strcasecmp(MatName.c_str(),"Elastic1DLinear") == 0){
        double E, nu, rho;
        InputFile >> E >> nu >> rho;

        //Instantiate the material object.
        theMaterial = std::make_unique<Elastic1DLinear>(E, nu, rho);
    }
    else if(strcasecmp(MatName.c_str(),"Hertzian1DLinear") == 0){
        double k1, k2, k3, rho;
        InputFile >> k1 >> k2 >> k3 >> rho;

        //Instantiate the material object.
        theMaterial = std::make_unique<Hertzian1DLinear>(k1, k2, k3, rho);
    }
    else if(strcasecmp(MatName.c_str(),"Viscous1DLinear") == 0){
        double eta;
        InputFile >> eta;

        //Instantiate the material object.
        theMaterial =std::make_unique<Viscous1DLinear>(eta);
    }
    else if(strcasecmp(MatName.c_str(),"Elastic2DPlaneStrain") == 0){
        double E, nu, rho;
        InputFile >> E >> nu >> rho;

        //Instantiate the material object.
        theMaterial = std::make_unique<Elastic2DPlaneStrain>(E, nu, rho);
    }
    else if(strcasecmp(MatName.c_str(),"Elastic2DPlaneStress") == 0){
        double E, nu, rho;
        InputFile >> E >> nu >> rho;

        //Instantiate the material object.
        theMaterial = std::make_unique<Elastic2DPlaneStress>(E, nu, rho);
    }
    else if(strcasecmp(MatName.c_str(),"Elastic3DLinear") == 0){
        double E, nu, rho;
        InputFile >> E >> nu >> rho;

        //Instantiate the material object.
        theMaterial = std::make_unique<Elastic3DLinear>(E, nu, rho);
    }
    else if(strcasecmp(MatName.c_str(),"Plastic1DJ2") == 0){
        double E, nu, rho, K,H, SigmaY;
        InputFile >> E >> nu >> rho >> K >> H >> SigmaY;
                
        //Instantiate the material object.
        theMaterial = std::make_unique<Plastic1DJ2>(E, nu, rho, K, H, SigmaY);
    }
    else if(strcasecmp(MatName.c_str(),"PlasticPlaneStrainJ2") == 0){
        double K, G, H, beta, rho, SigmaY;
        InputFile >> K >> G >> rho >> H >> beta >> SigmaY;

        //Instantiate the material object.
        theMaterial = std::make_unique<PlasticPlaneStrainJ2>(K, G, rho, H, beta, SigmaY);
    }
    else if(strcasecmp(MatName.c_str(),"PlasticPlaneStrainBA") == 0){
        double K, G, rho, H0, h, m, Su, beta;
        InputFile >> K >> G >> rho >> H0 >> h >> m >> Su >> beta;

        //Instantiate the material object.
        theMaterial = std::make_unique<PlasticPlaneStrainBA>(K, G, rho, H0, h, m, Su, beta);
    }
    else if(strcasecmp(MatName.c_str(),"Plastic3DJ2") == 0){
        double K, G, H, beta, rho, SigmaY;
        InputFile >> K >> G >> rho >> H >> beta >> SigmaY;

        //Instantiate the material object.
        theMaterial = std::make_unique<Plastic3DJ2>(K, G, rho, H, beta, SigmaY);
    }
    else if(strcasecmp(MatName.c_str(),"Plastic3DBA") == 0){
        double K, G, rho, H0, h, m, Su, beta;
        InputFile >> K >> G >> rho >> H0 >> h >> m >> Su >> beta;

        //Instantiate the material object.
        theMaterial = std::make_unique<Plastic3DBA>(K, G, rho, H0, h, m, Su, beta);
    }

    //TODO: Add more material models here.

    return Tag;
}

//Obtains a Section from Input File.
unsigned int 
Parser::CreateSection(std::ifstream& InputFile, std::shared_ptr<Mesh> &theMesh, std::unique_ptr<Section> &theSection){
    //Auxiliar variables.
    unsigned int Tag, matTag;
    std::string secName, secType, secShape;

    //Parser the section model.
    InputFile >> secType >> Tag;

    if(strcasecmp(secType.c_str(),"Plain") == 0){
        InputFile >> secType >> matTag >> secShape;
        if(strcasecmp(secShape.c_str(),"Rectangular") == 0){
            unsigned int ip=10;
            double h, b, theta=0.0;
            InputFile >> h >> b >> theta >> ip; 

            //Instantiate the section object.
            if(nDimensions == 2)
                theSection = std::make_unique<Lin2DRectangular>(h, b, theMesh->GetMaterial(matTag), theta, ip);
            else if(nDimensions == 3)
                theSection = std::make_unique<Lin3DRectangular>(h, b, theMesh->GetMaterial(matTag), theta, ip);
        }
        else if(strcasecmp(secShape.c_str(),"RectangularTube") == 0){
            unsigned int ip=10;
            double h, b, tw, tf, theta=0.0;
            InputFile >> h >> b >> tw >> tf >> theta >> ip;

            //Instantiate the section object.
            if(nDimensions == 2)
                theSection = std::make_unique<Lin2DRectangularTube>(h, b, tw, tf, theMesh->GetMaterial(matTag), theta, ip);
            else if(nDimensions == 3)
                theSection = std::make_unique<Lin3DRectangularTube>(h, b, tw, tf, theMesh->GetMaterial(matTag), theta, ip);
        }
        else if(strcasecmp(secShape.c_str(),"Circular") == 0){
            unsigned int ip=10;
            double r, theta=0.0;
            InputFile >> r >> theta >> ip;

            //Instantiate the section object.
            if(nDimensions == 2)
                theSection = std::make_unique<Lin2DCircular>(r, theMesh->GetMaterial(matTag), theta, ip);
            else if(nDimensions == 3)
                theSection = std::make_unique<Lin3DCircular>(r, theMesh->GetMaterial(matTag), theta, ip);
        }

        else if(strcasecmp(secShape.c_str(),"CircularTube") == 0){
            unsigned int ip=10;
            double re, ri, theta=0.0; 
            InputFile >> re >> ri >> theta >> ip;

            //Instantiate the section object.
            if(nDimensions == 2)
                theSection = std::make_unique<Lin2DCircularTube>(re, ri, theMesh->GetMaterial(matTag), theta, ip);
            else if(nDimensions == 3)
                theSection = std::make_unique<Lin3DCircularTube>(re, ri, theMesh->GetMaterial(matTag), theta, ip);
        }
        else if(strcasecmp(secShape.c_str(),"Angle") == 0){
            unsigned int ip=10;
            double h, b, tw, tf, theta=0.0; 
            InputFile >> h >> b >> tw >> tf >> theta >> ip;

            //Instantiate the section object.
            if(nDimensions == 2)
                theSection = std::make_unique<Lin2DAngle>(h, b, tw, tf, theMesh->GetMaterial(matTag), theta, ip);
            else if(nDimensions == 3)
                theSection = std::make_unique<Lin3DAngle>(h, b, tw, tf, theMesh->GetMaterial(matTag), theta, ip);
        }
        else if(strcasecmp(secShape.c_str(),"Channel") == 0){
            unsigned int ip=10;
            double h, b, tw, tf, theta=0.0; 
            InputFile >> h >> b >> tw >> tf >> theta >> ip;

            //Instantiate the section object.
            if(nDimensions == 2)
                theSection = std::make_unique<Lin2DChannel>(h, b, tw, tf, theMesh->GetMaterial(matTag), theta, ip);
            else if(nDimensions == 3)
                theSection = std::make_unique<Lin3DChannel>(h, b, tw, tf, theMesh->GetMaterial(matTag), theta, ip);
        }
        if(strcasecmp(secShape.c_str(),"Tee") == 0){
            unsigned int ip=10;
            double h, b, tw, tf, theta=0.0; 
            InputFile >> h >> b >> tw >> tf >> theta >> ip;

            //Instantiate the section object.
            if(nDimensions == 2)
                theSection = std::make_unique<Lin2DTee>(h, b, tw, tf, theMesh->GetMaterial(matTag), theta, ip);
            else if(nDimensions == 3)
                theSection = std::make_unique<Lin3DTee>(h, b, tw, tf, theMesh->GetMaterial(matTag), theta, ip);
        }
        else if(strcasecmp(secShape.c_str(),"WideFlange") == 0){
            unsigned int ip=10;
            double h, b, tw, tf, theta=0.0; 
            InputFile >> h >> b >> tw >> tf >> theta >> ip;

            //Instantiate the section object.
            if(nDimensions == 2)
                theSection = std::make_unique<Lin2DWideFlange>(h, b, tw, tf, theMesh->GetMaterial(matTag), theta, ip);
            else if(nDimensions == 3)
                theSection = std::make_unique<Lin3DWideFlange>(h, b, tw, tf, theMesh->GetMaterial(matTag), theta, ip);
        }
        else if(strcasecmp(secShape.c_str(),"ThinArea") == 0){
            double t;
            InputFile >> t;

            //Instantiate the section object.
            theSection = std::make_unique<Lin3DThinArea>(t, theMesh->GetMaterial(matTag));
        }
    }
    else if(strcasecmp(secType.c_str(),"Fiber") == 0){
        //TODO: Parse fiber sections components.
    }
    else if(strcasecmp(secType.c_str(),"General") == 0){
        //TODO: Parse general sections components.
    }

    //TODO: Add more section models here.

    return Tag;
}

//Obtains a Constraint from Input File.
int 
Parser::CreateConstraint(std::ifstream& InputFile, std::unique_ptr<Constraint> &theConstraint){
    //Auxiliar variables.
    int Tag;
    double value;
    unsigned int nCombs, Slave, dof;

    //Parser the constraint model.
    InputFile >> Tag >> Slave >> nCombs;

    //Add the loads that shares this combination.
    std::vector<unsigned int> Master(nCombs);
    std::vector<double> Coefficients(nCombs);

    for(unsigned int k = 0; k < nCombs; k++){
        InputFile >> dof >> value;
        
        Master[k] = dof;
        Coefficients[k] = value;    
    }

    //Creates a material object.
    theConstraint = std::make_unique<Constraint>(Slave, Master, Coefficients);

    return Tag;
}

//Obtains an Element from Input File.
unsigned int 
Parser::CreateElement(std::ifstream& InputFile, std::shared_ptr<Mesh> &theMesh, std::shared_ptr<Element> &theElement){
    //Auxiliar variables.
    std::string elemName;
    unsigned int dir, dim;
    std::vector<double> parameters;
    std::vector<unsigned int> nodes;

    //Element parse variables.
    std::string Quadrature;
    unsigned int Tag, matID, secID, nGauss;

    //Parser the element component.
    InputFile >> Tag >> elemName;

    if(strcasecmp(elemName.c_str(),"lin2DTruss2") == 0){        
        nodes.resize(2);
        parameters.resize(1);
        InputFile >> nodes[0] >> nodes[1] >> matID >> parameters[0];

        //Instantiate the lin2DTruss2 element.
        theElement = std::make_shared<lin2DTruss2>(nodes, theMesh->GetMaterial(matID), parameters[0], MassForm);
    }
    else if(strcasecmp(elemName.c_str(),"kin2DTruss2") == 0){        
        nodes.resize(2);
        parameters.resize(1);
        InputFile >> nodes[0] >> nodes[1] >> matID >> parameters[0];

        //Instantiate the kin2DTruss2 element.
        theElement = std::make_shared<kin2DTruss2>(nodes, theMesh->GetMaterial(matID), parameters[0], MassForm);
    }
    else if(strcasecmp(elemName.c_str(),"lin3DTruss2") == 0){        
        nodes.resize(2);
        parameters.resize(1);
        InputFile >> nodes[0] >> nodes[1] >> matID >> parameters[0];

        //Instantiate the lin3DTruss2 element.
        theElement = std::make_shared<lin3DTruss2>(nodes, theMesh->GetMaterial(matID), parameters[0], MassForm);
    }
    else if(strcasecmp(elemName.c_str(),"kin3DTruss2") == 0){
        nodes.resize(2);
        parameters.resize(1);
        InputFile >> nodes[0] >> nodes[1] >> matID >> parameters[0];

        //Instantiate the kin3DTruss2 element.
        theElement = std::make_shared<kin3DTruss2>(nodes, theMesh->GetMaterial(matID), parameters[0], MassForm);
    }
    else if(strcasecmp(elemName.c_str(),"lin2DTruss3") == 0){        
        nodes.resize(3);
        parameters.resize(1);
        InputFile >> nodes[0] >> nodes[1] >> nodes[2] >> matID >> parameters[0] >> Quadrature >> nGauss;

        //Instantiate the lin2DTruss3 element.
        theElement = std::make_shared<lin2DTruss3>(nodes, theMesh->GetMaterial(matID), parameters[0], Quadrature, nGauss, MassForm);
    }
    else if(strcasecmp(elemName.c_str(),"lin3DTruss3") == 0){        
        nodes.resize(3);
        parameters.resize(1);
        InputFile >> nodes[0] >> nodes[1] >> nodes[2] >> matID >> parameters[0] >> Quadrature >> nGauss;

        //Instantiate the lin3DTruss3 element.
        theElement = std::make_shared<lin3DTruss3>(nodes, theMesh->GetMaterial(matID), parameters[0], Quadrature, nGauss, MassForm);
    }
    else if(strcasecmp(elemName.c_str(),"ZeroLength1D") == 0){
        nodes.resize(2);
        InputFile >> nodes[0] >> nodes[1] >> matID >> dim >> dir;

        //Instantiate the zerolength1D element.
        theElement = std::make_shared<ZeroLength1D>(nodes, theMesh->GetMaterial(matID), dim, dir);
    }
    else if(strcasecmp(elemName.c_str(),"UnxBoucWen2DLink") == 0){
        nodes.resize(2);
        std::vector<double> vars(4);
        std::vector<double> params(6);

        InputFile >> nodes[0] >> nodes[1] >> params[0] >> params[1] >> params[2] >> params[3] >> params[4] >> params[5] >> vars[0] >> vars[1] >> vars[2] >> vars[3] >> dim >> dir;

        //Instantiate the zerolength1D element.
        theElement = std::make_shared<UnxBoucWen2DLink>(nodes, params, vars, dim, dir);
    }
    else if(strcasecmp(elemName.c_str(),"UnxBoucWen3DLink") == 0){
        nodes.resize(2);
        std::vector<double> vars(4);
        std::vector<double> params(6);

        InputFile >> nodes[0] >> nodes[1] >> params[0] >> params[1] >> params[2] >> params[3] >> params[4] >> params[5] >> vars[0] >> vars[1] >> vars[2] >> vars[3] >> dim >> dir;

        //Instantiate the zerolength1D element.
        theElement = std::make_shared<UnxBoucWen3DLink>(nodes, params, vars, dim, dir);
    }
    else if(strcasecmp(elemName.c_str(),"HDRBYamamoto2DLink") == 0){
        double De, Di, Hr;

        nodes.resize(2);
        InputFile >> nodes[0] >> nodes[1] >> De >> Di >> Hr >> dim;

        //Instantiate the zerolength1D element.
        theElement = std::make_shared<HDRBYamamoto2DLink>(nodes, De, Di, Hr, dim);
    }
    else if(strcasecmp(elemName.c_str(),"HDRBYamamoto3DLink") == 0){
        double De, Di, Hr;

        nodes.resize(2);
        InputFile >> nodes[0] >> nodes[1] >> De >> Di >> Hr >> dim;

        //Instantiate the zerolength1D element.
        theElement = std::make_shared<HDRBYamamoto3DLink>(nodes, De, Di, Hr, dim);
    }
    else if(strcasecmp(elemName.c_str(),"lin2DFrame2") == 0){
        nodes.resize(2);
        std::string Formulation;
        InputFile >> nodes[0] >> nodes[1] >> secID >> Formulation >> Quadrature >> nGauss;

        //Frame element formulation.
        bool Condition = false;
        if(strcasecmp(Formulation.c_str(),"Timoshenko") == 0)
            Condition = true;
                
        //Instantiate the lin2DFrame2 element.
        theElement = std::make_shared<lin2DFrame2>(nodes, theMesh->GetSection(secID), MassForm, Condition, Quadrature, nGauss);
    }
    else if(strcasecmp(elemName.c_str(),"kin2DFrame2") == 0){
        nodes.resize(2);
        std::string Formulation;
        InputFile >> nodes[0] >> nodes[1] >> secID;
                
        //Instantiate the kin2DFrame2 element.
        theElement = std::make_shared<kin2DFrame2>(nodes, theMesh->GetSection(secID), MassForm);
    }
    else if(strcasecmp(elemName.c_str(),"lin3DFrame2") == 0){
        nodes.resize(2);
        std::string Formulation;
        InputFile >> nodes[0] >> nodes[1] >> secID >> Formulation >> Quadrature >> nGauss;

        //Frame element formulation.
        bool Condition = false;
        if(strcasecmp(Formulation.c_str(),"Timoshenko") == 0)
            Condition = true;

        //Instantiate the lin3DFrame2 element.
        theElement = std::make_shared<lin3DFrame2>(nodes, theMesh->GetSection(secID), MassForm, Condition, Quadrature, nGauss);
    }
    else if(strcasecmp(elemName.c_str(),"lin2DQuad4") == 0){
        nodes.resize(4);
        parameters.resize(1);
        InputFile >> nodes[0] >> nodes[1] >> nodes[2] >> nodes[3] >> matID >> parameters[0] >> Quadrature  >> nGauss;

        //Instantiate the lin2DQuad4 element.
        theElement = std::make_shared<lin2DQuad4>(nodes, theMesh->GetMaterial(matID), parameters[0], Quadrature, nGauss, MassForm);
    }
    else if(strcasecmp(elemName.c_str(),"lin2DQuad8") == 0){
        nodes.resize(8);
        parameters.resize(1);
        InputFile >> nodes[0] >> nodes[1] >> nodes[2] >> nodes[3] >> nodes[4] >> nodes[5] >> nodes[6] >> nodes[7] >> matID >> parameters[0] >> Quadrature  >> nGauss;

        //Instantiate the lin2DQuad8 element.
        theElement = std::make_shared<lin2DQuad8>(nodes, theMesh->GetMaterial(matID), parameters[0], Quadrature, nGauss, MassForm);
    }
    else if(strcasecmp(elemName.c_str(),"PML2DQuad4") == 0){
        nodes.resize(4);
        parameters.resize(8);
        InputFile >> nodes[0] >> nodes[1] >> nodes[2] >> nodes[3] >> matID >> parameters[0] >> parameters[1] >> parameters[2] >> parameters[3] >> parameters[4] >> parameters[5] >> parameters[6] >> parameters[7] >> Quadrature >> nGauss;

        //Instantiate the PML2DQuad4 element.
        theElement = std::make_shared<PML2DQuad4>(nodes, theMesh->GetMaterial(matID), parameters, Quadrature, nGauss, MassForm);
    }
    else if(strcasecmp(elemName.c_str(),"PML2DQuad8") == 0){
        nodes.resize(8);
        parameters.resize(8);
        InputFile >> nodes[0] >> nodes[1] >> nodes[2] >> nodes[3] >> nodes[4] >> nodes[5] >> nodes[6] >> nodes[7] >> matID >> parameters[0] >> parameters[1] >> parameters[2] >> parameters[3] >> parameters[4] >> parameters[5] >> parameters[6] >> parameters[7] >> Quadrature >> nGauss;

        //Instantiate the PML2DQuad4 element.
        theElement = std::make_shared<PML2DQuad8>(nodes, theMesh->GetMaterial(matID), parameters, Quadrature, nGauss, MassForm);
    }
    else if(strcasecmp(elemName.c_str(),"kin2DQuad4") == 0){
        nodes.resize(4);
        parameters.resize(1);
        InputFile >> nodes[0] >> nodes[1] >> nodes[2] >> nodes[3] >> matID >> parameters[0] >> Quadrature  >> nGauss;

        //Instantiate the kin2DQuad4 element.
        theElement = std::make_shared<kin2DQuad4>(nodes, theMesh->GetMaterial(matID), parameters[0], Quadrature, nGauss, MassForm);
    }
    else if(strcasecmp(elemName.c_str(),"lin3DShell4") == 0){
        nodes.resize(4);
        InputFile >> nodes[0] >> nodes[1] >> nodes[2] >> nodes[3] >> secID >> Quadrature >> nGauss;
            
        //Instantiate the lin3DShell4 element.
        theElement = std::make_shared<lin3DShell4>(nodes, theMesh->GetSection(secID), Quadrature, nGauss, MassForm);
    }
    else if(strcasecmp(elemName.c_str(),"lin3DHexa8") == 0){
        nodes.resize(8);
        InputFile >> nodes[0] >> nodes[1] >> nodes[2] >> nodes[3] >> nodes[4] >> nodes[5] >> nodes[6] >> nodes[7] >> matID >> Quadrature >> nGauss;

        //Instantiate the lin3DHexa8 element.
        theElement = std::make_shared<lin3DHexa8>(nodes, theMesh->GetMaterial(matID), Quadrature, nGauss, MassForm);
    }
    else if(strcasecmp(elemName.c_str(),"kin3DHexa8") == 0){
        nodes.resize(8);
        InputFile >> nodes[0] >> nodes[1] >> nodes[2] >> nodes[3] >> nodes[4] >> nodes[5] >> nodes[6] >> nodes[7] >> matID >> Quadrature >> nGauss;

        //Instantiate the kin3DHexa8 element.
        theElement = std::make_shared<kin3DHexa8>(nodes, theMesh->GetMaterial(matID), Quadrature, nGauss, MassForm);
    }
    else if(strcasecmp(elemName.c_str(),"PML3DHexa8") == 0){
        nodes.resize(8);
        parameters.resize(9);
        InputFile >> nodes[0] >> nodes[1] >> nodes[2] >> nodes[3] >> nodes[4] >> nodes[5] >> nodes[6] >> nodes[7] >> matID >> parameters[0] >> parameters[1] >> parameters[2] >> parameters[3] >> parameters[4] >> parameters[5] >> parameters[6] >> parameters[7] >> parameters[8] >> Quadrature >> nGauss;

        //Instantiate the PML2DQuad4 element.
        theElement = std::make_shared<PML3DHexa8>(nodes, theMesh->GetMaterial(matID), parameters, Quadrature, nGauss, MassForm);
    }
    else if(strcasecmp(elemName.c_str(),"lin3DHexa20") == 0){
        nodes.resize(20);
        InputFile >> nodes[0] >> nodes[1] >> nodes[2] >> nodes[3] >> nodes[4] >> nodes[5] >> nodes[6] >> nodes[7] >> nodes[8] >> nodes[9] >> nodes[10] >> nodes[11] >> nodes[12] >> nodes[13] >> nodes[14] >> nodes[15] >> nodes[16] >> nodes[17] >> nodes[18] >> nodes[19] >> matID >> Quadrature >> nGauss;

        //Instantiate the lin3DHexa20 element.
        theElement = std::make_shared<lin3DHexa20>(nodes, theMesh->GetMaterial(matID), Quadrature, nGauss, MassForm);
    }
    else if(strcasecmp(elemName.c_str(),"TIEQlin2DQuad4") == 0){
        std::string eType;
        double zref, cf1, cf2, eref;

        nodes.resize(4);
        parameters.resize(1);
        InputFile >> nodes[0] >> nodes[1] >> nodes[2] >> nodes[3] >> matID >> parameters[0] >> Quadrature >> nGauss >> eType >> zref >> cf1 >> cf2 >> eref;

        //Instantiate the TIEQlin2DQuad4 element.
        theElement = std::make_shared<TIEQlin2DQuad4>(nodes, theMesh->GetMaterial(matID), parameters[0], Quadrature, nGauss, MassForm, eType, zref, cf1, cf2, eref);
    }
    else if(strcasecmp(elemName.c_str(),"EQlin2DQuad4") == 0){
        std::string eType;
        double zref, cf1, cf2;

        nodes.resize(4);
        parameters.resize(1);
        InputFile >> nodes[0] >> nodes[1] >> nodes[2] >> nodes[3] >> matID >> parameters[0] >> Quadrature  >> nGauss >> eType >> zref >> cf1 >> cf2;

        //Instantiate the EQlin2DQuad4 element.
        theElement = std::make_shared<EQlin2DQuad4>(nodes, theMesh->GetMaterial(matID), parameters[0], Quadrature, nGauss, MassForm, eType, zref, cf1, cf2);
    }
    else if(strcasecmp(elemName.c_str(),"null2DFrame2") == 0){
        nodes.resize(2);
        InputFile >> nodes[0] >> nodes[1];

        //Instantiate the lin2DTruss2 element.
        theElement = std::make_shared<null2DFrame2>(nodes);
    }
    else if(strcasecmp(elemName.c_str(),"null3DFrame2") == 0){
        nodes.resize(2);
        InputFile >> nodes[0] >> nodes[1];

        //Instantiate the lin2DTruss2 element.
        theElement = std::make_shared<null3DFrame2>(nodes);
    }

    //TODO: Add more elements models here.

    return Tag;
}

//Obtains the Support Motion from Input File.
unsigned int
Parser::CreateSupportMotion(std::ifstream& InputFile, std::vector<double>& Xo, unsigned int &dof){
    //Auxiliar variables.
    unsigned int Tag;
    std::string loadType;

    InputFile >> Tag >> loadType;

    if(strcasecmp(loadType.c_str(),"STATIC") == 0){
        Xo.resize(1);
        InputFile >> Xo[0] >> dof;
    }
    else if(strcasecmp(loadType.c_str(),"DYNAMIC") == 0){
        std::string LoadFile;

        InputFile >> LoadFile >> dof;
        LoadFile = GetSpacedName(LoadFile, " ");

        //Loads the time history values into memory.
        std::ifstream motion(LoadFile.c_str());

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

    return Tag;
}

//Obtains a Point Load from Input File.
unsigned int 
Parser::CreatePointLoad(std::ifstream& InputFile, std::shared_ptr<Load> &theLoad){
    //Auxiliar variables.
    unsigned int Tag, node, ndof, nnodes;
    std::vector<double> value;
    Eigen::VectorXd direction;
    std::vector<unsigned int> nodeTags;

    std::string loadType, LoadFile;

    //Parser the point load component.
    InputFile >> Tag >> loadType;

    if (strcasecmp(loadType.c_str(),"STATIC") == 0){
        value.resize(1);

        InputFile >> value[0] >> ndof;        
        direction.resize(ndof);

        for(unsigned int k = 0; k < ndof ; k++)
            InputFile >> direction(k);

        //Creates the point load.
        theLoad = std::make_shared<Load>(direction, value, 1);

        //Add the nodes that shares this load.
        InputFile >> nnodes;
        nodeTags.resize(nnodes);
        for(unsigned int j = 0; j < nnodes; j++){
            InputFile >> node;
            nodeTags[j] = node;
        }
        theLoad->AddNodes(nodeTags);
    }
    else if(strcasecmp(loadType.c_str(),"DYNAMIC") == 0){
        InputFile >> LoadFile >> ndof;
        LoadFile = GetSpacedName(LoadFile, " ");

        //Loads the time history values into memory.
        std::ifstream load(LoadFile.c_str());

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

        //Reads the load direction.
        direction.resize(ndof);
        for(unsigned int k = 0; k < ndof ; k++)
            InputFile >> direction(k);

        //Creates the point load.
        theLoad = std::make_shared<Load>(direction, value, 2);

        //Add the nodes that shares this load.
        InputFile >> nnodes;
        nodeTags.resize(nnodes);
        for(unsigned int j = 0; j < nnodes; j++){
            InputFile >> node;
            nodeTags[j] = node;
        }
        theLoad->AddNodes(nodeTags);
    }
    else if (strcasecmp(loadType.c_str(),"SUPPORTMOTION") == 0){
        //Creates the support motion as point load.
        theLoad = std::make_shared<Load>(8);

        //Add the nodes that shares this load.
        InputFile >> nnodes;
        nodeTags.resize(nnodes);
        for(unsigned int j = 0; j < nnodes; j++){
            InputFile >> node;
            nodeTags[j] = node;
        }

        theLoad->AddNodes(nodeTags);
    }

    return Tag;
}

//Obtains an Element Load from Input File.
unsigned int 
Parser::CreateElementLoad(std::ifstream& InputFile, std::shared_ptr<Mesh> &theMesh, std::shared_ptr<Load> &theLoad){
    //Auxiliar variables.
    unsigned int Tag, ndof, elem, nelems;
    std::string loadType, loadModel, LoadFile;
    std::vector<unsigned int> elemTags;
    std::vector<unsigned int> faceTags;
    
    //Direction/Values of load.
    Eigen::VectorXd direction;
    std::vector<double> value;

    //Parser the element load component.
    InputFile >> Tag >> loadModel >> loadType;

    if(strcasecmp(loadType.c_str(),"SURFACE") == 0){
        //Implement Element Surface Load.
        if (strcasecmp(loadModel.c_str(),"STATIC") == 0){
            value.resize(1);

            InputFile >> value[0] >> ndof;        
            direction.resize(ndof);

            for(unsigned int k = 0; k < ndof ; k++)
                InputFile >> direction(k);

            //Creates the point load.
            theLoad = std::make_shared<Load>(direction, value, 3);

            //Add the nodes that shares this load.
            InputFile >> nelems;
            elemTags.resize(nelems);
            faceTags.resize(nelems);
            for(unsigned int j = 0; j < nelems; j++){
                InputFile >> elem;
                elemTags[j] = elem;
            }
            for(unsigned int j = 0; j < nelems; j++){
                InputFile >> elem;
                faceTags[j] = elem;
            }
            theLoad->AddFaces(faceTags);
            theLoad->AddElements(elemTags);
        }
        else if(strcasecmp(loadModel.c_str(),"DYNAMIC") == 0){
            //TODO: Implement Dynamic Element Surface Load.
        }
    }
    else if(strcasecmp(loadType.c_str(),"BODY") == 0){
        if (strcasecmp(loadModel.c_str(),"STATIC") == 0){
            value.resize(1);

            InputFile >> value[0] >> ndof;        
            direction.resize(ndof);

            for(unsigned int k = 0; k < ndof ; k++)
                InputFile >> direction(k);

            //Creates the point load.
            theLoad = std::make_shared<Load>(direction, value, 4);

            //Add the nodes that shares this load.
            InputFile >> nelems;
            elemTags.resize(nelems);
            for(unsigned int j = 0; j < nelems; j++){
                InputFile >> elem;
                elemTags[j] = elem;
            }
            theLoad->AddElements(elemTags);
        }
        else if(strcasecmp(loadModel.c_str(),"DYNAMIC") == 0){
            InputFile >> LoadFile >> ndof;
            LoadFile = GetSpacedName(LoadFile, " ");

            //Loads the time history values into memory.
            std::ifstream load(LoadFile.c_str());

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

            //Reads the load direction.
            direction.resize(ndof);
            for(unsigned int k = 0; k < ndof ; k++)
                InputFile >> direction(k);

            //Creates the point load.
            theLoad = std::make_shared<Load>(direction, value, 5);

            //Add the nodes that shares this load.
            InputFile >> nelems;
            elemTags.resize(nelems);
            for(unsigned int j = 0; j < nelems; j++){
                InputFile >> elem;
                elemTags[j] = elem;
            }
            theLoad->AddElements(elemTags);
        }
    }
    else if(strcasecmp(loadType.c_str(),"GENERALWAVE") == 0){
        std::string theFile;
        InputFile >> theFile >> nelems;
        LoadFile = GetSpacedName(LoadFile, " ");

        //Creates the domain reduction load.
        theLoad = std::make_shared<Load>(7);

        //Add the elements that shares this load.
        elemTags.resize(nelems);
        for(unsigned int j = 0; j < nelems; j++){
            InputFile >> elem;
            elemTags[j] = elem;
        }
        theLoad->AddElements(elemTags);

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
        for(auto it : nodes){
            auto &ind = it.first;
            LoadFile  = GetPartitionName(theFile, ind, false);
                                
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
        //Clears node vector of pointers.
        nodes.clear();
    }

    return Tag;
}

//Obtains a Damping model from Input File.
unsigned int 
Parser::CreateDamping(std::ifstream& InputFile, std::vector<unsigned int>& eList, std::shared_ptr<Damping> &theDamping){
    //Auxiliar variables.
    std::string dampName;
    unsigned int Tag, nElems, elemID;
    std::vector<double> parameters;
        
    //Parser the damping model.
    InputFile >> Tag >> dampName;
        
    if(strcasecmp(dampName.c_str(),"Free") == 0){
    }
    else if(strcasecmp(dampName.c_str(),"Rayleigh") == 0){
        parameters.resize(2);
        InputFile >> parameters[0] >> parameters[1];
    }
    else if(strcasecmp(dampName.c_str(),"Caughey") == 0){
        //TODO: Caughey-viscous damping is not implemented yet.
        //parameters.resize(4);
        //InputFile >> parameters[0] >> parameters[1] >> parameters[2] >> parameters[3];
    }
    else if(strcasecmp(dampName.c_str(),"Capped") == 0){
        //TODO: Capped-viscous damping is not implemented yet.
    }

    //List of elements to apply damping.
    InputFile >> nElems;
    eList.resize(nElems);

    for(unsigned int k = 0; k < nElems; k++){
        InputFile >> elemID;
        eList[k] = elemID;    
    }                

    //Stores the damping in the mesh.
    theDamping = std::make_shared<Damping>(dampName, parameters);

    return Tag;
}

//Obtains a Combination object from Input File.
unsigned int 
Parser::CreateCombination(std::ifstream& InputFile, std::shared_ptr<LoadCombo> &theCombo){
    //Auxiliar variables.
    int Trial;
    double factor;
    std::string name;
    unsigned int Tag, IDs;
    std::vector<double> factors;
    std::vector<unsigned int> loads;

    //Parser the load combination component.
    InputFile >> Tag >> name >> Trial;

    //Check if the load combination has associated loads.
    if(Trial != -1){
        loads.resize(Trial);
        factors.resize(Trial);

        for(unsigned int k = 0; k < loads.size(); k++){
            InputFile >> factor >> IDs;
            loads[k]   = IDs;
            factors[k] = factor;
        }
    }

    //Creates the load combination.
    theCombo = std::make_unique<LoadCombo>(name, loads, factors); 

    return Tag;
}

//Obtains a Recorder object from Input File.
void 
Parser::CreateRecorder(std::ifstream& InputFile, std::shared_ptr<Recorder> &theRecorder){
    //Auxiliar variables.
    std::vector<unsigned int> IDs;
    std::string Path, type, file, response;
    unsigned int nlist, nsample, presicion, Index;

    InputFile >> type >> file >> nsample >> presicion;

    if(strcasecmp(type.c_str(),"PARAVIEW") == 0){
        InputFile >> nlist;
        
        //Creates the Paraview Recorder object.
        theRecorder = std::make_shared<Recorder>(Folder, file, type, nlist, nDimensions, nsample, presicion);
    }
    else if(strcasecmp(type.c_str(),"SECTION") == 0){
        unsigned int ndims;

        InputFile >> response >> ndims;
        std::vector<double> xcoord(ndims);

        for(unsigned int k = 0; k < ndims; k++)
            InputFile >> xcoord[k];

        //Gets the list of Sections for Elements for this recorder.
        InputFile >> nlist;
        IDs.resize(nlist);
        for(unsigned int k = 0; k < IDs.size(); k++){
            InputFile >> Index;
            IDs[k] = Index;
        }
        
        //Creates the Paraview Recorder object.
        theRecorder = std::make_shared<Recorder>(Folder, file, type, response, xcoord, IDs, nsample, presicion);
    }
    else{
        InputFile >> response >> nlist;

        //Gets the list of Points/Elements for this recorder.
        IDs.resize(nlist);
        for(unsigned int k = 0; k < IDs.size(); k++){
            InputFile >> Index;
            IDs[k] = Index;
        }

        //Creates the Point/Element Recorder object.
        theRecorder = std::make_shared<Recorder>(Folder, file, type, response, IDs, nsample, presicion);
    }
}

//Obtains an Analysis object from Input File.
void 
Parser::CreateAnalysis(std::ifstream& InputFile, std::shared_ptr<Mesh> &theMesh, std::unique_ptr<Analysis> &theAnalysis, std::map<unsigned int, std::shared_ptr<LoadCombo> > &LoadCombos){
    //Auxiliar variables.
    bool condition;
    std::string aux;
    std::string analysis;
    std::string algorithm;
    std::string integrator;
    std::string solver;
    unsigned int Tag, nSteps, option;    
    double dt, algtol, mtol, ktol, ftol, stol;
    unsigned int nMaxIteration, ConvergenceTest;

    //The pointes to define analysis.
    std::unique_ptr<LinearSystem> theSolver;
    std::shared_ptr<Integrator> theIntegrator;
    std::shared_ptr<Algorithm> theAlgorithm;

    //Parser the analysis component.
    InputFile >> Tag >> analysis >> nSteps >> aux >> algorithm >> algtol >> nMaxIteration >> ConvergenceTest >> aux >> integrator >> mtol >> ktol >> ftol >> dt >> aux >> solver;

    double Factor = (double)nSteps;
    Factor = 1.00/Factor; 

    //Creates the solver.
    if(strcasecmp(solver.c_str(),"EIGEN") == 0){
        InputFile >> condition;
        theSolver = std::make_unique<EigenSolver>(condition);
    }
    else if(strcasecmp(solver.c_str(),"MUMPS") == 0){
        InputFile >> option >> condition;
        theSolver = std::make_unique<MumpsSolver>(option, condition);
    }
    else if(strcasecmp(solver.c_str(),"PETSC") == 0){
        unsigned int kspnum, dnz, onz;
        InputFile >> kspnum >> stol >> dnz >> onz;
        
        theSolver = std::make_unique<PetscSolver>(dnz, onz, stol, kspnum);
    }

    //Creates the algorithm.
    if(strcasecmp(algorithm.c_str(),"LINEAR") == 0){
        theAlgorithm = std::make_shared<Linear>(theSolver, theMesh);
    }
    else if(strcasecmp(algorithm.c_str(),"NEWTON") == 0){
        theAlgorithm = std::make_shared<NewtonRaphson>(theSolver, theMesh, algtol, nMaxIteration, ConvergenceTest);
    }
    //else if(strcasecmp(algorithm.c_str(),"CONTROLPATH") == 0){
    //    theAlgorithm = std::make_shared<GeneralDisplacementPath>(theSolver, theIntegrator, theMesh); 
    //}

    //Creates the integrator.
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

    //Creates Circular dependency between algotithm and integrator.
    theIntegrator->SetAlgorithm(theAlgorithm);
    theAlgorithm->SetIntegrator(theIntegrator);

    //Creates the analysis case.
    if(strcasecmp(analysis.c_str(),"STATIC") == 0){
        theAlgorithm->SetLoadFactor(Factor);
        //Creates the static analysis.            
        theAnalysis = std::make_unique<StaticAnalysis>(theMesh, theAlgorithm, theIntegrator, LoadCombos[Tag], nSteps);
    }
    else if(strcasecmp(analysis.c_str(),"DYNAMIC") == 0){
        theAlgorithm->SetLoadFactor(1.00);
        //Creates the dynamic analysis.
        theAnalysis = std::make_unique<DynamicAnalysis>(theMesh, theAlgorithm, theIntegrator, LoadCombos[Tag], nSteps);
    }
}

//Parse the Input Mesh File.
bool 
Parser::ParseSVLFile(std::shared_ptr<Mesh> &theMesh, std::unique_ptr<Analysis> &theAnalysis, std::string FileName){
    //Parse the input file.
    std::ifstream InputFile(FileName.c_str());

    //Keyword from Input File
    std::string keyword;

    if(InputFile.is_open()){
        //Initializa the mesh subdomain.
        if(Flag)
            theMesh = std::make_shared<Mesh>();
        
        //Reads the first line keyword.
        InputFile >> keyword;
        while(!InputFile.eof()){
            
            if(strcasecmp(keyword.c_str(),"MODEL") == 0){
                //Parse the MODEL Input file line.
                std::string MassFormulation;
                InputFile >> nDimensions >> numberOfTotalDofs >> numberOfFreeDofs >> MassFormulation;
                
                //Mass formulation for elements
                if(strcasecmp(MassFormulation.c_str(),"LUMPED") == 0){
                    MassForm = true;
                }
                if(strcasecmp(MassFormulation.c_str(),"CONSISTENT") == 0){
                    MassForm = false;
                }
            }
            else if(strcasecmp(keyword.c_str(),"MASS") == 0){
                //Assings Nodal Mass.
                unsigned int Tag, nDofs;
                InputFile >> Tag >> nDofs;

                Eigen::VectorXd Mass(nDofs);
                for(unsigned int k = 0; k < nDofs; k++)
                    InputFile >> Mass(k);

                //Adds the nodal mass to the node.
                theMesh->AddMass(Tag, Mass);
            }
            else if(strcasecmp(keyword.c_str(),"NODE") == 0){
                //Creates a Point Object.
                std::shared_ptr<Node> theNode;
                unsigned int Tag = CreatePoint(InputFile, theNode);

                //Stores the information in node container.
                theMesh->AddNode(Tag, theNode);
            }
            else if(strcasecmp(keyword.c_str(),"CONSTRAINT") == 0){
                //Creates a constaint object.
                std::unique_ptr<Constraint> theConstraint;
                int Tag = CreateConstraint(InputFile, theConstraint);

                //Stores the information in material container.
                theMesh->AddConstraint(Tag, theConstraint);
            }
            else if(strcasecmp(keyword.c_str(),"ELEMENT") == 0){
                //Creates an element object.
                std::shared_ptr<Element> theElement;
                unsigned int Tag = CreateElement(InputFile, theMesh, theElement);

                //Stores the information in element container.
                theMesh->AddElement(Tag, theElement);
            }
            else if(strcasecmp(keyword.c_str(),"MATERIAL") == 0){
                //Creates a material object.
                std::unique_ptr<Material> theMaterial;
                unsigned int Tag = CreateMaterial(InputFile, theMaterial);

                //Stores the information in material container.
                theMesh->AddMaterial(Tag, theMaterial);
            }
            else if(strcasecmp(keyword.c_str(),"SECTION") == 0){
                //Creates a section object.
                std::unique_ptr<Section> theSection;
                unsigned int Tag = CreateSection(InputFile, theMesh, theSection);

                //Stores the section in the mesh.
                theMesh->AddSection(Tag, theSection);
            }
            else if(strcasecmp(keyword.c_str(),"DAMPING") == 0){
                //Creates a damping object.
                std::shared_ptr<Damping> theDamping;
                std::vector<unsigned int> eList;
                unsigned int Tag = CreateDamping(InputFile, eList, theDamping);

                //Stores the damping in the mesh.
                theMesh->AddDamping(Tag, theDamping);

                //Assign damping to group of elements.
                theMesh->SetDamping(Tag, eList);
            }
            else if(strcasecmp(keyword.c_str(),"INITIALSTATE") == 0){
                //Parse the INITIALSTATE Input file line.
                int Condition;
                unsigned int Tag, nDofs;
                InputFile >> Tag >> Condition >> nDofs;

                //Inittial State vector values.
                Eigen::VectorXd Xo(nDofs);
                for(unsigned int j = 0; j < nDofs; j++)
                    InputFile >> Xo(j);

                //Assign damping to group of elements.
                theMesh->SetInitialCondition(Tag, Condition, Xo);
            }
            else if(strcasecmp(keyword.c_str(),"SUPPORTMOTION") == 0){
                //Parse the INITIALSTATE Input file line.
                unsigned int Tag, dof;
                std::vector<double> Xo;

                Tag = CreateSupportMotion(InputFile, Xo, dof);

                //Assign damping to group of elements.
                theMesh->SetSupportMotion(Tag, dof, Xo);
            }
            else if(strcasecmp(keyword.c_str(),"POINTLOAD") == 0){
                //Creates a load object.
                std::shared_ptr<Load> theLoad;
                unsigned int Tag = CreatePointLoad(InputFile, theLoad);

                //Stores the load in the mesh.
                theMesh->AddLoad(Tag, theLoad); 
            }
            else if(strcasecmp(keyword.c_str(),"ELEMENTLOAD") == 0){
                //Creates a load object.
                std::shared_ptr<Load> theLoad;
                unsigned int Tag = CreateElementLoad(InputFile, theMesh, theLoad);

                //Stores the load in the mesh.
                theMesh->AddLoad(Tag, theLoad); 
            }
            else if(strcasecmp(keyword.c_str(),"COMBINATION") == 0){
                //Creates a Load Combination Object.
                std::shared_ptr<LoadCombo> theCombo;
                unsigned int Tag = CreateCombination(InputFile, theCombo);

                //Adds the combination to the analysis.
                LoadCombos[Tag] = theCombo;
            }
            else if(strcasecmp(keyword.c_str(),"RECORDER") == 0){
                //Creates a Recorder Object.
                std::shared_ptr<Recorder> theRecorder;
                CreateRecorder(InputFile, theRecorder);

                //Adds the recorder to the analysis.
                Recorders.push_back(theRecorder);
            }
            else if(strcasecmp(keyword.c_str(),"ANALYSIS") == 0){
                //Matrices storage is computed.
                theMesh->Initialize();

                //Creates an Analysis Object.
                CreateAnalysis(InputFile, theMesh, theAnalysis, LoadCombos);

                //Sets the Recorders for the analysis.
                for(unsigned int k = 0; k < Recorders.size(); k++)
                    theAnalysis->SetRecorder(Recorders[k]);
            }
            else if(strcasecmp(keyword.c_str(),"INCLUDE") == 0){
                Flag = false;

                //Sets the include file path.
                std::string IncludeFile;
                InputFile >> IncludeFile;
                IncludeFile = Folder + "/" + IncludeFile;

                bool condition = ParseSVLFile(theMesh, theAnalysis, IncludeFile);

                //Stops the analysis.
                if(condition){
                    return true;
                }
            }
            else{
                //Prints a warning on command line. This should not happens while running. 
                std::cout << "\x1B[33m ALERT: \x1B[0mThe KEYWORD='" << keyword << "' in Processor [" << rank <<  "] is not recognized. \n";
            }

            //Reads the next line keyword.
            InputFile >> keyword;
        }

        InputFile.close();
    }
    else{
        return true;    
    }

    return false;
}

//Parse the User's Input Model/Simulation.
bool 
Parser::GetFromFile(std::shared_ptr<Mesh> &theMesh, std::unique_ptr<Analysis> &theAnalysis){
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

    //Flag to check if something went wrong.
    bool stop;

    //Gets the partitioned mesh file name.
    std::string MeshFile = GetPartitionName(File, rank, true);
    MeshFile = GetSpacedName(MeshFile, " ");

    //Parse the corresponding mesh partition.
    stop = ParseSVLFile(theMesh, theAnalysis, MeshFile);
    if(stop){ 
        std::cout << "\x1B[31mERROR: \x1B[0mThe Mesh file in Processor [" << rank << "] couldn't be opened. \n";
        return true;
    }

    return false;
}
