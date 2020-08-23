#include <cmath>
#include "Lin3DCircularTube.hpp"
#include "Definitions.hpp"

//Define constant value PI:
const double PI = 3.1415926535897932;

//Overload Constructor.
Lin3DCircularTube::Lin3DCircularTube(double re, double ri, std::unique_ptr<Material> &material, double theta, unsigned int ip) : 
Section("Lin3DCircularTube"), re(re), ri(ri), Theta(theta), InsertPoint(ip){
    //Initialize material strain.
    Strain.resize(6);
    Strain << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

    //The section material:
    theMaterial = material->CopyMaterial();

    //Transform the rotation Anlgle into radians.
    Theta = PI*theta/180.0;
}

//Destructor.
Lin3DCircularTube::~Lin3DCircularTube(){
    //Does nothing.
}

//Clone the 'Lin3DCircularTube' section.
std::unique_ptr<Section>
Lin3DCircularTube::CopySection(){
    return std::make_unique<Lin3DCircularTube>(re, ri, theMaterial, 180.0*Theta/PI, InsertPoint);
}

//Returns the section generalized strain.
Eigen::VectorXd
Lin3DCircularTube::GetStrain(){
    return Strain;
}

//Returns the section generalized stress.
Eigen::VectorXd
Lin3DCircularTube::GetStress(){
    Eigen::VectorXd Stress = GetTangentStiffness()*Strain;
    return Stress;
}

//Returns the section axial stiffness.
Eigen::MatrixXd
Lin3DCircularTube::GetDensity(){
    //Area Properties.
    double A   = GetArea();
    double J   = GetInertiaAxis1();

    //Material density.
    double Rho = theMaterial->GetDensity(); 

    //Returns the section stiffness for 3-dimensions.
    Eigen::MatrixXd SectionDensity(4,4);
    SectionDensity << Rho*A, 0.0 ,  0.0 ,  0.0 ,
                       0.0 , J/A ,  0.0 ,  0.0 ,
                       0.0 , 0.0 , Rho*A,  0.0 ,
                       0.0 , 0.0 ,  0.0 , Rho*A;

    return SectionDensity;
}

//Returns the section axial stiffness.
Eigen::MatrixXd
Lin3DCircularTube::GetTangentStiffness(){
    //Area Properties.
    double A   = GetArea();
    double J   = GetInertiaAxis1();
    double As2 = GetShearArea2();
    double As3 = GetShearArea3();
    double I22 = GetInertiaAxis2();
    double I33 = GetInertiaAxis3();

    //Gets the section centroid.
    double zc, yc;
    ComputeSectionCenter(zc, yc);

    //Material properties.
    double G = theMaterial->GetShearModulus();
    Eigen::MatrixXd E = theMaterial->GetInitialTangentStiffness(); 

    //Section Shear areas at centroid.
    //TODO: Check this transformation is not right.
    double Asy = As2 + 2.0/PI*(As3 - As2)*Theta; 
    double Asz = As3 + 2.0/PI*(As2 - As3)*Theta; 

    //Section Stiffness at center of mass.
    Eigen::MatrixXd Cs(3,3);
    Cs << E(0,0)*A, 0.0       , 0.0       ,
          0.0     , E(0,0)*I22, 0.0       ,
          0.0     , 0.0       , E(0,0)*I33;

    //Compute the rotation matrix.
    Eigen::MatrixXd T = GetLineRotationMatrix(Theta);

    //Compute the translation matrix.
    Eigen::MatrixXd L = GetLineTranslationMatrix(re, ri, zc, yc, InsertPoint);

    //Computes the global section stiffness matrix.
    Cs = T.transpose()*L.transpose()*Cs*L*T;

    //Returns the section stiffness for 3-dimensions.
    Eigen::MatrixXd SectionStiffness(6,6);
    SectionStiffness << Cs(0,0), 0.0, Cs(0,1), Cs(0,2),  0.0 ,  0.0 ,
                          0.0  , G*J,   0.0  ,   0.0  ,  0.0 ,  0.0 ,
                        Cs(0,1), 0.0, Cs(1,1), Cs(2,1),  0.0 ,  0.0 ,
                        Cs(0,2), 0.0, Cs(2,1), Cs(2,2),  0.0 ,  0.0 ,
                          0.0  , 0.0,   0.0  ,   0.0  , G*Asy,  0.0 ,
                          0.0  , 0.0,   0.0  ,   0.0  ,  0.0 , G*Asz;

    return SectionStiffness;
}

//Returns the section initial tangent stiffness matrix.
Eigen::MatrixXd 
Lin3DCircularTube::GetInitialTangentStiffness(){
    return GetTangentStiffness();
}

//Returns the section strain at given position.
Eigen::VectorXd 
Lin3DCircularTube::GetStrainAt(double x3, double x2){
    //TODO: Compute the strain at point, needs shear value. 
    Eigen::VectorXd theStrain(6);
    theStrain.fill(0.0);

    return theStrain;
}

//Returns the section stress at given position.
Eigen::VectorXd 
Lin3DCircularTube::GetStressAt(double x3, double x2){
    //TODO: Compute the strain at point, needs shear value. 
    // Sigma = [Sxx, 0.0, 0.0, txy, 0.0, txz]
    Eigen::VectorXd theStress(6);
    theStress.fill(0.0);

    return theStress;
}

//Perform converged section state update.
void 
Lin3DCircularTube::CommitState(){
    theMaterial->CommitState();
}

//Update the section state for this iteration.
void
Lin3DCircularTube::UpdateState(Eigen::VectorXd strain, unsigned int cond){
    //Update the matrial behavior.
    Eigen::VectorXd matStrain(1);
    matStrain << strain(0);
    theMaterial->UpdateState(matStrain, cond);

    //Update the section strain.
    Strain = strain;
}

//Computes the Lin3DCircularTube area.
double 
Lin3DCircularTube::GetArea(){
    return PI*(re*re - ri*ri);
}

//Computes the Lin3DCircularTube shear area along axis 2.
double 
Lin3DCircularTube::GetShearArea2(){
    return PI*(re + ri)*(re - ri)/2.0;
}

//Computes the Lin3DCircularTube shear area along axis 3.
double 
Lin3DCircularTube::GetShearArea3(){
    return PI*(re + ri)*(re - ri)/2.0;
}

//Computes the Lin3DCircularTube torsional inertia.
double 
Lin3DCircularTube::GetInertiaAxis1(){
    return PI*(re*re*re*re - ri*ri*ri*ri)/2.0;
}

//Computes the Lin3DCircularTube flexural inertia.
double 
Lin3DCircularTube::GetInertiaAxis2(){
    return PI*(re*re*re*re - ri*ri*ri*ri)/4.0;
}

//Computes the Lin3DCircularTube flexural inertia.
double 
Lin3DCircularTube::GetInertiaAxis3(){
    return PI*(re*re*re*re - ri*ri*ri*ri)/4.0;
}

//Gets the section centroid.
void
Lin3DCircularTube::ComputeSectionCenter(double &zcm, double &ycm){
    //Cross Section-Centroid.
    ycm = re;
    zcm = re;
}
