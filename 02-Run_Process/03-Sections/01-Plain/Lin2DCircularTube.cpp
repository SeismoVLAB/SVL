#include <cmath>
#include "Lin2DCircularTube.hpp"
#include "Definitions.hpp"

//Define constant value PI:
const double PI = 3.1415926535897932;

//Overload Constructor.
Lin2DCircularTube::Lin2DCircularTube(double re, double ri, std::unique_ptr<Material> &material, double theta, unsigned int ip) : 
Section("Lin2DCircularTube"), re(re), ri(ri), Theta(theta), InsertPoint(ip){
    //Initialize material strain.
    Strain.resize(3);
    Strain << 0.0, 0.0, 0.0;

    //The section material:
    theMaterial = material->CopyMaterial();

    //Transform the rotation Anlgle into radians.
    Theta = PI*theta/180.0;
}

//Destructor.
Lin2DCircularTube::~Lin2DCircularTube(){
    //Does nothing.
}

//Clone the 'Lin2DCircularTube' section.
std::unique_ptr<Section>
Lin2DCircularTube::CopySection(){
    return std::make_unique<Lin2DCircularTube>(re, ri, theMaterial, 180.0*Theta/PI, InsertPoint);
}

//Returns the section generalized strain.
Eigen::VectorXd
Lin2DCircularTube::GetStrain(){
    return Strain;
}

//Returns the section generalized stress.
Eigen::VectorXd
Lin2DCircularTube::GetStress(){
    Eigen::VectorXd Stress = GetTangentStiffness()*Strain;
    return Stress;
}

//Returns the section axial stiffness.
Eigen::MatrixXd
Lin2DCircularTube::GetDensity(){
    //Area Properties.
    double A   = GetArea();

    //Material density.
    double Rho = theMaterial->GetDensity(); 

    //Returns the section stiffness for 2-dimensions.
    Eigen::MatrixXd SectionDensity(2,2);
    SectionDensity << Rho*A,  0.0 ,
                       0.0 , Rho*A;

    return SectionDensity;
}

//Returns the section axial stiffness.
Eigen::MatrixXd
Lin2DCircularTube::GetTangentStiffness(){
    //Area Properties.
    double A   = GetArea();
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
    double As = As2 + 2.0/PI*(As3 - As2)*Theta; 

    //Section Stiffness at centroid.
    Eigen::MatrixXd Cs(3,3);
    Cs << E(0,0)*A, 0.0       , 0.0       ,
          0.0     , E(0,0)*I22, 0.0       ,
          0.0     , 0.0       , E(0,0)*I33;

    //Compute the rotation matrix.
    Eigen::MatrixXd T = GetLineRotationMatrix(Theta);

    //Compute the translation matrix.
    Eigen::MatrixXd L = GetLineTranslationMatrix(re, ri, zc, yc, InsertPoint);

    //Computes the global section stiffness matrix.
    Cs = (T.transpose()*L.transpose())*Cs*(L*T);

    //Returns the section stiffness for 2-dimensions.
    Eigen::MatrixXd SectionStiffness(3,3);
    SectionStiffness << Cs(0,0),   0.0  , 0.0 ,
                          0.0  , Cs(2,2), 0.0 ,
                          0.0  ,   0.0  , G*As;

    return SectionStiffness;
}

//Returns the section initial tangent stiffness matrix.
Eigen::MatrixXd 
Lin2DCircularTube::GetInitialTangentStiffness(){
    return GetTangentStiffness();
}

//Returns the section strain at given position.
Eigen::VectorXd 
Lin2DCircularTube::GetStrainAt(double UNUSED(x3), double UNUSED(x2)){
    //TODO: Compute the strain at point, needs shear value. 
    Eigen::VectorXd theStrain(3);
    theStrain.fill(0.0);

    return theStrain;
}

//Returns the section stress at given position.
Eigen::VectorXd 
Lin2DCircularTube::GetStressAt(double UNUSED(x3), double UNUSED(x2)){
    //TODO: Compute the strain at point, needs shear value. 
    // Sigma = [Sxx, 0.0, txy]
    Eigen::VectorXd theStress(3);
    theStress.fill(0.0);

    return theStress;
}

//Perform converged section state update.
void 
Lin2DCircularTube::CommitState(){
    theMaterial->CommitState();
}

//Update the section state for this iteration.
void
Lin2DCircularTube::UpdateState(Eigen::VectorXd strain, unsigned int cond){
    //Update the matrial behavior.
    Eigen::VectorXd mStrain(1);
    mStrain << strain(0);
    theMaterial->UpdateState(mStrain, cond);

    //Update the section strain.
    Strain = strain;
}

//Computes the Lin2DCircularTube area.
double 
Lin2DCircularTube::GetArea(){
    return PI*(re*re - ri*ri);
}

//Computes the Lin2DCircularTube shear area along axis 2.
double 
Lin2DCircularTube::GetShearArea2(){
    return PI*(re + ri)*(re - ri);
}

//Computes the Lin2DCircularTube shear area along axis 3.
double 
Lin2DCircularTube::GetShearArea3(){
    return PI*(re + ri)*(re - ri);
}

//Computes the Lin2DCircularTube flexural inertia.
double 
Lin2DCircularTube::GetInertiaAxis2(){
    return PI*(re*re*re*re - ri*ri*ri*ri)/4.0;
}

//Computes the Lin2DCircularTube flexural inertia.
double 
Lin2DCircularTube::GetInertiaAxis3(){
    return PI*(re*re*re*re - ri*ri*ri*ri)/4.0;
}

//Gets the section centroid.
void
Lin2DCircularTube::ComputeSectionCenter(double &zcm, double &ycm){
    //Cross Section-Centroid.
    ycm = re;
    zcm = re;
}
