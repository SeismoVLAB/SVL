#include <cmath>
#include "Lin2DCircular.hpp"
#include "Definitions.hpp"

//Define constant value:
const double PI = 3.1415926535897932;

//Overload Constructor.
Lin2DCircular::Lin2DCircular(double r, std::unique_ptr<Material> &material, double theta, unsigned int ip) : 
Section("Lin2DCircular"), r(r), Theta(theta), InsertPoint(ip){
    //Initialize material strain.
    Strain.resize(3);
    Strain << 0.0, 0.0, 0.0;

    //The section material:
    theMaterial = material->CopyMaterial();

    //Transform the rotation Anlgle into radians.
    Theta = PI*theta/180.0;
}

//Destructor.
Lin2DCircular::~Lin2DCircular(){
    //Does nothing.
}

//Clone the 'Lin2DCircular' section.
std::unique_ptr<Section>
Lin2DCircular::CopySection(){
    return std::make_unique<Lin2DCircular>(r, theMaterial, 180.0*Theta/PI, InsertPoint);
}

//Returns the section generalized strain.
Eigen::VectorXd
Lin2DCircular::GetStrain(){
    return Strain;
}

//Returns the section generalized stress.
Eigen::VectorXd
Lin2DCircular::GetStress(){
    Eigen::VectorXd Stress = GetTangentStiffness()*Strain;
    return Stress;
}

//Returns the section axial stiffness.
Eigen::MatrixXd
Lin2DCircular::GetDensity(){
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
Lin2DCircular::GetTangentStiffness(){
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

    //Section Stiffness at center of centroid.
    Eigen::MatrixXd Cs(3,3);
    Cs << E(0,0)*A, 0.0       , 0.0       ,
          0.0     , E(0,0)*I22, 0.0       ,
          0.0     , 0.0       , E(0,0)*I33;

    //Compute the rotation matrix.
    Eigen::MatrixXd T = GetLineRotationMatrix(Theta);

    //Compute the translation matrix.
    Eigen::MatrixXd L = GetLineTranslationMatrix(r, r, zc, yc, InsertPoint);

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
Lin2DCircular::GetInitialTangentStiffness(){
    return GetTangentStiffness();
}

//Returns the section strain at given position.
Eigen::VectorXd 
Lin2DCircular::GetStrainAt(double x3, double x2){
    //Checks the coordinate is inside the section
    if (sqrt(x2*x2 + x3*x3) > r) {
        x3 = 0.0; x2 = r;
    }

    // Epsilon = [exx, 0.0, exy]
    Eigen::VectorXd theStrain(3);
    theStrain << Strain(0) - x2*Strain(1), 
                 0.0, 
                 Strain(2);

    return theStrain;
}

//Returns the section stress at given position.
Eigen::VectorXd 
Lin2DCircular::GetStressAt(double x3, double x2){
    double A   = GetArea();
    double I33 = GetInertiaAxis3();

    //Checks the coordinate is inside the section
    if (sqrt(x2*x2 + x3*x3) > r) {
        x3 = 0.0; x2 = r;
    }

    //Get the generalised stresses or section forces.
    Eigen::VectorXd Forces = GetStress();

    // Sigma = [Sxx, 0.0, txy]
    Eigen::VectorXd theStress(3);
    theStress << Forces(0)/A - Forces(1)*x2/I33, 
                 0.0, 
                 Forces(2)*(r*r - x2*x2)/I33/3.0;

    return theStress;
}

//Perform converged section state update.
void 
Lin2DCircular::CommitState(){
    theMaterial->CommitState();
}

//Update the section state for this iteration.
void
Lin2DCircular::UpdateState(Eigen::VectorXd strain, unsigned int cond){
    //Update the matrial behavior.
    Eigen::VectorXd mStrain(1);
    mStrain << strain(0);
    theMaterial->UpdateState(mStrain, cond);

    //Update the section strain.
    Strain = strain;
}

//Computes the Lin2DCircular area.
double 
Lin2DCircular::GetArea(){
    return PI*r*r;
}

//Computes the Lin2DCircular shear area along axis 2.
double 
Lin2DCircular::GetShearArea2(){
    double A = GetArea();
    return 0.9*A;
}

//Computes the Lin2DCircular shear area along axis 3.
double 
Lin2DCircular::GetShearArea3(){
    double A = GetArea();
    return 0.9*A;
}

//Computes the Lin2DCircular flexural inertia.
double 
Lin2DCircular::GetInertiaAxis2(){
    return PI*r*r*r*r/4.0;
}

//Computes the Lin2DCircular flexural inertia.
double 
Lin2DCircular::GetInertiaAxis3(){
    return PI*r*r*r*r/4.0;
}

//Gets the section centroid.
void
Lin2DCircular::ComputeSectionCenter(double &zcm, double &ycm){
    //Cross Section-Centroid.
    ycm = r;
    zcm = r;
}
