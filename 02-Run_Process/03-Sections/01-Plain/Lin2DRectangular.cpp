#include <cmath>
#include "Lin2DRectangular.hpp"
#include "Definitions.hpp"

//Define constant value PI:
const double PI = 3.1415926535897932;

//Overload Constructor.
Lin2DRectangular::Lin2DRectangular(double h, double b, std::unique_ptr<Material> &material, double theta, unsigned int ip) : 
Section("Lin2DRectangular"), h(h), b(b), Theta(theta), InsertPoint(ip){
    //Initialize material strain.
    Strain.resize(3);
    Strain << 0.0, 0.0, 0.0;

    //The section material:
    theMaterial = material->CopyMaterial();

    //Transform the rotation Anlgle into radians.
    Theta = PI*theta/180.0;
}

//Destructor.
Lin2DRectangular::~Lin2DRectangular(){
    //Does nothing.
}

//Clone the 'Lin2DRectangular' section.
std::unique_ptr<Section>
Lin2DRectangular::CopySection(){
    return std::make_unique<Lin2DRectangular>(h, b, theMaterial, 180.0*Theta/PI, InsertPoint);
}

//Returns the section generalized strain.
Eigen::VectorXd
Lin2DRectangular::GetStrain(){
    return Strain;
}

//Returns the section generalized stress.
Eigen::VectorXd
Lin2DRectangular::GetStress(){
    Eigen::VectorXd Stress = GetTangentStiffness()*Strain;
    return Stress;
}

//Returns the section axial stiffness.
Eigen::MatrixXd
Lin2DRectangular::GetDensity(){
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
Lin2DRectangular::GetTangentStiffness(){
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
    Eigen::MatrixXd L = GetLineTranslationMatrix(h, b, zc, yc, InsertPoint);

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
Lin2DRectangular::GetInitialTangentStiffness(){
    return GetTangentStiffness();
}

//Returns the section strain at given position.
Eigen::VectorXd 
Lin2DRectangular::GetStrainAt(double x3, double x2){
    //Checks the coordinate is inside the section
    if (!((x2 >= -h/2.0) & (x2 <= h/2.0) & (x3 >= -b/2.0) & (x3 <= b/2.0))) {
        x3 = -b/2.0; x2 = h/2.0;
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
Lin2DRectangular::GetStressAt(double x3, double x2){
    double A   = GetArea();
    double I33 = GetInertiaAxis3();

    //Computes statical moment of area for Zhuravskii shear stress formula.
    double Qs2, tb;
    if((x2 >= -h/2.0) & (x2 <= h/2.0) & (x3 >= -b/2.0) & (x3 <= b/2.0)){
        Qs2 = b/2.0*(h*h/4.0 - x2*x2); tb = b;
    }
    else{
        x3 = -b/2.0; x2 = h/2.0; Qs2 = 0.0; tb = b;
    }

    //Get the generalised stresses or section forces.
    Eigen::VectorXd Forces = GetStress();

    // Sigma = [Sxx, 0.0, txy]
    Eigen::VectorXd theStress(3);
    theStress << Forces(0)/A - Forces(1)*x2/I33, 
                 0.0, 
                 Forces(2)*Qs2/I33/tb;

    return theStress;
}

//Perform converged section state update.
void 
Lin2DRectangular::CommitState(){
    theMaterial->CommitState();
}

//Update the section state for this iteration.
void
Lin2DRectangular::UpdateState(Eigen::VectorXd strain, unsigned int cond){
    //Update the matrial behavior.
    Eigen::VectorXd mStrain(1);
    mStrain << strain(0);
    theMaterial->UpdateState(mStrain, cond);

    //Update the section strain.
    Strain = strain;
}

//Computes the Lin2DRectangular area.
double 
Lin2DRectangular::GetArea(){
    return h*b;
}

//Computes the Lin2DRectangular shear area along axis 2.
double 
Lin2DRectangular::GetShearArea2(){
    double A = GetArea();
    return 5.0/6.0*A;
}

//Computes the Lin2DRectangular shear area along axis 3.
double 
Lin2DRectangular::GetShearArea3(){
    double A = GetArea();
    return 5.0/6.0*A;
}

//Computes the Lin2DRectangular flexural inertia.
double 
Lin2DRectangular::GetInertiaAxis2(){
    return 1.0/12.0*h*b*b*b;
}

//Computes the Lin2DRectangular flexural inertia.
double 
Lin2DRectangular::GetInertiaAxis3(){
    return 1.0/12.0*b*h*h*h;
}

//Gets the section centroid.
void
Lin2DRectangular::ComputeSectionCenter(double &zcm, double &ycm){
    //Cross Section-Centroid.
    ycm = h/2.0;
    zcm = b/2.0;
}
