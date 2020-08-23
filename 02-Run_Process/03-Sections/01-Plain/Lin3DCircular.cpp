#include <cmath>
#include "Lin3DCircular.hpp"
#include "Definitions.hpp"

//Define constant value:
const double PI = 3.1415926535897932;

//Overload Constructor.
Lin3DCircular::Lin3DCircular(double r, std::unique_ptr<Material> &material, double theta, unsigned int ip) : 
Section("Lin3DCircular"), r(r), Theta(theta), InsertPoint(ip){
    //Initialize material strain.
    Strain.resize(6);
    Strain << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

    //The section material:
    theMaterial = material->CopyMaterial();

    //Transform the rotation Anlgle into radians.
    Theta = PI*theta/180.0;
}

//Destructor.
Lin3DCircular::~Lin3DCircular(){
    //Does nothing.
}

//Clone the 'Lin3DCircular' section.
std::unique_ptr<Section>
Lin3DCircular::CopySection(){
    return std::make_unique<Lin3DCircular>(r, theMaterial, 180.0*Theta/PI, InsertPoint);
}

//Returns the section generalized strain.
Eigen::VectorXd
Lin3DCircular::GetStrain(){
    return Strain;
}

//Returns the section generalized stress.
Eigen::VectorXd
Lin3DCircular::GetStress(){
    Eigen::VectorXd Stress = GetTangentStiffness()*Strain;
    return Stress;
}

//Returns the section axial stiffness.
Eigen::MatrixXd
Lin3DCircular::GetDensity(){
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
Lin3DCircular::GetTangentStiffness(){
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
    Eigen::MatrixXd L = GetLineTranslationMatrix(r, r, zc, yc, InsertPoint);

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
Lin3DCircular::GetInitialTangentStiffness(){
    return GetTangentStiffness();
}

//Returns the section strain at given position.
Eigen::VectorXd 
Lin3DCircular::GetStrainAt(double x3, double x2){
    //Checks the coordinate is inside the section
    if (sqrt(x2*x2 + x3*x3) > r) {
        x3 = 0.0; x2 = r;
    }

    // Epsilon = [exx, 0.0, 0.0, exy, 0.0, exz]
    Eigen::VectorXd theStrain(6);
    theStrain << Strain(0) - x2*Strain(2) + x3*Strain(3), 
                 0.0, 
                 0.0, 
                 Strain(4), 
                 0.0,
                 Strain(5);

    return theStrain;
}

//Returns the section stress at given position.
Eigen::VectorXd 
Lin3DCircular::GetStressAt(double x3, double x2){
    double A   = GetArea();
    double I22 = GetInertiaAxis2();
    double I33 = GetInertiaAxis3();

    //Checks the coordinate is inside the section
    if (sqrt(x2*x2 + x3*x3) > r) {
        x3 = 0.0; x2 = r;
    }

    //Get the generalised stresses or section forces.
    Eigen::VectorXd Forces = GetStress();

    // Sigma = [Sxx, 0.0, 0.0, txy, 0.0, txz]
    Eigen::VectorXd theStress(6);
    theStress << Forces(0)/A - Forces(2)*x2/I33 + Forces(3)*x3/I22, 
                 0.0, 
                 0.0, 
                 Forces(5)*(r*r -x2*x2)/I33/3.0, 
                 0.0,
                 Forces(6)*(r*r - x3*x3)/I22/3.0;

    return theStress;
}

//Perform converged section state update.
void 
Lin3DCircular::CommitState(){
    theMaterial->CommitState();
}

//Update the section state for this iteration.
void
Lin3DCircular::UpdateState(Eigen::VectorXd strain, unsigned int cond){
    theMaterial->UpdateState(strain, cond);
}

//Computes the Lin3DCircular area.
double 
Lin3DCircular::GetArea(){
    return PI*r*r;
}

//Computes the Lin3DCircular shear area along axis 2.
double 
Lin3DCircular::GetShearArea2(){
    double A = GetArea();
    return 0.9*A;
}

//Computes the Lin3DCircular shear area along axis 3.
double 
Lin3DCircular::GetShearArea3(){
    double A = GetArea();
    return 0.9*A;
}

//Computes the Lin3DCircular torsional inertia.
double 
Lin3DCircular::GetInertiaAxis1(){
    return PI*r*r*r*r/2.0;
}

//Computes the Lin3DCircular flexural inertia.
double 
Lin3DCircular::GetInertiaAxis2(){
    return PI*r*r*r*r/4.0;
}

//Computes the Lin3DCircular flexural inertia.
double 
Lin3DCircular::GetInertiaAxis3(){
    return PI*r*r*r*r/4.0;
}

//Gets the section centroid.
void
Lin3DCircular::ComputeSectionCenter(double &zcm, double &ycm){
    //Cross Section-Centroid.
    ycm = r;
    zcm = r;
}
