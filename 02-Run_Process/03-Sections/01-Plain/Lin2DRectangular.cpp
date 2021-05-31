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

    //Transform the rotation Angle into radians.
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

//Returns the section generalized strain in element coordinates.
Eigen::VectorXd
Lin2DRectangular::GetStrain(){
    return Strain;
}

//Returns the section generalized stress in element coordinates.
Eigen::VectorXd
Lin2DRectangular::GetStress(){
    Eigen::VectorXd Stress = GetTangentStiffness()*Strain;
    return Stress;
}

//Returns the section axial stiffness.
Eigen::MatrixXd
Lin2DRectangular::GetDensity(){
    //Area Properties.
    double A = GetArea();

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
    //Gets the section centroid.
    double x2, x3, zcm, ycm;
    ComputeSectionCenter(zcm, ycm);
    InsertionPointCoordinates(x3, x2, h, b, zcm, ycm, InsertPoint);

    //Area Properties.
    double A   = GetArea();
    double As2 = GetShearArea2();
    double As3 = GetShearArea3();
    double I22 = GetInertiaAxis2();
    double I33 = GetInertiaAxis3();

    //Section shear Stiffness at center of mass.
    double G = theMaterial->GetShearModulus();

    //Section flexural Stiffness at center of mass.
    Eigen::MatrixXd E = theMaterial->GetInitialTangentStiffness(); 

    //Computes the global section stiffness matrix.
    double EA = E(0,0)*A;
    double GAs2 = G*As2*cos(Theta)*cos(Theta) + G*As3*sin(Theta)*sin(Theta);
    double EI33 = E(0,0)*I33*cos(Theta)*cos(Theta) + E(0,0)*I22*sin(Theta)*sin(Theta) + E(0,0)*A*(x2*x2*cos(Theta)*cos(Theta) + x2*x3*sin(2.0*Theta) + x3*x3*sin(Theta)*sin(Theta));

    //Returns the section stiffness for 2-dimensions.
    Eigen::MatrixXd SectionStiffness(3,3);
    SectionStiffness <<  EA,  0.0,  0.0,
                        0.0, EI33,  0.0,
                        0.0,  0.0, GAs2;

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
    //The strain vector in local coordinates
    Eigen::VectorXd theStrain(3);
    theStrain.fill(0.0);

    //Gets the section centroid.
    double zcm, ycm;
    ComputeSectionCenter(zcm, ycm);

    //Transform coordinates into section local axis.
    x3 = zcm - x3;
    x2 = x2 - ycm;

    //Checks the coordinate is inside the section
    if ((x2 >= -h/2.0) & (x2 <= h/2.0) & (x3 >= -b/2.0) & (x3 <= b/2.0)) {
        // Epsilon = [exx, 0.0, exy]
        theStrain << Strain(0) - Strain(1)*cos(Theta)*x2 + Strain(1)*sin(Theta)*x3, 
                     0.0, 
                     Strain(2);
    }

    return theStrain;
}

//Returns the section stress at given position.
Eigen::VectorXd 
Lin2DRectangular::GetStressAt(double x3, double x2){
    //The stress vector in local coordinates
    Eigen::VectorXd theStress(3);
    theStress.fill(0.0);

    //Gets the section centroid.
    double zcm, ycm;
    ComputeSectionCenter(zcm, ycm);

    //Transform coordinates into section local axis.
    x3 = zcm - x3;
    x2 = x2 - ycm;

    //Computes statical moment of area for Zhuravskii shear stress formula.
    if((x2 >= -h/2.0) & (x2 <= h/2.0) & (x3 >= -b/2.0) & (x3 <= b/2.0)){
        //Section geometry properties.
        double A   = GetArea();
        double I22 = GetInertiaAxis2();
        double I33 = GetInertiaAxis3();

        //First moment area
        double tb = b;
        double th = h; 
        double Qs2 = b/2.0*(h*h/4.0 - x2*x2);
        double Qs3 = h/2.0*(b*b/4.0 - x3*x3);

        //Get the generalised stresses or section forces.
        Eigen::VectorXd Forces = GetStress();

        // Sigma = [Sxx, 0.0, txy]
        theStress << Forces(0)/A + Forces(1)*x3*sin(Theta)/I22 - Forces(1)*x2*cos(Theta)/I33, 
                     0.0, 
                     Forces(2)*cos(Theta)*Qs2/I33/tb - Forces(2)*sin(Theta)*Qs3/I22/th;
    }

    return theStress;
}

//Perform converged section state update.
void 
Lin2DRectangular::CommitState(){
    theMaterial->CommitState();
}

//Reverse the section states to previous converged state.
void 
Lin2DRectangular::ReverseState(){
    theMaterial->ReverseState();
}

//Brings the section states to its initial state.
void 
Lin2DRectangular::InitialState(){
    Strain.fill(0.0);
    theMaterial->InitialState();
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