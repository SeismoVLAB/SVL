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

    //Transform the rotation Angle into radians.
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

//Returns the section generalized strain in element coordinates.
Eigen::VectorXd
Lin2DCircularTube::GetStrain(){
    return Strain;
}

//Returns the section generalized stress in element coordinates.
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
    //Gets the section centroid.
    double x2, x3, zcm, ycm;
    ComputeSectionCenter(zcm, ycm);
    InsertionPointCoordinates(x3, x2, re, re, zcm, ycm, InsertPoint);

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
Lin2DCircularTube::GetInitialTangentStiffness(){
    return GetTangentStiffness();
}

//Returns the section strain at given position.
Eigen::VectorXd 
Lin2DCircularTube::GetStrainAt(double x3, double x2){
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
    if ((ri <= sqrt(x2*x2 + x3*x3)) & (sqrt(x2*x2 + x3*x3) <= re)) {
        // Epsilon = [exx, 0.0, exy]
        theStrain << Strain(0) - Strain(1)*cos(Theta)*x2 + Strain(1)*sin(Theta)*x3, 
                     0.0, 
                     Strain(2);
    }

    return theStrain;
}

//Returns the section stress at given position.
Eigen::VectorXd 
Lin2DCircularTube::GetStressAt(double x3, double x2){
    //The stress vector in local coordinates
    Eigen::VectorXd theStress(3);
    theStress.fill(0.0);

    //Gets the section centroid.
    double zcm, ycm;
    ComputeSectionCenter(zcm, ycm);

    //Transform coordinates into section local axis.
    x3 = zcm - x3;
    x2 = x2 - ycm;

    //Checks the coordinate is inside the section
    if ((ri <= sqrt(x2*x2 + x3*x3)) & (sqrt(x2*x2 + x3*x3) <= re)) {
        //Gets the section centroid.
        double A   = GetArea();
        double I22 = GetInertiaAxis2();
        double I33 = GetInertiaAxis3();

        //First moment area
        double rex2 = fabs(re*re - x2*x2);
        double rex3 = fabs(re*re - x3*x3);
        double rix2 = fabs(ri*ri - x2*x2);
        double rix3 = fabs(ri*ri - x3*x3);
        double Qs2 = ri < fabs(x2) ? 2.0/3.0*pow(rex2, 1.5) : 2.0/3.0*(pow(rex2, 1.5) - pow(rix2, 1.5)); 
        double Qs3 = ri < fabs(x3) ? 2.0/3.0*pow(rex3, 1.5) : 2.0/3.0*(pow(rex3, 1.5) - pow(rix3, 1.5)); 
        double tb = ri <= fabs(x2) ? 2.0*sqrt(rex2) : 2.0*(sqrt(rex2) - sqrt(rix2));
        double th = ri <= fabs(x3) ? 2.0*sqrt(rex3) : 2.0*(sqrt(rex3) - sqrt(rix3));

        //Get the generalised stresses or section forces.
        Eigen::VectorXd Forces = GetStress();

        // Sigma = [Sxx, 0.0, txy]
        theStress << Forces(0)/A + Forces(1)*x3*sin(Theta)/I22 - Forces(1)*x2*cos(Theta)/I33, 
                     0.0, 
                     Forces(2)*cos(Theta)*Qs2/I33/tb + Forces(2)*sin(Theta)*Qs3/I22/th;
    }

    return theStress;
}

//Perform converged section state update.
void 
Lin2DCircularTube::CommitState(){
    theMaterial->CommitState();
}

//Reverse the section states to previous converged state.
void 
Lin2DCircularTube::ReverseState(){
    theMaterial->ReverseState();
}

//Brings the section states to its initial state.
void 
Lin2DCircularTube::InitialState(){
    Strain.fill(0.0);
    theMaterial->InitialState();
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