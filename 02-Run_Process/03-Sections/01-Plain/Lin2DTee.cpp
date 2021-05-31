#include <cmath>
#include "Lin2DTee.hpp"
#include "Definitions.hpp"

//Define constant value PI:
const double PI = 3.1415926535897932;

//Overload Constructor.
Lin2DTee::Lin2DTee(double h, double b, double tw, double tf, std::unique_ptr<Material> &material, double theta, unsigned int ip) : 
Section("Lin2DTee"), h(h), b(b), tw(tw), tf(tf), Theta(theta), InsertPoint(ip){
    //Initialize material strain.
    Strain.resize(3);
    Strain << 0.0, 0.0, 0.0;

    //The section material:
    theMaterial = material->CopyMaterial();

    //Transform the rotation Angle into radians.
    Theta = PI*theta/180.0;
}

//Destructor.
Lin2DTee::~Lin2DTee(){
    //Does nothing.
}

//Clone the 'Lin2DTee' section.
std::unique_ptr<Section>
Lin2DTee::CopySection(){
    return std::make_unique<Lin2DTee>(h, b, tw, tf, theMaterial, 180.0*Theta/PI, InsertPoint);
}

//Returns the section generalized strain in element coordinates.
Eigen::VectorXd
Lin2DTee::GetStrain(){
    return Strain;
}

//Returns the section generalized stress in element coordinates.
Eigen::VectorXd
Lin2DTee::GetStress(){
    Eigen::VectorXd Stress = GetTangentStiffness()*Strain;
    return Stress;
}

//Returns the section axial stiffness.
Eigen::MatrixXd
Lin2DTee::GetDensity(){
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
Lin2DTee::GetTangentStiffness(){
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
Lin2DTee::GetInitialTangentStiffness(){
    return GetTangentStiffness();
}

//Returns the section strain at given position.
Eigen::VectorXd 
Lin2DTee::GetStrainAt(double x3, double x2){
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
    if (((x2 >= h - ycm - tf) & (x2 <= h - ycm) & (x3 >= -b/2.0) & (x3 <= b/2.0)) | ((x2 >= -ycm) & (x2 <= h - ycm - tf) & (x3 >= -tw/2.0) & (x3 <= tw/2.0))) {
        // Epsilon = [exx, 0.0, exy]
        theStrain << Strain(0) - Strain(1)*cos(Theta)*x2 + Strain(1)*sin(Theta)*x3, 
                     0.0, 
                     Strain(2);
    }

    return theStrain;
}

//Returns the section stress at given position.
Eigen::VectorXd 
Lin2DTee::GetStressAt(double x3, double x2){
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
    bool isInside = false;
    double Qs2 , Qs3, th, tb;
    if((x2 >= h - ycm - tf) & (x2 <= h - ycm) & (x3 >= -b/2.0) & (x3 <= b/2.0)){
        isInside = true;
        Qs2 =  b/2.0*(h*h/4.0 - x2*x2); tb = b;
        Qs3 = tf/2.0*(b*b/4.0 - x3*x3); th = tf; 
    }
    else if ((x2 >= -ycm) & (x2 <= h - ycm - tf) & (x3 >= -tw/2.0) & (x3 <= tw/2.0)){
        isInside = true;
        Qs2 = tw/2.0*(ycm*ycm - x2*x2); tb = tw;
        Qs3 = tf/2.0*(b*b/4.0 - x3*x3) + (h - tf)/2.0*(tw*tw/4.0 - x3*x3); th = h; 
    }

    //Computes the generalized section force
    if(isInside){
        //Section geometry properties.
        double A   = GetArea();
        double I22 = GetInertiaAxis2();
        double I33 = GetInertiaAxis3();

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
Lin2DTee::CommitState(){
    theMaterial->CommitState();
}

//Reverse the section states to previous converged state.
void 
Lin2DTee::ReverseState(){
    theMaterial->ReverseState();
}

//Brings the section states to its initial state.
void 
Lin2DTee::InitialState(){
    Strain.fill(0.0);
    theMaterial->InitialState();
}

//Update the section state for this iteration.
void
Lin2DTee::UpdateState(Eigen::VectorXd strain, unsigned int cond){
    //Update the matrial behavior.
    Eigen::VectorXd mStrain(1);
    mStrain << strain(0);
    theMaterial->UpdateState(mStrain, cond);

    //Update the section strain.
    Strain = strain;
}

//Computes the Lin2DTee area.
double 
Lin2DTee::GetArea(){
    return tw*h + (b - tw)*tf;
}

//Computes the Lin2DTee shear area along axis 2.
double 
Lin2DTee::GetShearArea2(){
    return h*tw;
}

//Computes the Lin2DTee shear area along axis 3.
double 
Lin2DTee::GetShearArea3(){
    return 5.0/6.0*b*tf;
}

//Computes the Lin2DTee flexural inertia.
double 
Lin2DTee::GetInertiaAxis2(){
    return 1.0/12.0*(tf*b*b*b + (h - tf)*tw*tw*tw);
}

//Computes the Lin2DTee flexural inertia.
double 
Lin2DTee::GetInertiaAxis3(){
    double Area = GetArea();
    double ycm  = 0.5*(tw*h*h + (b - tw)*tf*tf)/Area;
    double hc   = h - ycm;
    return 1.0/3.0*(b*ycm*ycm*ycm + tw*hc*hc*hc - (b - tw)*(ycm - tf)*(ycm - tf)*(ycm - tf));
}

//Gets the section centroid.
void
Lin2DTee::ComputeSectionCenter(double &zcm, double &ycm){
    //Cross Section-Centroid.
    double Area = GetArea();
    ycm = (tf*b*(h - tf/2.0) + tw*(h - tf)*(h - tf)/2.0)/Area;
    zcm = b/2.0;
}