#include <cmath>
#include "Lin2DAngle.hpp"
#include "Definitions.hpp"

//Define constant value PI:
const double PI = 3.1415926535897932;

//Overload Constructor.
Lin2DAngle::Lin2DAngle(double h, double b, double tw, double tf, std::unique_ptr<Material> &material, double theta, unsigned int ip) : 
Section("Lin2DAngle"), h(h), b(b), tw(tw), tf(tf), Theta(theta), InsertPoint(ip){
    //Initialize material strain.
    Strain.resize(3);
    Strain << 0.0, 0.0, 0.0;

    //The section material:
    theMaterial = material->CopyMaterial();

    //Transform the rotation Angle into radians.
    Theta = PI*theta/180.0;
}

//Destructor.
Lin2DAngle::~Lin2DAngle(){
    //Does nothing.
}

//Clone the 'Lin2DAngle' section.
std::unique_ptr<Section>
Lin2DAngle::CopySection(){
    return std::make_unique<Lin2DAngle>(h, b, tw, tf, theMaterial, 180.0*Theta/PI, InsertPoint);
}

//Returns the section generalized strain in element coordinates.
Eigen::VectorXd
Lin2DAngle::GetStrain(){
    return Strain;
}

//Returns the section generalized stress in element coordinates.
Eigen::VectorXd
Lin2DAngle::GetStress(){
    Eigen::VectorXd Stress = GetTangentStiffness()*Strain;
    return Stress;
}

//Returns the section axial stiffness.
Eigen::MatrixXd
Lin2DAngle::GetDensity(){
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
Lin2DAngle::GetTangentStiffness(){
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
    double I23 = GetInertiaAxis23();

    //Section shear Stiffness at center of mass.
    double G = theMaterial->GetShearModulus();

    //Section flexural Stiffness at center of mass.
    Eigen::MatrixXd E = theMaterial->GetInitialTangentStiffness(); 

    //Computes the global section stiffness matrix.
    double EA = E(0,0)*A;
    double GAs2 = G*As2*cos(Theta)*cos(Theta) + G*As3*sin(Theta)*sin(Theta);
    double EI33 = E(0,0)*I33*cos(Theta)*cos(Theta) + E(0,0)*I22*sin(Theta)*sin(Theta) + E(0,0)*I23*sin(2.0*Theta) + E(0,0)*A*(x2*x2*cos(Theta)*cos(Theta) + x2*x3*sin(2.0*Theta) + x3*x3*sin(Theta)*sin(Theta));

    //Returns the section stiffness for 2-dimensions.
    Eigen::MatrixXd SectionStiffness(3,3);
    SectionStiffness <<  EA,  0.0,  0.0,
                        0.0, EI33,  0.0,
                        0.0,  0.0, GAs2;

    return SectionStiffness;
}

//Returns the section initial tangent stiffness matrix.
Eigen::MatrixXd 
Lin2DAngle::GetInitialTangentStiffness(){
    return GetTangentStiffness();
}

//Returns the section strain at given position.
Eigen::VectorXd 
Lin2DAngle::GetStrainAt(double x3, double x2){
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
    if (((x2 >= h - ycm - tf) & (x2 <= h - ycm) & (x3 >= zcm - b + tw) & (x3 <= zcm)) | ((x2 >= -ycm) & (x2 <= h - ycm) & (x3 >= zcm - b) & (x3 <= zcm - b + tw))) {
        // Epsilon = [exx, 0.0, exy]
        theStrain << Strain(0) - Strain(1)*cos(Theta)*x2 + Strain(1)*sin(Theta)*x3, 
                     0.0, 
                     Strain(2);
    }

    return theStrain;
}

//Returns the section stress at given position.
Eigen::VectorXd 
Lin2DAngle::GetStressAt(double x3, double x2){
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
    double Qs2, Qs3, th, tb;
    if((x2 >= h - ycm - tf) & (x2 <= h - ycm) & (x3 >= zcm + tw - b) & (x3 <= zcm)){
        isInside = true;
        Qs2 =  b/2.0*((h - ycm)*(h - ycm) - x2*x2); tb = b;
        Qs3 = tf/2.0*(zcm*zcm - x3*x3); th = tf; 
    }
    else if((x2 >= -ycm) & (x2 <= h - ycm) & (x3 >= zcm - b) & (x3 <= zcm - b + tw)){
        isInside = true;
        Qs2 = tf*(b - tw)*(h - ycm - tf/2.0) + tw/2.0*((h - ycm)*(h - ycm) - x2*x2); tb = tw;
        Qs3 = h/2.0*((b - zcm)*(b - zcm) - x3*x3); th = h; 
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
Lin2DAngle::CommitState(){
    theMaterial->CommitState();
}

//Reverse the section states to previous converged state.
void 
Lin2DAngle::ReverseState(){
    theMaterial->ReverseState();
}

//Brings the section states to its initial state.
void 
Lin2DAngle::InitialState(){
    Strain.fill(0.0);
    theMaterial->InitialState();
}

//Update the section state for this iteration.
void
Lin2DAngle::UpdateState(Eigen::VectorXd strain, unsigned int cond){
    //Update the matrial behavior.
    Eigen::VectorXd mStrain(1);
    mStrain << strain(0);
    theMaterial->UpdateState(mStrain, cond);

    //Update the section strain.
    Strain = strain;
}

//Computes the Lin2DAngle area.
double 
Lin2DAngle::GetArea(){
    return tw*h + tf*(b - tw);
}

//Computes the Lin2DAngle shear area along axis 2.
double 
Lin2DAngle::GetShearArea2(){
    return 5.0/6.0*h*tw;
}

//Computes the Lin2DAngle shear area along axis 3.
double 
Lin2DAngle::GetShearArea3(){
    return 5.0/6.0*b*tf;
}

//Computes the Lin2DAngle flexural inertia.
double 
Lin2DAngle::GetInertiaAxis2(){
    double Area = GetArea();
    double xcm  = 0.5*(tf*b*b + tw*tw*(h - tf))/Area;
    double dc   = xcm - tw;
    return 1.0/3.0*(h*xcm*xcm*xcm + tf*(b - xcm)*(b - xcm)*(b - xcm) - (h - tf)*dc*dc*dc);
}

//Computes the Lin2DAngle flexural inertia.
double 
Lin2DAngle::GetInertiaAxis3(){
    double Area = GetArea();
    double ycm  = 0.5*(tw*h*h + tf*tf*(b - tw))/Area;
    double hc   = ycm - tf;
    return 1.0/3.0*(b*ycm*ycm*ycm + tw*(h - ycm)*(h - ycm)*(h - ycm) - (b - tw)*hc*hc*hc);
}

//Computes the Lin2DAngle product of inertia.
double 
Lin2DAngle::GetInertiaAxis23(){
    double Area = GetArea();
    double xcm  = 0.5*(tf*b*b + tw*tw*(h - tf))/Area;
    double ycm  = 0.5*(tw*h*h + tf*tf*(b - tw))/Area;
    return h*tw*(xcm - tw/2.0)*(h/2.0 - ycm) + tf*(b - tw)*((b + tw)/2.0 - xcm)*(ycm - tf/2.0);
}

//Gets the section centroid.
void
Lin2DAngle::ComputeSectionCenter(double &zcm, double &ycm){
    //Cross Section-Centroid.
    double Area = GetArea();
    zcm = (tf*b*b/2.0 + tw*(h - tf)*(b - tw/2.0))/Area;
    ycm = (tw*h*h/2.0 + tf*(b - tw)*(h - tf/2.0))/Area;
}