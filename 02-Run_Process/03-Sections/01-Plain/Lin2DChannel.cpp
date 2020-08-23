#include <cmath>
#include "Lin2DChannel.hpp"
#include "Definitions.hpp"

//Define constant value PI:
const double PI = 3.1415926535897932;

//Overload Constructor.
Lin2DChannel::Lin2DChannel(double h, double b, double tw, double tf, std::unique_ptr<Material> &material, double theta, unsigned int ip) : 
Section("Lin2DChannel"), h(h), b(b), tw(tw), tf(tf), Theta(theta), InsertPoint(ip){
    //Initialize material strain.
    Strain.resize(3);
    Strain << 0.0, 0.0, 0.0;

    //The section material:
    theMaterial = material->CopyMaterial();

    //Transform the rotation Anlgle into radians.
    Theta = PI*theta/180.0;
}

//Destructor.
Lin2DChannel::~Lin2DChannel(){
    //Does nothing.
}

//Clone the 'Lin2DChannel' section.
std::unique_ptr<Section>
Lin2DChannel::CopySection(){
    return std::make_unique<Lin2DChannel>(h, b, tw, tf, theMaterial, 180.0*Theta/PI, InsertPoint);
}

//Returns the section generalized strain.
Eigen::VectorXd
Lin2DChannel::GetStrain(){
    return Strain;
}

//Returns the section generalized stress.
Eigen::VectorXd
Lin2DChannel::GetStress(){
    Eigen::VectorXd Stress = GetTangentStiffness()*Strain;
    return Stress;
}

//Returns the section axial stiffness.
Eigen::MatrixXd
Lin2DChannel::GetDensity(){
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
Lin2DChannel::GetTangentStiffness(){
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
Lin2DChannel::GetInitialTangentStiffness(){
    return GetTangentStiffness();
}

//Returns the section strain at given position.
Eigen::VectorXd 
Lin2DChannel::GetStrainAt(double x3, double x2){
    //Checks the coordinate is inside the section
    double zcm, ycm;
    ComputeSectionCenter(zcm, ycm);
    if(!((fabs(x2) >= h/2.0 - tf) & (fabs(x2) <= h/2.0) & (x3 >= zcm - b) & (x3 <= zcm)) | !((x2 >= -h/2.0 + tf) & (x2 <= h/2.0 - tf) & (x3 >= zcm - b) & (x3 <= zcm - b + tw))){
        x3 = zcm - b; x2 = h/2.0;
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
Lin2DChannel::GetStressAt(double x3, double x2){
    double A   = GetArea();
    double I33 = GetInertiaAxis3();

    //Computes statical moment of area for Zhuravskii shear stress formula.
    double zcm, ycm;
    ComputeSectionCenter(zcm, ycm);

    double Qs2, tb;
    if((fabs(x2) >= h/2.0 - tf) & (fabs(x2) <= h/2.0) & (x3 >= zcm - b) & (x3 <= zcm)){
        Qs2 =  b/2.0*(h*h/4.0 - x2*x2); tb = b; 
    }
    else if((x2 >= -h/2.0 + tf) & (x2 <= h/2.0 - tf) & (x3 >= zcm - b) & (x3 <= zcm - b + tw)){
        Qs2 =  tf/2.0*(h - tf)*(b - tw) + tw/2.0*(h*h/4.0 - x2*x2); tb = tw; 
    }
    else{
        x2 = h/2.0; Qs2 = 0.0; tb = b;
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
Lin2DChannel::CommitState(){
    theMaterial->CommitState();
}

//Update the section state for this iteration.
void
Lin2DChannel::UpdateState(Eigen::VectorXd strain, unsigned int cond){
    //Update the matrial behavior.
    Eigen::VectorXd matStrain(1);
    matStrain << strain(0);
    theMaterial->UpdateState(matStrain, cond);

    //Update the section strain.
    Strain = strain;
}

//Computes the Lin2DChannel area.
double 
Lin2DChannel::GetArea(){
    return 2.0*tf*b + tw*(h - 2.0*tf);
}

//Computes the Lin2DChannel shear area along axis 2.
double 
Lin2DChannel::GetShearArea2(){
    return h*tw;
}

//Computes the Lin2DChannel shear area along axis 3.
double 
Lin2DChannel::GetShearArea3(){
    return 5.0/6.0*(2.0*b*tf);
}

//Computes the Lin2DChannel flexural inertia.
double 
Lin2DChannel::GetInertiaAxis2(){
    double Area = GetArea();
    double zcm  = (tf*b*b + 0.5*tw*tw*(h - 2.0*tf))/Area;
    double hc   = zcm - tw;
    double dc   = b - zcm;
    return 1.0/3.0*(h*zcm*zcm*zcm + 2.0*tf*dc*dc*dc - (h - 2.0*tf)*hc*hc*hc);
}

//Computes the Lin2DChannel flexural inertia.
double 
Lin2DChannel::GetInertiaAxis3(){
    return 1.0/12.0*(b*h*h*h - (b - tw)*(h - 2.0*tf)*(h - 2.0*tf)*(h - 2.0*tf));
}

//Gets the section centroid.
void
Lin2DChannel::ComputeSectionCenter(double &zcm, double &ycm){
    //Cross Section-Centroid.
    double Area = GetArea();
    ycm = h/2.0;
    zcm = (tf*b*b + tw*(h - 2.0*tf)*(b - tw/2.0))/Area;
}
