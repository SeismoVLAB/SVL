#include <cmath>
#include "Lin3DRectangular.hpp"
#include "Definitions.hpp"

//Define constant value PI:
const double PI = 3.1415926535897932;

//Overload Constructor.
Lin3DRectangular::Lin3DRectangular(double h, double b, std::unique_ptr<Material> &material, double theta, unsigned int ip) : 
Section("Lin3DRectangular"), h(h), b(b), Theta(theta), InsertPoint(ip){
    //Initialize material strain.
    Strain.resize(6);
    Strain << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

    //The section material:
    theMaterial = material->CopyMaterial();

    //Transform the rotation Angle into radians.
    Theta = PI*theta/180.0;
}

//Destructor.
Lin3DRectangular::~Lin3DRectangular(){
    //Does nothing.
}

//Clone the 'Lin3DRectangular' section.
std::unique_ptr<Section>
Lin3DRectangular::CopySection(){
    return std::make_unique<Lin3DRectangular>(h, b, theMaterial, 180.0*Theta/PI, InsertPoint);
}

//Returns the section generalized strain.
Eigen::VectorXd
Lin3DRectangular::GetStrain(){
    return Strain;
}

//Returns the section generalized stress.
Eigen::VectorXd
Lin3DRectangular::GetStress(){
    Eigen::VectorXd Stress = GetTangentStiffness()*Strain;
    return Stress;
}

//Returns the section axial stiffness.
Eigen::MatrixXd
Lin3DRectangular::GetDensity(){
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
Lin3DRectangular::GetTangentStiffness(){
    //Gets the section centroid.
    double zcm, ycm;
    ComputeSectionCenter(zcm, ycm);

    //Area Properties.
    double A   = GetArea();
    double J   = GetInertiaAxis1();
    double As2 = GetShearArea2();
    double As3 = GetShearArea3();
    double I22 = GetInertiaAxis2();
    double I33 = GetInertiaAxis3();

    //Section shear Stiffness at center of mass.
    double G = theMaterial->GetShearModulus();

    Eigen::MatrixXd Cs(3,3);
    Cs << G*J,  0.0 ,  0.0 ,
          0.0, G*As2,  0.0 ,
          0.0,  0.0 , G*As3;

    //Section flexural Stiffness at center of mass.
    Eigen::MatrixXd E = theMaterial->GetInitialTangentStiffness(); 

    Eigen::MatrixXd Cm(3,3);
    Cm << E(0,0)*A, 0.0       , 0.0       ,
          0.0     , E(0,0)*I22, 0.0       ,
          0.0     , 0.0       , E(0,0)*I33;

    //Compute the rotation matrix.
    Eigen::MatrixXd T = GetLineRotationMatrix(Theta);

    //Compute the translation matrix.
    Eigen::MatrixXd L = GetLineTranslationMatrix(h, b, zcm, ycm, InsertPoint);

    //Computes the global section stiffness matrix.
    Cs = T.transpose()*Cs*T;
    Cm = T.transpose()*L.transpose()*Cm*L*T;

    //Returns the section stiffness for 3-dimensions.
    Eigen::MatrixXd SectionStiffness(6,6);
    SectionStiffness << Cm(0,0),   0.0  , Cm(0,1), Cm(0,2),   0.0  ,   0.0  ,
                          0.0  , Cs(0,0),   0.0  ,   0.0  , Cs(0,1), Cs(0,2),
                        Cm(0,1),   0.0  , Cm(1,1), Cm(2,1),   0.0  ,   0.0  ,
                        Cm(0,2),   0.0  , Cm(2,1), Cm(2,2),   0.0  ,   0.0  ,
                          0.0  , Cs(0,1),   0.0  ,   0.0  , Cs(1,1), Cs(2,1),
                          0.0  , Cs(0,2),   0.0  ,   0.0  , Cs(2,1), Cs(2,2);

    return SectionStiffness;
}

//Returns the section initial tangent stiffness matrix.
Eigen::MatrixXd 
Lin3DRectangular::GetInitialTangentStiffness(){
    return GetTangentStiffness();
}

//Returns the section strain at given position.
Eigen::VectorXd 
Lin3DRectangular::GetStrainAt(double x3, double x2){
    //The strain vector in local coordinates
    Eigen::VectorXd theStrain(6);
    theStrain.fill(0.0);

    //Gets the section centroid.
    double zcm, ycm;
    ComputeSectionCenter(zcm, ycm);

    //Transform coordinates into section local axis.
    x3 = zcm - x3;
    x2 = x2 - ycm;

    //Checks the coordinate is inside the section
    if ( ((x2 >= -h/2.0) & (x2 <= h/2.0) & (x3 >= -b/2.0) & (x3 <= b/2.0)) ) {
        //Transforms generalised strains from Element to Section local coordinate
        Eigen::VectorXd strain = ComputeLineLocalAxes(h, b, zcm, ycm, Theta, InsertPoint)*Strain;

        // Epsilon = [exx, 0.0, 0.0, exy, 0.0, exz]
        theStrain << strain(0) + x3*strain(2) - x2*strain(3), 
                     0.0, 
                     0.0, 
                     strain(4),
                     0.0, 
                     strain(5);
    }

    return theStrain;
}

//Returns the section stress at given position.
Eigen::VectorXd 
Lin3DRectangular::GetStressAt(double x3, double x2){
    //The stress vector in local coordinates
    Eigen::VectorXd theStress(6);
    theStress.fill(0.0);

    //Gets the section centroid.
    double zcm, ycm;
    ComputeSectionCenter(zcm, ycm);

    //Transform coordinates into section local axis.
    x3 = zcm - x3;
    x2 = x2 - ycm;

    //Checks the coordinate is inside the section
    if((x2 >= -h/2.0) & (x2 <= h/2.0) & (x3 >= -b/2.0) & (x3 <= b/2.0)){
        //Section geometry properties.
        double A   = GetArea();
        double I22 = GetInertiaAxis2();
        double I33 = GetInertiaAxis3();

        //Statical moment of area for Zhuravskii shear stress formula.
        double tb = b;
        double th = h; 
        double Qs2 = b/2.0*(h*h/4.0 - x2*x2);
        double Qs3 = h/2.0*(b*b/4.0 - x3*x3);

        //Transforms generalised stresses from Element to Section local coordinate
        Eigen::VectorXd Forces = ComputeLineLocalAxes(h, b, zcm, ycm, Theta, InsertPoint)*GetStress();

        // Sigma = [Sxx, 0.0, 0.0, txy, 0.0, txz]
        theStress << Forces(0)/A + Forces(2)*x3/I22 - Forces(3)*x2/I33, 
                     0.0, 
                     0.0, 
                     Forces(4)*Qs2/I33/tb, 
                     0.0,
                     Forces(5)*Qs3/I22/th;        
    }

    return theStress;
}

//Perform converged section state update.
void 
Lin3DRectangular::CommitState(){
    theMaterial->CommitState();
}

//Reverse the section states to previous converged state.
void 
Lin3DRectangular::ReverseState(){
    theMaterial->ReverseState();
}

//Brings the section states to its initial state.
void 
Lin3DRectangular::InitialState(){
    Strain.fill(0.0);
    theMaterial->InitialState();
}

//Update the section state for this iteration.
void
Lin3DRectangular::UpdateState(Eigen::VectorXd strain, unsigned int cond){
    //Update the matrial behavior.
    Eigen::VectorXd mStrain(1);
    mStrain << strain(0);
    theMaterial->UpdateState(mStrain, cond);

    //Update the section strain.
    Strain = strain;
}

//Computes the Lin3DRectangular area.
double 
Lin3DRectangular::GetArea(){
    return h*b;
}

//Computes the Lin3DRectangular shear area along axis 2.
double 
Lin3DRectangular::GetShearArea2(){
    double A = GetArea();
    return 5.0/6.0*A;
}

//Computes the Lin3DRectangular shear area along axis 3.
double 
Lin3DRectangular::GetShearArea3(){
    double A = GetArea();
    return 5.0/6.0*A;
}

//Computes the Lin3DRectangular torsional inertia.
double 
Lin3DRectangular::GetInertiaAxis1(){
    return 1.0/12.0*h*b*(b*b + h*h);
}

//Computes the Lin3DRectangular flexural inertia.
double 
Lin3DRectangular::GetInertiaAxis2(){
    return 1.0/12.0*h*b*b*b;
}

//Computes the Lin3DRectangular flexural inertia.
double 
Lin3DRectangular::GetInertiaAxis3(){
    return 1.0/12.0*b*h*h*h;
}

//Gets the section centroid.
void
Lin3DRectangular::ComputeSectionCenter(double &zcm, double &ycm){
    //Cross Section-Centroid.
    ycm = h/2.0;
    zcm = b/2.0;
}
