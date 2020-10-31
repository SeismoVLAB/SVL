#include <cmath>
#include "Lin2DRectangularTube.hpp"
#include "Definitions.hpp"

//Define constant value PI:
const double PI = 3.1415926535897932;

//Overload Constructor.
Lin2DRectangularTube::Lin2DRectangularTube(double h, double b, double tw, double tf, std::unique_ptr<Material> &material, double theta, unsigned int ip) : 
Section("Lin2DRectangularTube"), h(h), b(b), tw(tw), tf(tf), Theta(theta), InsertPoint(ip){
    //Initialize material strain.
    Strain.resize(3);
    Strain << 0.0, 0.0, 0.0;

    //The section material:
    theMaterial = material->CopyMaterial();

    //Transform the rotation Anlgle into radians.
    Theta = PI*theta/180.0;
}

//Destructor.
Lin2DRectangularTube::~Lin2DRectangularTube(){
    //Does nothing.
}

//Clone the 'Lin2DRectangularTube' section.
std::unique_ptr<Section>
Lin2DRectangularTube::CopySection(){
    return std::make_unique<Lin2DRectangularTube>(h, b, tw, tf, theMaterial, 180.0*Theta/PI, InsertPoint);
}

//Returns the section generalized strain.
Eigen::VectorXd
Lin2DRectangularTube::GetStrain(){
    return Strain;
}

//Returns the section generalized stress.
Eigen::VectorXd
Lin2DRectangularTube::GetStress(){
    Eigen::VectorXd Stress = GetTangentStiffness()*Strain;
    return Stress;
}

//Returns the section axial stiffness.
Eigen::MatrixXd
Lin2DRectangularTube::GetDensity(){
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
Lin2DRectangularTube::GetTangentStiffness(){
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
Lin2DRectangularTube::GetInitialTangentStiffness(){
    return GetTangentStiffness();
}

//Returns the section strain at given position.
Eigen::VectorXd 
Lin2DRectangularTube::GetStrainAt(double x3, double x2){
    //Checks the coordinate is inside the section
    if (!((fabs(x2) >= h/2.0 - tf) & (fabs(x2) <= h/2.0) & (x3 >= -b/2.0) & (x3 <= b/2.0)) | !((x2 >= -h/2.0 + tf) & (x2 <= h/2.0 - tf) & (fabs(x3) >= b/2.0 - tw) & (fabs(x3) <= b/2.0))) {
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
Lin2DRectangularTube::GetStressAt(double x3, double x2){
    double A   = GetArea();
    double I33 = GetInertiaAxis3();

    //Computes statical moment of area for Zhuravskii shear stress formula.
    double Qs2, tb;
    if((fabs(x2) >= h/2.0 - tf) & (fabs(x2) <= h/2.0) & (x3 >= -b/2.0) & (x3 <= b/2.0)){
        Qs2 =  b/2.0*(h*h/4.0 - x2*x2); tb = b;
    }
    else if((x2 >= -h/2.0 + tf) & (x2 <= h/2.0 - tf) & (fabs(x3) >= b/2.0 - tw) & (fabs(x3) <= b/2.0)){
        Qs2 = tf/2.0*(b - 2.0*tw)*(h - tf) + tw*(h*h/4.0 - x2*x2); tb = 2.0*tw;
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
Lin2DRectangularTube::CommitState(){
    theMaterial->CommitState();
}

//Update the section state for this iteration.
void
Lin2DRectangularTube::UpdateState(Eigen::VectorXd strain, unsigned int cond){
    //Update the matrial behavior.
    Eigen::VectorXd mStrain(1);
    mStrain << strain(0);
    theMaterial->UpdateState(mStrain, cond);

    //Update the section strain.
    Strain = strain;
}

//Computes the Lin2DRectangularTube area.
double 
Lin2DRectangularTube::GetArea(){
    return h*b - (h - 2.0*tf)*(b - 2.0*tw);
}

//Computes the Lin2DRectangularTube shear area along axis 2.
double 
Lin2DRectangularTube::GetShearArea2(){
    return 2*h*tw;
}

//Computes the Lin2DRectangularTube shear area along axis 3.
double 
Lin2DRectangularTube::GetShearArea3(){
    return 2*b*tf;
}

//Computes the Lin2DRectangularTube flexural inertia.
double 
Lin2DRectangularTube::GetInertiaAxis2(){
    return 1.0/12.0*(h*b*b*b - (h - 2.0*tf)*(b - 2.0*tw)*(b - 2.0*tw)*(b - 2.0*tw));
}

//Computes the Lin2DRectangularTube flexural inertia.
double 
Lin2DRectangularTube::GetInertiaAxis3(){
    return 1.0/12.0*(b*h*h*h - (b - 2.0*tw)*(h - 2.0*tf)*(h - 2.0*tf)*(h - 2.0*tf));
}

//Gets the section centroid.
void
Lin2DRectangularTube::ComputeSectionCenter(double &zcm, double &ycm){
    //Cross Section-Centroid.
    ycm = h/2.0;
    zcm = b/2.0;
}
