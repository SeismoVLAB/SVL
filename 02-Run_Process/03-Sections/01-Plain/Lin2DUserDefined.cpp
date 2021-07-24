#include <cmath>
#include "Lin2DUserDefined.hpp"
#include "Definitions.hpp"

//Define constant value PI:
const double PI = 3.1415926535897932;

//Overload Constructor.
Lin2DUserDefined::Lin2DUserDefined(std::vector<double> properties, std::unique_ptr<Material> &material, double theta) : 
Section("Lin2DUserDefined"), Theta(theta){
    //Assign section properties:
    A   = properties[0];
    As2 = properties[1];
    I33 = properties[2];

    //Initialize material strain.
    Strain.resize(3);
    Strain << 0.0, 0.0, 0.0;

    //The section material:
    theMaterial = material->CopyMaterial();

    //Transform the rotation Angle into radians.
    Theta = PI*theta/180.0;
}

//Destructor.
Lin2DUserDefined::~Lin2DUserDefined(){
    //Does nothing.
}

//Clone the 'Lin2DUserDefined' section.
std::unique_ptr<Section>
Lin2DUserDefined::CopySection(){
    std::vector<double> properties{A, As2, I33};
    return std::make_unique<Lin2DUserDefined>(properties, theMaterial, 180.0*Theta/PI);
}

//Returns the section generalized strain.
Eigen::VectorXd
Lin2DUserDefined::GetStrain(){
    return Strain;
}

//Returns the section generalized stress.
Eigen::VectorXd
Lin2DUserDefined::GetStress(){
    Eigen::VectorXd Stress = GetTangentStiffness()*Strain;
    return Stress;
}

//Returns the section axial stiffness.
Eigen::MatrixXd
Lin2DUserDefined::GetDensity(){
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
Lin2DUserDefined::GetTangentStiffness(){
    //Section shear Stiffness at center of mass.
    double G = theMaterial->GetShearModulus();

    //Section flexural Stiffness at center of mass.
    Eigen::MatrixXd E = theMaterial->GetInitialTangentStiffness(); 

    //Returns the section stiffness for 2-dimensions.
    Eigen::MatrixXd SectionStiffness(3,3);
    SectionStiffness << E(0,0)*A, 0.0       , 0.0,
                        0.0     , E(0,0)*I33, 0.0,
                        0.0     , 0.0       , G*As2;

    return SectionStiffness;
}

//Returns the section initial tangent stiffness matrix.
Eigen::MatrixXd 
Lin2DUserDefined::GetInitialTangentStiffness(){
    return GetTangentStiffness();
}

//Returns the section strain at given position.
Eigen::VectorXd 
Lin2DUserDefined::GetStrainAt(double UNUSED(x3), double x2){
    //The strain vector in local coordinates
    Eigen::VectorXd theStrain(3);
    theStrain.fill(0.0);

    // Epsilon = [exx, 0.0, exy]
    theStrain << Strain(0) - Strain(1)*x2, 
                 0.0, 
                 Strain(2);

    return theStrain;
}

//Returns the section stress at given position.
Eigen::VectorXd 
Lin2DUserDefined::GetStressAt(double UNUSED(x3), double x2){
    //Section shear Stiffness at center of mass.
    double G = theMaterial->GetShearModulus();

    //Get the generalised stresses or section forces.
    Eigen::VectorXd Forces = GetStress();

    //The stress vector in local coordinates
    Eigen::VectorXd theStress(3);
    theStress.fill(0.0);

    // Sigma = [Sxx, 0.0, txy]
    theStress << Forces(0)/A - Forces(1)*x2/I33, 
                 0.0, 
                 G*Strain(2); 

    return theStress;
}

//Perform converged section state update.
void 
Lin2DUserDefined::CommitState(){
    theMaterial->CommitState();
}

//Reverse the section states to previous converged state.
void 
Lin2DUserDefined::ReverseState(){
    theMaterial->ReverseState();
}

//Brings the section states to its initial state.
void 
Lin2DUserDefined::InitialState(){
    Strain.fill(0.0);
    theMaterial->InitialState();
}

//Update the section state for this iteration.
void
Lin2DUserDefined::UpdateState(Eigen::VectorXd strain, unsigned int cond){
    //Update the matrial behavior.
    Eigen::VectorXd mStrain(1);
    mStrain << strain(0);
    theMaterial->UpdateState(mStrain, cond);

    //Update the section strain.
    Strain = strain;
}
