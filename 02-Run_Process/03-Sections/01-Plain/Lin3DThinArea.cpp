#include <cmath>
#include "Lin3DThinArea.hpp"
#include "Definitions.hpp"

//Overload Constructor.
Lin3DThinArea::Lin3DThinArea(double t, std::unique_ptr<Material> &material) : 
Section("Lin3DThinArea"), t(t) {
    //Initialize material strain.
    Strain.resize(6);
    Strain << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

    //The section material:
    theMaterial = material->CopyMaterial();
}

//Destructor.
Lin3DThinArea::~Lin3DThinArea(){
    //Does nothing.
}

//Clone the 'Lin3DThinArea' section.
std::unique_ptr<Section>
Lin3DThinArea::CopySection(){
    return std::make_unique<Lin3DThinArea>(t, theMaterial);
}

//Returns the section generalized strain.
Eigen::VectorXd
Lin3DThinArea::GetStrain(){
    return Strain;
}

//Returns the section generalized stress.
Eigen::VectorXd
Lin3DThinArea::GetStress(){
    //Forces = [Nx, Ny, Nxy, Mx, My, Mxy]
    Eigen::VectorXd Stress = GetTangentStiffness()*Strain;
    return Stress;
}

//Access the material rotational density.
Eigen::MatrixXd 
Lin3DThinArea::GetDensity(){
    //Gets the material density.
    double Rho = theMaterial->GetDensity();

    //Cross-Section Area.
    double I = t*t*t/12.0;

    //Density matrix for membrane/plate effect.
    Eigen::MatrixXd Density(6,6);
    Density << t*Rho,  0.0 ,  0.0,  0.0 ,  0.0 ,  0.0 ,
               0.0 , t*Rho,  0.0 ,  0.0 ,  0.0 ,  0.0 ,
               0.0 ,  0.0 , t*Rho,  0.0 ,  0.0 ,  0.0 ,
               0.0 ,  0.0 ,  0.0 , t*Rho,  0.0 ,  0.0 ,
               0.0 ,  0.0 ,  0.0 ,  0.0 , I*Rho,  0.0 ,
               0.0 ,  0.0 ,  0.0 ,  0.0 ,  0.0 , I*Rho;
    
    return Density;
}


//Returns the section axial stiffness.
Eigen::MatrixXd
Lin3DThinArea::GetTangentStiffness(){
    //Section Inertia.
    double I = t*t*t/12.0;

    //Material tangent stiffness matrix.
    Eigen::MatrixXd C = theMaterial->GetTangentStiffness(); 

    //Stiffness matrix for membrane/plate effect.
    Eigen::MatrixXd SectionStiffness(6,6);
    SectionStiffness << t*C(0,0), t*C(0,1), t*C(0,2),   0.0   ,   0.0   ,   0.0   ,
                        t*C(1,0), t*C(1,1), t*C(1,2),   0.0   ,   0.0   ,   0.0   ,
                        t*C(2,0), t*C(2,1), t*C(2,2),   0.0   ,   0.0   ,   0.0   ,
                          0.0   ,   0.0   ,   0.0   , I*C(0,0), I*C(0,1), I*C(0,2),
                          0.0   ,   0.0   ,   0.0   , I*C(1,0), I*C(1,1), I*C(1,2),
                          0.0   ,   0.0   ,   0.0   , I*C(2,0), I*C(2,1), I*C(2,2);

    return SectionStiffness;
}

//Returns the section initial tangent stiffness matrix.
Eigen::MatrixXd 
Lin3DThinArea::GetInitialTangentStiffness(){
    return GetTangentStiffness();
}

//Returns the section strain at given position.
Eigen::VectorXd 
Lin3DThinArea::GetStrainAt(double x3, double UNUSED(x2)){
    //The stress vector in local coordinates
    Eigen::VectorXd theStrain(6);
    theStrain.fill(0.0);

    //Checks if coordinate is inside the section.
    if(fabs(x3) <= t/2.0){ 
        x3 = x3 - t/2.0;

        // Epsilon = [exx, eyy, 0.0, exy, 0.0, 0.0
        theStrain << Strain(0) + x3*Strain(3), 
                     Strain(1) + x3*Strain(4), 
                     0.0, 
                     Strain(2) + 2.0*x3*Strain(5), 
                     0.0, 
                     0.0; 
    }

    return theStrain;
}

//Returns the section stress at given position.
Eigen::VectorXd 
Lin3DThinArea::GetStressAt(double x3, double UNUSED(x2)){
    //The stress vector in local coordinates
    Eigen::VectorXd theStress(6);
    theStress.fill(0.0);

    //Checks if coordinate is inside the section.
    if(fabs(x3) > t/2.0){
        //Transform coordinates into section local axis.
        x3 = x3 - t/2.0;

        //Gets material properties.s
        double E  = theMaterial->GetElasticityModulus();
        double nu = theMaterial->GetPoissonRatio();

        //TODO: Compute the strain at point, needs shear value. 
        // Sigma = [Sxx, Syy, 0.0, txy, 0.0, 0.0]
        Eigen::VectorXd theStress(6);
        theStress << E/(1.0 -nu*nu)*(Strain(0) + nu*Strain(1) + x3*(Strain(3) + nu*Strain(4))),
                     E/(1.0 -nu*nu)*(Strain(1) + nu*Strain(0) + x3*(Strain(4) + nu*Strain(3))),
                     0.0,
                     E/(2.0*(1.0 + nu))*(Strain(2) + 2.0*x3*Strain(5)),
                     0.0,
                     0.0;
    }

    return theStress;
}

//Perform converged section state update.
void 
Lin3DThinArea::CommitState(){
    theMaterial->CommitState();
}

//Reverse the section states to previous converged state.
void 
Lin3DThinArea::ReverseState(){
    theMaterial->ReverseState();
}

//Brings the section states to its initial state.
void 
Lin3DThinArea::InitialState(){
    Strain.fill(0.0);
    theMaterial->InitialState();
}

//Update the section state for this iteration.
void
Lin3DThinArea::UpdateState(Eigen::VectorXd strain, unsigned int cond){
    //Update the matrial behaviour.
    Eigen::VectorXd mStrain(3);
    mStrain << strain(0), strain(1), strain(2);
    theMaterial->UpdateState(mStrain, cond);

    //Update the strain.
    Strain = strain;
}