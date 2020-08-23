#include <iostream>
#include "Viscous1DLinear.hpp" 
#include "Definitions.hpp"

//Overload constructor.
Viscous1DLinear::Viscous1DLinear(const double eta) : 
Material("Viscous1DLinear", true), eta(eta){
    //Initialize material strain.
    StrainRate.resize(1);
    StrainRate << 0.0;

    //Initialize material stress.
    Stress.resize(1);
    Stress << 0.0;

    //Material model coefficients:
    Damping.resize(1,1);
    Damping << eta; 
}

//Destructor.
Viscous1DLinear::~Viscous1DLinear(){
    //Does nothing.
}

//Clone the 'Viscous1DLinear' material.
std::unique_ptr<Material>
Viscous1DLinear::CopyMaterial(){
    return std::make_unique<Viscous1DLinear>(eta);
}

//Access material density.
double 
Viscous1DLinear::GetDensity() const{
    return 0.0;
}

//Returns the Poisson's ratio.
double 
Viscous1DLinear::GetPoissonRatio() const{
    return 0.0;
}

//Access bulk modulus.
double 
Viscous1DLinear::GetBulkModulus() const{
    return 0.0;
}

//Access shear modulus.
double 
Viscous1DLinear::GetShearModulus() const{
    return 0.0;
}

//Access modulus of elasticity.
double 
Viscous1DLinear::GetElasticityModulus() const{
    return 0.0;
}

//Returns consistent material matrix.
Eigen::MatrixXd
Viscous1DLinear::GetDamping() const{
    return Damping;
}

//Returns material strain vector.
Eigen::VectorXd
Viscous1DLinear::GetStrain() const{
    //Compute the strain rate.
    Eigen::VectorXd Strain(1);
    Strain << 0.0; 

    return Strain;
}

//Returns material strain rate vector.
Eigen::VectorXd
Viscous1DLinear::GetStrainRate() const{
    return StrainRate;
}

//Returns material stress vector.
Eigen::VectorXd
Viscous1DLinear::GetStress() const{
    //Compute the elastic stress.
    Eigen::VectorXd theStress(1);
    theStress << 0.0; 

    return theStress;
}

//Computes the material total stress.
Eigen::VectorXd 
Viscous1DLinear::GetTotalStress() const{
    return Stress;
}

//Returns consistent material matrix.
Eigen::MatrixXd
Viscous1DLinear::GetTangentStiffness() const{
    //Creates the stiffness.
    Eigen::MatrixXd TangentStiffness(1,1);
    TangentStiffness << 0.0;

    return TangentStiffness;
}

//Returns the initial material stiffness.
Eigen::MatrixXd
Viscous1DLinear::GetInitialTangentStiffness() const{
    //Creates the stiffness.
    Eigen::MatrixXd TangentStiffness(1,1);
    TangentStiffness << 0.0; 

    return TangentStiffness;
}

//Perform converged material state update.
void 
Viscous1DLinear::CommitState(){
}

//Update the material state for this iteration.
void
Viscous1DLinear::UpdateState(const Eigen::VectorXd strainrate, const unsigned int cond){
    //Updates the viscous material components.    
    if(cond == 2){
        //Update the strain.
        StrainRate = strainrate;

        //Update the stress.
        Stress = eta*strainrate;
    }
}
