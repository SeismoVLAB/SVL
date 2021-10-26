#include "Elastic1DFiber.hpp" 
#include "Definitions.hpp"

//Overload constructor.
Elastic1DFiber::Elastic1DFiber(const double E, const double nu, const double rho) : 
Material("Elastic1DFiber", false), E(E), nu(nu), Rho(rho){
    //Initialize material strain.
    newStrain.resize(1);
    newStrain << 0.0;

    oldStrain.resize(1);
    oldStrain << 0.0;
}

//Destructor.
Elastic1DFiber::~Elastic1DFiber(){
    //Does nothing.
}

//Clone the 'Elastic1DFiber' material.
std::unique_ptr<Material>
Elastic1DFiber::CopyMaterial(){
    return std::make_unique<Elastic1DFiber>(E, nu, Rho);
}

//Access material density.
double 
Elastic1DFiber::GetDensity() const{
    return Rho;
}

//Returns the Poisson's ratio.
double 
Elastic1DFiber::GetPoissonRatio() const{
    return nu;
}

//Access bulk modulus.
double 
Elastic1DFiber::GetBulkModulus() const{
    return E/3.0/(1.0 - 2.0*nu);
}

//Gets the linear Shear Modulus.
double 
Elastic1DFiber:: GetShearModulus() const{
    return E/2.0/(1 + nu);
}

//Access modulus of elasticity.
double 
Elastic1DFiber::GetElasticityModulus() const{
    return E;
}

//Access the material's energy at current strain.
double 
Elastic1DFiber::GetEnergy() const{
    double W = 1.0/2.0*oldStrain(0)*E*oldStrain(0);
    return W;
}

//Returns the material viscous damping.
Eigen::MatrixXd 
Elastic1DFiber::GetDamping() const{
    //Compute the damping.
    Eigen::MatrixXd Damping(1,1);
    Damping.fill(0.0); 

    return Damping;
}

//Returns material strain vector.
Eigen::VectorXd
Elastic1DFiber::GetStrain() const{
    return oldStrain;
}

//Returns material stress vector.
Eigen::VectorXd
Elastic1DFiber::GetStress() const{
	Eigen::VectorXd Stress = E*oldStrain;

    return Stress;
}

//Returns material strain rate vector.
Eigen::VectorXd 
Elastic1DFiber::GetStrainRate()  const{
    //Compute the strain rate.
    Eigen::VectorXd StrainRate(1);
    StrainRate.fill(0.0); 

    return StrainRate;
}

//Computes the material total stress.
Eigen::VectorXd 
Elastic1DFiber::GetTotalStress() const{
    return GetStress();
}

//Computes consistent material matrix.
Eigen::MatrixXd
Elastic1DFiber::GetTangentStiffness() const{
	Eigen::MatrixXd TangentStiffness(1,1);
    TangentStiffness << E;

    return TangentStiffness;
}

//Returns the initial material stiffness.
Eigen::MatrixXd 
Elastic1DFiber::GetInitialTangentStiffness() const{
	Eigen::MatrixXd TangentStiffness(1,1);
    TangentStiffness << E;

    return TangentStiffness;
}

//Perform converged material state update.
void 
Elastic1DFiber::CommitState(){
    newStrain = oldStrain;
}

//Reverse the material states to previous converged state.
void 
Elastic1DFiber::ReverseState(){
    oldStrain = newStrain;
}

//Brings the material states to its initial state in the element.
void 
Elastic1DFiber::InitialState(){
    oldStrain << 0.0;
    newStrain << 0.0;
}

//Update the material state for this iteration.
void
Elastic1DFiber::UpdateState(const Eigen::VectorXd strain, unsigned int UNUSED(cond)){
    //Update the strain.
    oldStrain << strain(0);
}
