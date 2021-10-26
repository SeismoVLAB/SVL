#include "Elastic1DLinear.hpp" 
#include "Definitions.hpp"

//Overload constructor.
Elastic1DLinear::Elastic1DLinear(const double E, const double nu, const double rho) : 
Material("Elastic1DLinear", false), E(E), nu(nu), Rho(rho){
    //Initialize material strain.
    Strain.resize(1);
    Strain << 0.0;

    //Initialize material commited strain.
    newStrain.resize(1);
    newStrain << 0.0;

    //Initialize material stiffness.
    TangentStiffness.resize(1,1);
    TangentStiffness << E;
}

//Destructor.
Elastic1DLinear::~Elastic1DLinear(){
    //Does nothing.
}

//Clone the 'Elastic1DLinear' material.
std::unique_ptr<Material>
Elastic1DLinear::CopyMaterial(){
    return std::make_unique<Elastic1DLinear>(E, nu, Rho);
}

//Access material density.
double 
Elastic1DLinear::GetDensity() const{
    return Rho;
}

//Returns the Poisson's ratio.
double 
Elastic1DLinear::GetPoissonRatio() const{
    return nu;
}

//Access bulk modulus.
double 
Elastic1DLinear::GetBulkModulus() const{
    return E/3.0/(1.0 - 2.0*nu);
}

//Access shear modulus.
double 
Elastic1DLinear::GetShearModulus() const{
    return E/2.0/(1 + nu);
}

//Access modulus of elasticity.
double 
Elastic1DLinear::GetElasticityModulus() const{
    return E;
}

//Access the material's energy at current strain.
double 
Elastic1DLinear::GetEnergy() const{
    double W = 1.0/2.0*Strain(0)*TangentStiffness(0,0)*Strain(0);
    return W;
}

//Returns the material viscous damping.
Eigen::MatrixXd 
Elastic1DLinear::GetDamping() const{
    //Compute the damping.
    Eigen::MatrixXd Damping(1,1);
    Damping.fill(0.0); 

    return Damping;
}

//Returns material strain vector.
Eigen::VectorXd
Elastic1DLinear::GetStrain() const{
    return Strain;
}

//Returns material stress vector.
Eigen::VectorXd
Elastic1DLinear::GetStress() const{
    Eigen::VectorXd Stress = TangentStiffness*Strain;
    return Stress;
}

//Returns material strain rate vector.
Eigen::VectorXd 
Elastic1DLinear::GetStrainRate()  const{
    //Compute the strain rate.
    Eigen::VectorXd StrainRate(1);
    StrainRate.fill(0.0); 

    return StrainRate;
}

//Computes the material total stress.
Eigen::VectorXd 
Elastic1DLinear::GetTotalStress() const{
    Eigen::VectorXd Stress = TangentStiffness*Strain;
    return Stress;
}

//Computes consistent material matrix.
Eigen::MatrixXd
Elastic1DLinear::GetTangentStiffness() const{
    return TangentStiffness;
}

//Returns the initial material stiffness.
Eigen::MatrixXd 
Elastic1DLinear::GetInitialTangentStiffness() const{
    return TangentStiffness;
}

//Perform converged material state update.
void 
Elastic1DLinear::CommitState(){
    newStrain = Strain;
}

//Reverse the material states to previous converged state.
void 
Elastic1DLinear::ReverseState(){
    Strain = newStrain;
}

//Brings the material states to its initial state in the element.
void 
Elastic1DLinear::InitialState(){
    Strain << 0.0;
    newStrain << 0.0;
}

//Update the material state for this iteration.
void
Elastic1DLinear::UpdateState(const Eigen::VectorXd strain, const unsigned int cond){
    //Updates the elastic/plastic material components.    
    if(cond == 1){
        //Update the strain.
        Strain = strain;
    }
}
