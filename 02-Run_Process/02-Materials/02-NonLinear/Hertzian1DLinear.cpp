#include "Hertzian1DLinear.hpp"
#include "Definitions.hpp"

//Overload constructor.
Hertzian1DLinear::Hertzian1DLinear(const double k1, const double k2, const double k3, const double rho) : 
Material("Hertzian1DLinear", false), k1(k1), k2(k2), k3(k3), Rho(rho){
    //Initialize material strain.
    Strain.resize(1);
    Strain << 0.0;

    //Initialize material stress.
    Stress.resize(1);
    Stress << 0.0;

    //Initialize material stiffness.
    TangentStiffness.resize(1,1);
    TangentStiffness << k1;
}

//Destructor.
Hertzian1DLinear::~Hertzian1DLinear(){
    //Does nothing.
}

//Clone the 'Hertzian1DLinear' material.
std::unique_ptr<Material>
Hertzian1DLinear::CopyMaterial(){
    return std::make_unique<Hertzian1DLinear>(k1, k2, k3, Rho);
}

//Access material density.
double 
Hertzian1DLinear::GetDensity() const{
    return Rho;
}

//Returns the Poisson's ratio.
double 
Hertzian1DLinear::GetPoissonRatio() const{
    return 0.0;
}

//Access bulk modulus.
double 
Hertzian1DLinear::GetBulkModulus() const{
    return 0.0;
}

//Access shear modulus.
double 
Hertzian1DLinear::GetShearModulus() const{
    return 0.0;
}

//Access modulus of elasticity.
double 
Hertzian1DLinear::GetElasticityModulus() const{
    return 0.0;
}

//Returns the material viscous damping.
Eigen::MatrixXd 
Hertzian1DLinear::GetDamping() const{
    //Compute the damping.
    Eigen::MatrixXd Damping(1,1);
    Damping.fill(0.0); 

    return Damping;
}

//Returns material strain vector.
Eigen::VectorXd
Hertzian1DLinear::GetStrain() const{
    return Strain;
}

//Returns material stress vector.
Eigen::VectorXd
Hertzian1DLinear::GetStress() const{
    return Stress;
}

//Returns material strain rate vector.
Eigen::VectorXd 
Hertzian1DLinear::GetStrainRate() const{
    //Compute the strain rate.
    Eigen::VectorXd StrainRate(1);
    StrainRate.fill(0.0); 

    return StrainRate;
}

//Computes the material total stress.
Eigen::VectorXd 
Hertzian1DLinear::GetTotalStress() const{
    return Stress;
}

//Computes consistent material matrix.
Eigen::MatrixXd
Hertzian1DLinear::GetTangentStiffness() const{
    return TangentStiffness;
}

//Returns the initial material stiffness.
Eigen::MatrixXd
Hertzian1DLinear::GetInitialTangentStiffness() const{
    return TangentStiffness;
}

//Perform converged material state update.
void 
Hertzian1DLinear::CommitState(){
}

//Update the material state for this iteration.
void
Hertzian1DLinear::UpdateState(const Eigen::VectorXd strain, const unsigned int cond){
    //Updates the elatic/platic material components.    
    if(cond == 1){
        //Update the strain.
        Strain = strain;

        //Update the stress.
        Stress = k1*strain + k2*strain*strain + k3*strain*strain*strain;
    }
}
