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

//Access the material's energy at current strain.
double 
Hertzian1DLinear::GetEnergy() const{
    double e2 = Strain(0)*Strain(0);
    double W  = 1.0/2.0*k1*e2 + 1.0/3.0*k2*e2*Strain(0) + 1.0/4.0*k3*e2*e2;
    return W;
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
    //The material initial stiffness matrix.
    Eigen::MatrixXd InitialTangentStiffness(1,1);
    InitialTangentStiffness << k1;

    return InitialTangentStiffness;
}

//Perform converged material state update.
void 
Hertzian1DLinear::CommitState(){
}

//Reverse the material states to previous converged state.
void 
Hertzian1DLinear::ReverseState(){
    //TODO: Get back to previous commited state
}

//Brings the material states to its initial state in the element.
void 
Hertzian1DLinear::InitialState(){
    Strain.fill(0.0);
    Stress.fill(0.0);
}

//Update the material state for this iteration.
void
Hertzian1DLinear::UpdateState(const Eigen::VectorXd strain, const unsigned int cond){
    //Updates the elastic/plastic material components.    
    if(cond == 1){
        //Update the strain.
        Strain = strain;

        //Update the stress.
        Stress << k1*strain(0) + k2*strain(0)*strain(0) + k3*strain(0)*strain(0)*strain(0);

        //Consistent Tangent Operator.
        TangentStiffness << k1 + 2.0*k2*strain(0) + 3.0*strain(0)*strain(0);
    }
}