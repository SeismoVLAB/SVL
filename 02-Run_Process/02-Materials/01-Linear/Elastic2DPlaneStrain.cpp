#include "Elastic2DPlaneStrain.hpp"
#include "Definitions.hpp"

//Overload Constructor.
Elastic2DPlaneStrain::Elastic2DPlaneStrain(const double E, const double nu, const double rho) : 
Material("Elastic2DPlaneStrain", false), E(E), nu(nu), Rho(rho){
    //Initialize material strain.
    Strain.resize(3);
    Strain.fill(0.0);

    //Initialize material commited strain.
    newStrain.resize(3);
    newStrain.fill(0.0);

    //Material model coefficients:
    TangentStiffness.resize(3,3);

    //Material model coefficients.
    double c1 = E*(1.0 - nu)/(1.0 - 2.0*nu)/(1.0 + nu);
    double c2 = E*nu/(1.0 - 2.0*nu)/(1.0 + nu); // lame parameters, lambda
    double c3 = E/(2.0*(1.0 + nu));             // lame parameters, mu

    TangentStiffness.resize(3,3);
    TangentStiffness << c1 , c2 , 0.0,
                        c2 , c1 , 0.0,
                        0.0, 0.0,  c3;
}

//Destructor.
Elastic2DPlaneStrain::~Elastic2DPlaneStrain(){
    //Does nothing.
}

//Clone the 'Elastic2DPlaneStrain' material.
std::unique_ptr<Material>
Elastic2DPlaneStrain::CopyMaterial(){
    return std::make_unique<Elastic2DPlaneStrain>(E, nu, Rho);
}


//Access material density.
double 
Elastic2DPlaneStrain::GetDensity() const{
    return Rho;
}

//Returns the Poisson's ratio.
double 
Elastic2DPlaneStrain::GetPoissonRatio() const{
    return nu;
}

//Access bulk modulus.
double 
Elastic2DPlaneStrain::GetBulkModulus() const{
    return E/3.0/(1.0 - 2.0*nu);
}

//Access shear modulus.
double 
Elastic2DPlaneStrain::GetShearModulus() const{
    return E/2.0/(1 + nu);
}

//Access modulus of elasticity.
double 
Elastic2DPlaneStrain::GetElasticityModulus() const{
    return E;
}

//Access the material's energy at current strain.
double 
Elastic2DPlaneStrain::GetEnergy() const{
    double W = 1.0/2.0*Strain.transpose()*TangentStiffness*Strain;
    return W;
}

//Returns the material viscous damping.
Eigen::MatrixXd 
Elastic2DPlaneStrain::GetDamping() const{
    //Compute the damping.
    Eigen::MatrixXd Damping(3,3);
    Damping.fill(0.0); 

    return Damping;
}

//Returns the material strain.
Eigen::VectorXd
Elastic2DPlaneStrain::GetStrain() const{
    return Strain;
}

//Returns the material stress.
Eigen::VectorXd
Elastic2DPlaneStrain::GetStress() const{
    Eigen::VectorXd Stress = TangentStiffness*Strain;
    return Stress;
}

//Returns material strain rate vector.
Eigen::VectorXd 
Elastic2DPlaneStrain::GetStrainRate() const{
    //Compute the strain rate.
    Eigen::VectorXd StrainRate(3);
    StrainRate.fill(0.0); 

    return StrainRate;
}

//Computes the material total stress.
Eigen::VectorXd 
Elastic2DPlaneStrain::GetTotalStress() const{
    Eigen::VectorXd Stress = TangentStiffness*Strain;
    return Stress;
}

//Returns the material stiffness.
Eigen::MatrixXd
Elastic2DPlaneStrain::GetTangentStiffness() const{
    return TangentStiffness;
}

//Returns the initial material stiffness.
Eigen::MatrixXd 
Elastic2DPlaneStrain::GetInitialTangentStiffness() const{
    return TangentStiffness;
}

//Perform converged material state update.
void 
Elastic2DPlaneStrain::CommitState(){
    newStrain = Strain;
}

//Reverse the material states to previous converged state.
void 
Elastic2DPlaneStrain::ReverseState(){
    Strain = newStrain;
}

//Brings the material states to its initial state in the element.
void 
Elastic2DPlaneStrain::InitialState(){
    Strain.fill(0.0);
    newStrain.fill(0.0);
}

//Update the material state for this iteration.
void
Elastic2DPlaneStrain::UpdateState(const Eigen::VectorXd strain, const unsigned int cond){
    //Updates the elastic/plastic material components.    
    if(cond == 1){
        //Update the strain.
        Strain = strain;
    }
}
