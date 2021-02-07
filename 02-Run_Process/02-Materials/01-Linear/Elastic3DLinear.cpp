#include "Elastic3DLinear.hpp"
#include "Definitions.hpp"

//Overload Constructor.
Elastic3DLinear::Elastic3DLinear(const double E, const double nu, const double rho) : 
Material("Elastic3DLinear", false), E(E), nu(nu), Rho(rho){
    //Initialize material strain.
    Strain.resize(6);
    Strain.fill(0.0);

    //Initialize material stress.
    Stress.resize(6);
    Stress.fill(0.0);

    //Material model coefficients:
    TangentStiffness.resize(6,6);

    //Material model coefficients:
    double c1 = E*(1.0 - nu)/(1.0 - 2.0*nu)/(1.0 + nu);
    double c2 = E*nu/(1.0 - 2.0*nu)/(1.0 + nu);
    double c3 = E/(2.0*(1.0 + nu));

    TangentStiffness << c1 ,  c2,  c2, 0.0, 0.0, 0.0,
                        c2 ,  c1,  c2, 0.0, 0.0, 0.0,
                        c2 ,  c2,  c1, 0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0,  c3, 0.0, 0.0,
                        0.0, 0.0, 0.0, 0.0,  c3, 0.0,
                        0.0, 0.0, 0.0, 0.0, 0.0,  c3;
}

//Destructor.
Elastic3DLinear::~Elastic3DLinear(){
    //Does nothing.
}

//Clone the 'Elastic3DLinear' material.
std::unique_ptr<Material>
Elastic3DLinear::CopyMaterial(){
    return std::make_unique<Elastic3DLinear>(E, nu, Rho);
}

//Access material density.
double 
Elastic3DLinear::GetDensity() const{
    return Rho;
}

//Returns the Poisson's ratio.
double 
Elastic3DLinear::GetPoissonRatio() const{
    return nu;
}

//Access bulk modulus.
double 
Elastic3DLinear::GetBulkModulus() const{
    return E/3.0/(1.0 - 2.0*nu);
}

//Access shear modulus.
double 
Elastic3DLinear::GetShearModulus() const{
    return E/2.0/(1 + nu);
}

//Access modulus of elasticity.
double 
Elastic3DLinear::GetElasticityModulus() const{
    return E;
}

//Access the material's energy at current strain.
double 
Elastic3DLinear::GetEnergy() const{
    double W = 1.0/2.0*Strain.transpose()*TangentStiffness*Strain;
    return W;
}

//Returns the material viscous damping.
Eigen::MatrixXd 
Elastic3DLinear::GetDamping() const{
    //Compute the damping.
    Eigen::MatrixXd Damping(6,6);
    Damping.fill(0.0); 

    return Damping;
}

//Returns the material strain.
Eigen::VectorXd
Elastic3DLinear::GetStrain() const{
    return Strain;
}

//Returns the material stress.
Eigen::VectorXd
Elastic3DLinear::GetStress() const{
    return Stress;
}

//Returns material strain rate vector.
Eigen::VectorXd 
Elastic3DLinear::GetStrainRate() const{
    //Compute the strain rate.
    Eigen::VectorXd StrainRate(6);
    StrainRate.fill(0.0); 

    return StrainRate;
}

//Computes the material total stress.
Eigen::VectorXd 
Elastic3DLinear::GetTotalStress() const{
    return Stress;
}

//Returns the material stiffness.
Eigen::MatrixXd
Elastic3DLinear::GetTangentStiffness() const{
    return TangentStiffness;
}

//Returns the initial material stiffness.
Eigen::MatrixXd
Elastic3DLinear::GetInitialTangentStiffness() const{
    return TangentStiffness;
}

//Perform converged material state update.
void 
Elastic3DLinear::CommitState(){
}

//Update the material state for this iteration.
void
Elastic3DLinear::UpdateState(const Eigen::VectorXd strain, const unsigned int cond){
    //Updates the elatic/platic material components.    
    if(cond == 1){
        //Update the strain.
        Strain = strain;

        //Update the stress.
        Stress = TangentStiffness*strain;
    }
}
