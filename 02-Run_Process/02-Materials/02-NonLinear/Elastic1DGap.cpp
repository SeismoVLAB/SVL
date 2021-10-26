#include <cfloat>
#include "Elastic1DGap.hpp" 
#include "Definitions.hpp"

//Define constant tolerance value:
const double TOL = 1.0E-06;

//Overload constructor.
Elastic1DGap::Elastic1DGap(double E, double gap, bool behavior) : 
Material("Elastic1DGap", false), E(E), Behavior(behavior){  
    //Initialize fiber internal variables.
    Gap = fabs(gap),
	oldStrain = 0.0;
	newStrain = 0.0;
}

//Destructor.
Elastic1DGap::~Elastic1DGap(){
    //Does nothing.
}

//Clone the 'Elastic1DGap' material.
std::unique_ptr<Material>
Elastic1DGap::CopyMaterial(){
    return std::make_unique<Elastic1DGap>(E, Gap, Behavior);
}

//Access material density.
double 
Elastic1DGap::GetDensity() const{
    return 0.0;
}

//Returns the Poisson's ratio.
double 
Elastic1DGap::GetPoissonRatio() const{
    return 0.0;
}

//Access bulk modulus.
double 
Elastic1DGap::GetBulkModulus() const{
    return 0.0;
}

//Gets the linear Shear Modulus.
double 
Elastic1DGap:: GetShearModulus() const{
    return 0.0;
}

//Access modulus of elasticity.
double 
Elastic1DGap::GetElasticityModulus() const{
    return E;
}

//Access the material's energy at current strain.
double 
Elastic1DGap::GetEnergy() const{
	//TODO: Compute the energy
    double W = 0.0;
    return W;
}

//Returns the material viscous damping.
Eigen::MatrixXd 
Elastic1DGap::GetDamping() const{
    //Compute the damping.
    Eigen::MatrixXd Damping(1,1);
    Damping.fill(0.0); 

    return Damping;
}

//Returns material strain vector.
Eigen::VectorXd
Elastic1DGap::GetStrain() const{
    Eigen::VectorXd Strain(1);
    Strain << oldStrain;
    return Strain;
}

//Returns material stress vector.
Eigen::VectorXd
Elastic1DGap::GetStress() const{
    //Material stress vector
	Eigen::VectorXd Stress(1);
    Stress << E*(oldStrain - Gap)*(oldStrain > Gap)*(Behavior) + E*(oldStrain + Gap)*(oldStrain < -Gap)*(!Behavior);

    return Stress;
}

//Returns material strain rate vector.
Eigen::VectorXd 
Elastic1DGap::GetStrainRate()  const{
    //Compute the strain rate.
    Eigen::VectorXd StrainRate(1);
    StrainRate.fill(0.0); 

    return StrainRate;
}

//Computes the material total stress.
Eigen::VectorXd 
Elastic1DGap::GetTotalStress() const{
    return GetStress();
}

//Computes consistent material matrix.
Eigen::MatrixXd
Elastic1DGap::GetTangentStiffness() const{
    //Consistent material stiffness
	Eigen::MatrixXd TangentStiffness(1,1);
    TangentStiffness << E*(oldStrain > Gap)*(Behavior) + E*(oldStrain < -Gap)*(!Behavior);

    return TangentStiffness;
}

//Returns the initial material stiffness.
Eigen::MatrixXd 
Elastic1DGap::GetInitialTangentStiffness() const{
	Eigen::MatrixXd TangentStiffness(1,1);
    TangentStiffness << E*(Gap < TOL);

    return TangentStiffness;
}

//Perform converged material state update.
void 
Elastic1DGap::CommitState(){
    newStrain = oldStrain;
}

//Reverse the material states to previous converged state.
void 
Elastic1DGap::ReverseState(){
    oldStrain = newStrain;
}

//Brings the material states to its initial state in the element.
void 
Elastic1DGap::InitialState(){
	oldStrain = 0.0;
    newStrain = 0.0;
}

//Update the material state for this iteration.
void
Elastic1DGap::UpdateState(const Eigen::VectorXd strain, unsigned int UNUSED(cond)){
    //Update the strain and internal variables.
	oldStrain = strain(0);
}
