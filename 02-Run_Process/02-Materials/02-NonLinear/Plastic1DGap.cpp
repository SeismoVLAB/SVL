#include <cfloat>
#include "Plastic1DGap.hpp" 
#include "Definitions.hpp"

//Define constant tolerance value:
const double TOL = 1.0E-06;

//Overload constructor.
Plastic1DGap::Plastic1DGap(double E, double Sy, double gap, double eta, bool behavior) : 
Material("Plastic1DGap", false), E(E), Ratio(eta), Behavior(behavior){
    //Sets the proper signs according to behavior.
    fy  = 1.0*(Behavior)*fabs(Sy)  - 1.0*(!Behavior)*fabs(Sy);
    Gap = 1.0*(Behavior)*fabs(gap) - 1.0*(!Behavior)*fabs(gap);

    //Initialize Plastic internal variables.
    minYieldStrain = Gap;
    maxYieldStrain = Gap + fy/E;

	//Elastic Gap history variables:
    newStrain = 0.0;
    oldStrain = 0.0;
}

//Destructor.
Plastic1DGap::~Plastic1DGap(){
    //Does nothing.
}

//Clone the 'Plastic1DGap' material.
std::unique_ptr<Material>
Plastic1DGap::CopyMaterial(){
    return std::make_unique<Plastic1DGap>(E, fy, Gap, Ratio, Behavior);
}

//Access material density.
double 
Plastic1DGap::GetDensity() const{
    return 0.0;
}

//Returns the Poisson's ratio.
double 
Plastic1DGap::GetPoissonRatio() const{
    return 0.0;
}

//Access bulk modulus.
double 
Plastic1DGap::GetBulkModulus() const{
    return 0.0;
}

//Gets the linear Shear Modulus.
double 
Plastic1DGap:: GetShearModulus() const{
    return 0.0;
}

//Access modulus of elasticity.
double 
Plastic1DGap::GetElasticityModulus() const{
    return E;
}

//Access the material's energy at current strain.
double 
Plastic1DGap::GetEnergy() const{
	//TODO: Compute the energy
    double W = 0.0;
    return W;
}

//Returns the material viscous damping.
Eigen::MatrixXd 
Plastic1DGap::GetDamping() const{
    //Compute the damping.
    Eigen::MatrixXd Damping(1,1);
    Damping.fill(0.0); 

    return Damping;
}

//Returns material strain vector.
Eigen::VectorXd
Plastic1DGap::GetStrain() const{
    Eigen::VectorXd Strain(1);
    Strain << oldStrain;

    return Strain;
}

//Returns material stress vector.
Eigen::VectorXd
Plastic1DGap::GetStress() const{
    //Elastic and Plastic Stresses
    double Se = E*(oldStrain - minYieldStrain);
    double Sp = fy + (oldStrain - Gap - fy/E)*Ratio*E;

    //Material stress vector
	Eigen::VectorXd Stress(1);

    if(Behavior){
        Stress << Se*(oldStrain <= maxYieldStrain)*(minYieldStrain < oldStrain) + Sp*(oldStrain > maxYieldStrain);
    }
    else{
        Stress << Se*(oldStrain < minYieldStrain)*(maxYieldStrain <= oldStrain) + Sp*(oldStrain < maxYieldStrain);
    }

    return Stress;
}

//Returns material strain rate vector.
Eigen::VectorXd 
Plastic1DGap::GetStrainRate()  const{
    //Compute the strain rate.
    Eigen::VectorXd StrainRate(1);
    StrainRate.fill(0.0); 

    return StrainRate;
}

//Computes the material total stress.
Eigen::VectorXd 
Plastic1DGap::GetTotalStress() const{
    return GetStress();
}

//Computes consistent material matrix.
Eigen::MatrixXd
Plastic1DGap::GetTangentStiffness() const{
    double Ee = E;
    double Ep = Ratio*E;

    //Consistent material stiffness
    Eigen::MatrixXd TangentStiffness(1,1);

    if(Behavior){
        TangentStiffness << Ee*(oldStrain <= maxYieldStrain)*(minYieldStrain < oldStrain) + Ep*(oldStrain > maxYieldStrain);
    }
    else{
        TangentStiffness << Ee*(oldStrain < minYieldStrain)*(maxYieldStrain <= oldStrain) + Ep*(oldStrain < maxYieldStrain);
    }

    return TangentStiffness;
}

//Returns the initial material stiffness.
Eigen::MatrixXd 
Plastic1DGap::GetInitialTangentStiffness() const{
	Eigen::MatrixXd TangentStiffness(1,1);
    TangentStiffness << E*(fabs(Gap) < TOL);

    return TangentStiffness;
}

//Perform converged material state update.
void 
Plastic1DGap::CommitState(){
    //Computes the associated stress
    Eigen::VectorXd oldStress = GetStress();

    if(Behavior){
        if(oldStrain > maxYieldStrain){
            maxYieldStrain = oldStrain;
            minYieldStrain = oldStrain - oldStress(0)/E; 
        }
        else if(oldStrain < minYieldStrain && Gap < oldStrain){
            maxYieldStrain = fy/E + (oldStrain - Ratio*Gap)/(1.0 - Ratio);
            minYieldStrain = oldStrain;
        }
    }
    else{
        if(oldStrain < maxYieldStrain){
            maxYieldStrain = oldStrain;
            minYieldStrain = oldStrain - oldStress(0)/E; 
        }
        else if(minYieldStrain < oldStrain  && oldStrain < Gap){
            maxYieldStrain = fy/E + (oldStrain - Ratio*Gap)/(1.0 - Ratio);
            minYieldStrain = oldStrain;
        }
    }

    newStrain = oldStrain;
}

//Reverse the material states to previous converged state.
void 
Plastic1DGap::ReverseState(){
    oldStrain = newStrain;
}

//Brings the material states to its initial state in the element.
void 
Plastic1DGap::InitialState(){
    //Initialize fiber internal variables.
    minYieldStrain = Gap;
    maxYieldStrain = Gap + fy/E;

	//Elastic Gap history variables:
    newStrain = 0.0;
    oldStrain = 0.0;
}

//Update the material state for this iteration.
void
Plastic1DGap::UpdateState(const Eigen::VectorXd strain, unsigned int UNUSED(cond)){
    //Update the strain and internal variables.
	oldStrain = strain(0);
}
