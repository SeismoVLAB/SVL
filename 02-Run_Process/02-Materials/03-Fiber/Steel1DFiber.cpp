#include <cfloat>
#include "Steel1DFiber.hpp" 
#include "Definitions.hpp"

//Overload constructor.
Steel1DFiber::Steel1DFiber(double fy, double E, double b, double R0, double cR1, double cR2, double a1, double a2, double a3, double a4, double nu, double Rho) : 
Material("Steel1DFiber", false), nu(nu), fy(fy), E(E), b(b), R0(R0), cR1(cR1), cR2(cR2), a1(a1), a2(a2), a3(a3), a4(a4), Rho(Rho) {
    //Initialize fiber internal variables.
	newLoadUnloadIndex = 0;
	oldLoadUnloadIndex = 0;

	oldStrain = 0.0;
	oldStress = 0.0;
	oldTangentStiffness = E;

	newStrain = 0.0;
	newStress = 0.0;
	newTangentStiffness = E;

	newMinStrain = -fy/E;
	newMaxStrain =  fy/E;
	newPlasticStrain = 0.0;
	newStrainInterception = 0.0;
	newStressInterception = 0.0;
	newStrainLastInversion = 0.0;
	newStressLastInversion = 0.0;
}

//Destructor.
Steel1DFiber::~Steel1DFiber(){
    //Does nothing.
}

//Clone the 'Steel1DFiber' material.
std::unique_ptr<Material>
Steel1DFiber::CopyMaterial(){
    return std::make_unique<Steel1DFiber>(fy, E, b, R0, cR1, cR2, a1, a2, a3, a4, nu, Rho);
}

//Access material density.
double 
Steel1DFiber::GetDensity() const{
    return Rho;
}

//Returns the Poisson's ratio.
double 
Steel1DFiber::GetPoissonRatio() const{
    return nu;
}

//Access bulk modulus.
double 
Steel1DFiber::GetBulkModulus() const{
    return E/3.0/(1.0 - 2.0*nu);
}

//Gets the linear Shear Modulus.
double 
Steel1DFiber:: GetShearModulus() const{
    return E/(2.0*(1.0 + nu));
}

//Access modulus of elasticity.
double 
Steel1DFiber::GetElasticityModulus() const{
    return E;
}

//Access the material's energy at current strain.
double 
Steel1DFiber::GetEnergy() const{
	//TODO: Compute the energy
    double W = 0.0;
    return W;
}

//Returns the material viscous damping.
Eigen::MatrixXd 
Steel1DFiber::GetDamping() const{
    //Compute the damping.
    Eigen::MatrixXd Damping(1,1);
    Damping.fill(0.0); 

    return Damping;
}

//Returns material strain vector.
Eigen::VectorXd
Steel1DFiber::GetStrain() const{
    Eigen::VectorXd Strain(1);
    Strain << oldStrain;
    return Strain;
}

//Returns material stress vector.
Eigen::VectorXd
Steel1DFiber::GetStress() const{
	Eigen::VectorXd Stress(1);
    Stress << oldStress;

    return Stress;
}

//Returns material strain rate vector.
Eigen::VectorXd 
Steel1DFiber::GetStrainRate()  const{
    //Compute the strain rate.
    Eigen::VectorXd StrainRate(1);
    StrainRate.fill(0.0); 

    return StrainRate;
}

//Computes the material total stress.
Eigen::VectorXd 
Steel1DFiber::GetTotalStress() const{
    return GetStress();
}

//Computes consistent material matrix.
Eigen::MatrixXd
Steel1DFiber::GetTangentStiffness() const{
	Eigen::MatrixXd TangentStiffness(1,1);
    TangentStiffness << oldTangentStiffness;

    return TangentStiffness;
}

//Returns the initial material stiffness.
Eigen::MatrixXd 
Steel1DFiber::GetInitialTangentStiffness() const{
	Eigen::MatrixXd TangentStiffness(1,1);
    TangentStiffness << E;

    return TangentStiffness;
}

//Perform converged material state update.
void 
Steel1DFiber::CommitState(){
	newStrain = oldStrain;
	newStress = oldStress;
	newMinStrain = oldMinStrain;
	newMaxStrain = oldMaxStrain;
	newPlasticStrain = oldPlasticStrain;
	newTangentStiffness = oldTangentStiffness;
	newStrainInterception  = oldStrainInterception;
	newStressInterception  = oldStressInterception;
	newStrainLastInversion = oldStrainLastInversion;
	newStressLastInversion = oldStressLastInversion;

	newLoadUnloadIndex = oldLoadUnloadIndex;
}

//Reverse the material states to previous converged state.
void 
Steel1DFiber::ReverseState(){
	oldStrain = newStrain;
	oldStress = newStress;
	oldMinStrain = newMinStrain;
	oldMaxStrain = newMaxStrain;
	oldPlasticStrain = newPlasticStrain;
	oldTangentStiffness = newTangentStiffness;
	oldStrainInterception  = newStrainInterception;
	oldStressInterception  = newStressInterception;
	oldStrainLastInversion = newStrainLastInversion;
	oldStressLastInversion = newStressLastInversion;

	oldLoadUnloadIndex = newLoadUnloadIndex;
}

//Brings the material states to its initial state in the element.
void 
Steel1DFiber::InitialState(){
	oldStrain = 0.0;
	oldStress = 0.0;
	oldTangentStiffness = E;

	newStrain = 0.0;
	newStress = 0.0;
	newMaxStrain = fy/E;
	newMinStrain = -newMaxStrain;
	newPlasticStrain = 0.0;
	newTangentStiffness = E;
	newStrainInterception  = 0.0;
	newStressInterception  = 0.0;
	newStrainLastInversion = 0.0;
	newStressLastInversion = 0.0;

	newLoadUnloadIndex = 0;
}

//Update the material state for this iteration.
void
Steel1DFiber::UpdateState(const Eigen::VectorXd strain, unsigned int UNUSED(cond)){
    //Update the strain and internal variables.
	double Esh  = b*E;
	double epsy = fy/E;

	oldStrain = strain(0);

	double deps = oldStrain - newStrain;

	oldMinStrain = newMinStrain;
	oldMaxStrain = newMaxStrain;
	oldPlasticStrain = newPlasticStrain;
	oldStrainInterception  = newStrainInterception;
	oldStressInterception  = newStressInterception;
	oldStrainLastInversion = newStrainLastInversion;
	oldStressLastInversion = newStressLastInversion;
	oldLoadUnloadIndex = newLoadUnloadIndex;

	if(oldLoadUnloadIndex == 0 || oldLoadUnloadIndex == 3){

		if(fabs(deps) < 10.0*DBL_EPSILON){
			oldTangentStiffness = E;
			oldStress = 0.0;
			oldLoadUnloadIndex = 3;
			return; 
		}
		else{

			oldMaxStrain =  epsy;
			oldMinStrain = -epsy;

			if(deps < 0.0){
				oldLoadUnloadIndex = 2;
				oldStrainInterception = oldMinStrain;
				oldStressInterception = -fy;
				oldPlasticStrain = oldMinStrain;
			}
			else{
				oldLoadUnloadIndex = 1;
				oldStrainInterception = oldMaxStrain;
				oldStressInterception = fy;
				oldPlasticStrain = oldMaxStrain;
			}
		}
	}

	//Load reversal from negative to positive strain increment
	if(oldLoadUnloadIndex == 2 && deps > 0.0){

		oldLoadUnloadIndex = 1;
		oldStrainLastInversion = newStrain;
		oldStressLastInversion = newStress;

		if(newStrain < oldMinStrain)
			oldMinStrain = newStrain;

		double d1 = (oldMaxStrain - oldMinStrain)/(2.0*a4*epsy);
		double shift = 1.0 + a3*pow(d1, 0.8);

		oldStrainInterception = (fy*shift - Esh*epsy*shift - oldStressLastInversion + E*oldStrainLastInversion)/(E - Esh);
		oldStressInterception = fy*shift + Esh*(oldStrainInterception - epsy*shift);
		oldPlasticStrain = oldMaxStrain;
	}
	else if(oldLoadUnloadIndex == 1 && deps < 0.0){

		oldLoadUnloadIndex = 2;
		oldStrainLastInversion = newStrain;
		oldStressLastInversion = newStress;

		if(newStrain > oldMaxStrain)
			oldMaxStrain = newStrain;

		double d1 = (oldMaxStrain - oldMinStrain)/(2.0*a2*epsy);
		double shift = 1.0 + a1*pow(d1, 0.8);

		oldStrainInterception = (-fy*shift + Esh*epsy*shift - oldStressLastInversion + E*oldStrainLastInversion)/(E - Esh);
		oldStressInterception = -fy*shift + Esh*(oldStrainInterception + epsy*shift);
		oldPlasticStrain = oldMinStrain;
	}

	//Compute current Strain, Stress, and TangetStiffness
	double xi   = fabs((oldPlasticStrain - oldStrainInterception)/epsy);
	double R    = R0*(1.0 - (cR1*xi)/(cR2 + xi));
	double epsr = (oldStrain - oldStrainLastInversion)/(oldStrainInterception - oldStrainLastInversion);
	double fac1 = 1.0 + pow(fabs(epsr),R);
	double fac2 = pow(fac1, 1.0/R);

	oldStress = b*epsr + (1.0 - b)*epsr/fac2;
	oldStress = oldStress*(oldStressInterception - oldStressLastInversion) + oldStressLastInversion;

	oldTangentStiffness = b + (1.0 - b)/(fac1*fac2);
	oldTangentStiffness = oldTangentStiffness*(oldStressInterception - oldStressLastInversion)/(oldStrainInterception - oldStrainLastInversion);
}
