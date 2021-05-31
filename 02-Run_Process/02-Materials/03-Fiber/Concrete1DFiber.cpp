#include <cfloat>
#include "Concrete1DFiber.hpp" 
#include "Definitions.hpp"

//Overload constructor.
Concrete1DFiber::Concrete1DFiber(double fc, double ecc, double fcu, double ecu, double ratio, double ft, double Et, double nu, double Rho) : 
Material("Concrete1DFiber", false), nu(nu), Et(Et), fc(fc), ft(ft), fcu(fcu), ecu(ecu), ecc(ecc), Rho(Rho), lambda(ratio){
    //Initialize fiber internal variables.
	newStrain = 0.0;
	newStress = 0.0;
	newMinStrain = 0.0; 
	newMaxStrain = 0.0;
	newTangentStiffness = 2.0*fc/ecc;

	oldStrain = 0.0;
	oldStress = 0.0;
	oldTangentStiffness = 2.0*fc/ecc;
}

//Destructor.
Concrete1DFiber::~Concrete1DFiber(){
    //Does nothing.
}

//Clone the 'Concrete1DFiber' material.
std::unique_ptr<Material>
Concrete1DFiber::CopyMaterial(){
    return std::make_unique<Concrete1DFiber>(fc, ecc, fcu, ecu, lambda, ft, Et, nu, Rho);
}

//Access material density.
double 
Concrete1DFiber::GetDensity() const{
    return Rho;
}

//Returns the Poisson's ratio.
double 
Concrete1DFiber::GetPoissonRatio() const{
    return nu;
}

//Access bulk modulus.
double 
Concrete1DFiber::GetBulkModulus() const{
	double E = 2.0*fc/ecc;
    return E/3.0/(1.0 - 2.0*nu);
}

//Gets the linear Shear Modulus.
double 
Concrete1DFiber:: GetShearModulus() const{
	double Ec = 2.0*fc/ecc;
    return Ec/(2.0*(1.0 + nu));
}

//Access modulus of elasticity.
double 
Concrete1DFiber::GetElasticityModulus() const{
    return 2.0*fc/ecc;
}

//Access the material's energy at current strain.
double 
Concrete1DFiber::GetEnergy() const{
	//TODO: Compute material energy 
    double W = 0.0;
    return W;
}

//Returns the material viscous damping.
Eigen::MatrixXd 
Concrete1DFiber::GetDamping() const{
    //Compute the damping.
    Eigen::MatrixXd Damping(1,1);
    Damping.fill(0.0); 

    return Damping;
}

//Returns material strain vector.
Eigen::VectorXd
Concrete1DFiber::GetStrain() const{
    Eigen::VectorXd Strain(1);
    Strain << oldStrain;
    return Strain;
}

//Returns material stress vector.
Eigen::VectorXd
Concrete1DFiber::GetStress() const{
	Eigen::VectorXd Stress(1);
    Stress << oldStress;

    return Stress;
}

//Returns material strain rate vector.
Eigen::VectorXd 
Concrete1DFiber::GetStrainRate()  const{
    //Compute the strain rate.
    Eigen::VectorXd StrainRate(1);
    StrainRate.fill(0.0); 

    return StrainRate;
}

//Computes the material total stress.
Eigen::VectorXd 
Concrete1DFiber::GetTotalStress() const{
    return GetStress();
}

//Computes consistent material matrix.
Eigen::MatrixXd
Concrete1DFiber::GetTangentStiffness() const{
	Eigen::MatrixXd TangentStiffness(1,1);
    TangentStiffness << oldTangentStiffness;

    return TangentStiffness;
}

//Returns the initial material stiffness.
Eigen::MatrixXd 
Concrete1DFiber::GetInitialTangentStiffness() const{
	Eigen::MatrixXd TangentStiffness(1,1);
    TangentStiffness << 2.0*fc/ecc;

    return TangentStiffness;
}

//Perform converged material state update.
void 
Concrete1DFiber::CommitState(){
    newStrain = oldStrain;
	newStress = oldStress;
	newTangentStiffness = oldTangentStiffness;
	newMinStrain = oldMinStrain;
	newMaxStrain = oldMaxStrain;
}

//Reverse the material states to previous converged state.
void 
Concrete1DFiber::ReverseState(){
    oldStrain = newStrain;
	oldStress = newStress;
	oldTangentStiffness = newTangentStiffness;
	oldMinStrain = newMinStrain;
	oldMaxStrain = newMaxStrain;
}

//Brings the material states to its initial state in the element.
void 
Concrete1DFiber::InitialState(){
    newStrain = 0.0;
	newStress = 0.0;
	newMinStrain = 0.0; 
	newMaxStrain = 0.0;
	newTangentStiffness = 2.0*fc/ecc;

	oldStrain = 0.0;
	oldStress = 0.0;
	oldTangentStiffness = 2.0*fc/ecc;
}

//Update the material state for this iteration.
void
Concrete1DFiber::UpdateState(const Eigen::VectorXd strain, unsigned int UNUSED(cond)){
    //Update the strain and internal variables.
	double ec0 = 2.0*fc/ecc;

	oldMinStrain = newMinStrain;
	oldMaxStrain = newMaxStrain;

	oldStrain = strain(0);
	double deps = oldStrain - newStrain;

	//Computes current strain:
	if(fabs(deps) < 10.0*DBL_EPSILON)
		return;

	if(oldStrain < oldMinStrain){
		//Reset the minimum current strain 
		Compressive_Envelope(oldStrain, oldStress, oldTangentStiffness); 
		oldMinStrain = oldStrain;

	}
	else{

		//unloading-reloading branch
		double reloadStrain = (fcu - lambda*ec0*ecu)/(ec0*(1.0 - lambda)); 
		double reloadStress = ec0*reloadStrain;

		//The previous minimum stress: 
		double minStress;
		double tanSlope;
		Compressive_Envelope(oldMinStrain, minStress, tanSlope); 

		//Computes the reloading slope:
		double Er = (minStress - reloadStress)/(oldMinStrain - reloadStrain);
		double ept  = oldMinStrain - minStress/Er;

		if(oldStrain <= ept){

			double sigmin = minStress + Er*(oldStrain - oldMinStrain);
			double sigmax = 0.5*Er*(oldStrain - ept);

			oldStress =  newStress + ec0*deps;
			oldTangentStiffness = ec0;

			if(oldStress <= sigmin){
				oldStress = sigmin;
				oldTangentStiffness = Er;
			}
			if(oldStress >= sigmax){
				oldStress = sigmax;
				oldTangentStiffness = 0.5*Er;
			}
		}
		else{

			//Maximum remaining tensile strength
			double epn = ept + oldMaxStrain;
			
			if(oldStrain <= epn){
                double remainTensileStress;
				Tensile_Envelope(oldMaxStrain, remainTensileStress, oldTangentStiffness); 

				if(oldMaxStrain > 0.0){
					oldTangentStiffness = remainTensileStress/oldMaxStrain;
				}
				else{
					oldTangentStiffness = ec0;
				}

				oldStress = oldTangentStiffness*(oldStrain - ept);
			}
			else{

				//Tensile envelope curve shifted by ept:
				double tempStrain = oldStrain - ept;
				Tensile_Envelope(tempStrain, oldStress, oldTangentStiffness); 

				oldMaxStrain = oldStrain - ept;
			}
		}
	}
}

//Compute back-bone tensile envelope
void 
Concrete1DFiber::Tensile_Envelope(double eccn, double &fccn, double &Etn){
	double Ec = 2.0*fc/ecc;
	double eps0 = ft/Ec;
	double epst = 0.0;

	if(Et > 0.0)
		epst = ft/Et;

	double epsu = eps0 + epst; 
	if(eccn <= eps0){
		//Linear Positive Branch
		Etn  = Ec;
		fccn = eccn*Ec;
	}
	else{
		//Linear Negative Branch
		if(eccn <= epsu){
			Etn  = -Et;
			fccn = ft - Et*(eccn - eps0);
		}
		else{
			Etn  = 0.0;
			fccn = 0.0;
		}
	}
}

//Compute back-bone compressive envelope
void 
Concrete1DFiber::Compressive_Envelope(double eccn, double &fccn, double &Ecn){
	double Ec = 2.0*fc/ecc;
	double Ratio = eccn/ecc;

	if(eccn >= ecc){
		//Parabolic Branch
		Ecn  = Ec*(1.0 - Ratio);
		fccn = fc*Ratio*(2.0 - Ratio);
	}
	else{
		if(eccn > ecu){
			//Linear Descending Branch
			Ecn  = (fcu - fc)/(ecu - ecc);
			fccn = fc + (fcu - fc)*(eccn - ecc)/(ecu - ecc);
		}
		else{
			//Constant Flat Friction Branch
			Ecn  = 0.0;
			fccn = fcu;
        }
	}
}
