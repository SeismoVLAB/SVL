#include "Plastic1DJ2.hpp"
#include "Definitions.hpp"

//Overload constructor.
Plastic1DJ2::Plastic1DJ2(const double E, const double nu, const double rho, const double K, const double H, const double SigmaY) :
Material("Plastic1DJ2", false), E(E), nu(nu), Rho(rho), K(K), H(H), SigmaY(SigmaY){
    //Initialize internal hardening variable.
    alpha = 0.0;

    //Initialize strain.
    Strain.resize(1);
    Strain.fill(0.0);

    //Initialize strain.
    Stress.resize(1);
    Stress.fill(0.0);

    //Initialize plastic strain.
    PlasticStrain.resize(1);
    PlasticStrain.fill(0.0);
    
    //Initialize back stress.
    BackStress.resize(1);
    BackStress.fill(0.0);
    
    //Initialize consistent tangent stiffness.
    TangentStiffness.resize(1,1);
    TangentStiffness << E;
}

//Destructor.
Plastic1DJ2::~Plastic1DJ2(){
    //Does nothing.
}

//Clone the 'Plastic1DJ2' material.
std::unique_ptr<Material>
Plastic1DJ2::CopyMaterial(){
    return std::make_unique<Plastic1DJ2>(E, nu, Rho, K, H, SigmaY);
}

//Access material density.
double 
Plastic1DJ2::GetDensity() const{
    return Rho;
}

//Returns the Poisson's ratio.
double 
Plastic1DJ2::GetPoissonRatio() const{
    return nu;
}

//Access bulk modulus.
double 
Plastic1DJ2::GetBulkModulus() const{
    return E/3.0/(1.0 - 2.0*nu);
}

//Access shear modulus.
double 
Plastic1DJ2::GetShearModulus() const{
    return E/2.0/(1 + nu);
}

//Access modulus of elasticity.
double 
Plastic1DJ2::GetElasticityModulus() const{
    return E;
}

//Returns the material viscous damping.
Eigen::MatrixXd 
Plastic1DJ2::GetDamping() const{
    //Compute the damping.
    Eigen::MatrixXd Damping(1,1);
    Damping.fill(0.0); 

    return Damping;
}

//Returns material strain vector.
Eigen::VectorXd
Plastic1DJ2::GetStrain() const{
    return Strain;
}

//Returns material stress vector.
Eigen::VectorXd
Plastic1DJ2::GetStress() const{
    return Stress;
}

//Returns material strain rate vector.
Eigen::VectorXd 
Plastic1DJ2::GetStrainRate() const{
    //Compute the strain rate.
    Eigen::VectorXd StrainRate(1);
    StrainRate.fill(0.0); 

    return StrainRate;
}

//Computes the material total stress.
Eigen::VectorXd 
Plastic1DJ2::GetTotalStress() const{
    return Stress;
}

//Returns consistent tangent stiffness matrix.
Eigen::MatrixXd
Plastic1DJ2::GetTangentStiffness() const{
    return TangentStiffness;
}

//Returns the initial material stiffness.
Eigen::MatrixXd
Plastic1DJ2::GetInitialTangentStiffness() const{
    //Creates the stiffness.
    Eigen::MatrixXd InitialTangentStiffness(1,1);
    InitialTangentStiffness << E;

    return InitialTangentStiffness;
}

//Perform converged material state update.
void 
Plastic1DJ2::CommitState(){
}

//Update the material state for this iteration.
void
Plastic1DJ2::UpdateState(const Eigen::VectorXd strain, const unsigned int cond){
    //Updates the elatic/platic material components.    
    if(cond == 1){
        //Trial stress tensor.
        Eigen::VectorXd TrialStress = E*(strain - PlasticStrain);

        //Trial relative stress tensor.
        Eigen::VectorXd TrialRelStress = TrialStress - BackStress;

        //Direction of the trial relative stress tensor.
        Eigen::VectorXd TrialRelStressSign = TrialRelStress/TrialRelStress.norm();

        //Trial flow condition
        double TrialF = TrialRelStress.norm() - SigmaY - K*alpha;
    
        if (TrialF <= 0){
            //Elastic regime: Consistent tangent stiffness.
            TangentStiffness << E;

        }
        else{
            //Consistency parameter.
            double DeltaGamma;

            //Plastic regime: Return mapping.
            DeltaGamma    =   TrialF/(E + K + H);
            TrialStress   += -DeltaGamma*E*TrialRelStressSign;
            PlasticStrain +=  DeltaGamma*TrialRelStressSign;
            BackStress    +=  DeltaGamma*H*TrialRelStressSign;
            alpha         +=  DeltaGamma;

            //Consistent tangent stiffness.
            TangentStiffness << E*(K + H)/(E + K + H);
        }

        //Update total strain.
        Strain = strain;
        Stress = TrialStress;
    }
}
