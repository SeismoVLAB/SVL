#include <cmath>
#include "Plastic3DJ2.hpp"
#include "Definitions.hpp"

//Overload constructor.
Plastic3DJ2::Plastic3DJ2(const double K, const double G, const double rho, const double Hbar, const double beta, const double SigmaY) :
Material("Plastic3DJ2", false), K(K), G(G), Rho(rho), H(Hbar), beta(beta), SigmaY(SigmaY){
    //Initialize internal hardening variable.
    alpha = 0.0;

    //Initialize strain.
    Strain.resize(6);
    Strain.fill(0.0);

    //Initialize strain.
    Stress.resize(6);
    Stress.fill(0.0);

    //Initialize plastic strain.
    PlasticStrain.resize(6);
    PlasticStrain.fill(0.0);
    
    //Initialize back stress.
    BackStress.resize(6);
    BackStress.fill(0.0);

    //Fourth-rank identity tensor.
    Eigen::MatrixXd D = ComputeIdentityTensor();

    //Rank-four symmetric identity operator.
    Eigen::MatrixXd I = ComputeIdentityOperator();

    //Initialize stiffness matrix.
    TangentStiffness.resize(6,6);
    TangentStiffness = K*D + 2.0*G*(I - 1.0/3.0*D);
}

//Destructor.
Plastic3DJ2::~Plastic3DJ2(){
    //Does nothing.
}

//Clone the 'Plastic3DJ2' material.
std::unique_ptr<Material>
Plastic3DJ2::CopyMaterial(){
    return std::make_unique<Plastic3DJ2>(K, G, Rho, H, beta, SigmaY);
}

//Access material density.
double 
Plastic3DJ2::GetDensity() const{
    return Rho;
}

//Returns the Poisson's ratio.
double 
Plastic3DJ2::GetPoissonRatio() const{
    return (3.0*K - 2.0*G)/(2.0*(3.0*K + G));
}

//Access bulk modulus.
double 
Plastic3DJ2::GetBulkModulus() const{
    return K;
}

//Access shear modulus.
double 
Plastic3DJ2::GetShearModulus() const{
    return G;
}

//Access modulus of elasticity.
double 
Plastic3DJ2::GetElasticityModulus() const{
    return 9.0*K*G/(3.0*K + G);
}

//Access the material's energy at current strain.
double 
Plastic3DJ2::GetEnergy() const{
    //TODO: Compute/write the energy density for this material 
    return 0.0;
}

//Returns the material viscous damping.
Eigen::MatrixXd 
Plastic3DJ2::GetDamping() const{
    //Compute the damping.
    Eigen::MatrixXd Damping(6,6);
    Damping.fill(0.0); 

    return Damping;
}

//Returns material strain vector.
Eigen::VectorXd
Plastic3DJ2::GetStrain() const{
    return Strain;
}

//Returns material stress vector.
Eigen::VectorXd
Plastic3DJ2::GetStress() const{
    return Stress;
}

//Returns material strain rate vector.
Eigen::VectorXd 
Plastic3DJ2::GetStrainRate() const{
    //Compute the strain rate.
    Eigen::VectorXd StrainRate(6);
    StrainRate.fill(0.0); 

    return StrainRate;
}

//Computes the material total stress.
Eigen::VectorXd 
Plastic3DJ2::GetTotalStress() const{
    return Stress;
}

//Returns consistent tangent stiffness matrix.
Eigen::MatrixXd
Plastic3DJ2::GetTangentStiffness() const{
    return TangentStiffness;
}

//Returns the initial material stiffness.
Eigen::MatrixXd
Plastic3DJ2::GetInitialTangentStiffness() const{
    //Fourth-rank identity tensor.
    Eigen::MatrixXd D = ComputeIdentityTensor();

    //Rank-four symmetric identity operator.
    Eigen::MatrixXd I = ComputeIdentityOperator();

    //Initialize stiffness matrix.
    Eigen::MatrixXd InitialTangentStiffness = K*D + 2.0*G*(I - 1.0/3.0*D);

    return InitialTangentStiffness;
}

//Perform converged material state update.
void 
Plastic3DJ2::CommitState(){
}

//Update the material state for this iteration.
void
Plastic3DJ2::UpdateState(const Eigen::VectorXd strain, const unsigned int cond){
    //Updates the elatic/plastic material components.    
    if(cond == 1){
        //Auxiliar tensors.
        Eigen::MatrixXd D   = ComputeIdentityTensor();
        Eigen::VectorXd One = ComputeIdentityVector();        
        Eigen::MatrixXd Id  = ComputeDeviatoricTensor();
        
        //Engineering strain tensor.
        Strain << strain(0), strain(1), strain(2), 0.5*strain(3), 0.5*strain(4), 0.5*strain(5); 

        //Strain trace.
        double StrainTrace = ComputeTensorTrace(Strain); 

        //Deviatoric strain.
        Eigen::VectorXd DeviatoricStrain = Strain - 1.0/3.0*StrainTrace*One;

        //Trial deviatoric stress tensor.
        Eigen::VectorXd TrialDeviatoricStress = 2.0*G*(DeviatoricStrain - PlasticStrain);

        //Trial relative stress tensor.
        Eigen::VectorXd TrialRelativeStress = TrialDeviatoricStress - BackStress;

        //Direction of the trial relative stress tensor.
        double TrialStressNorm = ComputeTensorNorm(TrialRelativeStress);  

        //Trial flow condition.
        double TrialF = TrialStressNorm - sqrt(2.0/3.0)*(SigmaY + alpha*beta*H); 
    
        if (TrialF <= 0){
            //Elastic regime.
            TangentStiffness = K*D + 2.0*G*Id;
            Stress = K*StrainTrace*One + TrialDeviatoricStress;
        }
        else{
            //Consistency parameter.
            double DeltaGamma;

            //Unit normal vector.
            Eigen::VectorXd n = TrialRelativeStress/TrialStressNorm; 

            //Plastic regime. Return mapping [1].
            DeltaGamma     =  TrialF/(2.0*G + 2.0/3.0*H);
            alpha         +=  sqrt(2.0/3.0)*DeltaGamma;
            BackStress    +=  2.0/3.0*(1.0 - beta)*H*DeltaGamma*n;
            PlasticStrain +=  DeltaGamma*n;
            Stress         =  K*StrainTrace*One + TrialDeviatoricStress - 2.0*G*DeltaGamma*n; 

            //Consistent tangent stiffness [2].
            TangentStiffness = K*D + 2.0*G*(Id - 1.0/(1.0 + H/G/3.0)*n*n.transpose()) - 4.0*G*G*DeltaGamma/TrialStressNorm*(Id - n*n.transpose());
        }
    }
}

//Second Order Identity Tensor.
Eigen::VectorXd 
Plastic3DJ2::ComputeIdentityVector() const{
    //Second Order Identity Tensor.
    Eigen::VectorXd One(6);
    One << 1.0, 1.0, 1.0, 0.0, 0.0, 0.0;

    return One;
}

//Constructs deviatoric tensor.
Eigen::MatrixXd 
Plastic3DJ2::ComputeDeviatoricTensor() const{
    //Fourth-rank identity tensor.
    Eigen::MatrixXd D = ComputeIdentityTensor();

    //Rank-four symmetric identity operator.
    Eigen::MatrixXd I = ComputeIdentityOperator();

    //Deviatoric tensor.
    Eigen::MatrixXd P = I - 1.0/3.0*D;

    return P;
}

//Constructs fourth-rank identity tensor.
Eigen::MatrixXd 
Plastic3DJ2::ComputeIdentityTensor() const{
    //Fourth-rank identity tensor.
    Eigen::MatrixXd D(6,6);

    D << 1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
         1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
         1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

    return D;
}

//Constructs rank-four symmetric identity operator.
Eigen::MatrixXd 
Plastic3DJ2::ComputeIdentityOperator() const{
    //Rank-four symmetric identity operator.
    Eigen::MatrixXd I(6,6);

    I << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.5, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.5, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.5;

    return I;
}

//Computes trace of tensor.
double 
Plastic3DJ2::ComputeTensorTrace(const Eigen::VectorXd &T){
    double Tr = T(0) + T(1) + T(2); 
    return Tr;
}

//Computes tensor norm.
double 
Plastic3DJ2::ComputeTensorNorm(const Eigen::VectorXd &T){
    double Norm = sqrt(T(0)*T(0) + T(1)*T(1) + T(2)*T(2) + 2.0*(T(3)*T(3) + T(4)*T(4) + T(5)*T(5)));
    return Norm;
}
