#include <math.h>
#include <iostream>
#include <stdlib.h>
#include "PlasticPlaneStrainBA.hpp"
#include "Definitions.hpp"

//Overload constructor.
PlasticPlaneStrainBA::PlasticPlaneStrainBA(const double K, const double G, const double rho, const double H0, const double h, const double m, const double Su, const double beta) :
Material("PlasticPlaneStrainBA", false), K(K), G(G), Rho(rho), H0(H0), h(h), m(m), Su(Su), beta(beta){
    //Initialize strain.
    Strain.resize(6);
    Strain.fill(0.0);

    Strain_n.resize(6);
    Strain_n.fill(0.0);

    //Initialize trial stress;
    Stress.resize(6);
    Stress.fill(0.0);

    //Initialize stress_n
    Stress_n.resize(6);
    Stress_n.fill(0.0);

    //Initialize deviatoric stress at F0.
    DeviatoricStress0.resize(6);
    DeviatoricStress0.fill(0.0);
    
    //Fourth-rank identity tensor.
    Eigen::MatrixXd D = ComputeIdentityTensor();

    //Rank-four deviatoric identity tensor.
    Eigen::MatrixXd Id = ComputeDeviatoricTensor();

    //Initialize stiffness matrix.
    TangentStiffness.resize(6,6); 
    TangentStiffness = K*D + 2.0*G*Id;

    //Bounding surface radius
    R = sqrt(8.0/3.0)*Su;

    //initialize hardening internal variables
    psi = 2.0*G; 
    kappa = 1.0E12;

    FirstLoadFlag = 0;

    rootFlag = 0;
}

//Destructor.
PlasticPlaneStrainBA::~PlasticPlaneStrainBA(){
    //Does nothing.
}

//Clone the 'PlasticPlaneStrainBA' material.
std::unique_ptr<Material>
PlasticPlaneStrainBA::CopyMaterial(){
    return std::make_unique<PlasticPlaneStrainBA>(K, G, Rho, H0, h, m, Su, beta);
}

//Access material density.
double 
PlasticPlaneStrainBA::GetDensity() const{
    return Rho;
}

//Returns the Poisson's ratio.
double 
PlasticPlaneStrainBA::GetPoissonRatio() const{
    return (3.0*K - 2.0*G)/(2.0*(3.0*K + G));
}

//Access bulk modulus.
double 
PlasticPlaneStrainBA::GetBulkModulus() const{
    return K;
}

//Access shear modulus.
double 
PlasticPlaneStrainBA::GetShearModulus() const{
    return G;
}

//Access modulus of elasticity.
double 
PlasticPlaneStrainBA::GetElasticityModulus() const{
    return 9.0*K*G/(3.0*K + G);
}

//Access the material's energy at current strain.
double 
PlasticPlaneStrainBA::GetEnergy() const{
    //TODO: Compute/write the energy density for this material 
    return 0.0;
}

//Returns the material viscous damping.
Eigen::MatrixXd 
PlasticPlaneStrainBA::GetDamping() const{
    //Compute the damping.
    Eigen::MatrixXd Damping(3,3);
    Damping.fill(0.0); 

    return Damping;
}

//Returns material strain vector.
Eigen::VectorXd
PlasticPlaneStrainBA::GetStrain() const{
    Eigen::VectorXd StrainPlaneStrain(3);
    StrainPlaneStrain << Strain(1), Strain(2), Strain(4);

    return StrainPlaneStrain;
}

//Returns material stress vector.
Eigen::VectorXd
PlasticPlaneStrainBA::GetStress() const{
    Eigen::VectorXd StressPlaneStrain(3);
    StressPlaneStrain << Stress(1), Stress(2), Stress(4);

    return StressPlaneStrain;
}

//Returns material strain rate vector.
Eigen::VectorXd 
PlasticPlaneStrainBA::GetStrainRate() const{
    //Compute the strain rate.
    Eigen::VectorXd StrainRate(3);
    StrainRate.fill(0.0); 

    return StrainRate;
}

//Computes the material total stress.
Eigen::VectorXd 
PlasticPlaneStrainBA::GetTotalStress() const{
    Eigen::VectorXd StressPlaneStrain(3);
    StressPlaneStrain << Stress(1), Stress(2), Stress(4);

    return StressPlaneStrain;
}

//Returns consistent tangent stiffness matrix.
Eigen::MatrixXd
PlasticPlaneStrainBA::GetTangentStiffness() const{
    Eigen::MatrixXd TangentStiffnessPlaneStrain(3,3);
    TangentStiffnessPlaneStrain << TangentStiffness(1,1), TangentStiffness(1,2), TangentStiffness(1,4),
                                   TangentStiffness(2,1), TangentStiffness(2,2), TangentStiffness(2,4),
                                   TangentStiffness(4,1), TangentStiffness(4,2), TangentStiffness(4,4); 

    return TangentStiffnessPlaneStrain;
}

//Returns the initial material stiffness.
Eigen::MatrixXd
PlasticPlaneStrainBA::GetInitialTangentStiffness() const{
    //Fourth-rank identity tensor.
    Eigen::MatrixXd D = ComputeIdentityTensor();

    //Rank-four deviatoric identity tensor.
    Eigen::MatrixXd Id = ComputeDeviatoricTensor();

    //Initialize stiffness matrix.
    Eigen::MatrixXd Initial3DStiffness = K*D + 2.0*G*Id;

    //The material plain strain matrix.
    Eigen::MatrixXd InitialTangentStiffness(3,3);
    InitialTangentStiffness <<  Initial3DStiffness(1,1), Initial3DStiffness(1,2), Initial3DStiffness(1,4),
                                Initial3DStiffness(2,1), Initial3DStiffness(2,2), Initial3DStiffness(2,4),
                                Initial3DStiffness(4,1), Initial3DStiffness(4,2), Initial3DStiffness(4,4); 

    return InitialTangentStiffness;
}

//Perform converged material state update.
void 
PlasticPlaneStrainBA::CommitState(){
    Stress_n = Stress;
    Strain_n = Strain;
}

//Reverse the material states to previous converged state.
void 
PlasticPlaneStrainBA::ReverseState(){
    //TODO: Get back to previous commited state
}

//Brings the material states to its initial state in the element.
void 
PlasticPlaneStrainBA::InitialState(){
    psi   = 2.0*G; 
    kappa = 1.0E12;

    rootFlag = 0;
    FirstLoadFlag = 0;

    Strain.fill(0.0);
    Stress.fill(0.0);
    Strain_n.fill(0.0);
    Stress_n.fill(0.0);
    DeviatoricStress0.fill(0.0);

    Eigen::MatrixXd D = ComputeIdentityTensor();
    Eigen::MatrixXd Id = ComputeDeviatoricTensor();
    TangentStiffness = K*D + 2.0*G*Id;
}

//Update the material state for this iteration.
void
PlasticPlaneStrainBA::UpdateState(const Eigen::VectorXd strain, const unsigned int cond) {
//Updates the elastic/plastic material components.    
    if(cond == 1) {

        //Engineering strain tensor.
        Strain << 0.0, strain(0), strain(1), 0.0, 0.5*strain(2), 0.0; 
        
        //Second-rank identity tensor.
        Eigen::VectorXd One = ComputeIdentityVector();
        //Fourth-rank identity tensor.
        Eigen::MatrixXd D   = ComputeIdentityTensor();
        //Rank-four deviatoric identity tensor.
        Eigen::MatrixXd Id  = ComputeDeviatoricTensor();

        //Incremental strain
        Eigen::VectorXd IncrStrain = Strain - Strain_n;

        Eigen::VectorXd DeviatoricIncrStrain = IncrStrain - 1.0/3.0*ComputeTensorTrace(IncrStrain)*One;
        
        //Deviatoric stress tensor.
        Eigen::VectorXd DeviatoricStress = Stress - 1.0/3.0*ComputeTensorTrace(Stress)*One;
        
        Eigen::VectorXd DeviatoricStrain = Strain - 1.0/3.0*ComputeTensorTrace(Strain)*One;

        double StrainNorm = ComputeTensorNorm(DeviatoricIncrStrain);
        
        double infty = 1.0E12 ;
        
        double DSTol = 1.0E-12;

        Hn = h*pow(kappa,m);

        if (FirstLoadFlag == 0) {// to check if this is first loading
            FirstLoadFlag = 1;
            Stress = Stress_n + K*ComputeTensorTrace(IncrStrain)*One + psi*DeviatoricIncrStrain;
            TangentStiffness = K*D+psi*Id;
            return;
        }

        Eigen::VectorXd a = DeviatoricStress - DeviatoricStress0*kappa/(1.0+kappa);
        double LoadCond = -1.0*ComputeInnerProduct(a,DeviatoricIncrStrain);
        if (LoadCond > 0.0) {
            DeviatoricStress0 = DeviatoricStress; kappa = infty;
        }

        if (StrainNorm < DSTol) {
            Stress           = Stress_n + K*ComputeTensorTrace(IncrStrain)*One + psi*DeviatoricIncrStrain;
            TangentStiffness = K*D+psi*Id;
            return;
        }
        
        Eigen::VectorXd X  = ComputeHardening(DeviatoricIncrStrain,DeviatoricStress);
        
        if (rootFlag == 1) {
            Stress           = Stress_n + K*ComputeTensorTrace(IncrStrain)*One + psi*DeviatoricIncrStrain;
            TangentStiffness = K*D+psi*Id;
            return;
        }

        psi   = X(0);
        kappa = X(1);

        Stress           = Stress_n + K*ComputeTensorTrace(IncrStrain)*One + psi*DeviatoricIncrStrain;
        TangentStiffness = K*D+psi*Id;

    }
}

//Second Order Identity Tensor.
Eigen::VectorXd 
PlasticPlaneStrainBA::ComputeIdentityVector() const{
    //Second Order Identity Tensor.
    Eigen::VectorXd One(6);
    One << 1.0, 1.0, 1.0, 0.0, 0.0, 0.0;

    return One;
}

//Constructs deviatoric tensor.
Eigen::MatrixXd 
PlasticPlaneStrainBA::ComputeDeviatoricTensor() const{
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
PlasticPlaneStrainBA::ComputeIdentityTensor() const{
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
PlasticPlaneStrainBA::ComputeIdentityOperator() const{
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
PlasticPlaneStrainBA::ComputeTensorTrace(const Eigen::VectorXd &T){
    double Tr = T(0) + T(1) + T(2); 
    return Tr;
}

//Computes tensor norm.
double 
PlasticPlaneStrainBA::ComputeTensorNorm(const Eigen::VectorXd &T){
    double Norm = sqrt(T(0)*T(0) + T(1)*T(1) + T(2)*T(2) + 2.0*(T(3)*T(3) + T(4)*T(4) + T(5)*T(5)));
    return Norm;
}

//Computes inner product.
double 
PlasticPlaneStrainBA::ComputeInnerProduct(const Eigen::VectorXd &V1, const Eigen::VectorXd &V2){
    double V = V1(0)*V2(0)+V1(1)*V2(1)+V1(2)*V2(2)+2.0*V1(3)*V2(3)+2.0*V1(4)*V2(4)+2.0*V1(5)*V2(5);
    return V;
}

//Compute kappa & psi using bisection method
Eigen::VectorXd
PlasticPlaneStrainBA::ComputeHardeningBisection(const Eigen::VectorXd &de, const Eigen::VectorXd &s) {

    Eigen::VectorXd Y(2);

    int iter = 0, maxiter = 100; double TOL  = 1E-4;
    
    rootFlag = 0;

    double x2L, x2R, x2C, x1, x2, fL, fR, fC, Hnp1;

    x2L = 1.0E12; Hnp1 = h*pow(x2L,m);
    x1  = 2.0*G/(1.0+3.0*G*((1.0-beta)/Hn+beta/Hnp1));
    Eigen::VectorXd a = s + x1*de + x2L*(s + x1*de - DeviatoricStress0);
    fL = ComputeTensorNorm(a)-R;

    x2R = 0.0;
    x1 = 2.0*G/(1.0+3.0*G*((1.0-beta)/Hn));
    a = s + x1*de + x2R*(s + x1*de - DeviatoricStress0);
    fR = ComputeTensorNorm(a)-R;

    if (fL*fR > 0.0) {
        rootFlag = 1;
    }
    else {
        x2C = (x2L+x2R)/2.0;
        Hnp1 = h*pow(x2C,m) ;
        x1   = 2.0*G/(1.0+3.0*G*((1.0-beta)/Hn+beta/Hnp1));
        a  = s + x1*de + x2C*(s + x1*de - DeviatoricStress0);
        fC  = ComputeTensorNorm(a)-R;
    
        while (fabs(fC/R) > TOL && iter < maxiter) {
            if (fC*fL < 0.0) {
                x2R = x2C;
            }
            else if (fC*fR < 0.0) {
                x2L = x2C;
            }
            x2C = (x2L+x2R)/2.0;
            Hnp1 = h*pow(x2C,m);
            x1   = 2.0*G/(1.0+3.0*G*((1.0-beta)/Hn+beta/Hnp1));
            a  = s + x1*de + x2C*(s + x1*de - DeviatoricStress0);
            fC  = ComputeTensorNorm(a)-R;
            iter ++;
        }
    }
    
    x2 = x2C;

    Y  << x1, x2;

    return Y;

}

//Compute kappa and psi using newton method
Eigen::VectorXd
PlasticPlaneStrainBA::ComputeHardening(const Eigen::VectorXd &de, const Eigen::VectorXd &s) {

    Eigen::VectorXd Y(2); 

    int iter = 0, maxiter = 100; double TOL  = 1E-6;
    
    rootFlag = 0;

    double x1, x2, Hnp1, f1, f2, err, J11, J12, J21, J22, detJ;

    x1 = 2.0*G; x2 = kappa;

    Hnp1 = h*pow(x2, m);

    Eigen::VectorXd a = s + x1*de + x2*(s + x1*de - DeviatoricStress0);
    
    Eigen::VectorXd dadx1, dadx2;

    f1   = x1 + 3.0*G*x1*((1.0-beta)/Hn+beta/Hnp1) - 2.0*G;
    f2   = ComputeInnerProduct(a,a)-R*R;
    err  = sqrt(f1*f1/G/G/4.0+f2*f2/R/R/R/R);

    while (err > TOL && iter < maxiter && x2 > 0.0) {
    
        dadx1 = (1.0 + x2)*de ;
        dadx2 = s + x1*de - DeviatoricStress0;
    
        J11  =  1.0+3.0*G*((1.0-beta)/Hn+beta/Hnp1);
        J12  = -3.0*G*x1*beta*m/h/pow(x2,m+1.0);
        J21  =  2.0*ComputeInnerProduct(a,dadx1);
        J22  =  2.0*ComputeInnerProduct(a,dadx2);
    
        detJ = J11*J22-J12*J21;
    
        x1 += (-J22*f1+J12*f2)/detJ;
        x2 += ( J21*f1-J11*f2)/detJ;

        Hnp1 = h*pow(x2,m);
    
        a    = s + x1*de + x2*(s + x1*de - DeviatoricStress0);

        f1   = x1 + 3.0*G*x1*((1.0-beta)/Hn+beta/Hnp1) - 2.0*G;
        f2   = ComputeInnerProduct(a,a)-R*R;
    
        err  = sqrt(f1*f1/G/G/4.0+f2*f2/R/R/R/R);
    
        iter ++;
    }

    if (x2 <= 0) rootFlag = 1;

    if (iter == maxiter) std::cout << "BA newton reached max iteration.\n";

    Y  << x1, x2;

    return Y;
}