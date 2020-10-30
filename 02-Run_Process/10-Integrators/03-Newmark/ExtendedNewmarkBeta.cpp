#include "ExtendedNewmarkBeta.hpp"
#include "Definitions.hpp"
#include "Profiler.hpp"

//Overload constructor.
ExtendedNewmarkBeta::ExtendedNewmarkBeta(std::shared_ptr<Mesh> &mesh, double TimeStep, double mtol, double ktol, double ftol) : 
Integrator(mesh), dt(TimeStep){
    //Allocate memory for total state vector. 
    U.resize(numberOfTotalDofs); U.fill(0.0);
    V.resize(numberOfTotalDofs); V.fill(0.0);
    A.resize(numberOfTotalDofs); A.fill(0.0);
    Ubar.resize(numberOfTotalDofs); Ubar.fill(0.0);

    //Allocate memory for total model matrices. 
    M.resize(numberOfTotalDofs, numberOfTotalDofs); 
    C.resize(numberOfTotalDofs, numberOfTotalDofs);
    K.resize(numberOfTotalDofs, numberOfTotalDofs);
    G.resize(numberOfTotalDofs, numberOfTotalDofs);

    //Creates the dynamic assembler for this integrator.
    theAssembler = std::make_unique<Assembler>();
    theAssembler->SetMassTolerance(mtol);
    theAssembler->SetForceTolerance(ftol);
    theAssembler->SetStiffnessTolerance(ktol);

    Initialize(mesh);
}

//Default destructor.
ExtendedNewmarkBeta::~ExtendedNewmarkBeta(){
    //Does nothing.
}

//Initialize model matrices.
void 
ExtendedNewmarkBeta::Initialize(std::shared_ptr<Mesh> &mesh){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Computes the mass matrix of the model.
    M = theAssembler->ComputeMassMatrix(mesh);

    //Computes the total damping matrix.
    C = theAssembler->ComputeDampingMatrix(mesh);

    // Compute the PML history matrix
    G = theAssembler->ComputePMLHistoryMatrix(mesh);
}

//Sets the load combination to be used.
void 
ExtendedNewmarkBeta::SetLoadCombination(std::shared_ptr<LoadCombo> &combo){
    theAssembler->SetLoadCombination(combo);
}

//Sets the algorithm  to be used.
void 
ExtendedNewmarkBeta::SetAlgorithm(std::shared_ptr<Algorithm> &algorithm){
    theAlgorithm = algorithm;
}

//Gets the PML history vector.
Eigen::VectorXd& 
ExtendedNewmarkBeta::GetPMLHistoryVector(){
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

    return Ubar;
} 

//Gets the displacement vector.
Eigen::VectorXd& 
ExtendedNewmarkBeta::GetDisplacements(){
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

    return U;
}    

//Gets the velocity vector.
Eigen::VectorXd& 
ExtendedNewmarkBeta::GetVelocities(){
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

    return V;
}

//Gets the acceleration vector.
Eigen::VectorXd& 
ExtendedNewmarkBeta::GetAccelerations(){
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

    return A;
}

//Computes a new time step.
bool 
ExtendedNewmarkBeta::ComputeNewStep(std::shared_ptr<Mesh> &mesh, unsigned int k){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Gets the shared_ptr information from the weak_ptr object.
    std::shared_ptr<Algorithm> p = theAlgorithm.lock();

    //Computes a new displacement increment.
    bool stop = p->ComputeNewIncrement(mesh, k);

    //Checks the solution has issues.
    if(stop) return stop;

    //Obtains the displacement increment from algorithm.
    Eigen::VectorXd dU = Total2FreeMatrix*(p->GetDisplacementIncrement());

    //Update displacement states.
    U += (dU + SupportMotion);

    //Update acceleration states.
    A = 4.0/dt/dt*dU - 4.0/dt*V - A;

    //Update velocity states.
    V = 2.0/dt*dU - V;

    //Return the integrator status.
    return false;
}

//Gets the reaction force ins this step.
Eigen::VectorXd 
ExtendedNewmarkBeta::ComputeReactionForce(std::shared_ptr<Mesh> &mesh, unsigned int k){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Assemble the total external force vector.
    Eigen::VectorXd Fint = theAssembler->ComputeDynamicInternalForceVector(mesh);

    //Assemble the total external force vector.
    Eigen::VectorXd Fext = theAssembler->ComputeExternalForceVector(mesh, k);

    //Computes the reaction forces.
    Eigen::VectorXd Reaction = Fint - Fext;

    return Reaction;
}

//Gets the incremental nodal support motion vector.
Eigen::VectorXd 
ExtendedNewmarkBeta::ComputeSupportMotionVector(std::shared_ptr<Mesh> &mesh, double factor, unsigned int k){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Assembles the incremental support motion displacements.
    SupportMotion = theAssembler->ComputeSupportMotionIncrement(mesh, k);

    //Computes the required forces to impose these displacements.
    Eigen::VectorXd Lg = Total2FreeMatrix.transpose()*(K*SupportMotion);

    return Lg;
}

//Gets the effective force associated to this integrator.
void
ExtendedNewmarkBeta::ComputeEffectiveForce(std::shared_ptr<Mesh> &mesh, Eigen::VectorXd &Feff, double factor, unsigned int k){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Gets the shared_ptr information from the weak_ptr object.
    std::shared_ptr<Algorithm> p = theAlgorithm.lock();

    //Obtains the displacement increment from algorithm.
    Eigen::VectorXd dU = Total2FreeMatrix*(p->GetDisplacementIncrement());

    //Assemble the total external force vector.
    Eigen::VectorXd Fint = theAssembler->ComputeInternalForceVector(mesh);

    //Assemble the total internal force vector.
    Eigen::VectorXd Fext = theAssembler->ComputeExternalForceVector(mesh, k);

    //Computes the effective force vector.
    Fext = Fext - Fint + M*(4.0/dt*V + A - 4.0/dt/dt*dU) + C*(V - 2.0/dt*dU);

    //Impose boundary conditions on effective force vector.
    Feff = Total2FreeMatrix.transpose()*Fext;
}

//Gets the effective stiffness assiciated to this integrator.
void
ExtendedNewmarkBeta::ComputeEffectiveStiffness(std::shared_ptr<Mesh> &mesh, Eigen::SparseMatrix<double> &Keff){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Computes the stiffness matrix of the model.
    K = theAssembler->ComputeStiffnessMatrix(mesh);

    //Assemble the effective stiffness matrix.
    K = K + 4.0/dt/dt*M + 2.0/dt*C;

    //Impose boundary conditions on effective stiffness matrix.
    Keff = Eigen::SparseMatrix<double>(Total2FreeMatrix.transpose())*K*Total2FreeMatrix;
}