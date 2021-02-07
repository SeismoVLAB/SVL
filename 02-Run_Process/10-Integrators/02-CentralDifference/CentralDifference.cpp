#include "CentralDifference.hpp"
#include "Definitions.hpp"
#include "Profiler.hpp"

//Overload constructor.
CentralDifference::CentralDifference(std::shared_ptr<Mesh> &mesh, double TimeStep, double mtol, double ktol, double ftol) : 
Integrator(mesh), dt(TimeStep){
    //Allocate memory for total state vector. 
    U.resize(numberOfTotalDofs); U.fill(0.0);
    V.resize(numberOfTotalDofs); V.fill(0.0);
    A.resize(numberOfTotalDofs); A.fill(0.0);

    //Allocate memory for previous displacement state vector. 
    Up.resize(numberOfTotalDofs); 
    Up = U - dt*V + dt*dt/2.0*A; 

    //Allocate memory for total model matrices.
    M.resize(numberOfTotalDofs, numberOfTotalDofs);     
    C.resize(numberOfTotalDofs, numberOfTotalDofs); 
    K.resize(numberOfTotalDofs, numberOfTotalDofs);

    //Creates the dynamic assembler for this integrator.
    theAssembler = std::make_unique<Assembler>();
    theAssembler->SetMassTolerance(mtol);
    theAssembler->SetForceTolerance(ftol);
    theAssembler->SetStiffnessTolerance(ktol);

    Initialize(mesh);
}

//Default destructor.
CentralDifference::~CentralDifference(){
    //Does nothing.
}

//Initialize model matrices.
void 
CentralDifference::Initialize(std::shared_ptr<Mesh> &mesh){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Computes the mass matrix of the model.
    M = theAssembler->ComputeMassMatrix(mesh);

    //Computes the total damping matrix.
    C = theAssembler->ComputeDampingMatrix(mesh);

    //Model effective stiffness matrix.
    K = 1.0/dt/dt*M + 1.0/2.0/dt*C;
}

//Sets the load combination to be used.
void 
CentralDifference::SetLoadCombination(std::shared_ptr<LoadCombo> &combo){
    theAssembler->SetLoadCombination(combo);
}

//Sets the algorithm  to be used.
void 
CentralDifference::SetAlgorithm(std::shared_ptr<Algorithm> &algorithm){
    theAlgorithm = algorithm;
}

//Gets the displacement vector.
Eigen::VectorXd& 
CentralDifference::GetDisplacements(){
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

    return U;
}    

//Gets the velocity vector.
Eigen::VectorXd& 
CentralDifference::GetVelocities(){
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

    return V;
}

//Gets the acceleration vector.
Eigen::VectorXd& 
CentralDifference::GetAccelerations(){
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

    return A;
}

//Gets the perfectly-matched layer history vector.
Eigen::VectorXd& 
CentralDifference::GetPMLHistoryVector(){
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

    //Empty PML history vector (not used).
    return Ubar;
}

//Computes a new time step.
bool 
CentralDifference::ComputeNewStep(std::shared_ptr<Mesh> &mesh, unsigned int k){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Gets the shared_ptr information from the weak_ptr object.
    std::shared_ptr<Algorithm> p = theAlgorithm.lock();

    //Computes a new displacement increment.
    bool stop = p->ComputeNewIncrement(mesh, k);

    //Checks the solution has issues.
    if(stop) return stop;

    //Obtains the displacement from algorithm.
    Eigen::VectorXd dU = Total2FreeMatrix*(p->GetDisplacementIncrement()) + SupportMotion;

    //Update velocity states.
    V = 1.0/2.0/dt*(U + dU - Up); 

    //Update acceleration states.
    A = 1.0/dt/dt*(dU - U + Up);

    //Update displacements states.
    Up  = U;
    U  += dU; 

    //Return the integrator status.
    return false;
}

//Gets the reaction force ins this step.
Eigen::VectorXd 
CentralDifference::ComputeReactionForce(std::shared_ptr<Mesh> &mesh, unsigned int k){
    //TODO:Include some if statement that performs the addition only if there is support motion applied.
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
void
CentralDifference::ComputeSupportMotionVector(std::shared_ptr<Mesh> &mesh, Eigen::VectorXd &Feff, double UNUSED(factor), unsigned int k){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Assembles the incremental support motion displacements.
    SupportMotion = theAssembler->ComputeSupportMotionIncrement(mesh, k);

    //Computes the required forces to impose these displacements.
    Eigen::VectorXd Lg = Total2FreeMatrix.transpose()*(K*SupportMotion);

    //Add the contribution to the current effective force vector.
    Feff -= Lg;
}

//Gets the effective force associated to this integrator.
void
CentralDifference::ComputeEffectiveForce(std::shared_ptr<Mesh> &mesh, Eigen::VectorXd &Feff, double UNUSED(factor), unsigned int k){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Assemble the total external force vector.
    Eigen::VectorXd Fint = theAssembler->ComputeInternalForceVector(mesh);

    //Assemble the total internal force vector.
    Eigen::VectorXd Fext = theAssembler->ComputeExternalForceVector(mesh, k);

    //Computes the effective force vector.
    Fext = Fext - Fint + (1.0/dt/dt*M - 1.0/2.0/dt*C)*(U - Up);

    //Impose boundary conditions on effective force vector.
    Feff = Total2FreeMatrix.transpose()*Fext;
}

//Gets the effective stiffness assiciated to this integrator.
void
CentralDifference::ComputeEffectiveStiffness(std::shared_ptr<Mesh>& UNUSED(mesh), Eigen::SparseMatrix<double> &Keff){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Transform the effective stiffness matrix into free degree-of-freedom.
    Keff = Eigen::SparseMatrix<double>(Total2FreeMatrix.transpose())*K*Total2FreeMatrix;
}
