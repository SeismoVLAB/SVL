#include "CompositeBathe.hpp"
#include "Definitions.hpp"
#include "Profiler.hpp"

//Overload constructor.
CompositeBathe::CompositeBathe(std::shared_ptr<Mesh> &mesh, double TimeStep, double mtol, double ktol, double ftol) : 
Integrator(mesh), dt(TimeStep){
    //Allocate memory for total state vector. 
    U.resize(numberOfTotalDofs); U.fill(0.0);
    V.resize(numberOfTotalDofs); V.fill(0.0);
    A.resize(numberOfTotalDofs); A.fill(0.0);

    //Allocate memory for total model matrices. 
    M.resize(numberOfTotalDofs, numberOfTotalDofs); 
    C.resize(numberOfTotalDofs, numberOfTotalDofs);
    K.resize(numberOfTotalDofs, numberOfTotalDofs);

    //Creates the dynamic assembler for this integrator.
    theAssembler = std::make_unique<Assembler>();
    theAssembler->SetMassTolerance(mtol);
    theAssembler->SetForceTolerance(ftol);
    theAssembler->SetStiffnessTolerance(ktol);

    //Assemble the external force vector from previous analysis.
    Fbar = theAssembler->ComputeProgressiveForceVector(mesh);
}

//Default destructor.
CompositeBathe::~CompositeBathe(){
    //Does nothing.
}

//Initialize model matrices.
void 
CompositeBathe::Initialize(std::shared_ptr<Mesh> &mesh){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Sets up Initial Condition from previous simulation
    std::map<unsigned int, std::shared_ptr<Node> > Nodes = mesh->GetNodes();  
    for(auto it : Nodes){
        auto &Tag = it.first;

        //Gets the associated nodal degree-of-freedom.
        std::vector<int> TotalDofs = Nodes[Tag]->GetTotalDegreeOfFreedom();

        //Creates the nodal/incremental vector state.
        Eigen::VectorXd Uij = Nodes[Tag]->GetDisplacements();
        Eigen::VectorXd Vij = Nodes[Tag]->GetVelocities();
        Eigen::VectorXd Aij = Nodes[Tag]->GetAccelerations();
        
        for(unsigned int j = 0; j < TotalDofs.size(); j++){
            U(TotalDofs[j]) = Uij(j);
            V(TotalDofs[j]) = Vij(j);
            A(TotalDofs[j]) = Aij(j);
        }
    }

    //Computes the mass matrix of the model.
    M = theAssembler->ComputeMassMatrix(mesh);

    //Computes the total damping matrix.
    C = theAssembler->ComputeDampingMatrix(mesh);
}

//Sets the load combination to be used.
void 
CompositeBathe::SetLoadCombination(std::shared_ptr<LoadCombo> &combo){
    theAssembler->SetLoadCombination(combo);
}

//Sets the algorithm  to be used.
void 
CompositeBathe::SetAlgorithm(std::shared_ptr<Algorithm> &algorithm){
    theAlgorithm = algorithm;
}

//Gets the displacement vector.
const Eigen::VectorXd& 
CompositeBathe::GetDisplacements(){
    return U;
}    

//Gets the velocity vector.
const Eigen::VectorXd& 
CompositeBathe::GetVelocities(){
    return V;
}

//Gets the acceleration vector.
const Eigen::VectorXd& 
CompositeBathe::GetAccelerations(){
    return A;
}

//Gets the perfectly-matched layer history vector.
const Eigen::VectorXd& 
CompositeBathe::GetPMLHistoryVector(){
    //Empty PML history vector (not used).
    return Ubar;
}

//Computes a new time step.
bool 
CompositeBathe::ComputeNewStep(std::shared_ptr<Mesh> &mesh, unsigned int k){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Gets the shared_ptr information from the weak_ptr object.
    std::shared_ptr<Algorithm> p = theAlgorithm.lock();

    //If linear system solver fails.
    bool stop;

    //Incremental displacement vector.
    Eigen::VectorXd dU(numberOfTotalDofs);

    //Stores displacement/velocities states.
    Eigen::VectorXd Up = U; 
    Eigen::VectorXd Vp = V;

    //FIRST SUB-STEP SOLUTION: Trapezoidal rule.
    flag = true;

    //Computes a new displacement increment.
    stop = p->ComputeNewIncrement(mesh, k);

    //Checks the solution has issues.
    if(stop) return stop;

    //Obtains the displacement increment from algorithm.
    dU = Total2FreeMatrix*(p->GetDisplacementIncrement());

    //Update displacement states.
    Um = U + dU + SupportMotion;

    //Update velocity states.
    Vm = 4.0/dt*dU - V;

    //Reverse Node and Element Internal variables to last commited state
    p->ReverseStatesIncrements(mesh);

    //SECOND SUB-STEP SOLUTION: 3-point backward Euler rule.
    flag = false; 

    //Computes a new displacement increment.
    stop = p->ComputeNewIncrement(mesh, k);

    //Checks the solution has issues.
    if(stop) return stop;

    //Obtains the displacement increment from algorithm.
    dU = Total2FreeMatrix*(p->GetDisplacementIncrement());

    //Update displacement states.
    U += (dU + SupportMotion);

    //Update velocity states.
    V = 1.0/dt*Up - 4.0/dt*Um + 3.0/dt*U;

    //Update velocity states.
    A = 1.0/dt*Vp - 4.0/dt*Vm + 3.0/dt*V;

    //Return the integrator status.
    return false;
}

//Gets the reaction force ins this step.
Eigen::VectorXd 
CompositeBathe::ComputeReactionForce(std::shared_ptr<Mesh> &mesh, unsigned int k){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Assemble the total external force vector.
    Eigen::VectorXd Fint = theAssembler->ComputeDynamicInternalForceVector(mesh);

    //Assemble the total external force vector.
    Eigen::VectorXd Fext = theAssembler->ComputeExternalForceVector(mesh, k);

    //Computes the reaction forces.
    Eigen::VectorXd Reaction = Fint - Fext - Fbar;

    return Reaction;
}

//Gets the external force vector for next phase analysis.
Eigen::VectorXd 
CompositeBathe::ComputeProgressiveForce(std::shared_ptr<Mesh> &mesh, unsigned int k){
//Starts profiling this function.
    PROFILE_FUNCTION();

    //Assemble the total internal force vector.
    Eigen::VectorXd Fext = theAssembler->ComputeExternalForceVector(mesh, k);

    //Update the stage force vector.
    Eigen::VectorXd Force = Fext + Fbar;

    return Force;
}

//Gets the incremental nodal support motion vector.
void
CompositeBathe::ComputeSupportMotionVector(std::shared_ptr<Mesh>& mesh, Eigen::VectorXd &Feff, double UNUSED(factor), unsigned int k){
    //TODO:Include some if statement that performs the addition only if there is support motion applied.
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Assembles the incremental support motion displacements.
    SupportMotion = theAssembler->ComputeSupportMotionIncrement(mesh, k-1);

    if(flag){
        //half-time step for trapezoidal rule.
        SupportMotion *= 0.5;
    }

    //Computes the required forces to impose these displacements.
    Eigen::VectorXd Lg = Total2FreeMatrix.transpose()*(K*SupportMotion);

    //Add the contribution to the current effective force vector.
    Feff -= Lg;
}

//Gets the effective force associated to this integrator.
void
CompositeBathe::ComputeEffectiveForce(std::shared_ptr<Mesh> &mesh, Eigen::VectorXd &Feff, double UNUSED(factor), unsigned int k){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Gets the shared_ptr information from the weak_ptr object.
    std::shared_ptr<Algorithm> p = theAlgorithm.lock();

    //Obtains the displacement increment from algorithm.
    Eigen::VectorXd dU = Total2FreeMatrix*(p->GetDisplacementIncrement());

    //Assemble the total internal force vector.
    Eigen::VectorXd Fint = theAssembler->ComputeInternalForceVector(mesh);

    //Computes the effective force vector.
    Eigen::VectorXd Fext = theAssembler->ComputeExternalForceVector(mesh, k);

    //Computes the effective force depending on the step.
    if(flag){
        //Force using trapezoidal rule.
        Eigen::VectorXd Faux = theAssembler->ComputeExternalForceVector(mesh, k-1);
        Fext = 1.0/2.0*(Fext + Faux) + Fbar - Fint - M*(16.0/dt/dt*dU - 8.0/dt*V - A) - C*(4.0/dt*dU - V);
    }
    else{
        //Force using 3-point backward Euler rule.
        Fext = Fext + Fbar - Fint - M*(9.0/dt/dt*(U + dU) - 12.0/dt/dt*Um + 3.0/dt/dt*U - 4.0/dt*Vm + 1.0/dt*V) - C*(3.0/dt*(U + dU) - 4.0/dt*Um + 1.0/dt*U);
    }

    //Impose boundary conditions on effective force vector.
    Feff = Total2FreeMatrix.transpose()*Fext;
}

//Gets the effective stiffness associated to this integrator.
void
CompositeBathe::ComputeEffectiveStiffness(std::shared_ptr<Mesh> &mesh, Eigen::SparseMatrix<double> &Keff){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Computes the stiffness matrix of the model.
    K = theAssembler->ComputeStiffnessMatrix(mesh);

    //Computes the effective stiffness matrix depending on the step.
    if(flag){
        //Stiffness using trapezoidal rule.
        K = 16.0/dt/dt*M + 4.0/dt*C + K;
    }
    else{
        //Stiffness using 3-point backward Euler rule.
        K = 9.0/dt/dt*M + 3.0/dt*C + K;
    }

    //Impose boundary conditions on effective stiffness matrix.
    Keff = Eigen::SparseMatrix<double>(Total2FreeMatrix.transpose())*K*Total2FreeMatrix;
}
