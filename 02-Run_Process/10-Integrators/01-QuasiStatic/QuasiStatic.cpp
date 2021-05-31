#include "QuasiStatic.hpp"
#include "Definitions.hpp"
#include "Profiler.hpp"

//Overload constructor.
QuasiStatic::QuasiStatic(std::shared_ptr<Mesh> &mesh, double mtol, double ktol, double ftol) : 
Integrator(mesh), TimeStep(0.0){
    //Allocate memory for total state vector. 
    U.resize(numberOfTotalDofs); U.fill(0.0);

    //Allocate memory for total model matrices. 
    K.resize(numberOfTotalDofs, numberOfTotalDofs);

    //Creates the static assembler for this integrator.
    theAssembler = std::make_unique<Assembler>();
    theAssembler->SetMassTolerance(mtol);
    theAssembler->SetForceTolerance(ftol);
    theAssembler->SetStiffnessTolerance(ktol);
}

//Default destructor.
QuasiStatic::~QuasiStatic(){
    //Does nothing.
}

//Initialize model matrices.
void 
QuasiStatic::Initialize(std::shared_ptr<Mesh> &mesh){
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

        for(unsigned int j = 0; j < TotalDofs.size(); j++)
            U(TotalDofs[j]) = Uij(j);
    } 

    //Assemble the effective stiffness matrix.
    K = theAssembler->ComputeStiffnessMatrix(mesh);
}

//Sets the load combination to be used.
void 
QuasiStatic::SetLoadCombination(std::shared_ptr<LoadCombo> &combo){
    theAssembler->SetLoadCombination(combo);
}

//Sets the algorithm  to be used.
void 
QuasiStatic::SetAlgorithm(std::shared_ptr<Algorithm> &algorithm){
    theAlgorithm = algorithm;
}

//Gets the displacement vector.
Eigen::VectorXd& 
QuasiStatic::GetDisplacements(){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    return U; 
}    

//Gets the velocity vector.
Eigen::VectorXd&
QuasiStatic::GetVelocities(){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Velocity vector is no used.
    return U;
}

//Gets the acceleration vector.
Eigen::VectorXd& 
QuasiStatic::GetAccelerations(){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Acceleration vector is no used.
    return U;
}

//Gets the perfectly-matched layer history vector.
Eigen::VectorXd& 
QuasiStatic::GetPMLHistoryVector(){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Empty PML history vector (not used).
    return Ubar;
}

//Computes a new time step.
bool 
QuasiStatic::ComputeNewStep(std::shared_ptr<Mesh> &mesh, unsigned int k){
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

    //Return the integrator status.
    return false;
}

//Gets the reaction force ins this step.
Eigen::VectorXd 
QuasiStatic::ComputeReactionForce(std::shared_ptr<Mesh> &mesh, unsigned int k){
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
QuasiStatic::ComputeSupportMotionVector(std::shared_ptr<Mesh> &mesh, Eigen::VectorXd &Feff, double factor, unsigned int UNUSED(k)){
    //TODO:Include some if statement that performs the addition only if there is support motion applied.
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Assembles the incremental support motion displacements.
    SupportMotion = factor*(theAssembler->ComputeSupportMotionIncrement(mesh, 0));

    //Computes the required forces to impose these displacements.
    Eigen::VectorXd Lg = Total2FreeMatrix.transpose()*(K*SupportMotion);

    //Add the contribution to the current effective force vector.
    Feff -= Lg;
}

//Gets the effective force associated to this integrator.
void
QuasiStatic::ComputeEffectiveForce(std::shared_ptr<Mesh> &mesh, Eigen::VectorXd &Feff, double factor, unsigned int k){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Assemble the total external force vector.
    Eigen::VectorXd Fint = theAssembler->ComputeInternalForceVector(mesh);

    //Assemble the total internal force vector.
    Eigen::VectorXd Fext = theAssembler->ComputeExternalForceVector(mesh, k);

    //Computes the effective-right-hand-side force vector.
    Fext = (1.00 + k)*factor*Fext - Fint;

    //Impose boundary conditions on effective force vector.
    Feff = Total2FreeMatrix.transpose()*Fext;
}

//Gets the effective stiffness associated to this integrator.
void
QuasiStatic::ComputeEffectiveStiffness(std::shared_ptr<Mesh> &mesh, Eigen::SparseMatrix<double> &Keff){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Assemble the effective stiffness matrix.
    K = theAssembler->ComputeStiffnessMatrix(mesh);

    //Impose boundary conditions on effective stiffness matrix.
    Keff = Eigen::SparseMatrix<double>(Total2FreeMatrix.transpose())*K*Total2FreeMatrix;
}