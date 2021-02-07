#include "Linear.hpp"
#include "Definitions.hpp"
#include "Profiler.hpp"

//Overload constructor.
Linear::Linear(std::unique_ptr<LinearSystem> &solver, std::shared_ptr<Mesh> &mesh, double tol, double ldFactor, unsigned int nIters) : 
Algorithm(mesh), Tolerance(tol), LoadFactor(ldFactor), nMaxIterations(nIters){
    //Move solver pointer to algorithm.
    theSolver = std::move(solver);

    //Allocate memory for the incremental displacement vector.
    dU.resize(numberOfFreeDofs);
}

//Default destructor.
Linear::~Linear(){
    //Does nothing.
}

//Compute a new solution.
bool 
Linear::ComputeNewIncrement(std::shared_ptr<Mesh> &mesh, unsigned int i){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Initialize the displacement increment.
    dU.fill(0.0);

    //Initialize Stiffness Matrix and Force Vector.
    Eigen::VectorXd Feff(numberOfFreeDofs);
    Eigen::SparseMatrix<double> Keff(numberOfFreeDofs, numberOfFreeDofs);

    //Assemble the total force vector.
    theIntegrator->ComputeEffectiveForce(mesh, Feff, LoadFactor, i);

    //Assemble the total stiffness matrix.
    theIntegrator->ComputeEffectiveStiffness(mesh, Keff);

    //Incorporates the effects of support motion.
    theIntegrator->ComputeSupportMotionVector(mesh, Feff, LoadFactor, i);

    //Solve the linear system.
    bool stop = theSolver->SolveSystem(Keff, Feff);

    //Checks the solution has issues.
    if(stop) return stop;

    //Gets the total incremental displacement vector.
    dU = theSolver->GetSolution();

    //Update the incremental state variables in the mesh.
    UpdateStatesIncrements(mesh, dU);

    //Return the algorithm status.
    return false;
}

//Gets the displacement increment.
Eigen::VectorXd 
Linear::GetDisplacementIncrement(){
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

    return dU;
}

//Set the load factor.
void 
Linear::SetLoadFactor(double factor){
    LoadFactor = factor;
}

//Sets the integrator for this algorithm.
void 
Linear::SetIntegrator(std::shared_ptr<Integrator> &integrator){
    theIntegrator = integrator;
}
