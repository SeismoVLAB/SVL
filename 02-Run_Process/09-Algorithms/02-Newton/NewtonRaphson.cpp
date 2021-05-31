#include <mpi.h>
#include <iostream>
#include "NewtonRaphson.hpp"
#include "Definitions.hpp"
#include "Profiler.hpp"

//Overload constructor.
NewtonRaphson::NewtonRaphson(std::unique_ptr<LinearSystem> &solver, std::shared_ptr<Mesh> &mesh, double tol, unsigned int nIters, unsigned int flag, double ldFactor) : 
Algorithm(mesh, flag), Tolerance(tol), LoadFactor(ldFactor), nMaxIterations(nIters){    
    //Move solver pointer to algorithm.
    theSolver = std::move(solver);

    //Allocate memory for the incremental displacement vector.
    dU.resize(numberOfFreeDofs);
}

//Default destructor.
NewtonRaphson::~NewtonRaphson(){
    //Does nothing.
}

//Compute a new solution.
bool 
NewtonRaphson::ComputeNewIncrement(std::shared_ptr<Mesh> &mesh, unsigned int i){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Auxiliar variables.
    double Residual;
    unsigned int k = 0;

    //Initialize the displacement increment.
    dU.fill(0.0);

    //Initialize Stiffness Matrix and Force Vector.
    Eigen::VectorXd Feff(numberOfFreeDofs); 
    Eigen::SparseMatrix<double> Keff(numberOfFreeDofs, numberOfFreeDofs); 

    //Assemble the total force vector.
    theIntegrator->ComputeEffectiveForce(mesh, Feff, LoadFactor, i);

    //Computes the residual.
    Residual = ComputeConvergence(Feff, 0.0, true);

    do{
        //Assemble the total stiffness matrix.
        theIntegrator->ComputeEffectiveStiffness(mesh, Keff);

        //Incorporates the effects of support motion.
        theIntegrator->ComputeSupportMotionVector(mesh, Feff, LoadFactor, i);

        //Solve the linear system.
        bool stop = theSolver->SolveSystem(Keff, Feff);

        //Checks the solution has issues.
        if(stop) return stop;

        //Gets the total incremental displacement vector.
        dU += theSolver->GetSolution();

        //Update the incremental state variables in the mesh.
        UpdateStatesIncrements(mesh, dU);

        //Assemble the total force vector.
        theIntegrator->ComputeEffectiveForce(mesh, Feff, LoadFactor, i);

        //Updates the residual.
        Residual = ComputeConvergence(Feff, dU.norm(), false);

        k++;
    }
    while((Residual > Tolerance) & (k < nMaxIterations));

    //Alert that maximun number of iterations was reached.
    if((k == nMaxIterations) & (rank == 0))
        std::cout << "\n Newton Raphson algorithm reached maximun number of iterations. The residual is : " << Residual << "\n";

    return false;
}

//Gets the displacement increment.
Eigen::VectorXd 
NewtonRaphson::GetDisplacementIncrement(){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    return dU;
}

//Set the load factor.
void 
NewtonRaphson::SetLoadFactor(double factor){
    LoadFactor = factor;
}

//Sets the integrator for this algorithm.
void 
NewtonRaphson::SetIntegrator(std::shared_ptr<Integrator> &integrator){
    theIntegrator = integrator;
}