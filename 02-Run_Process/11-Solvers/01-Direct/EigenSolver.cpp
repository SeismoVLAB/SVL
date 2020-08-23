#include <ctime>
#include <iostream>
#include "EigenSolver.hpp"
#include "Definitions.hpp"
#include "Profiler.hpp"

//Overload constructor.
EigenSolver::EigenSolver(bool flag) : Flag(flag), Factored(false), Initialized(false){
    //Does nothing
}

//Destructor:
EigenSolver::~EigenSolver(){
    //Clears full vectors.    
    x.resize(0);
}

//Solve the linear system.
bool 
EigenSolver::SolveSystem(Eigen::SparseMatrix<double> &A, Eigen::VectorXd &b){
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

    //LINEAR ANALYSIS.    
    if(Flag){
        if(!Initialized){
            Solver.analyzePattern(A);
            Initialized = true;
        }

        //Factorize the matrix just once.
        if(!Factored){
            Solver.factorize(A);
            Factored = true;
        }

        //Solve the linear system.
        x = Solver.solve(b);
    }
    //NONLINEAR ANALYSIS.
    else{
        //Re-order the nonzero elements of A.
        Solver.analyzePattern(A);

        //Factorize the matrix.    
        Solver.factorize(A);

        //Solve the linear system.
        x = Solver.solve(b);
    }

    //Check solution status.
    if(Solver.info() != 0){
        std::cout << "\n \x1B[31mERROR: \x1B[0mEigenSolver::SolveSystem() : Solution for the linear system failed.";
        return true;
    }

    //Return the solver status.
    return false;
}

//Gets the soultion vector.    
Eigen::VectorXd 
EigenSolver::GetSolution(){
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

    return x;
}
