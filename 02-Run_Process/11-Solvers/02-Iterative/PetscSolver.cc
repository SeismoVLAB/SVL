#include <ctime>
#include <sstream> 
#include <iostream>
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>

#include "PetscSolver.hpp"
#include "Definitions.hpp"
#include "Profiler.hpp"

//The KSP Solver Type.
KSPType KSPSolverName;

//Overload constructor.
PetscSolver::PetscSolver(unsigned int dnz, unsigned int onz, double tol, unsigned int kspnum) : 
d_nz(dnz), o_nz(onz), Tolerance(tol) {
    //Assigns the iterative solver.    
    switch (kspnum){
        case 0: {    //The Preconditioned Conjugate Gradient (PCG) iterative method 
            KSPSolverName = KSPCG; 
            break; 
        }
        case 1: {     //The BiCGStab (Stabilized version of BiConjugate Gradient) method.
            KSPSolverName = KSPBCGS;
            break;
        }
        case 2: {    //The CGS (Conjugate Gradient Squared) method.
            KSPSolverName = KSPCGS; 
            break;
        }
        case 3:    {    //The Biconjugate gradient iterative method.
            KSPSolverName = KSPBICG;
            break;
        }
        default: {     //Does nothing.
            KSPSolverName = KSPBICG;
        }
    }

    //Sets vector sizes.
    x.resize(numberOfFreeDofs);
    x.fill(0.0);
}

//Destructor:
PetscSolver::~PetscSolver(){
    //Clears full vectors.
    x.resize(0);
}

//Solve the linear system.
bool 
PetscSolver::SolveSystem(Eigen::SparseMatrix<double> &K, Eigen::VectorXd &b){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Creates and Assembles the right-hand side vector.
    Vec B;
    VecCreate(PETSC_COMM_WORLD, &B);
    VecSetSizes(B, PETSC_DECIDE, numberOfFreeDofs);
    VecSetFromOptions(B); 
    VecSet(B, 0.0); 

    for (unsigned int k = 0; k < numberOfFreeDofs; k++) {
        double value = b(k);
        if(value != 0.0){
            PetscInt Jc = k;
            PetscScalar val = value;
            VecSetValues(B, 1, &Jc, &val, ADD_VALUES);
        }
    }

    VecAssemblyBegin(B);
    VecAssemblyEnd(B);

    //Creates and Assembles the left-hand side matrix.
    Mat A;
    MatCreate(PETSC_COMM_WORLD, &A); 
    MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, numberOfFreeDofs, numberOfFreeDofs);     
    MatSetFromOptions(A);
    MatMPIAIJSetPreallocation(A, d_nz, PETSC_NULL, o_nz, PETSC_NULL);

    //Sequential Matrix for single core execution
    if(size == 1)
        MatSeqAIJSetPreallocation(A, d_nz, PETSC_NULL);

    for(unsigned int k = 0; k < K.outerSize(); ++k){
        for(Eigen::SparseMatrix<double>::InnerIterator it(K,k); it; ++it) {
            unsigned int ir = it.row();
            unsigned int jc = it.col();
            double  value = it.value();
            MatSetValue(A, ir, jc, value, ADD_VALUES);
        }
    }

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);     

    //Creates the Solution Vector.
    Vec X;
    VecCreate(PETSC_COMM_WORLD, &X);
    VecSetSizes(X, PETSC_DECIDE, numberOfFreeDofs); 
    VecSetFromOptions(X);
    VecSet(X, 0.0);

    //Solves the Linear System.
    KSP ksp;
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, A, A);
    KSPSetType(ksp, KSPSolverName); 
    KSPSetTolerances(ksp, Tolerance, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    KSPSetFromOptions(ksp);
    KSPSolve(ksp, B, X);

    //Distribute solution to each processor.
    Vec V_SEQ;
    VecCreate(PETSC_COMM_WORLD, &V_SEQ);
    VecSetSizes(V_SEQ, PETSC_DECIDE, numberOfFreeDofs);
    VecSetFromOptions(V_SEQ);

    VecScatter ctx;
    VecScatterCreateToAll(X, &ctx, &V_SEQ);
    VecScatterBegin(ctx, X, V_SEQ, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctx, X, V_SEQ, INSERT_VALUES, SCATTER_FORWARD);

    //Transform the solution into Eigen's format.
    for(unsigned int k = 0; k < numberOfFreeDofs; k++) {
        PetscInt Jc = k;
        PetscScalar value;
        VecGetValues(V_SEQ, 1, &Jc, &value);
        double newValue = value;
        x(k) = newValue;
    }

    //PetscInt its;
    //KSPGetIterationNumber(ksp, &its);
    //PetscPrintf(PETSC_COMM_WORLD,"Iterations %D\n", its);

    //Cleans the PETSC's variables.
    MatDestroy(&A);
    VecDestroy(&B);
    VecDestroy(&X);

    KSPDestroy(&ksp);
    VecDestroy(&V_SEQ);
    VecScatterDestroy(&ctx);

    //Return the solver status.
    return false;
}

//Gets the soultion vector.    
Eigen::VectorXd 
PetscSolver::GetSolution(){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    return x;
}
