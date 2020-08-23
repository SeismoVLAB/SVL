#include <mpi.h>
#include <cstdlib>
#include <iostream>
#include "MumpsSolver.hpp"
#include "Definitions.hpp"
#include "Profiler.hpp"

//Overload constructor.
MumpsSolver::MumpsSolver(unsigned int option, bool flag) : Option(option), Flag(flag), Factored(false), Initialized(false){ 
    //Initialize a MUMPS instance.
    id.job = JOB_INIT; 
    id.par = 1;

    switch(Option){
        //Symmetric-Positive Definite.
        case  0: id.sym = 1; break;
        //Symmetric matrix pattern.
        case  1: id.sym = 2; break;
        //UnSymmetric matrix pattern.
        default: id.sym = 0; break;
    }

    id.comm_fortran = USE_COMM_WORLD;

    dmumps_c(&id);

    //Display Outputs Information.
    id.ICNTL(1) = 0; 
    id.ICNTL(2) = 0; 
    id.ICNTL(3) = 0; 
    id.ICNTL(4) = 0; 

    id.ICNTL(5)  = 0; 
    id.ICNTL(6)  = 0;
    id.ICNTL(18) = 3;

    //Dynamic Memory Allocation.
    sol = new double[numberOfFreeDofs];
}

//Destructor:
MumpsSolver::~MumpsSolver(){
    //Terminates MUMPS Instance.
    id.job = JOB_END; 
    dmumps_c(&id);

    //Delete Dynamic Pointers.
    delete[] sol;
}

//Solve the linear system.
bool 
MumpsSolver::SolveSystem(Eigen::SparseMatrix<double> &A, Eigen::VectorXd &b){
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

    //Transform Eigen vector to pointer.
    double *rhs = new double[numberOfFreeDofs];

    for(unsigned int k = 0; k < numberOfFreeDofs; k++)
        rhs[k] = b(k);

    //Assembly Force Vector.
    MPI_Reduce(rhs, sol, numberOfFreeDofs, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    //Vector Coefficients.
    id.rhs = sol;
     
    //Number of non-zeros (Full-Matrix or Upper-Triangular).
    if(Option == 2){
        nz = A.nonZeros();
    }
    else{
        nz = 0;
        for(int i = 0; i < A.outerSize(); ++i){
          for(Eigen::SparseMatrix<double>::InnerIterator it(A,i); it; ++it){
            if(it.row() >= i)
                nz++;
          }
        }
    }

    //Pointer to Upper-Triangular Matrix indeces.
    int *irn = new int[nz];
    int *jcn = new int[nz];
    double *ptrs = new double[nz];

    if(Option == 2){
        //UnSymmetric matrix pattern: Full matrix Values
        int j = 0;
        for(int i = 0; i < A.outerSize(); ++i){
            for(Eigen::SparseMatrix<double>::InnerIterator it(A,i); it; ++it){
                irn[j]  = it.row() + 1;   
                jcn[j]  = it.col() + 1;
                ptrs[j] = it.value();
                j++;
            }
        }
    }
    else{
        //Symmetric matrix pattern: Upper-Triangular Values
        int j = 0;
        for(int i = 0; i < A.outerSize(); ++i){
            for(Eigen::SparseMatrix<double>::InnerIterator it(A,i); it; ++it){
                if(it.row() >= i){
                    irn[j]  = it.row() + 1;   
                    jcn[j]  = it.col() + 1; 
                    ptrs[j] = it.value();
                    j++;
                }
            }
        }
    }

    //Linear System Size.
    id.n       = numberOfFreeDofs; 
    id.nz_loc  = nz; 
    id.irn_loc = irn; 
    id.jcn_loc = jcn;
    id.a_loc   = ptrs; 

    //LINEAR ANALYSIS
    if(Flag){ 
        if(!Initialized){
            //Call the MUMPS package to analyze the system
            id.job = 1;
            dmumps_c(&id);
            Initialized = true;
        }

        if(!Factored){
            //Call the MUMPS package to factor and solve the system
            id.job = 2;
            dmumps_c(&id);
            Factored = true;
        }

        //Call the MUMPS package to only solve the system
        id.job = 3;
        dmumps_c(&id);
    }
    //For NONLINEAR ANALYSIS
    else{
        //Call the MUMPS package to analyze, factor and solve the system
        id.job = 6;
        dmumps_c(&id);
    }

    //Check solution status.
    int info = id.infog[0];
    if (info != 0){
        if (rank == 0){    
            std::cout << "\n \x1B[31mError: \x1B[0mMumpsSolver::SolveSystem() : ";
            switch(info) {
                case  -5: std::cout << "OUT OF MEMORY allocation error.\n" ; break;
                case  -6: std::cout << "MATRIX IS SINGULAR in Structure.\n"; break;
                case  -7: std::cout << "OUT OF MEMORY allocation error.\n" ; break;
                case  -8: std::cout << "SMALL WORK ARRAY use -ICNTL14 option, the default is -ICNTL 20 make 20 larger.\n"; break;
                case  -9: std::cout << "SMALL WORK ARRAY use -ICNTL14 option, the default is -ICNTL 20 make 20 larger.\n"; break;
                case -10: std::cout << "MATRIX IS SINGULAR Numerically.\n"; break;
                default : std::cout << "Solution for the linear system failed.\n"; break;
            }
        }
        return true;
    }

    //Copy Solution to All Processors.
    MPI_Bcast(sol, numberOfFreeDofs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //Delete Dynamic Pointers.
    delete[] irn;
    delete[] jcn;
    delete[] ptrs;
    delete[] rhs;

    //Return the integrator status.
    return false;
}

//Gets the soultion vector.    
Eigen::VectorXd 
MumpsSolver::GetSolution(){
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

    //Transform Mumps solution into Eigen Vector.
    Eigen::VectorXd x(numberOfFreeDofs);
    for(unsigned int k = 0; k < numberOfFreeDofs; k++)
        x(k) = sol[k];

    return x;
}
