//==============================================================================
//
//                    (Seismo)s (V)irtual (Lab)oratory
//             Module for Serial and Parallel Analysis of seismic 
//         wave propagation and soil-structure interaction simulation
//         Copyright (C) 2018, The California Institute of Technology
//                           All Rights Reserved.
//
// Commercial use of this program without express permission of the California
// Institute of Technology, is strictly  prohibited. See  file "COPYRIGHT"  in
// main  directory  for  information on  usage  and  redistribution, and for a
// DISCLAIMER OF ALL WARRANTIES.
//
//==============================================================================
//
// Written by:
//   Danilo S. Kusanovic (dkusanov@caltech.edu)
//   Elnaz E. Seylabi    (elnaze@unr.edu)
//
// Supervised by:
//   Domniki M. Asimaki  (domniki@caltech.edu)
//
// References : 
//   [1] A Fully Asynchronous Multifrontal Solver Using Distributed Dynamic 
//       Scheduling, P. R. Amestoy and I. S. Duff and J. Koster and J.-Y. 
//       L'Excellent, SIAM Journal on Matrix Analysis and Applications, 
//       volume 23 (1) 15-41, 2003. 
//   [2] Hybrid scheduling for the parallel solution of linear systems, P. R. 
//       Amestoy and A. Guermouche and J.-Y. L'Excellent and S. Pralet, Journal 
//       of Parallel Computing, volume 32(2), 136-156, 2006.
//
// Description:
///This file contains the "MumpsSolver" solver declaration; this object 
///performs a (MU)ltifrontal (M)assively (P)arallel (S)parse direct Solver on 
///the matrix A, assuming that such matrix is symmetric positive definite (SPD), 
///general symmetric (SYM) or unsymmetric (USYM) matrix
//------------------------------------------------------------------------------

#ifndef _MUMPSSOLVER_HPP_
#define _MUMPSSOLVER_HPP_

#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654
#define ICNTL(I) icntl[(I)-1]

#include <dmumps_c.h>
#include <Eigen/Dense>
#include <Eigen/SparseCore>

#include "LinearSystem.hpp"

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      June 8, 2020
/// @version   1.0
/// @file      MumpsSolver.hpp
/// @class     MumpsSolver
/// @see       LinearSystem.hpp Algorithm.hpp
/// @brief     Class for solving a linear system using the (MU)ltifrontal (M)assively (P)arallel (S)parse direct Solver assuming the matrix is symmetric positive definite (SPD), general symmetric (SYM) or unsymmetric (USYM).
class MumpsSolver : public LinearSystem {

    public:
        ///Creates a MumpsSolver object.
        ///@param option The type of left-hand side matrix: SPD, SYM, or USYM.
        ///@param flag Whether the analysis is linear or non-linear.
        ///@note More details can be found at @ref linkMumpsSolver.
        ///@see MumpsSolver::Flag, MumpsSolver::Option, MumpsSolver::n.
        MumpsSolver(unsigned int option=1, bool flag = false);

        ///Destroys this MumpsSolver object.
        ~MumpsSolver();

        ///Solve the linear system.
        ///@param A The matrix with the left-hand side.
        ///@param b The vector with the right-hand side.
        ///@return Whether or not the system was successfully solved.
        ///@note More details can be found at @ref linkMumpsSolver.
        ///@see Integrator::ComputeEffectiveForce(), Integrator::ComputeEffectiveStiffness().
        bool SolveSystem(Eigen::SparseMatrix<double> &A, Eigen::VectorXd &b);

        ///Gets the soultion vector.
        ///@return The solution vector.
        ///@see MumpsSolver::sol.
        Eigen::VectorXd GetSolution();

    protected:
        ///Number of non-zero values.
        unsigned int nz;

        ///Matrix structure 0: SPD, 1:Symmetric, 0: Unsymmetric.
        unsigned int Option;

        ///Analysis flag.
        bool Flag;

        ///MUMPS Factorization Flag.
        bool Factored;

        ///MUMPS Initialization Flag.
        bool Initialized;

        ///Solution vector.
        double *sol;

        ///MUMPS Structure of parameters.
        DMUMPS_STRUC_C id;
};

#endif
