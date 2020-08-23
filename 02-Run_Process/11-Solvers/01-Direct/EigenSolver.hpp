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
//   [1] https://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html
//
// Description:
///This file contains the "EigenSolver" solver declaration; the object 
///performs Cholesky-decomposition on the matrix A, and then solves the 
///upper triangular system. This function is mean to be amployed to 
///symmetric and positive-definite structure of A.
//------------------------------------------------------------------------------

#ifndef _EIGENSOLVER_HPP_
#define _EIGENSOLVER_HPP_

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>

#include "LinearSystem.hpp"

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      March 8, 2018
/// @version   1.0
/// @file      EigenSolver.hpp
/// @class     EigenSolver
/// @see       LinearSystem.hpp
/// @brief     Class for solving a linear system using a Cholesky-decomposition on the effective stiffness matrix
class EigenSolver : public LinearSystem {

    public:
        ///Creates a EigenSolver object.
        ///@param flag Whether the analysis is linear or non-linear.
        ///@note More details can be found at @ref linkEigenSolver.
        ///@see EigenSolver::Flag.
        EigenSolver(bool flag = false);

        ///Destroys this EigenSolver object.
        ~EigenSolver();

        ///Solve the linear system.
        ///@param A The matrix with the left-hand side.
        ///@param b The vector with the right-hand side.
        ///@return Whether or not the system was successfully solved.
        ///@note More details can be found at @ref linkEigenSolver.
        ///@see Integrator::ComputeEffectiveForce(), Integrator::ComputeEffectiveStiffness().
        bool SolveSystem(Eigen::SparseMatrix<double> &A, Eigen::VectorXd &b);

        ///Gets the soultion vector.
        ///@return The solution vector.
        ///@see EigenSolver::x.
        Eigen::VectorXd GetSolution();

    protected:
        ///Analysis flag.
        bool Flag;

        ///Eigen Factorization Flag.
        bool Factored;

        ///Eigen Initialization Flag.
        bool Initialized;

        ///Vector of unknowns. 
        Eigen::VectorXd x;

        ///Direct LDLT factorization solver. 
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > Solver;
};

#endif
