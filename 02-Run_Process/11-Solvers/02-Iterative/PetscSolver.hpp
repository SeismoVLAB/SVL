//==============================================================================
//
//                       Seismo Virtual Laboratory
//             Module for Serial and Parallel Analysis of seismic 
//         wave propagation and soil-structure interaction simulation
//         Copyright (C) 2018-2021, The California Institute of Technology
//                         All Rights Reserved.
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
//   [1] Kris Buschelman, Lisandro Dalcin, Alp Dener, Victor Eijkhout, William D.
//       Gropp, Dmitry Karpeyev, Dinesh Kaushik, Matthew G. Knepley, Dave A.
//       May, Lois Curfman McInnes, Richard Tran Mills, Todd Munson, Karl Rupp,
//       Patrick Sanan, Barry F. Smith, Stefano Zampini, Hong Zhang, and Hong
//       Zhang. PETSc Web page. https://www.mcs.anl.gov/petsc, 2019.
//       URL https://www.mcs.anl.gov/petsc
//
//   [2] Satish Balay, Shrirang Abhyankar, Mark F. Adams, Jed Brown, Peter Brune,
//       Kris Buschelman, Lisandro Dalcin, Alp Dener, Victor Eijkhout, William D.
//       Gropp, Dmitry Karpeyev, Dinesh Kaushik, Matthew G. Knepley, Dave A.
//       May, Lois Curfman McInnes, Richard Tran Mills, Todd Munson, Karl Rupp,
//       Patrick Sanan, Barry F. Smith, Stefano Zampini, Hong Zhang, and Hong
//       Zhang. PETSc users manual. Technical Report ANL-95/11 - Revision 3.12,
//       Argonne National Laboratory, 2019. URL https://www.mcs.anl.gov/petsc
//
// Description:
///This file contains the "PetscSolver" solver declaration; this object is an 
///iterative solver, and it can be employed to symmetric positive-definite (SPD),  
///general symmetric (SYM) or unsymmetric (USYM) structure of A.
//------------------------------------------------------------------------------

#ifndef _PETSCSOLVER_HPP_
#define _PETSCSOLVER_HPP_

#include <Eigen/Dense>
#include <Eigen/SparseCore>

#include "LinearSystem.hpp"

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu) 
/// @date      June 16, 2020 
/// @version   1.0 
/// @file      PetscSolver.hpp 
/// @class     PetscSolver 
/// @see       LinearSystem.hpp Algorithm.hpp 
/// @brief     Class for solving a linear system using an iterative solver assuming the matrix can be symmetric positive-definite, general symmetric or asymmetric 
class PetscSolver : public LinearSystem {

    public:
        ///Creates a PetscSolver object.
        ///@param dnz Number of nonzero in diagonal per row.
        ///@param onz Number of nonzero in off-diagonal per row.
        ///@param tol The tolerance at which convergence is reached.
        ///@param kspnum The iterative algorithm to be employed.
        ///@note More details can be found at @ref linkPetscSolver.
        ///@see PetscSolver::n, PetscSolver::d_nz, PetscSolver::o_nz, PetscSolver::Tolerance.
        PetscSolver(unsigned int dnz, unsigned int onz, double tol=1E-12, unsigned int kspnum=0);

        ///Destroys this PetscSolver object.
        ~PetscSolver();

        ///Solve the linear system.
        ///@param A The matrix with the left-hand side.
        ///@param b The vector with the right-hand side.
        ///@return Whether or not the system was successfully solved.
        ///@note More details can be found at @ref linkPetscSolver.
        ///@see Integrator::ComputeEffectiveForce(), Integrator::ComputeEffectiveStiffness().
        bool SolveSystem(Eigen::SparseMatrix<double> &A, Eigen::VectorXd &b);

        ///Gets the soultion vector.
        ///@return The solution vector.
        ///@see PetscSolver::x.
        const Eigen::VectorXd& GetSolution();

    protected:
        ///Number of nonzero in diagonal per row. 
        unsigned int d_nz;

        ///Number of nonzero in off-diagonal per row. 
        unsigned int o_nz;

        ///Residual tolerance error.
        double Tolerance;

        ///Vector of unknowns. 
        Eigen::VectorXd x;
};

#endif
