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
//   [1]
//
// Description:
///This file contains the abstract "LinearSystem object" declarations, which 
///defines the library to solve the linear system of equations Ax = b.
//------------------------------------------------------------------------------

#ifndef _LINEARSYSTEM_HHP_
#define _LINEARSYSTEM_HHP_

#include <Eigen/Dense>
#include <Eigen/SparseCore>

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      February 21, 2018
/// @version   1.0
/// @file      LinearSystem.hpp
/// @class     LinearSystem
/// @see       EigenSolver.hpp MumpsSolver.hpp PetscSolver.hpp
/// @brief     Virtual class for solving a linear system of the form Ax = b
class LinearSystem {

    public:
        ///Creates a LinearSystem object.
        LinearSystem();

        ///Destroys this LinearSystem object.
        virtual ~LinearSystem() = 0;    

        ///Gets the solution vector.
        ///@return The solution vector.
        virtual Eigen::VectorXd GetSolution() = 0;

        ///Solve the linear system.
        ///@param A The matrix with the left-hand side.
        ///@param b The vector with the right-hand side.
        ///@return Whether or not the system was successfully solved.
        ///@see Integrator::ComputeEffectiveForce(), Integrator::ComputeEffectiveStiffness().
        virtual bool SolveSystem(Eigen::SparseMatrix<double> &A, Eigen::VectorXd &b) = 0;

    private:
        //No private varibles.
};

#endif
