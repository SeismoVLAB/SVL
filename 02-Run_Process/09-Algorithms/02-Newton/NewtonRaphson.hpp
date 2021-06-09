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
//   [1] Finite Element Procedures, Bathe, K.J., Chapter 8: pages 755-759. 
//       Prentice-Hall, 1996. 
//
// Description: 
///This file contains the "NewtonRaphson" algorithm declarations, which solves 
///the non-linear system of equations in a few iterations until the norm of the 
///residuals reaches a certain tolerance or a maximun number of iterations is 
///reached; therefore, this routine is meant to be used for non-linear analysis.
//------------------------------------------------------------------------------

#ifndef _NEWTONRAPHSON_HPP_
#define _NEWTONRAPHSON_HPP_

#include <memory>
#include <Eigen/Dense>
#include <Eigen/SparseCore>

#include "Mesh.hpp"
#include "LoadCombo.hpp"
#include "Algorithm.hpp"
#include "Integrator.hpp"
#include "LinearSystem.hpp"

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      September 26, 2018
/// @version   1.0
/// @file      NewtonRaphson.hpp
/// @class     NewtonRaphson
/// @see       Algorithm.hpp Integrator.hpp
/// @brief     Class for solving a non-linear system in a few iterations until the norm of the residuals reaches a certain tolerance or a maximun number of iterations is reached
class NewtonRaphson : public Algorithm {

    public:
        ///Creates a Linear object.
        ///@param solver Pointer to the LinearSystem solver to be employed.
        ///@param mesh Pointer to the Mesh container to extract Node and Element.
        ///@param tol The tolerance for which the solution is assumed to converge.
        ///@param nIters The number of maximum itertions for convergence.
        ///@param flag The convergence test to be performeds.
        ///@param ldFactor The incremental load factor for external loads.
        ///@note More details can be found at @ref linkNewtonRaphson.
        ///@see NewtonRaphson::Tolerance, NewtonRaphson::LoadFactor, NewtonRaphson::nMaxIterations.
        NewtonRaphson(std::unique_ptr<LinearSystem> &solver, std::shared_ptr<Mesh> &mesh, double tol=1E-6, unsigned int nIters=50, unsigned int flag=1,  double ldFactor=0.02);

        ///Destroys this NewtonRaphson object.
        ~NewtonRaphson();

        ///Computes a new incremental solution.
        ///@param mesh The finite element Mesh object.
        ///@param i The time step number to be solved.
        ///@return Whether or not the algorithm was successful.
        ///@note More details can be found at @ref linkNewtonRaphson.
        bool ComputeNewIncrement(std::shared_ptr<Mesh> &mesh, unsigned int i=0);

        ///Gets the displacement increment.
        ///@return The incremental displacement vector.
        ///@note More details can be found at @ref linkNewtonRaphson.
        Eigen::VectorXd GetDisplacementIncrement();    

        ///Set the load factor.
        ///@param factor The incremental load factor.
        void SetLoadFactor(double factor);

        ///Sets the integrator for this algorithm.
        ///@param integrator Pointer to the Integrator to obtain the effective stiffness and force.
        void SetIntegrator(std::shared_ptr<Integrator> &integrator);

    protected:
        ///Convergence tolerance.
        double Tolerance;

        ///Incremental load factor.
        double LoadFactor;

        ///Maximum allowed number of iterations.
        unsigned int nMaxIterations;

        ///Incremental displacement.
        Eigen::VectorXd dU;
        
        ///The linear system solver.
        std::unique_ptr<LinearSystem> theSolver;

        ///The time-domain integrator.
        std::shared_ptr<Integrator> theIntegrator;
};

#endif
