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
//   [1]
//
// Description:
///This file contains the "Linear" algorithm declarations, which solves the 
///linear system of equations in just one-step; therefore, no-iteration or 
///convergence error is performed, and is meant to be used for linear analysis. 
//------------------------------------------------------------------------------

#ifndef _LINEAR_HPP_
#define _LINEAR_HPP_

#include <memory>
#include <Eigen/Dense>
#include <Eigen/SparseCore>

#include "Mesh.hpp"
#include "LoadCombo.hpp"
#include "Algorithm.hpp"
#include "Integrator.hpp"
#include "LinearSystem.hpp"

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      March 21, 2018
/// @version   1.0
/// @file      Linear.hpp
/// @class     Linear
/// @see       Algorithm.hpp Integrator.hpp
/// @brief     Class for solving a linear system in just one-step; therefore, no-iteration or convergence test is performed
class Linear : public Algorithm{

    public:
        ///Creates a Linear object.
        ///@param solver Pointer to the LinearSystem solver to be employed.
        ///@param mesh Pointer to the Mesh container to extract Node and Element.
        ///@param tol The tolerance for which the solution is assumed to converge.
        ///@param ldFactor The incremental load factor for external loads.
        ///@param nIters The number of maximum itertions for convergence.
        ///@note More details can be found at @ref linkLinear.
        ///@see Linear::Tolerance, Linear::LoadFactor, Linear::nMaxIterations.
        Linear(std::unique_ptr<LinearSystem> &solver, std::shared_ptr<Mesh> &mesh, double tol=1E-6, double ldFactor=1.00, unsigned int nIters=1);

        ///Destroys this Linear object.
        ~Linear();

        ///Computes a new incremental solution.
        ///@param mesh The finite element Mesh object.
        ///@param i The time step number to be solved.
        ///@return Whether or not the algorithm was successful.
        ///@note More details can be found at @ref linkLinear.
        bool ComputeNewIncrement(std::shared_ptr<Mesh> &mesh, unsigned int i=0);

        ///Gets the displacement increment.
        ///@return The incremental displacement vector.
        ///@note More details can be found at @ref linkLinear.
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
