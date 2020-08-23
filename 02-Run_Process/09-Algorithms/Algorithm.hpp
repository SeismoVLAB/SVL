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
///This file contains the "Algorithm object" declarations, and defines how the 
///soluction between two steps is going to be carried out.
//------------------------------------------------------------------------------

#ifndef _ALGORITHM_HPP_
#define _ALGORITHM_HPP_

#include <map>
#include <mpi.h>
#include <memory>
#include <Eigen/Dense>
#include <Eigen/SparseCore>

#include "Mesh.hpp"
#include "LoadCombo.hpp"

class Integrator;

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      November 21, 2018
/// @version   1.0
/// @file      Algorithm.hpp
/// @class     Algorithm
/// @see       Linear.hpp NewtonRaphson.hpp
/// @brief     Virtual class for defining how the solution between two steps is going to be carried out
class Algorithm{

    public:
        ///Creates a Algorithm object.
        ///@param mesh Pointer to the Mesh container to extract Node and Element.
        ///@param flag The convergence test to be performed.
        ///@param NormFactor The relative unbalanced factor for convergence.
        ///@see Algorithm::flag, Algorithm::NormFactor.
        Algorithm(const std::shared_ptr<Mesh> &mesh, unsigned int flag=1, double NormFactor=1.0);

        ///Destroys this Algorithm object.
        virtual ~Algorithm() = 0;

        ///Computes a new incremental solution.
        ///@param mesh The finite element Mesh object.
        ///@param i The time step number to be solved.
        ///@return Whether or not the algorithm was successful.	
        virtual bool ComputeNewIncrement(std::shared_ptr<Mesh> &mesh, unsigned int i) = 0;

        ///Gets the displacement increment.
        ///@return The incremental displacement vector.	
        virtual Eigen::VectorXd GetDisplacementIncrement() = 0;

        ///Sets the integrator for this algorithm.
        ///@param integrator Pointer to the Integrator to obtain the effective stiffness and force.
        virtual void SetIntegrator(std::shared_ptr<Integrator> &integrator) = 0;

        ///Set the load factor.
        ///@param factor The incremental load factor.
        virtual void SetLoadFactor(double factor) = 0;

        ///Update the incremental state variables in the mesh.
        ///@param mesh Pointer to the Mesh container to extract Node and Element.
        ///@param dU Vector with the converged incremental displacements.
        ///note This function sets incremental displacements at each Node in Mesh.
        ///@see Node::SetIncrementalDisplacements().
        void UpdateStatesIncrements(std::shared_ptr<Mesh> &mesh, const Eigen::VectorXd &dU);

        ///Construct the residual vector force from each processor.
        ///@param Feff The residual vector in this partition.
        ///@param Residual The value of the residual vector norm.
        ///@note In parallel execution Feff is gathered across processors.
        void ReducedParallelResidual(const Eigen::VectorXd &Feff, double &Residual);

        ///Computes convergence tests for this algorithm.
        ///@param Force The residual force vector.
        ///@param Delta The incremental displacement norm.
        ///@param isFirstIteration Whether this step in Algorithm is the first iteration.
        ///@note If Algorithm::flag = 1: Unbalanced Force Norm. 
        ///@note If Algorithm::flag = 2: Increment Displacement Norm. 
        ///@note If Algorithm::flag = 3: Relative Unbalanced Force Norm. 
        ///@note If Algorithm::flag = 4: Relative Increment Displacement Norm.
        double ComputeConvergence(const Eigen::VectorXd &Force, double Delta, bool isFirstIteration);

    protected:
        ///Operator that enforced restrain/constraint. 
        ///@see Mesh::GetTotalToFreeMatrix().
        Eigen::SparseMatrix<double> Total2FreeMatrix;

    private:
        ///Convergence test to be performed.
        unsigned int flag;

        ///Relative Unbalanced Factor.
        double NormFactor;
};

#endif
