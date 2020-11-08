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
///This file contains the "Integrator object" declarations, and defines how the 
///dynamic solution between two steps is going to be performed.
//------------------------------------------------------------------------------

#ifndef _INTEGRATOR_HPP_
#define _INTEGRATOR_HPP_

#include <memory>
#include <Eigen/Dense>
#include <Eigen/SparseCore>

#include "Mesh.hpp"
#include "LoadCombo.hpp"

class Algorithm;

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      March 17, 2018
/// @version   1.0
/// @file      Integrator.hpp
/// @class     Integrator
/// @see       QuasiStatic.hpp CentralDifference.hpp NewmarkBeta.hpp
/// @brief     Class for defining how the static/dynamic snalysis between two steps is going to be performed
class Integrator{

    public:
        ///Creates a Integrator object.
        ///@param mesh Pointer to the Mesh container to extract Node and Element.
        Integrator(const std::shared_ptr<Mesh> &mesh);

        ///Destroys this Integrator object.
        virtual ~Integrator() = 0;

        ///Initialize model matrices.
        ///@param mesh Pointer to the Mesh container to extract Node and Element.
        ///@note This function computes matrices that are constant through the analysis.
        virtual void Initialize(std::shared_ptr<Mesh> &mesh) = 0;

        ///Set the load combination.
        ///@param combo Pointer to the LoadCombo to be simulated.
        virtual void SetLoadCombination(std::shared_ptr<LoadCombo> &combo) = 0;

        ///Sets the integrator for this algorithm.
        ///@param algorithm Pointer to the Algorithm to obtain the effective stiffness and force.
        virtual void SetAlgorithm(std::shared_ptr<Algorithm> &algorithm) = 0;

        ///Gets the displacement vector.
        ///@return Vector with the displacement states at current time step.
        virtual Eigen::VectorXd& GetDisplacements() = 0;    

        ///Gets the velocity vector.
        ///@return Vector with the velocity states at current time step.
        virtual Eigen::VectorXd& GetVelocities() = 0;

        ///Gets the acceleration vector.
        ///@return Vector with the acceleration states at current time step.
        virtual Eigen::VectorXd& GetAccelerations() = 0;

        ///Gets the perfectly-matched layer history vector.
        ///@return Vector with the PML history states at current time step.
        virtual Eigen::VectorXd& GetPMLHistoryVector() = 0;

        ///Computes a new time step.
        ///@param mesh The finite element Mesh object.
        ///@param k The time step number to be solved.
        ///@return Whether or not the Integrator was successfully applied.
        virtual bool ComputeNewStep(std::shared_ptr<Mesh> &mesh, unsigned int k=0) = 0;

        ///Gets the reaction force ins this step.
        ///@param mesh Pointer to the Mesh object where Node are stored.
        ///@param k The time step number to be solved.
        ///@return Vector with the reaction forces and external forces.
        ///@note More details can be found at @ref linkReaction.
        ///@see Node::GetReaction(), Assembler::ComputeDynamicInternalForceVector().
        virtual Eigen::VectorXd ComputeReactionForce(std::shared_ptr<Mesh> &mesh, unsigned int k=0) = 0;

        ///Gets the incremental nodal support motion vector.
        ///@param mesh Pointer to the Mesh object where Node are stored.
        ///@param factor The incremental load factor.
        ///@param k The time step number to be solved.
        ///@return Vector with the incremental support motion displacement.
        ///@see Node::GetSupportMotion(), Assembler::ComputeSupportMotionIncrement().
        virtual Eigen::VectorXd ComputeSupportMotionVector(std::shared_ptr<Mesh> &mesh, double factor=1.00, unsigned int k=0) = 0;

        ///Gets the effective force assiciated to the integrator.
        ///@param mesh Pointer to the Mesh object where Node and Element are stored.
        ///@param Feff Vector that stores the effective force.
        ///@param factor The incremental load factor.
        ///@param k The time step number to be solved.
        ///@see Assembler::ComputeInternalForceVector(), Assembler::ComputeExternalForceVector().
        virtual void ComputeEffectiveForce(std::shared_ptr<Mesh> &mesh, Eigen::VectorXd &Feff, double factor=1.00, unsigned int k=0) = 0;

        ///Gets the effective stiffness assiciated to the integrator.
        ///@param mesh Pointer to the Mesh object where Node and Element are stored.
        ///@param Keff Matrix that stores the effective stiffness.
        ///@see Assembler::ComputeMassMatrix(), Assembler::ComputeStiffnessMatrix(), Assembler::ComputeDampingMatrix().
        virtual void ComputeEffectiveStiffness(std::shared_ptr<Mesh> &mesh, Eigen::SparseMatrix<double> &Keff) = 0;

    protected:
        ///Nodal support motion vector.
        Eigen::VectorXd SupportMotion;

        ///Operator that enforced restrain/constraint. 
        ///@see Mesh::GetTotalToFreeMatrix().
        Eigen::SparseMatrix<double> Total2FreeMatrix;

    private:
        //No private variables.
};

#endif
