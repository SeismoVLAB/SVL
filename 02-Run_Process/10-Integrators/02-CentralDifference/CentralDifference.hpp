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
//   [1] Finite Element Procedures, Bathe, K.J., Chapter 9: pages 770-774. 
//       Prentice-Hall, 1996. 
//
// Description:
///This file contains the "CentralDifference" integration method, which 
///inncludes any damping and inertial forces present in the analysis, also 
///this integrator is meant to be used in dynamic analysis with linear or 
///non-linear algorithms. 
//------------------------------------------------------------------------------

#ifndef _CENTRALDIFFERENCE_HPP_
#define _CENTRALDIFFERENCE_HPP_

#include <Eigen/Dense>
#include <Eigen/SparseCore>

#include "Mesh.hpp"
#include "LoadCombo.hpp"
#include "Assembler.hpp"
#include "Algorithm.hpp"
#include "Integrator.hpp"

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      June 8, 2019
/// @version   1.0
/// @file      CentralDifference.hpp
/// @class     CentralDifference
/// @see       Algorithm.hpp Integrator.hpp Analysis.hpp
/// @brief     Class for integrating the equation of motion using an explicit central difference method
class CentralDifference : public Integrator{

    public:
        ///Creates a CentralDifference object.
        ///@param mesh Pointer to the Mesh container to extract Node and Element.
        ///@param mtol Threshold for which a mass value will be neglected.
        ///@param ktol Threshold for which a stiffness value will be neglected.
        ///@param ftol Threshold for which a force value will be neglected.
        ///@note More details can be found at @ref linkCentralDifference.
        ///@see Assembler::MassTolerance, Assembler::StiffnessTolerance, Assembler::ForceTolerance.
        CentralDifference(std::shared_ptr<Mesh> &mesh, double TimeStep, double mtol=1E-12, double ktol=1E-12, double ftol=1E-12);

        ///Destroys this CentralDifference object.
        ~CentralDifference();

        ///Initialize model matrices.
        ///@param mesh Pointer to the Mesh container to extract Node and Element.
        ///@note This function computes matrices that are constant through the analysis.
        void Initialize(std::shared_ptr<Mesh> &mesh);

        ///Set the load combination.
        ///@param combo Pointer to the LoadCombo to be simulated.
        void SetLoadCombination(std::shared_ptr<LoadCombo> &combo);

        ///Sets the integrator for this algorithm.
        ///@param algorithm Pointer to the Algorithm to obtain the effective stiffness and force.
        void SetAlgorithm(std::shared_ptr<Algorithm> &algorithm);

        ///Gets the displacement vector.
        ///@return Vector with the displacement states at current time step.
        ///@note More details can be found at @ref linkCentralDifference.
        const Eigen::VectorXd& GetDisplacements();    

        ///Gets the velocity vector.
        ///@return Vector with the velocity states at current time step.
        ///@note More details can be found at @ref linkCentralDifference.
        const Eigen::VectorXd& GetVelocities();

        ///Gets the acceleration vector.
        ///@return Vector with the acceleration states at current time step.
        ///@note More details can be found at @ref linkCentralDifference.
        const Eigen::VectorXd& GetAccelerations();

        ///Gets the PML history vector.
        ///@return Vector with the displacement states at current time step.
        ///@note More details can be found at @ref linkExtendedNewmarkBeta.
        const Eigen::VectorXd& GetPMLHistoryVector();  

        ///Computes a new time step.
        ///@param mesh The finite element Mesh object.
        ///@param k The time step number to be solved.
        ///@return Whether or not the Integrator was successfully applied.
        bool ComputeNewStep(std::shared_ptr<Mesh> &mesh, unsigned int k=0);

        ///Gets the reaction force ins this step.
        ///@param mesh Pointer to the Mesh object where Node are stored.
        ///@param k The time step number to be solved.
        ///@return Vector with the reaction forces and external forces.
        ///@note More details can be found at @ref linkReaction.
        ///@see Node::GetReaction(), Assembler::ComputeDynamicInternalForceVector().
        Eigen::VectorXd ComputeReactionForce(std::shared_ptr<Mesh> &mesh, unsigned int k=0);

        ///Gets the external force vector from previous analysis.
        ///@param mesh Pointer to the Mesh object where Node and Element are stored.
        ///@param k The time step number to be solved.
        ///@see Assembler::ComputeExternalForceVector().
        Eigen::VectorXd ComputeProgressiveForce(std::shared_ptr<Mesh> &mesh, unsigned int k=0);

        ///Gets the incremental nodal support motion vector.
        ///@param mesh Pointer to the Mesh object where Node are stored.
        ///@param Feff The effective force vector to incorporate support motion forces.
        ///@param factor The incremental load factor.
        ///@param k The time step number to be solved.
        ///@return Vector with the incremental support motion displacement.
        ///@see Node::GetSupportMotion(), Assembler::ComputeSupportMotionIncrement().
        void ComputeSupportMotionVector(std::shared_ptr<Mesh> &mesh, Eigen::VectorXd &Feff, double factor=1.00, unsigned int k=0);

        ///Gets the effective force associated to the CentralDifference integrator.
        ///@param mesh Pointer to the Mesh object where Node and Element are stored.
        ///@param Feff Vector that stores the effective force.
        ///@param factor The incremental load factor.
        ///@param k The time step number to be solved.
        ///@note More details can be found at @ref linkCentralDifference.
        ///@see Assembler::ComputeInternalForceVector(), Assembler::ComputeExternalForceVector().
        void ComputeEffectiveForce(std::shared_ptr<Mesh> &mesh, Eigen::VectorXd &Feff, double factor=1.00, unsigned int k=0);

        ///Gets the effective stiffness associated to the CentralDifference integrator.
        ///@param mesh Pointer to the Mesh object where Node and Element are stored.
        ///@param Keff Matrix that stores the effective stiffness.
        ///@note More details can be found at @ref linkCentralDifference.
        ///@see Assembler::ComputeMassMatrix(), Assembler::ComputeStiffnessMatrix(), Assembler::ComputeDampingMatrix().
        void ComputeEffectiveStiffness(std::shared_ptr<Mesh> &mesh, Eigen::SparseMatrix<double> &Keff);

    protected:
        ///Integration time step.
        double dt;

        ///Total current displacements.
        Eigen::VectorXd U;

        ///Total current velocity.
        Eigen::VectorXd V;

        ///Total current acceleration.
        Eigen::VectorXd A;

        ///Total previous displacements (t - dt).
        Eigen::VectorXd Up;

        ///Total previous pml history values.
        Eigen::VectorXd Ubar;

        ///The previous stage Force vector.
        Eigen::VectorXd Fbar;

        ///Model mass matrix.
        Eigen::SparseMatrix<double> M; 

        ///Model damping matrix.
        Eigen::SparseMatrix<double> C; 

        ///Model stiffness matrix.
        Eigen::SparseMatrix<double> K; 

        ///The static solver algorithm.
        std::weak_ptr<Algorithm> theAlgorithm;

        ///The finite element assembler.
        std::unique_ptr<Assembler> theAssembler;
};

#endif
