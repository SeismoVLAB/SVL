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
///This file contains the "QuasiStatic" integrator declaration, this integrator 
///neglects any damping and inertial forces present in the analysis, and is 
///meant to be used in Static analysis with linear or non-linear algorithms.
//------------------------------------------------------------------------------

#ifndef _QUASISTATIC_HPP_
#define _QUASISTATIC_HPP_

#include <Eigen/Dense>
#include <Eigen/SparseCore>

#include "Mesh.hpp"
#include "LoadCombo.hpp"
#include "Assembler.hpp"
#include "Algorithm.hpp"
#include "Integrator.hpp"

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      February 27, 2018
/// @version   1.0
/// @file      QuasiStatic.hpp
/// @class     QuasiStatic
/// @see       Algorithm.hpp Integrator.hpp Analysis.hpp
/// @brief     Class for defining an static integrator that neglects any damping and inertial forces present in the analysis
class QuasiStatic : public Integrator{

    public:
        ///Creates a QuasiStatic object.
        ///@param mesh Pointer to the Mesh container to extract Node and Element.
        ///@param mtol Threshold for which a mass value will be neglected.
        ///@param ktol Threshold for which a stiffness value will be neglected.
        ///@param ftol Threshold for which a force value will be neglected.
        ///@note More details can be found at @ref linkQuasiStatic.
        ///@see Assembler::MassTolerance, Assembler::StiffnessTolerance, Assembler::ForceTolerance.
        QuasiStatic(std::shared_ptr<Mesh> &mesh, double mtol=1E-12, double ktol=1E-12, double ftol=1E-12);

        ///Destroys this QuasiStatic object.
        ~QuasiStatic();

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
        ///@note More details can be found at @ref linkQuasiStatic.
        Eigen::VectorXd& GetDisplacements();    

        ///Gets the velocity vector.
        ///@return Vector with the velocity states at current time step.
        ///@note More details can be found at @ref linkQuasiStatic.
        Eigen::VectorXd& GetVelocities();

        ///Gets the acceleration vector.
        ///@return Vector with the acceleration states at current time step.
        ///@note More details can be found at @ref linkQuasiStatic.
        Eigen::VectorXd& GetAccelerations();

        ///Gets the PML history vector.
        ///@return Vector with the displacement states at current time step.
        ///@note More details can be found at @ref linkExtendedNewmarkBeta.
        Eigen::VectorXd& GetPMLHistoryVector();  

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

        ///Gets the incremental nodal support motion vector.
        ///@param mesh Pointer to the Mesh object where Node are stored.
        ///@param factor The incremental load factor.
        ///@param k The time step number to be solved.
        ///@return Vector with the incremental support motion displacement.
        ///@see Node::GetSupportMotion(), Assembler::ComputeSupportMotionIncrement().
        void ComputeSupportMotionVector(std::shared_ptr<Mesh> &mesh, Eigen::VectorXd &Feff, double factor=1.00, unsigned int k=0);

        ///Gets the effective force assiciated to the QuasiStatic integrator.
        ///@param mesh Pointer to the Mesh object where Node and Element are stored.
        ///@param Feff Vector that stores the effective force.
        ///@param factor The incremental load factor.
        ///@param k The time step number to be solved.
        ///@note More details can be found at @ref linkQuasiStatic.
        ///@see Assembler::ComputeInternalForceVector(), Assembler::ComputeExternalForceVector().
        void ComputeEffectiveForce(std::shared_ptr<Mesh> &mesh, Eigen::VectorXd &Feff, double factor=1.00, unsigned int k=0);

        ///Gets the effective stiffness assiciated to the QuasiStatic integrator.
        ///@param mesh Pointer to the Mesh object where Node and Element are stored.
        ///@param Keff Matrix that stores the effective stiffness.
        ///@note More details can be found at @ref linkQuasiStatic.
        ///@see Assembler::ComputeMassMatrix(), Assembler::ComputeStiffnessMatrix(), Assembler::ComputeDampingMatrix().
        void ComputeEffectiveStiffness(std::shared_ptr<Mesh> &mesh, Eigen::SparseMatrix<double> &Keff);

    protected:
        ///Integration time step.
        double TimeStep;

        ///Total initial/previous displacements.
        Eigen::VectorXd U;

        ///Total previous pml history values.
        Eigen::VectorXd Ubar;

        ///Model stiffness matrix.
        Eigen::SparseMatrix<double> K; 

        ///The static solver algorithm.
        std::weak_ptr<Algorithm> theAlgorithm;

        ///The finite element assembler.
        std::unique_ptr<Assembler> theAssembler;
};

#endif
