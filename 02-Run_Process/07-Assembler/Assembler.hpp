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
///This file contains the "Assembler object" declarations, this object computes 
///the system mass, damping and stiffness matrices as well as the internal and 
///external force vector. 
//------------------------------------------------------------------------------

#ifndef _ASSEMBLER_HPP_
#define _ASSEMBLER_HPP_

#include <vector>
#include <memory>
#include <Eigen/Dense>
#include <Eigen/SparseCore>

#include "Mesh.hpp"
#include "LoadCombo.hpp"

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      November 22, 2018
/// @version   1.0
/// @file      Assembler.hpp
/// @class     Assembler
/// @see       QuasiStatic.hpp CentralDifference.hpp NewmarkBeta.hpp
/// @brief     Class for assembling the finite element matrices and vector
class Assembler{

    public:
        ///Creates an Assembler object.
        ///@note More details can be found at @ref linkAssembler.
        Assembler();

        ///Destroys this Assembler object.
        ~Assembler();

        ///Set the mass assembly allowed tolerance.
        ///@param tol Threshold for which a mass value will be neglected.
        void SetMassTolerance(double tol);

        ///Set the mass assembly allowed tolerance.
        ///@param tol Threshold for which a force value will be neglected.
        void SetForceTolerance(double tol);

        ///Set the stiffness assembly allowed tolerance.
        ///@param tol Threshold for which a stiffness value will be neglected.
        void SetStiffnessTolerance(double tol);

        ///Sets the load combination to be used.
        ///@param combo Pointer to the load combination to be simulated.
        void SetLoadCombination(std::shared_ptr<LoadCombo> &combo);

        ///Assemble the mass matrix.
        ///@param mesh The finite element Mesh object.
        ///return A matrix with the model mass.
        ///@note More details can be found at @ref linkAssembler.
        ///@see Element::ComputeMassMatrix(), Integrator::ComputeEffectiveStiffness().
        Eigen::SparseMatrix<double> ComputeMassMatrix(std::shared_ptr<Mesh> &mesh); 

        ///Assemble the stiffness matrix.
        ///@param mesh The finite element Mesh object.  
        ///@return A matrix with the model stiffness.
        ///@note More details can be found at @ref linkAssembler.
        ///@see Element::ComputeStiffnessMatrix(), Integrator::ComputeEffectiveStiffness().
        Eigen::SparseMatrix<double> ComputeStiffnessMatrix(std::shared_ptr<Mesh> &mesh);

        ///Assemble the stiffness matrix.
        ///@param mesh The finite element Mesh object.
        ///@return A matrix with the model damping.
        ///@note More details can be found at @ref linkAssembler.
        ///@see Element::ComputeDampingMatrix(), Integrator::ComputeEffectiveStiffness().
        Eigen::SparseMatrix<double> ComputeDampingMatrix(std::shared_ptr<Mesh> &mesh);

        ///Assemble the integrated history matrix for Perfectly-Matched Layer (PML).
        ///@param mesh The finite element Mesh object.
        ///@return A vector with the PML history values.
        ///@note More details can be found at @ref linkAssembler.
        ///@see Element::ComputeDampingMatrix(), Integrator::ComputeEffectiveStiffness().
        Eigen::SparseMatrix<double> ComputePMLHistoryMatrix(std::shared_ptr<Mesh> &mesh);

        ///Assemble the internal force vector.
        ///@param mesh The finite element Mesh object.
        ///return A vector with the model total internal force.
        ///@see Element::ComputeInternalForces(), Integrator::ComputeEffectiveForce().   
        Eigen::VectorXd ComputeInternalForceVector(std::shared_ptr<Mesh> &mesh);

        ///Assemble the internal elastic, inertial, and viscous force vector.
        ///@param mesh The finite element Mesh object.
        ///return A vector with the model total (inertial, viscous, elastic) internal force.
        ///@see Element::ComputeInternalDynamicForces().   
        Eigen::VectorXd ComputeDynamicInternalForceVector(std::shared_ptr<Mesh> &mesh);

        ///Assemble the external force vector.
        ///@param mesh The finite element Mesh object.
        ///@param k The time step when the force is evaluated.
        ///@return A vector with the model total external force.
        ///@note More details can be found at @ref linkAssembler.
        ///@see Element::ComputeSurfaceForces(), Element::ComputeBodyForces(), Integrator::ComputeEffectiveForce().
        Eigen::VectorXd ComputeExternalForceVector(std::shared_ptr<Mesh> &mesh, unsigned int k);

        ///Assemble the support motion displacement vector.
        ///@param mesh The finite element Mesh object.
        ///@param k The time step when the support motion is evaluated.
        ///@return A vector with the model support motion displacements.
        ///@see Integrator::ComputeEffectiveForce(), Integrator::ComputeEffectiveStiffness() 
        Eigen::VectorXd ComputeSupportMotionIncrement(std::shared_ptr<Mesh> &mesh, unsigned int k);

        ///Assemble the integrated history vector for Perfectly-Matched Layer (PML).
        ///@param mesh The finite element Mesh object.
        ///@return A vector with the PML history values.
        ///@see Integrator::ComputeEffectiveForce(), Integrator::ComputeEffectiveStiffness() 
        Eigen::VectorXd ComputePMLHistoryVector(std::shared_ptr<Mesh> &mesh);

    protected:
        ///Adds the inertial forces contribution associated to the nodes.
        ///@param mesh The finite element Mesh object.
        ///@param DynamicForces The vector with the node inertial forces.
        ///@see Node, Node::GetInertialForces(). 
        void AddNodeInertiaForces(std::shared_ptr<Mesh> &mesh, Eigen::VectorXd &DynamicForces);

        ///Adds the elastic, inertial, and viscous forces associated to the elements.
        ///@param mesh The finite element Mesh object.
        ///@param DynamicForces The vector with the Element dynamic internal forces.
        ///@see Element, Element::ComputeInternalDynamicForces(). 
        void AddElementDynamicForces(std::shared_ptr<Mesh> &mesh, Eigen::VectorXd &DynamicForces);

        ///Assemble Nodal Mass matrix.  
        ///@param mesh The finite element Mesh object.
        ///@param nodemass A matrix that stores the Node mass contribution. 
        ///@see Node, Node::GetMass(). 
        void AssembleNodalMass(std::shared_ptr<Mesh> &mesh, Eigen::SparseMatrix<double>& nodemass);

        ///Assemble Element Mass matrix.
        ///@param mesh The finite element Mesh object.
        ///@param elemmass A matrix that stores the Element mass contribution. 
        ///@see Element, Element::ComputeMassMatrix().
        void AssembleElementMass(std::shared_ptr<Mesh> &mesh, Eigen::SparseMatrix<double>& elemmass);

    private:
        ///Define mass tolerance.
        double MassTolerance;

        ///Define stiffness tolerance.
        double StiffnessTolerance;

        ///Define stiffness tolerance.
        double ForceTolerance;

        ///Load combination number.
        std::shared_ptr<LoadCombo> LoadCombination;
};

#endif