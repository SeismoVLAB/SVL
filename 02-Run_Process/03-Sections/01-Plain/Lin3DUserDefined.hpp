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
///This file contains the "Lin3DUserDefined" section declarations, which defines 
///a general geometry for 3D analysis with linear elastic material for 
///frame elements.
//------------------------------------------------------------------------------

#ifndef _LIN3DUSERDEFINED_HPP_
#define _LIN3DUSERDEFINED_HPP_

#include <vector>
#include <string>
#include <memory>
#include <Eigen/Dense>

#include "Section.hpp"
#include "Material.hpp"

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      July 16, 2021
/// @version   1.0
/// @file      Lin3DUserDefined.hpp
/// @class     Lin3DUserDefined
/// @see       Section.hpp
/// @brief     Class for creating a general solid section for 3D analysis with linear elastic material for a frame elements
class Lin3DUserDefined : public Section{

    public:
        ///Creates a Section to be specified at a Gauss-point in an Element.
        ///@param properties The geometric properties of the section.
        ///@param material The Material that this Section is made out of.
        ///@param theta The angle of rotation about x3-axis.
        ///@note Details about local axis rotation can be found at @ref linkSectionLocalAxis.
        ///@see Lin3DUserDefined::theMaterial, Lin3DUserDefined::Theta.
        Lin3DUserDefined(std::vector<double> properties, std::unique_ptr<Material> &material, double theta=0.0);

        ///Destroys this Lin3DUserDefined object.
        ~Lin3DUserDefined();

        ///Clone the 'Lin3DUserDefined' section.
        ///@return A Section pointer to the derived class.
        ///@see Section.
        std::unique_ptr<Section> CopySection();

        ///Returns the resultant (generalised) strain vector over the section. 
        ///@return Vector with the resultant strain components. 
        Eigen::VectorXd GetStrain();

        ///Returns the resultant (generalised) stress vector over the section. 
        ///@return Vector with the resultant stress components.  
        Eigen::VectorXd GetStress();

        ///Access the section density matrix.
        ///@return The section density matrix.
        Eigen::MatrixXd GetDensity();

        ///Returns the section stiffness matrix.
        ///@return Matrix with the section consistent tangent stiffness matrix.
        Eigen::MatrixXd GetTangentStiffness();

        ///Returns the section initial stiffness matrix.
        ///@return Matrix with the initial section tangent stiffness matrix.
        ///@note The initial tangent stiffness matrix is computed when the generalized strain vector is zero.
        Eigen::MatrixXd GetInitialTangentStiffness();

        ///Returns the section strain at given position.
        ///@param x3 Position on the x3-axis where the strain is evaluated.
        ///@param x2 Position on the x2-axis where the strain is evaluated.
        ///@return Vector with the strain components at coordinate (x3,x2).
        Eigen::VectorXd GetStrainAt(double x3=0, double x2=0);

        ///Returns the section stress at given position.
        ///@param x3 Position on the x3-axis where the stress is evaluated.
        ///@param x2 Position on the x2-axis where the stress is evaluated.
        ///@return Vector with the stress components at coordinate (x3,x2).
        Eigen::VectorXd GetStressAt(double x3=0, double x2=0);

        ///Perform converged section state update.
        ///@note This function sets the trail stress and strain as converged.
        void CommitState();

        ///Reverse the section states to previous converged state.
        ///@note This function returns the material states to previous converged states.
        void ReverseState();

        ///Brings the section states to its initial state.
        ///@note This function returns the material states to the beginning.
        void InitialState();

        ///Update the section state for this iteration.
        ///@param strain Vector with the strain components at this Gauss-point.
        ///@param cond If the the elastic/plastic material components will be updated.
        ///@note This function computes the strain and tanget stiffness matrix once the trial strain converged.
        void UpdateState(const Eigen::VectorXd strain, const unsigned int cond);

    private:
        ///Area.
        double A;

        ///Shear area 2-axis.
        double As2;

        ///Shear area 3-axis.
        double As3;

        ///Polar Inertia 1-axis.
        double J;

        ///Inertia 2-axis.
        double I22;

        ///Inertia 3-axis.
        double I33;

        ///Product of Inertia.
        double I23;

        ///Rotation angle.
        double Theta;

        ///Generalized Strain vector.
        Eigen::VectorXd Strain;

        ///The section's material.
        std::unique_ptr<Material> theMaterial;
};

#endif
