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
///This file contains the "Lin3DThinArea" section declarations, which defines 
///an area geometry for 3D analysis with a material for Shell, Plate, and 
///Membrane elements.
//------------------------------------------------------------------------------

#ifndef _LIN3DTHINAREA_HPP_
#define _LIN3DTHINAREA_HPP_

#include <string>
#include <memory>
#include <Eigen/Dense>

#include "Section.hpp"
#include "Material.hpp"

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      May 8, 2020
/// @version   1.0
/// @file      Lin3DThinArea.hpp
/// @class     Lin3DThinArea
/// @see       Section.hpp
/// @brief     Class for creating an area section for 3D analysis with linear elastic material for a shell elements
class Lin3DThinArea : public Section{

    public:
        ///Creates a Section to be especified at a Gauss-point in an Element.
        ///@param t The thickness the section.
        ///@param material The Material that this Section is made out of.
        ///@note Details about local axis rotation can be found at @ref linkSectionLocalAxis.
        ///@see Lin3DThinArea::h, Lin3DThinArea::theMaterial.
        Lin3DThinArea(double t, std::unique_ptr<Material> &material);

        ///Destroys this Lin3DThinArea object.
        ~Lin3DThinArea();

        ///Clone the 'Lin3DThinArea' section.
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
        ///@note The initial tangent stiifness matrix is computed when the generalized strain vector is zero.
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

        ///Update the section state for this iteration.
        ///@param strain Vector with the strain components at this Gauss-point.
        ///@param cond If the the elatic/platic material components will be updated.
        ///@note This function computes the strain and tanget stiffness matrix once the trial strain converged.
        void UpdateState(const Eigen::VectorXd strain, const unsigned int cond);

    private:
        ///Thickness.
        double t;

        ///Generalized Strain vector.
        Eigen::VectorXd Strain;

        ///The section's material.
        std::unique_ptr<Material> theMaterial;
};

#endif
