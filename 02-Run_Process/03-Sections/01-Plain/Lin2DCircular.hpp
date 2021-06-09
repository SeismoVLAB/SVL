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
///This file contains the "Lin2DCircular" section declarations, which defines a 
///Circular geometry to be used in 2D analysis with linear elastic material for 
///frame elements.
//------------------------------------------------------------------------------

#ifndef _LIN2DCIRCULAR_HPP_
#define _LIN2DCIRCULAR_HPP_

#include <string>
#include <memory>
#include <Eigen/Dense>

#include "Section.hpp"
#include "Material.hpp"

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      May 7, 2020
/// @version   1.0
/// @file      Lin2DCircular.hpp
/// @class     Lin2DCircular
/// @see       Section.hpp
/// @brief     Class for creating a solid circular section for 2D analysis with linear elastic material for a frame elements
class Lin2DCircular : public Section{

    public:
        ///Creates a Section to be specified at a Gauss-point in an Element.
        ///@param r The radius.
        ///@param material The Material that this Section is made out of.
        ///@param theta The angle of rotation about x3-axis.
        ///@param ip The insertion point where the section axes are located.
        ///@note Details about insertion point and local axis rotation can be found at @ref linkSectionInsertionPoint and @ref linkSectionLocalAxis.
        ///@see Lin2DCircular::r, Lin2DCircular::theMaterial, Lin2DCircular::Theta, Lin2DCircular::InsertPoint.
        Lin2DCircular(double r, std::unique_ptr<Material> &material, double theta=0.0, unsigned int ip=10);

        ///Destroys this Lin2DCircular object.
        ~Lin2DCircular();

        ///Clone the 'Lin2DCircular' section.
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

    protected:
        ///Computes the Angle area.
        ///@return A constant with the cross-section area.
        ///@see The geometric properties can be revisited at @ref linkLin3DCircular.
        double GetArea();

        ///Computes the Angle shear area.
        ///@return A constant with the shear area along axis x2.
        ///@see The geometric properties can be revisited at @ref linkLin3DCircular.
        double GetShearArea2();

        ///Computes the Angle shear area.
        ///@return A constant with the shear area along axis x3.
        ///@see The geometric properties can be revisited at @ref linkLin3DCircular.
        double GetShearArea3();

        ///Computes the Angle flexural inertia.
        ///@return A constant with the moment of inertia about axis x2.
        ///@see The geometric properties can be revisited at @ref linkLin3DCircular.
        double GetInertiaAxis2();

        ///Computes the Angle flexural inertia.
        ///@return A constant with the moment of inertia about axis x3.
        ///@see The geometric properties can be revisited at @ref linkLin3DCircular.
        double GetInertiaAxis3();

        ///Gets the section centroid.
        ///@param zcm The centroid position along x3 axis.
        ///@param ycm The centroid position along x2 axis.
        ///@see The geometric properties can be revisited at @ref linkLin3DCircular.
        void ComputeSectionCenter(double &zcm, double &ycm);

    private:
        ///Circle Radius.
        double r;

        ///Rotation angle.
        double Theta;

        ///Local axis coordinate.
        unsigned int InsertPoint;

        ///Generalized Strain vector.
        Eigen::VectorXd Strain;

        ///The section's material.
        std::unique_ptr<Material> theMaterial;
};

#endif
