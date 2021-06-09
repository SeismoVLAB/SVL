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
///This file contains the abstract "Section object" declarations, which computes 
///the strain, stress, and section tangent stiffness at a Gauss-point.
//------------------------------------------------------------------------------

#ifndef _SECTION_HPP_
#define _SECTION_HPP_

#include <string>
#include <memory>
#include <Eigen/Dense>

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      May 8, 2020
/// @version   1.0
/// @file      Section.hpp
/// @class     Section
/// @see       Mesh.hpp
/// @brief     Virtual class for creating a section object
class Section {

    public:
        ///Creates a Section to be specified at a Gauss-point in an Element.
        ///@param name Name of the derived class.
        ///@see Section::Name.
        Section(std::string name);

        ///Destroys this Section object.
        virtual ~Section() = 0;

        ///Clone the selected Section.
        ///@return A Section pointer to the derived class.
        ///@see Section.
        virtual std::unique_ptr<Section> CopySection() = 0;

        ///Returns the resultant (generalised) strain vector over the section.
        ///@return Vector with the resultant strain components. 
        virtual Eigen::VectorXd GetStrain() = 0;

        ///Returns the resultant (generalised) stress vector over the section.
        ///@return Vector with the resultant stress components.  
        virtual Eigen::VectorXd GetStress() = 0;

        ///Access the section density matrix.
        ///@return The section density matrix.
        virtual Eigen::MatrixXd GetDensity() = 0;

        ///Returns the section stiffness matrix.
        ///@return Matrix with the section consistent tangent stiffness matrix.
        virtual Eigen::MatrixXd GetTangentStiffness() = 0;

        ///Returns the section initial stiffness matrix.
        ///@return Matrix with the initial section tangent stiffness matrix.
        ///@note The initial tangent stiffness matrix is computed when the generalized strain vector is zero.
        virtual Eigen::MatrixXd GetInitialTangentStiffness() = 0;

        ///Returns the section strain at given position.
        ///@param x3 Position on the x3-axis where the strain is evaluated.
        ///@param x2 Position on the x2-axis where the strain is evaluated.
        ///@return Vector with the strain components at coordinate (x3,x2).
        virtual Eigen::VectorXd GetStrainAt(double x3=0, double x2=0) = 0;

        ///Returns the section stress at given position.
        ///@param x3 Position on the x3-axis where the stress is evaluated.
        ///@param x2 Position on the x2-axis where the stress is evaluated.
        ///@return Vector with the stress components at coordinate (x3,x2).
        virtual Eigen::VectorXd GetStressAt(double x3=0, double x2=0) = 0;

        ///Perform converged section state update.
        ///@note This function sets the trail stress and strain as converged.
        virtual void CommitState() = 0;

        ///Reverse the section states to previous converged state.
        ///@note This function returns the material states to previous converged states.
        virtual void ReverseState() = 0;

        ///Brings the section states to its initial state.
        ///@note This function returns the material states to the beginning.
        virtual void InitialState() = 0;

        ///Update the section state for this iteration.
        ///@param strain Vector with the strain components at this Gauss-point.
        ///@param cond If the the elastic/plastic material components will be updated.
        ///@note This function computes the strain and tanget stiffness matrix once the trial strain converged.
        virtual void UpdateState(Eigen::VectorXd strain, unsigned int cond=0) = 0;

        ///Gets Section name.
        ///@see Section::Name.
        std::string GetName();

    protected:
        ///Gets Local Axis Rotation for Line Section according to provided angle.
        ///@param theta The angle of rotation about x1-axis.
        ///@return The section rotation matrix in theta degrees.
        ///@note More details can be found at @ref linkSectionLocalAxis.
        Eigen::MatrixXd GetLineRotationMatrix(double theta);

        ///Gets Local Axis Rotation for Area Section according to provided angle.
        ///@param theta The angle of rotation about x3-axis.
        ///@return The section rotation matrix in theta degrees.
        ///@note More details can be found at @ref linkSectionLocalAxis.
        Eigen::MatrixXd GetAreaRotationMatrix(double Theta);

        ///Gets centroid translation axis for Line Section according to insertion point.
        ///@param h The total height of the section.
        ///@param b The total width of the section.
        ///@param zc Position of the section axis along x3-axis.
        ///@param yc Position of the section axis along x2-axis.
        ///@param ip The insertion point where the section axes are located.
        ///@return The section translation matrix according to ip.
        ///@note More details can be found at @ref linkSectionInsertionPoint.
        Eigen::MatrixXd GetLineTranslationMatrix(double h, double b, double zc, double yc, unsigned int ip);

        ///Gets the coordinate according to insertion point.
        ///@param x3 The coordinate along the 3-axis.
        ///@param x2 The coordinate along the 2-axis. 
        ///@param h The total height of the section.
        ///@param b The total width of the section.
        ///@param zc Position of the section axis along x3-axis.
        ///@param yc Position of the section axis along x2-axis.
        ///@param ip The insertion point where the section axes are located.
        ///@return The section translation matrix according to ip.
        ///@note More details can be found at @ref linkSectionInsertionPoint.
        void InsertionPointCoordinates(double &x3, double &x2, double h, double b, double zc, double yc, unsigned int ip);

        ///Transforms generalised strain/stresses from Element to Section local coordinate
        ///@param h The total height of the section.
        ///@param b The total width of the section.
        ///@param zcm Position of the section axis along x3-axis.
        ///@param ycm Position of the section axis along x2-axis.
        ///@param ip The insertion point where the section axes are located.
        ///@return The section translation matrix according to ip.
        ///@note More details can be found at @ref linkSectionInsertionPoint.
        Eigen::MatrixXd ComputeLineLocalAxes(double h, double b, double zcm, double ycm, double angle, unsigned int ip);

    private:
        ///Section Name.
        std::string Name;
};

#endif
