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
///This file contains the "Fib3DLineSection" section declarations, which defines 
///a Rectangular geometry for 3D analysis with linear elastic material for 
///frame elements.
//------------------------------------------------------------------------------

#ifndef _FIB3DLINESECTION_HPP_
#define _FIB3DLINESECTION_HPP_

#include <string>
#include <memory>
#include <vector>
#include <Eigen/Dense>

#include "Material.hpp"
#include "Section.hpp"

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      March 24, 2021
/// @version   1.0
/// @file      Fib3DLineSection.hpp
/// @class     Fib3DLineSection
/// @see       Section.hpp
/// @brief     Class for creating a generic fiber section for 3D analysis with linear/nonlinear material for a frame elements
class Fib3DLineSection : public Section{

    public:
        ///Creates a Section to be specified at a Gauss-point in an Element.
        ///@param fibers A vector of fibers for each discretized area. 
        ///@param zi A vector with horizontal coordinate (x3) of each fiber.
        ///@param yi A vector with vertical coordinate (x2) of each fiber.
        ///@param Ai A vector with the area values of each fiber.
        ///@param k2 Shear area factor along the x2-axis.
        ///@param k3 Shear area factor along the x2-axis.
        ///@param theta The angle of rotation about x3-axis.
        ///@param ip The insertion point where the section axes are located.
        ///@note Details about insertion point and local axis rotation can be found at @ref linkSectionInsertionPoint and @ref linkSectionLocalAxis.
        ///@see Fib3DLineSection::h, Fib3DLineSection::b, Fib3DLineSection::theMaterial, Fib3DLineSection::Theta, Fib3DLineSection::InsertPoint.
        Fib3DLineSection(double h, double b, const std::vector<std::unique_ptr<Material> > &fibers, const std::vector<double> &zi, const std::vector<double> &yi, const std::vector<double> &Ai, double k2=1.0, double k3=1.0, double theta=0.0, unsigned int ip=10);

        ///Destroys this Fib3DLineSection object.
        ~Fib3DLineSection();

        ///Clone the 'Fib3DLineSection' section.
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
        ///@note This function returns the fiber states to previous converged states.
        void ReverseState();

        ///Brings the section states to its initial state.
        ///@note This function returns the fiber states to the beginning.
        void InitialState();

        ///Update the section state for this iteration.
        ///@param strain Vector with the strain components at this Gauss-point.
        ///@param cond If the the elastic/plastic fiber components will be updated.
        ///@note This function computes the strain and tanget stiffness matrix once the trial strain converged.
        void UpdateState(const Eigen::VectorXd strain, const unsigned int cond=0);

    protected:
        ///Gets the section centroid.
        ///@see The geometric properties can be revisited at @ref linkFib3DLineSection.
        void ComputeSectionCenter();

        ///Return the fiber strain closer to (x3, x2)
        ///@param x3 Position on the x3-axis where the stress is evaluated.
        ///@param x2 Position on the x2-axis where the stress is evaluated.
        ///@return Index of the fiber closer to coordinate (x3,x2).
        unsigned int FiberIndexAt(double x3=0, double x2=0);

    private:
        ///Height.
        double h;

        ///Width.
        double b;

        ///Rotation angle.
        double zcm;

        ///Rotation angle.
        double ycm;

        ///Shear area factor (As2/A) along the x2-axis.
        double kappa2;

        ///Shear area factor (As3/A) along the x3-axis.
        double kappa3;

        ///Rotation angle.
        double Theta;

        ///Local axis coordinate.
        unsigned int InsertPoint;

        ///Local axis coordinate.
        unsigned int NumberOfFibers;

        ///Generalized Strain vector.
        Eigen::VectorXd Strain;

        ///The horizontal (x3) local coordinates of each fiber.
        std::vector<double> zi;

        ///The vertical (x2) local coordinates of each fiber.
        std::vector<double> yi;

        ///The discretized area value of each fiber.
        std::vector<double> Ai;

        ///The section's fiber list.
        std::vector<std::unique_ptr<Material> > theFiber;
};

#endif
