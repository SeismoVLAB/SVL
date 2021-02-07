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
//  [1] Finite Element Procedures, Bathe, K.J., Chapter 5: pages 342-348. 
//      Prentice-Hall, 1996.
//
// Description:
///This file contains the "kin2DQuad4" kinematic four-node element declarations, 
///which defines an element in a finite element mesh.   
//------------------------------------------------------------------------------

#ifndef _KIN2DQUAD4_HPP_
#define _KIN2DQUAD4_HPP_

#include <map>
#include <memory>
#include <string>
#include <Eigen/Dense>

#include "Node.hpp"
#include "Load.hpp"
#include "Material.hpp"
#include "Element.hpp"
#include "QuadratureRule.hpp"
#include "Damping.hpp"

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      August 26, 2019
/// @version   1.0
/// @file      kin2DQuad4.hpp
/// @class     kin2DQuad4
/// @see       Element.hpp Mesh.hpp
/// @brief     Class for creating a 2D finite kinematic four-node quadrilateral element in a mesh
class kin2DQuad4 : public Element{

    public:
        ///Creates a kin2DQuad4 in a finite element Mesh.
        ///@param nodes The Node connectivity array of this Element.
        ///@param material Pointer to the Material that this Element is made out of.
        ///@param th The thickness of the lin2DQuad4 Element.
        ///@param quadrature The integration rule to be employed.
        ///@param nGauss Number of Gauss points for Element integration.
        ///@note More details can be found at @ref linkkin2DQuad4.
        ///@see kin2DQuad4::theNodes, kin2DQuad4::theMaterial, kin2DQuad4::QuadraturePoints.
        kin2DQuad4(const std::vector<unsigned int> nodes, std::unique_ptr<Material> &material, const double th, const std::string quadrature="GAUSS", const unsigned int nGauss=4);

        ///Destroys this kin2DQuad4 object.
        ~kin2DQuad4();

        ///Save the material states in the element.
        ///@note This funtion sets the trial states as converged ones in Material.
        void CommitState();

        ///Update the material states in the element.
        ///@note This funtion update the trial states at the Material level.
        void UpdateState();

        ///Sets the finite element dependance among objects.
        ///@param nodes The Node list of the Mesh object.
        ///@note This funtion sets the relation between Node and Element objects.
        ///@see lin2DQuad4::theNodes.
        void SetDomain(std::map<unsigned int, std::shared_ptr<Node> > &nodes);

        ///Sets the damping model.
        ///@param damping Pointer to the damping model.
        ///@note Several Element objects can share the same damping model.
        void SetDamping(const std::shared_ptr<Damping> &damping);

        ///Gets the list of total-degree of freedom of this Element.
        ///@return Vector with the list of degree-of-freedom of this Element.
        std::vector<unsigned int> GetTotalDegreeOfFreedom() const;

        ///Gets the material/section (generalised) strain.
        ///@return Matrix with the strain at each integration point.
        ///@note The index (i,j) are the strain and Gauss-point respectively. 
        Eigen::MatrixXd GetStrain() const;

        ///Gets the material/section (generalised) stress.
        ///@return Matrix with the stress at each integration point.
        ///@note The index (i,j) are the stress and Gauss-point respectively. 
        Eigen::MatrixXd GetStress() const;

        ///Gets the material/section (generalised) strain-rate.
        ///@return Matrix with the strain-rate at each integration point.
        ///@note The index (i,j) are the strain-rate and Gauss-point respectively.
        Eigen::MatrixXd GetStrainRate() const;

        ///Gets the material strain in section at  coordinate (x3,x2).
        ///@param x3 Local coordinate along the x3-axis.
        ///@param x2 Local coordinate along the x2-axis.
        ///@return Matrix with the strain at coordinate (x3,x2).
        ///@note The strains are interpolated at this coordinate.
        Eigen::MatrixXd GetStrainAt(double x3, double x2) const;

        ///Gets the material stress in section at  coordinate (x3,x2).
        ///@param x3 Local coordinate along the x3-axis.
        ///@param x2 Local coordinate along the x2-axis.
        ///@return Matrix with the stresses at coordinate (x3,x2).
        ///@note The stresses are interpolated at this coordinate.
        Eigen::MatrixXd GetStressAt(double x3, double x2) const;

        ///Gets the element internal response in VTK format for Paraview display.
        ///@param response The response to be display in Paraview.
        ///@return Vector with the response at the Element center.
        ///@note The current responses are: "Strain", "Stress".
        Eigen::VectorXd GetVTKResponse(std::string response) const;

        ///Computes the element energy for a given deformation.
        ///@return Scalar with the element deformation energy.
        double ComputeEnergy();

        ///Compute the lumped/consistent mass matrix of the element.
        ///@return Matrix with the Element mass matrix.
        ///@note The mass matrix can be revisited in @ref linkkin2DQuad4.
        ///@see Assembler::ComputeMassMatrix(), Integrator::ComputeEffectiveStiffness().
        Eigen::MatrixXd ComputeMassMatrix();

        ///Compute the stiffness matrix of the element using gauss-integration.
        ///@return Matrix with the Element stiffness matrix.
        ///@note The stiffness matrix can be revisited in @ref linkkin2DQuad4.
        ///@see Assembler::ComputeStiffnessMatrix(), Integrator::ComputeEffectiveStiffness().
        Eigen::MatrixXd ComputeStiffnessMatrix();

        ///Compute the damping matrix of the element using gauss-integration.
        ///@return Matrix with the Element damping matrix.
        ///@note The damping matrix can be revisited in @ref linkkin2DQuad4.
        ///@see Assembler::ComputeDampingMatrix(), Integrator::ComputeEffectiveStiffness().
        Eigen::MatrixXd ComputeDampingMatrix();

        ///Compute the PML history matrix using gauss-integration.
        ///@return Matrix with the PML Element matrix.
        ///@note The PML matrix is none existent for this element.
        ///@see Assembler::ComputePMLHistoryMatrix(), Integrator::ComputeEffectiveStiffness().
        Eigen::MatrixXd ComputePMLMatrix();

        ///Compute the internal (elastic) forces acting on the element.
        ///@return Vector with the Element internal force.
        ///@note The internal force vector can be revisited in @ref linkkin2DQuad4.
        ///@see Assembler::ComputeInternalForceVector(), Integrator::ComputeEffectiveForce().
        Eigen::VectorXd ComputeInternalForces();

        ///Compute the elastic, inertial, and vicous forces acting on the element.
        ///@return Vector with the Element dynamic internal force.
        ///@note The internal force vector can be revisited in @ref linkElement.
        ///@see Assembler::ComputeDynamicInternalForceVector().
        Eigen::VectorXd ComputeInternalDynamicForces();

        ///Compute the surface forces acting on the element.
        ///@param surface Pointer to the Load object that contains this information.
        ///@param k The time step at which the surface load is evaluated.
        ///@return Vector with the Element surface force.
        ///@note The surface force vector can be revisited in @ref linkZeroLength1D.
        ///@see Assembler::ComputeExternalForceVector(), Integrator::ComputeEffectiveForce().
        Eigen::VectorXd ComputeSurfaceForces(const std::shared_ptr<Load> &surface, unsigned int face);

        ///Compute the body forces acting on the element.
        ///@param body Pointer to the Load object that contains this information.
        ///@param k The time step at which the body load is evaluated.
        ///@return Vector with the Element surface force.
        ///@note The body force vector can be revisited in @ref linkZeroLength1D.
        ///@see Assembler::ComputeExternalForceVector(), Integrator::ComputeEffectiveForce().
        Eigen::VectorXd ComputeBodyForces(const std::shared_ptr<Load> &body, unsigned int k=0);

        ///Compute the domain reduction forces acting on the element.
        ///@param drm Pointer to the DRM Load object that contains this information.
        ///@param k The time step at which the body load is evaluated.
        ///@return Vector with the Element domain reduction forces.
        ///@note The DRM force vector can be revisited in @ref linkZeroLength1D.
        ///@see Assembler::ComputeExternalForceVector(), Integrator::ComputeEffectiveForce().
        Eigen::VectorXd ComputeDomainReductionForces(const std::shared_ptr<Load> &drm, unsigned int k);

    protected:
        ///Transform tensor components into vector.
        ///@param vector The tensor to be transformed into a vector.
        ///@return The equivalent vector object.
        Eigen::VectorXd TransformTensorToVector(const Eigen::MatrixXd &vector) const;

        ///Transform vector components into tensor.
        ///@param vector The vector to be transformed into a tensor.
        ///@return The equivalent tensor object.
        Eigen::MatrixXd TransformVectorToTensor(const Eigen::VectorXd &vector) const;

    private:
        ///Element thickness.
        double t;

        ///The Damping model.
        std::shared_ptr<Damping> theDamping;

        ///The Element's Nodes.
        std::vector<std::shared_ptr<Node> > theNodes;

        ///The Element's material.
        std::vector<std::unique_ptr<Material> > theMaterial;

        ///Coordinate of Gauss points.
        std::unique_ptr<QuadratureRule> QuadraturePoints;

        ///Update strain in the element.
        ///@param ri The i-th Gauss coordinate in the r-axis.
        ///@param si The i-th Gauss coordinate in the s-axis.
        ///@param Jij The Jacobian matrix evaluated at (ri,si).
        ///@return Vector with the strain at integration point.
        Eigen::VectorXd ComputeStrain(const double ri, const double si, const Eigen::MatrixXd &Jij) const;

        ///Update strain rate in the element.
        ///@param ri The i-th Gauss coordinate in the r-axis.
        ///@param si The i-th Gauss coordinate in the s-axis.
        ///@param Bij The strain-displacement matrix.
        ///@return Vector with the strain-rate at integration point.
        Eigen::VectorXd ComputeStrainRate(const double ri, const double si, const Eigen::MatrixXd &Bij) const;

        ///Computes the jacobian of the transformation. 
        ///@param ri The i-th Gauss coordinate in the r-axis.
        ///@param si The i-th Gauss coordinate in the s-axis.
        ///@return The Jacobian matrix evaluated at i-th Gauss point (ri,si).
        Eigen::MatrixXd ComputeJacobianMatrix(const double ri, const double si) const;

        ///Computes the deformation gradient. 
        ///@param ri The i-th Gauss coordinate in the r-axis.
        ///@param si The i-th Gauss coordinate in the s-axis.
        ///@param Jij The Jacobian matrix evaluated at (ri,si).
        ///@return Matrix with the deformation gradient.
        Eigen::MatrixXd ComputeDeformationGradientMatrix(const double ri, const double si, const Eigen::MatrixXd &Jij) const;

        ///Compute the Second Piola-Kirchhoff Stress Tensor.
        ///@param Stress The stress vector to be mapped to current configuration.
        Eigen::MatrixXd ComputeSecondPiolaKirchhoffMatrix(const Eigen::VectorXd &Stress) const;

        ///Evaluates the shape function matrix at a given Gauss point.
        ///@param ri The i-th Gauss coordinate in the r-axis.
        ///@param si The i-th Gauss coordinate in the s-axis.
        ///@return Matrix with the shape function at the Gauss coordinate.
        Eigen::MatrixXd ComputeShapeFunctionMatrix(const double ri, const double si) const;

        ///Evaluates the strain-displacement matrix at a given Gauss point.
        ///@param ri The i-th Gauss coordinate in the r-axis.
        ///@param si The i-th Gauss coordinate in the s-axis.
        ///@param Jij The Jacobian matrix evaluated at (ri,si).
        ///@return Matrix with the linear stran-displacement operator.
        Eigen::MatrixXd ComputeLinearStrainDisplacementMatrix(const double ri, const double si, const Eigen::MatrixXd &Jij) const;

        ///Evaluates the strain-displacement matrix at a given Gauss point.
        ///@param ri The i-th Gauss coordinate in the r-axis.
        ///@param si The i-th Gauss coordinate in the s-axis.
        ///@param Jij The Jacobian matrix evaluated at (ri,si).
        ///@return Matrix with the geometric stran-displacement operator.
        Eigen::MatrixXd ComputeNonLinearStrainDisplacementMatrix(const double ri, const double si, const Eigen::MatrixXd &Jij) const;

        ///Evaluates the linearized (small strain) strain-displacement matrix at a given Gauss point.
        ///@param ri The i-th Gauss coordinate in the r-axis.
        ///@param si The i-th Gauss coordinate in the s-axis.
        ///@param Jij The Jacobian matrix evaluated at (ri,si).
        ///@return Matrix with the linear stran-displacement operator.
        Eigen::MatrixXd ComputeLinearizedStrainDisplacementMatrix(const double ri, const double si, const Eigen::MatrixXd &Jij) const;

        ///Compute the initial stiffness matrix of the element using gauss-integration.
        ///@return Matrix with the initial Tanget stiffness matrix.
        ///@note This matrix is computed when the strain deformation is zero.
        Eigen::MatrixXd ComputeInitialStiffnessMatrix() const;
};

#endif
