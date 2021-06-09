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
//  [1] 2D Co-rotational Truss Formulation, Louie L. Yaw, E. F. Cross School of 
//      Engineering, Walla Walla University, April 23, 2009. 
//
// Description:
///This file contains the "kin2DTruss2" kinematic two-node element declarations, 
///which defines an element in a finite element mesh.
//------------------------------------------------------------------------------

#ifndef _KIN2DTRUSS2_HPP_
#define _KIN2DTRUSS2_HPP_

#include <map>
#include <memory>
#include <string>
#include <Eigen/Dense>

#include "Node.hpp"
#include "Load.hpp"
#include "Material.hpp"
#include "Element.hpp"
#include "Damping.hpp"

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      May 10, 2019
/// @version   1.0
/// @file      kin2DTruss2.hpp
/// @class     kin2DTruss2
/// @see       Element.hpp Mesh.hpp
/// @brief     Class for creating a 2D finite kinematic two-node truss element in a mesh
class kin2DTruss2 : public Element{

    public:
        ///Creates a kin2DTruss2 in a finite element Mesh.
        ///@param nodes The Node connectivity array of this Element.
        ///@param material Pointer to the Material that this Element is made out of.
        ///@param area The cross-section area of the kin2DTruss2 Element.
        ///@note More details can be found at @ref linkkin2DTruss2.
        ///@see kin2DTruss2::theNodes, kin2DTruss2::theMaterial.
        kin2DTruss2(const std::vector<unsigned int> nodes, std::unique_ptr<Material> &material, const double area);

        ///Destroys this kin3DTruss2 object.
        ~kin2DTruss2();

        ///Save the material states in the element.
        ///@note This function sets the trial states as converged ones in Material.
        void CommitState();

        ///Reverse the material/section states to previous converged state in this element.
        ///@note This function returns the trial states to previous converged states at the Material level.
        void ReverseState();

        ///Brings the material/section state to its initial state in this element.
        ///@note This function returns the meterial states to the beginning.
        void InitialState();

        ///Update the material states in the element.
        ///@note This function update the trial states at the Material level.
        void UpdateState();

        ///Sets the finite element dependance among objects.
        ///@param nodes The Node list of the Mesh object.
        ///@note This function sets the relation between Node and Element objects.
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
        ///@note The mass matrix can be revisited in @ref linkkin2DTruss2.
        ///@see Assembler::ComputeMassMatrix(), Integrator::ComputeEffectiveStiffness().
        Eigen::MatrixXd ComputeMassMatrix();

        ///Compute the stiffness matrix of the element using gauss-integration.
        ///@return Matrix with the Element stiffness matrix.
        ///@note The stiffness matrix can be revisited in @ref linkkin2DTruss2.
        ///@see Assembler::ComputeStiffnessMatrix(), Integrator::ComputeEffectiveStiffness().
        Eigen::MatrixXd ComputeStiffnessMatrix();

        ///Compute the damping matrix of the element using gauss-integration.
        ///@return Matrix with the Element damping matrix.
        ///@note The damping matrix can be revisited in @ref linkkin2DTruss2.
        ///@see Assembler::ComputeDampingMatrix(), Integrator::ComputeEffectiveStiffness().
        Eigen::MatrixXd ComputeDampingMatrix();

        ///Compute the PML history matrix using gauss-integration.
        ///@return Matrix with the PML Element matrix.
        ///@note The PML matrix is none existent for this element.
        ///@see Assembler::ComputePMLHistoryMatrix(), Integrator::ComputeEffectiveStiffness().
        Eigen::MatrixXd ComputePMLMatrix();

        ///Compute the internal (elastic) forces acting on the element.
        ///@return Vector with the Element internal force.
        ///@note The internal force vector can be revisited in @ref linkkin2DTruss2.
        ///@see Assembler::ComputeInternalForceVector(), Integrator::ComputeEffectiveForce().
        Eigen::VectorXd ComputeInternalForces();

        ///Compute the elastic, inertial, and viscous forces acting on the element.
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

    private:    
        ///Length of the element undeformed configuration.
        double Lo;

        ///Area of cross-section.
        double A;
    
        ///The Damping model.
        std::shared_ptr<Damping> theDamping;

        ///The Element's Nodes.
        std::vector<std::shared_ptr<Node> > theNodes;

        ///The Element's material.
        std::unique_ptr<Material> theMaterial;

        ///Compute the current length of the element.
        double ComputeLength() const;

        ///Compute the axial force in element.
        double ComputeAxialForce() const;

        ///Compute/update the local axis-1 of the element.
        ///@return Matrix with the local x1-axis transformation.
        ///@note Axes transformation are according to @ref linkkin2DTruss2.
        Eigen::VectorXd ComputeLocalAxes() const;

        ///Compute/update the geometric global axis of the element.
        ///@return Matrix with the geometric global axis of the element.
        ///@note Axes transformation are according to @ref linkkin2DTruss2.
        Eigen::MatrixXd ComputeGeometricAxes() const;

        ///Compute/update the tranformation matrix from local to global axis.
        ///@return Matrix with the local axes x1-x2-x3 transformation.
        ///@note Axes transformation are according to @ref linkkin2DTruss2.
        Eigen::MatrixXd ComputeTransformationAxes() const;

        ///Update strain in the element:
        ///@return Vector with the strain at Element center.
        Eigen::VectorXd ComputeStrain() const;

        ///Update strain rate in the element:
        ///@return Vector with the strain-rate at integration point.
        Eigen::VectorXd ComputeStrainRate() const;

        ///Compute the initial stiffness matrix of the element.
        ///@return Matrix with the initial Tanget stiffness matrix.
        ///@note This matrix is computed when the strain deformation is zero.
        Eigen::MatrixXd ComputeInitialStiffnessMatrix() const;
};

#endif
