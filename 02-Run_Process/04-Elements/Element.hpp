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
///This file contains the "Element" object declarations, which defines an 
///element in a finite element mesh.
//------------------------------------------------------------------------------

#ifndef _ELEMENT_HPP_
#define _ELEMENT_HPP_

#include <map>
#include <vector>
#include <memory>
#include <string>
#include <Eigen/Dense>

#include "Node.hpp"
#include "Load.hpp"
#include "Damping.hpp"

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      January 7, 2018
/// @version   1.0
/// @file      Element.hpp
/// @class     Element
/// @see       Mesh.hpp
/// @brief     Virtual class for creating an element in a mesh
class Element{

    public:
        ///Creates an Element in a finite element Mesh.
        ///@param name The name of this Element object.
        ///@param nodes The Node connectivity array.
        ///@param ndofs The number of degree-of-freedom of the Element object.
        ///@param VTKcell The VTK cell number for Paraview display.
        ///@note More details can be found at @ref linkElement.
        ///@see Element::Name, Element::VTKCell, Element::NumberOfNodes, Element::Nodes.
        Element(std::string name, const std::vector<unsigned int> nodes, unsigned int ndofs, unsigned int VTKcell);

        ///Destroys this Element object.
        virtual ~Element() = 0;

        ///Save the material/section states in the element.
        ///@note This funtion sets the trial states as converged ones in Material/Section.
        virtual void CommitState() = 0;

        ///Update the material/section states in the element.
        ///@note This funtion update the trial states at the Material/Section level.
        virtual void UpdateState() = 0;

        ///Reverse the material/section states to previous converged state in this element.
        ///@note This funtion returns the trial states to previous converged states at the Material/Section level.
        virtual void ReverseState() = 0;

        ///Brings the material/section state to its initial state in this element.
        ///@note This funtion returns the meterial states to the beginning.
        virtual void InitialState() = 0;

        ///Sets the finite element dependance among objects.
        ///@param nodes The Node list of the Mesh object.
        ///@note This funtion sets the relation between Node and Element objects.
        ///@see Element::Nodes.
        virtual void SetDomain(std::map<unsigned int, std::shared_ptr<Node> > &nodes) = 0;

        ///Sets the damping model.
        ///@param damping Pointer to the damping model.
        ///@note Several Element objects can share the same damping model.
        virtual void SetDamping(const std::shared_ptr<Damping> &damping) = 0;

        ///Gets the list of total-degree of freedom of this Element.
        ///@return Vector with the list of degree-of-freedom of this Element.
        virtual std::vector<unsigned int> GetTotalDegreeOfFreedom() const = 0; 

        ///Gets the material/section (generalised) strain.
        ///@return Matrix with the strain at each integration point.
        ///@note The index (i,j) are the strain and Gauss-point respectively. 
        virtual Eigen::MatrixXd GetStrain() const = 0;

        ///Gets the material/section (generalised) stress.
        ///@return Matrix with the stress at each integration point.
        ///@note The index (i,j) are the stress and Gauss-point respectively. 
        virtual Eigen::MatrixXd GetStress() const = 0;

        ///Gets the material/section (generalised) strain-rate.
        ///@return Matrix with the strain-rate at each integration point.
        ///@note The index (i,j) are the strain-rate and Gauss-point respectively.
        virtual Eigen::MatrixXd GetStrainRate() const = 0;

        ///Gets the material strain in section at  coordinate (x3,x2).
        ///@param x3 Local coordinate along the x3-axis.
        ///@param x2 Local coordinate along the x2-axis.
        ///@return Matrix with the strain at coordinate (x3,x2).
        ///@note The strains are interpolated at this coordinate.
        virtual Eigen::MatrixXd GetStrainAt(double x3=0.0, double x2=0.0) const = 0;

        ///Gets the material stress in section at  coordinate (x3,x2).
        ///@param x3 Local coordinate along the x3-axis.
        ///@param x2 Local coordinate along the x2-axis.
        ///@return Matrix with the stresses at coordinate (x3,x2).
        ///@note The stresses are interpolated at this coordinate.
        virtual Eigen::MatrixXd GetStressAt(double x3=0.0, double x2=0.0) const = 0;

        ///Gets the element internal response in VTK format for Paraview display.
        ///@param response The response to be display in Paraview.
        ///@return Vector with the response at the Element center.
        ///@note The current responses are: "Strain", "Stress".
        virtual Eigen::VectorXd GetVTKResponse(std::string response) const = 0;

        ///Computes the element energy for a given deformation.
        ///@return Scalar with the element deformation energy.
        virtual double ComputeEnergy() = 0;

        ///Compute the lumped/consistent mass matrix of the element or close form solutions.
        ///@return Matrix with the Element mass matrix.
        ///@note The mass matrix can be revisited in @ref linkElement.
        ///@see Assembler::ComputeMassMatrix(), Integrator::ComputeEffectiveStiffness().
        virtual Eigen::MatrixXd ComputeMassMatrix() = 0;

        ///Compute the stiffness matrix of the element using gauss-integration or close form solutions.
        ///@return Matrix with the Element stiffness matrix.
        ///@note The stiffness matrix can be revisited in @ref linkElement.
        ///@see Assembler::ComputeStiffnessMatrix(), Integrator::ComputeEffectiveStiffness().
        virtual Eigen::MatrixXd ComputeStiffnessMatrix() = 0;

        ///Compute the damping matrix of the element using gauss-integration or close form solutions.
        ///@return Matrix with the Element damping matrix.
        ///@note The damping matrix can be revisited in @ref linkElement.
        ///@see Assembler::ComputeDampingMatrix(), Integrator::ComputeEffectiveStiffness().
        virtual Eigen::MatrixXd ComputeDampingMatrix() = 0;

        ///Compute the PML history matrix for Perfectly-Matched Layer (PML).
        ///@return Matrix with the PML history variables.
        ///@note The history matrix only applies for PML in 3D, see @ref linkPML3DHexa8 linkPML3DHexa20.
        ///@see Assembler::ComputeDampingMatrix(), Integrator::ComputeEffectiveStiffness().
        virtual Eigen::MatrixXd ComputePMLMatrix() = 0;

        ///Compute the internal forces acting on the element.
        ///@return Vector with the Element internal force.
        ///@note The internal force vector can be revisited in @ref linkElement.
        ///@see Assembler::ComputeInternalForceVector(), Integrator::ComputeEffectiveForce().
        virtual Eigen::VectorXd ComputeInternalForces() = 0;

        ///Compute the elastic, inertial, and viscous forces acting on the element.
        ///@return Vector with the Element dynamic internal force.
        ///@note The internal force vector can be revisited in @ref linkElement.
        ///@see Assembler::ComputeDynamicInternalForceVector().
        virtual Eigen::VectorXd ComputeInternalDynamicForces() = 0;

        ///Compute the surface forces acting on the element.
        ///@param surface Pointer to the Load object that contains this information.
        ///@param k The time step at which the surface load is evaluated.
        ///@return Vector with the Element surface force.
        ///@note The surface force vector can be revisited in @ref linkElement.
        ///@see Assembler::ComputeExternalForceVector(), Integrator::ComputeEffectiveForce().
        virtual Eigen::VectorXd ComputeSurfaceForces(const std::shared_ptr<Load> &surface, unsigned int k) = 0;

        ///Compute the body forces acting on the element.
        ///@param body Pointer to the Load object that contains this information.
        ///@param k The time step at which the body load is evaluated.
        ///@return Vector with the Element surface force.
        ///@note The body force vector can be revisited in @ref linkElement.
        ///@see Assembler::ComputeExternalForceVector(), Integrator::ComputeEffectiveForce().
        virtual Eigen::VectorXd ComputeBodyForces(const std::shared_ptr<Load> &body, unsigned int k=0) = 0;

        ///Compute the domain reduction forces acting on the element.
        ///@param drm Pointer to the DRM Load object that contains this information.
        ///@param k The time step at which the body load is evaluated.
        ///@return Vector with the Element domain reduction forces.
        ///@note The DRM force vector can be revisited in @ref linkElement.
        ///@see Assembler::ComputeExternalForceVector(), Integrator::ComputeEffectiveForce().
        virtual Eigen::VectorXd ComputeDomainReductionForces(const std::shared_ptr<Load> &drm, unsigned int k) = 0;

        ///Gets the Element Name.
        ///@return String with the Element name.
        ///@see Element::Name.
        std::string GetName() const;

        ///Gets the Element VTK cell type.
        ///@return Value with the VTK cell type for Paraview.
        ///@see Element::VTKCell.
        unsigned int GetVTKCellType() const;

        ///Returns the number of nodes in element.
        ///@return The number of Node in this Element.
        ///@see Element::NumberOfNodes.
        unsigned int GetNumberOfNodes() const;

        ///Returns total number of degree of freedom in the element.
        ///@return The number of degree-of-freedom in this Element.
        ///@see Element::NumberOfDegreeOfFreedom.
        unsigned int GetNumberOfDegreeOfFreedom() const;

        ///Returns the Node Connectivity Indexes.
        ///@return The Node connectivity array in this Element.
        ///@see Element::Nodes.
        std::vector<unsigned int> GetNodes() const;

        ///Returns if the element has fixed nodes.
        ///@see Node::IsFixed.
        bool HasFixedNode(const std::vector<std::shared_ptr<Node> > &nodes) const;

    private:
        ///The element name.
        std::string Name;

        ///The VTK element number.
        unsigned int VTKCell;

        ///The number of nodes in element.
        unsigned int NumberOfNodes;

        ///The number of nodes in element.
        unsigned int NumberOfDegreeOfFreedom;

        ///The node connectivity indexes.
        std::vector<unsigned int> Nodes;
};

#endif
