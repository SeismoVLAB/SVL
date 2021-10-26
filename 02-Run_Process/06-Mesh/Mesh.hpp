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
///This file contains the "Mesh object" declarations, which stores nodes, 
///materials, elements, and loads in a finite element model, and also creates 
///the links between the formers.
//------------------------------------------------------------------------------

#ifndef _MESH_HPP_
#define _MESH_HPP_

#include <map>
#include <vector>
#include <memory>
#include <Eigen/SparseCore>

#include "Node.hpp"
#include "Load.hpp"
#include "Section.hpp"
#include "Material.hpp"
#include "Damping.hpp"
#include "Element.hpp"
#include "Constraint.hpp"

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      April 25, 2018
/// @version   1.0
/// @file      Mesh.hpp
/// @class     Mesh
/// @see       Node.hpp Element.hpp Material.hpp Section.hpp Load.hpp
/// @brief     Class for defining a mesh that stores the nodes, materials, elements, and loads
class Mesh{

    public:
        ///Creates a Mesh object (container) to store Node, Material, Element, Load.
        ///@note More details can be found at @ref linkMesh.
        Mesh(); 

        ///Destroys this Mesh object.
        ~Mesh();

        ///Initialize the Mesh and compute memory storage.
        void Initialize();

        ///Add a Node object to the Mesh.
        ///@param tag the Node identifier.
        ///@param node the Node object to be added.
        ///@see Node.
        void AddNode(unsigned int tag, std::shared_ptr<Node> &node);

        ///Add a Constraint to the Mesh.
        ///@param tag the Constraint identifier.
        ///@param constraint the Constraint object to be added.
        ///@see Constraint.
        void AddConstraint(unsigned int tag, std::shared_ptr<Constraint> &constraint);

        ///Add a Material to the Mesh.
        ///@param tag the Material identifier.
        ///@param fiber the Material object to be added in fiber.
        ///@see Material.
        void AddFiber(unsigned int tag, std::unique_ptr<Material> &fiber);

        ///Add a Material to the Mesh.
        ///@param tag the Material identifier.
        ///@param material the Material object to be added.
        ///@see Material.
        void AddMaterial(unsigned int tag, std::unique_ptr<Material> &material);

        ///Add a Section to the Mesh.
        ///@param tag the Section identifier.
        ///@param section the Section object to be added.
        ///@see Section.
        void AddSection(unsigned int tag, std::unique_ptr<Section> &section);
    
        ///Add a Damping model to the Mesh.
        ///@param tag the Damping identifier.
        ///@param damping the Damping object to be added.
        ///@see Damping.
        void AddDamping(unsigned int tag, std::shared_ptr<Damping> &damping);

        ///Add a Element to the Mesh.
        ///@param tag the Element identifier.
        ///@param element the Element object to be added.
        ///@see Element.
        void AddElement(unsigned int tag, std::shared_ptr<Element> &element);

        ///Add a Load to the Mesh.
        ///@param tag the Load identifier.
        ///@param load the Load object to be added.
        ///@see Load.
        void AddLoad(unsigned int tag, std::shared_ptr<Load> &load);

        ///Add a point mass to a Node.
        ///@param tag the Node identifier.
        ///@param mass the vector with masses at each degree-of-freedom.
        ///@see Node::SetMass(), Node::GetMass(), Node::Mass.
        void AddMass(unsigned int tag, Eigen::VectorXd& mass);

        ///Remove a Node from Mesh.
        ///@param tag the Node identifier to be removed.
        ///@see Node.
        void DelNode(unsigned int tag);

        ///Remove a Point Mass from Mesh.
        ///@param tag the Node identifier where the mass is removed.
        ///@see Node.
        void DelMass(unsigned int tag);

        ///Remove the support motions associated to Node in from Mesh.
        ///@param tag the Node identifier where the mass is removed.
        ///@see Node.
        void DelSupportMotion(unsigned int tag);

        ///Remove a Constraint from Mesh.
        ///@param tag the Constraint identifier to be removed.
        ///@see Constraint.
        void DelConstraint(int tag);

        ///Remove a Material from Mesh.
        ///@param tag the Material identifier to be removed.
        ///@see Material.
        void DelMaterial(unsigned int tag);

        ///Remove a Section from Mesh.
        ///@param tag the Section identifier to be removed.
        ///@see Section.
        void DelSection(unsigned int tag);

        ///Remove an Element from Mesh.
        ///@param tag the Element identifier to be removed.
        ///@see Element.
        void DelElement(unsigned int tag);

        //Remove a Damping from Mesh.
        ///@param tag the Damping identifier to be removed.
        ///@see Damping.
        void DelDamping(unsigned int tag);

        ///Remove a Load from Mesh.
        ///@param tag the Load identifier to be removed.
        ///@see Load.
        void DelLoad(unsigned int tag);

        ///Add a Damping model to the Mesh.
        ///@param tag the Damping identifier.
        ///@param group the Element object that shares this damping model.
        ///@see Damping.
        void SetDamping(unsigned int tag, std::vector<unsigned int> &group);

        ///Specifies the Node initial condition.
        ///@param tag The Node identifier.
        ///@param cond Option to identify displacement, velocity or acceleration.
        ///@param Xo The vector of initial conditions.
        ///@see Node::SetDisplacements, Node::SetVelocities, Node::SetAccelerations.
        void SetInitialCondition(unsigned int tag, int cond, Eigen::VectorXd& Xo);

        ///Specifies the support motion for a certain Node object.
        ///@param tag The Node identifier.
        ///@param dof The degree-of-freedom where the support motion is imposed.
        ///@param Xo The vector of support motion displacement.
        ///@see Node::SetSupportMotion, Node::GetSupportMotion, Node::SupportMotion.
        void SetSupportMotion(unsigned int tag, unsigned int dof, std::vector<double>& Xo);

        ///Gets a material from the mesh.
        ///@param tag The material tag to be obtained.
        ///@return A Material pointer.
        ///@note More details can be found at @ref linkMesh.
        ///@see Meterial.
        std::unique_ptr<Material>& GetMaterial(unsigned int tag);

        ///Gets a section from the mesh.
        ///@param tag The section tag to be obtained.
        ///@return A Section pointer.
        ///@note More details can be found at @ref linkMesh.
        ///@see Section.
        std::unique_ptr<Section>& GetSection(unsigned int tag);

        ///Gets damping from the mesh.
        ///@param tag The damping tag to be obtained.
        ///@return A Damping pointer.
        ///@note More details can be found at @ref linkMesh.
        ///@see Damping.
        std::shared_ptr<Damping>& GetDamping(unsigned int tag);

        ///Gets nodes from the mesh.
        ///@return A map with all Node in the Mesh object.
        ///@note More details can be found at @ref linkMesh.
        ///@see Node.
        std::map<unsigned int, std::shared_ptr<Node> >& GetNodes();

        ///Gets constraint from the mesh.
        ///@return A map with all Constraint in the Mesh object.
        ///@note More details can be found at @ref linkMesh.
        ///@see Constraint.
        std::map<int, std::shared_ptr<Constraint> >& GetConstraints();
    
        ///Gets elements from the mesh.
        ///@return A map with all Element in the Mesh object.
        ///@note More details can be found at @ref linkMesh.
        ///@see Element.
        std::map<unsigned int, std::shared_ptr<Element> >& GetElements();

        ///Gets all damping from the mesh.
        ///@return A map with all Damping in the Mesh object.
        ///@note More details can be found at @ref linkMesh.
        ///@see Damping.
        std::map<unsigned int, std::shared_ptr<Damping> >& GetDampings();

        ///Gets loads from the mesh.
        ///@return A map with all Load in the Mesh object.
        ///@note More details can be found at @ref linkMesh.
        ///@see Load.
        std::map<unsigned int, std::shared_ptr<Load> >& GetLoads();

        ///Returns a vector with identifier of the property requested
        ///@param Name The entity to get the identifiers or tags.
        ///@return A vector with the identifiers
        template<typename T> std::vector<T> GetVectorIDs(std::string Name);

        ///Gets operator that impose restrains/constraints on the model. 
        ///@return A matrix that enforce the linear kinematic constraints.
        ///@note More details can be found at @ref linkMesh, and @ref linkAssembler.
        ///@see Constraint, Integrator::ComputeEffectiveForce, Integrator::ComputeEffectiveStiffness.
        Eigen::SparseMatrix<double> GetTotalToFreeMatrix();

    protected:
        ///Container/map of nodes for this partition.
        std::map<unsigned int, std::shared_ptr<Node> > Nodes;

        ///Container/map of constraints for this partition.
        std::map<int, std::shared_ptr<Constraint> > Constraints;

        ///Container/map of materials for this partition.
        std::map<unsigned int, std::unique_ptr<Material> > Materials;

        ///Container/map of materials for this partition.
        std::map<unsigned int, std::unique_ptr<Section> > Sections;
    
        ///Container/map of damping for this partition.
        std::map<unsigned int, std::shared_ptr<Damping> > Dampings;

        ///Container/map of elements for this partition.
        std::map<unsigned int, std::shared_ptr<Element> > Elements;

        ///Container/map of loads for this partition.
        std::map<unsigned int, std::shared_ptr<Load> > Loads;
};

#endif