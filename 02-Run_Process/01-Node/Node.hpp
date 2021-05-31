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
///This file contains the "Node object" declarations, which is stores the 
///coordinates, state variables, degrees-of-freedom, and total-degree of 
///freedom lists of a node in a finite element mesh.
//------------------------------------------------------------------------------

#ifndef _NODE_HPP_
#define _NODE_HPP_

#include <map>
#include <vector>
#include <Eigen/Dense>

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      January 15, 2018
/// @version   1.0
/// @file      Node.hpp
/// @class     Node
/// @see       Mesh.hpp
/// @brief     Class for creating a node in space that handles coordinates, state variables, free and total degrees-of-freedom numbering
class Node{

    public:
        ///Creates a Node in a finite element Mesh.
        ///@param numberDofs Number of degree-of-freedom for this node.
        ///@param coordinates Matrix that specifies the coordinates of this node.
        ///@param isFixed Whether the node is fixed or free.
        ///@param isAbsorbent Whether the node is belongs is absorbent (PML 3D).
        ///@note More details can be found at @ref linkNode.
        ///@see Node::Coordinates Node::NumberOfDegreeOfFreedom.
        Node(unsigned int numberDofs, Eigen::VectorXd coordinates, bool isFixed);

        ///Destroys this Node object.
        ~Node();

        ///Whether the node is fixed or free.
        ///@see Node::Fixed.
        bool IsFixed() const;

        ///Set all node variables to be zero
        void InitialState();

        ///Set the node's mass at each degree of freedom.
        ///@param mass Vector of masses.
        ///@return Nothing.
        ///@see Node::Mass.
        void SetMass(Eigen::VectorXd &mass);

        ///Sets the current reaction force of this node.
        ///@param mass Vector of masses.
        ///@return Nothing.
        ///@see Node::Reaction.
        void SetReaction(Eigen::VectorXd &reaction);

        ///Set the node's position in a finite element Mesh.
        ///@param coordinates Vector of coordinates or position.
        ///@note More details can be found at @ref linkNode.
        ///@see Node::Coordinates.
        void SetCoordinates(Eigen::VectorXd &coordinates);

        ///Set the free degree of freedom number list for this Node.
        ///@param freeDOF Vector of integers with the 'free' degree-of-freedom numbering.
        ///@note More details can be found at @ref linkNode.
        ///@see Node::FreeDegreeOfFreedom.
        void SetFreeDegreeOfFreedom(std::vector<int> &freeDOF);

        ///Set the global degree of freedom number list for this Node.
        ///@param totalDOF Vector of integers with the 'total' degree-of-freedom numbering.
        ///@note More details can be found at @ref linkNode.
        ///@see Node::TotalDegreeOfFreedom.
        void SetTotalDegreeOfFreedom(std::vector<int> &totalDOF);

        ///Sets the vector of displacement associated with this Node.
        ///@param Uo Vector of displacement.
        ///@see Node::Displacements.
        void SetDisplacements(Eigen::VectorXd &Uo);

        ///Sets the current velocity state associated with this Node.
        ///@param Vo Vector of velocities.
        ///@see Node::Velocities.
        void SetVelocities(Eigen::VectorXd &Vo);

        ///Sets the current acceleration state associated with this Node.
        ///@param Ao Vector of accelertions.
        ///@see Node::Accelerations.
        void SetAccelerations(Eigen::VectorXd &Ao);

        ///Sets the current incremental displacement associated with this Node.
        ///@param dU Vector of incremental displacements.
        ///@return Nothing.
        ///@see Node::IncrementalDisplacements.
        void SetIncrementalDisplacements(Eigen::VectorXd &dU);

        ///Sets the domain-reduction nodal displacements.
        ///@param Uo Matrix of displacements, velocities, and accelerations.
        ///@see Node::DomainReductionMotion.
        void SetDomainReductionMotion(Eigen::MatrixXd &Uo);

        ///Sets the vector of absorbent integrated history values.
        ///@param Ub Vector of integrated (states) history values.
        ///@see Node::PMLIntegratedVector Node::GetPMLVector.
        void SetPMLVector(Eigen::VectorXd &Ub);

        ///Sets the displacements support motion associated with this Node.
        ///@param k Integer that represents the degree-of-freedom where the displacement is imposed.
        ///@param Uo Vector of displacements in time.
        ///@see Node::SupportMotion.
        void SetSupportMotion(unsigned int k, std::vector<double> &Uo);

        ///Gets the node's mass as a vector.
        ///@return Vector with the mass at each degree-of-freedom.
        ///@note More details can be found at @ref linkNode, and @ref linkMass.
        ///@see Node::Mass.
        Eigen::VectorXd GetMass() const;

        ///Gets the node reaction force
        ///@return Vector with the reaction at each degree-of-freedom.
        ///@see Node::Reaction.
        Eigen::VectorXd GetReaction() const;

        ///Gets the node's coordinate as a vector.
        ///@return Vector with the coordinates or position of this node.
        ///@note More details can be found at @ref linkNode.
        ///@see Node::Coordinates.
        Eigen::VectorXd GetCoordinates() const;

        ///Gets the current displacement state of this node.
        ///@return Vector with the displacements for each degree-of-freedom.
        ///@see Node::Displacements.
        Eigen::VectorXd GetDisplacements() const;

        ///Gets the current velocity state of this node.
        ///@return Vector with the velocities for each degree-of-freedom.
        ///@see Node::Velocities.
        Eigen::VectorXd GetVelocities() const;

        ///Gets the current acceleration state of this node.
        ///@return Vector with the accelerations for each degree-of-freedom.
        ///@see Node::Accelerations.
        Eigen::VectorXd GetAccelerations() const;

        ///Returns the inertial forces associated to this node.
        ///@return Vector with the inertial forces for each degree-of-freedom.
        ///@see Node::Mass.
        Eigen::VectorXd GetInertialForces() const;

        ///Gets the current PML history vector of this node.
        ///@return Vector with the integrated PML states at each degree-of-freedom.
        ///@see Node::PMLIntegratedVector Node::SetPMLVector.
        Eigen::VectorXd GetPMLVector() const;

        ///Gets the current incremental displacement of this node.
        ///@return Vector with the incremental displacements for each degree-of-freedom.
        ///@see Node::IncrementalDisplacements.
        Eigen::VectorXd GetIncrementalDisplacements() const;

        ///Gets the current domain reduction displacements, velocities and accelerations.
        ///@param k Time step at which the domain reduction information is retrieved.
        ///@return Vector with the domain reduction information for computing forces.
        ///@note More details can be found at @ref linkNode, and @ref linkLoad.
        ///@see Node::DomainReductionMotion.
        Eigen::VectorXd GetDomainReductionMotion(unsigned int k) const;

        ///Gets the nodal support motion displacement of this node.
        ///@param k Degree-of-freedom at which the displacement supports are retrieved.
        ///@return Vector with the associated support displacements.
        ///@see Node::SupportMotion.
        Eigen::VectorXd GetSupportMotion(unsigned int k);

        ///Gets the number of degree-of-freedom of this node.
        ///@return Integer that represents the number of degree-of-freedom.
        ///@see Node::TotalDegreeOfFreedom.
        unsigned int GetNumberOfDegreeOfFreedom() const;

        ///Gets the node's restrained degree of freedom as a vector.
        ///@return Array of integers with the free degree-of-freedom numbering.
        ///@note More details can be found at @ref linkNode.
        ///@see Node::FreeDegreeOfFreedom.
        std::vector<int> GetFreeDegreeOfFreedom() const;

        ///Gets the node's total degree of freedom as a vector.
        ///@return Array of integers with the total degree-of-freedom numbering.
        ///@note More details can be found at @ref linkNode.
        ///@see Node::TotalDegreeOfFreedom.
        std::vector<int> GetTotalDegreeOfFreedom() const;

    private:
        //The node condition.
        bool Fixed;

        ///Number of Degree of Freedom:
        unsigned int NumberOfDegreeOfFreedom;

        ///The location of this node within the Mesh..
        Eigen::VectorXd Coordinates;

        ///The current displacement vector of this Node.
        Eigen::VectorXd Displacements;

        ///The current velocity vector of this Node.
        Eigen::VectorXd Velocities;

        ///The current acceleration vectors of this Node.
        Eigen::VectorXd Accelerations;

        ///The current incremental displacement vector of Node.
        Eigen::VectorXd IncrementalDisplacements;

        ///The imposed domain-reduction nodal displacements, velocity and acceleration.
        Eigen::MatrixXd DomainReductionMotion;

        ///The applied mass of this Node.
        Eigen::VectorXd Mass;

        ///The reaction forces of this node.
        Eigen::VectorXd Reaction;

        ///The 3D perfectly-matched history verctor.
        Eigen::VectorXd PMLIntegratedVector;

        ///The free-degree-of-freedom list.
        std::vector<int> FreeDegreeOfFreedom;

        ///The total-degree-of-freedom list.
        std::vector<int> TotalDegreeOfFreedom;

        ///The imposed nodal support motion displacements.
        std::map<unsigned int, std::vector<double> > SupportMotion;
};

#endif