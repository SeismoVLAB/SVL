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
// Information:
//   The load formulation are given as,
//       [1] Static point-load.
//       [2] Dynamic point-load.
//       [3] Static surface-load.
//       [4] Static body-load.
//       [5] Dynamic body-load.
//       [6] Dynamic wave-front load.
//       [7] Dynamic user's define wave load.
//       [8] Static/Dynamic displacement support motion.
//
// Description:   
///This file contains the "Load" class declarations, which defines a load that 
///can act in a node or an element of a finite element mesh.
//------------------------------------------------------------------------------

#ifndef _LOAD_HPP_
#define _LOAD_HPP_

#include <map>
#include <vector>
#include <Eigen/Dense>

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      March 25, 2018
/// @version   1.0
/// @file      Load.hpp
/// @class     Load
/// @see       Node.hpp Element.hpp LoadCombo.hpp Mesh.hpp
/// @brief     Class for creating a load that can be applied to a node or to a element
class Load{

    public:
        ///Creates a Load in a finite element Mesh.
        ///@param type The load formulation.
        ///@note More details can be found at @ref linkLoad.
        ///@see Load::Classification
        Load(unsigned int type);

        ///Creates a Load in a finite element Mesh.
        ///@param dir The unit load direction.
        ///@param val The load magnitude in time.
        ///@param type The load formulation.
        ///@note More details can be found at @ref linkLoad.
        ///@see Load::Classification, Load::ForceDirection, Load::ForceAmplitude.
        Load(Eigen::VectorXd dir, std::vector<double> val, unsigned int type);

        ///Destroys this Load object.
        ~Load();

        ///Adds node that share this loaded.
        ///@param tags The list of Node identifier.
        ///@see Node, Load::Tags.
        void AddNodes(std::vector<unsigned int> tags);

        ///Adds face that share this loaded.
        ///@param tags The list of Element faces.
        ///@see Element::ComputeSurfaceForces(), Load::Faces.
        void AddFaces(std::vector<unsigned int> tags);

        ///Adds element that share this loaded.
        ///@param tags The list of Element identifiers.
        ///@see Element::ComputeBodyForces(), Element::ComputeSurfaceForces(), Load::Tags.
        void AddElements(std::vector<unsigned int> tags);

        ///Adds the exterior/interior condition for the domain reduction node.
        ///@param tag The Node identifier that belongs to domain reduction.
        ///@param cond If the Node identifier is interior (0) or exterior (1).
        ///@see Element::ComputeDomainReductionForces(), Load::DRMConditions.
        void AddDRMCondition(unsigned int tag, bool cond);

        ///Returns the load classification.
        ///@return The Load type.
        ///@note More details can be found at @ref linkLoad.
        ///@see Load::Classification.
        unsigned int GetClassification() const;

        ///Returns the nodes index that share this load.
        ///@return The Node list that this load is applied.
        ///@see Load::Tags.
        std::vector<unsigned int> GetNodes() const;

        ///Returns the faces index that share this load.
        ///@return The Element face list that this load is applied.
        ///@see Load::Faces.
        std::vector<unsigned int> GetFaces() const;

        ///Returns the nodes index that share this load.
        ///@return The Element list that this load is applied.
        ///@see Load::Tags.
        std::vector<unsigned int> GetElements() const;

        ///Returns the force vector to be applied.
        ///@param step The time step at which the load is evaluated.
        ///@return The force vector to be applied at step.
        ///@note More details can be found at @ref linkLoad.
        Eigen::VectorXd GetLoadVector(unsigned int step=0) const;

        ///Returns the exterior/interior condition for the domain reduction node.
        ///@param tag The domain reduction Node tag.
        ///@return whether the Node is interior or exterior.
        ///@note More details can be found at @ref linkLoad.
        ///@see Load::DRMConditions, Element::ComputeDomainReductionForces()
        bool GetDRMCondition(unsigned int tag);

    private:
        ///Load classification.
        unsigned int Classification;

        ///Direction of applied force.
        Eigen::VectorXd ForceDirection;

        ///Index of nodes/element which share this load.
        std::vector<unsigned int> Tags; 

        ///Index of element faces which share this load.
        std::vector<unsigned int> Faces;

        ///Time varing coefficients:
        std::vector<double> ForceAmplitude;

        ///Exterior/Interior condition for the domain reduction node.
        std::map<unsigned int, bool> DRMConditions;
};

#endif