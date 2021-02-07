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
// References   : 
//  [1] 
//
// Description:
///This file contains the "UnxBoucWen3DLink" two-node link declarations, which 
///defines a uniaxial Bouc-Wen model for 3D analyses.
//------------------------------------------------------------------------------

#ifndef _UNXBOUCWEN3DLINK_HPP_
#define _UNXBOUCWEN3DLINK_HPP_

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
/// @date      October 6, 2020
/// @version   1.0
/// @file      UnxBoucWen3DLink.hpp
/// @class     UnxBoucWen3DLink
/// @see       Element.hpp Mesh.hpp
/// @brief     Class for creating an uniaxial 3D two-node link element using Bouc-Wen formulation
class UnxBoucWen3DLink : public Element{

    public:
        ///Creates a UnxBoucWen3DLink in a finite element Mesh.
        ///@param nodes The Node connectivity array of this Element.
        ///@param parameters Vector that contains the Bouc-Wen model parameters.
        ///@param variables Vector that contains auxiliar model parameters.
        ///@param dim The model dimension.
        ///@param dir The local direction where this Element is acting.
        ///@param tol Tolerance to stop Newton-Raphson iterations.
        ///@param nmax Maximum number of iteration for Newton-Raphson.
        ///@note More details can be found at @ref linkUnxBoucWen3DLink.
        ///@see UnxBoucWen3DLink::theNodes, UnxBoucWen3DLink::theDirection, UnxBoucWen3DLink::theDimension.
        UnxBoucWen3DLink(const std::vector<unsigned int> nodes, std::vector<double> parameters, std::vector<double> variables, const unsigned int dim, const unsigned int dir, double tol=1E-06, unsigned int nmax=50);

        ///Destroys this UnxBoucWen3DLink object.
        ~UnxBoucWen3DLink();

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
        ///@note The mass matrix can be revisited in @ref linkUnxBoucWen3DLink.
        ///@see Assembler::ComputeMassMatrix(), Integrator::ComputeEffectiveStiffness().
        Eigen::MatrixXd ComputeMassMatrix();

        ///Compute the stiffness matrix of the element using gauss-integration.
        ///@return Matrix with the Element stiffness matrix.
        ///@note The stiffness matrix can be revisited in @ref linkUnxBoucWen3DLink.
        ///@see Assembler::ComputeStiffnessMatrix(), Integrator::ComputeEffectiveStiffness().
        Eigen::MatrixXd ComputeStiffnessMatrix();

        ///Compute the damping matrix of the element using gauss-integration.
        ///@return Matrix with the Element damping matrix.
        ///@note The damping matrix can be revisited in @ref linkUnxBoucWen3DLink.
        ///@see Assembler::ComputeDampingMatrix(), Integrator::ComputeEffectiveStiffness().
        Eigen::MatrixXd ComputeDampingMatrix();

        ///Compute the PML history matrix using gauss-integration.
        ///@return Matrix with the PML Element matrix.
        ///@note The PML matrix is none existent for this element.
        ///@see Assembler::ComputePMLHistoryMatrix(), Integrator::ComputeEffectiveStiffness().
        Eigen::MatrixXd ComputePMLMatrix();

        ///Compute the internal (elastic) forces acting on the element.
        ///@return Vector with the Element internal force.
        ///@note The internal force vector can be revisited in @ref linkUnxBoucWen3DLink.
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
        ///@note The surface force vector can be revisited in @ref linkUnxBoucWen3DLink.
        ///@see Assembler::ComputeExternalForceVector(), Integrator::ComputeEffectiveForce().
        Eigen::VectorXd ComputeSurfaceForces(const std::shared_ptr<Load> &surface, unsigned int face);

        ///Compute the body forces acting on the element.
        ///@param body Pointer to the Load object that contains this information.
        ///@param k The time step at which the body load is evaluated.
        ///@return Vector with the Element surface force.
        ///@note The body force vector can be revisited in @ref linkUnxBoucWen3DLink.
        ///@see Assembler::ComputeExternalForceVector(), Integrator::ComputeEffectiveForce().
        Eigen::VectorXd ComputeBodyForces(const std::shared_ptr<Load> &body, unsigned int k=0);

        ///Compute the domain reduction forces acting on the element.
        ///@param drm Pointer to the DRM Load object that contains this information.
        ///@param k The time step at which the body load is evaluated.
        ///@return Vector with the Element domain reduction forces.
        ///@note The DRM force vector can be revisited in @ref linkUnxBoucWen3DLink.
        ///@see Assembler::ComputeExternalForceVector(), Integrator::ComputeEffectiveForce().
        Eigen::VectorXd ComputeDomainReductionForces(const std::shared_ptr<Load> &drm, unsigned int k);

    private:
        ///Initial stiffness of hysteretic component
        double Ko;

        ///Yield force 
        double Fy;

        ///Yielding exponent (sharpness of hysteresis loop corners) 
        double alpha;

        ///Yielding exponent (sharpness of hysteresis loop corners) 
        double eta;

        ///First hysteretic shape parameter
        double beta;

        ///second hysteretic shape parameter
        double gamma;

        ///Newton-Raphson Tolerance.
        double Tol;

        ///Trial Hysteretic evolution parameters
        double z;

        ///Hysteretic evolution parameters
        double zn;

        ///Trial displacements in local coordinates
        double U;

        ///Displacements in local coordinates
        double Un;

        ///Non-linear Bouc-Wen internal force.
        double qbw;

        ///Consistent Bouc-Wen stiffness.
        double kbw;

        ///Newton-Raphson maximum number of iterations.
        unsigned int nMax;

        ///The element dimension in global coordinates.
        unsigned int Dimension;

        ///The element direction in global coordinates.
        unsigned int Direction;
    
        ///The Element's Nodes.
        std::vector<std::shared_ptr<Node> > theNodes;

        ///Sign function.
        double sign(double x) const;

        ///Update strain in the element:
        ///@return Vector with the total relative deformation.
        Eigen::VectorXd ComputeRelativeDeformation() const;
    
        ///Compute local axes.
        ///@return Matrix with the local axis transformation.
        ///@note Axes transformation are according to @ref linkUnxBoucWen3DLink.
        Eigen::MatrixXd ComputeLocalAxes() const;

        ///Compute rotation matrix.
        ///@return Matrix with the rotation.
        ///@note Axes transformation are according to @ref linkUnxBoucWen3DLink.
        Eigen::MatrixXd ComputeRotationMatrix() const;
};

#endif
