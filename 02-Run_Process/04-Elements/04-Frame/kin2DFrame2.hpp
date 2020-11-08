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
//  [1] Finite Element Procedures, Bathe, K.J., Chapter 5: pages 414-422. 
//      Prentice-Hall, 1996. 
//  [2] 2D Co-rotational Beam Formulation, Louie L. Yaw, E. F. Cross School of 
//      Engineering, Walla Walla University, November 30, 2009.
//
// Description:
///This file contains the "kin2DFrame2" finite kinematic two-node element 
///declarations, which defines a frame element in a finite element mesh. 
//------------------------------------------------------------------------------

#ifndef _KIN2DFRAME2_HPP_
#define _KIN2DFRAME2_HPP_

#include <map>
#include <memory>
#include <string>
#include <Eigen/Dense>

#include "Node.hpp"
#include "Load.hpp"
#include "Section.hpp"
#include "Element.hpp"
#include "Damping.hpp"

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      July 15, 2019
/// @version   1.0
/// @file      kin2DFrame2.hpp
/// @class     kin2DFrame2
/// @see       Element.hpp Mesh.hpp
/// @brief     Class for creating a 2D finite kinematic two-node frame element in a mesh
class kin2DFrame2 : public Element{

    public:
        ///Creates a kin2DFrame2 in a finite element Mesh.
        ///@param nodes The Node connectivity array of this Element.
        ///@param material Pointer to the Material that this Element is made out of.
        ///@param area The cross-section area of the kin3DTruss2 Element.
        ///@param massform The mass formulation to compute the mass matrix.
        ///@note More details can be found at @ref linkkin2DFrame2.
        ///@see kin3DTruss2::theNodes, kin3DTruss2::theMaterial.
        kin2DFrame2(const std::vector<unsigned int> nodes, std::unique_ptr<Section> &section, bool massform=true);

        ///Destroys this kin3DTruss2 object.
        ~kin2DFrame2();

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

        ///Compute the lumped/consistent mass matrix of the element.
        ///@return Matrix with the Element mass matrix.
        ///@note The mass matrix can be revisited in @ref linkkin2DFrame2.
        ///@see Assembler::ComputeMassMatrix(), Integrator::ComputeEffectiveStiffness().
        Eigen::MatrixXd ComputeMassMatrix();

        ///Compute the stiffness matrix of the element using gauss-integration.
        ///@return Matrix with the Element stiffness matrix.
        ///@note The stiffness matrix can be revisited in @ref linkkin2DFrame2.
        ///@see Assembler::ComputeStiffnessMatrix(), Integrator::ComputeEffectiveStiffness().
        Eigen::MatrixXd ComputeStiffnessMatrix();

        ///Compute the damping matrix of the element using gauss-integration.
        ///@return Matrix with the Element damping matrix.
        ///@note The damping matrix can be revisited in @ref linkkin2DFrame2.
        ///@see Assembler::ComputeDampingMatrix(), Integrator::ComputeEffectiveStiffness().
        Eigen::MatrixXd ComputeDampingMatrix();

        ///Compute the PML history matrix using gauss-integration.
        ///@return Matrix with the PML Element matrix.
        ///@note The PML matrix is none existent for this element.
        ///@see Assembler::ComputePMLHistoryMatrix(), Integrator::ComputeEffectiveStiffness().
        Eigen::MatrixXd ComputePMLMatrix();

        ///Compute the internal (elastic) forces acting on the element.
        ///@return Vector with the Element internal force.
        ///@note The internal force vector can be revisited in @ref linkkin2DFrame2.
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

    private:    
        ///Length of the element.
        double Lo;

        ///Angle of the element.
        double Alpha;

        ///Mass Formulation.
        bool MassForm;

        ///The Damping model.
        std::shared_ptr<Damping> theDamping;

        ///The Element's Nodes.
        std::vector<std::shared_ptr<Node> > theNodes;

        ///The Element's Section.
        std::unique_ptr<Section> theSection;

        ///Compute the current length of the element.
        double ComputeLength() const;

        ///Compute/update the local axis-1-2-3 of the element.
        ///@return Matrix with the local axes x1-x2-x3 transformation.
        ///@note Axes transformation are according to @ref linkkin2DFrame2.
        Eigen::MatrixXd ComputeLocalAxes() const;

        ///Update strain in the element.
        ///@return Vector with the strain at Element center.
        Eigen::VectorXd ComputeStrain() const;

        ///Update strain rate in the element.
        ///@return Vector with the strain-rate at Element center.
        Eigen::VectorXd ComputeStrainRate() const;

        ///Computes the current axial strain.
        ///@return A value with the constant axial strain.
        double ComputeAxialStrain() const;

        ///Computes the current axial force.
        ///@return A value with the constant axial force.
        double ComputeAxialForce() const;

        ///Computes the current bending moments.
        ///@return Vector with the bending moments at each Node.
        ///@note The bending momments are according to @ref linkkin2DFrame2.
        Eigen::VectorXd ComputeBendingMoment() const;

        ///Computes the current axial and bending moments.
        ///@return Vector with the local internal forces at each Node.
        ///@note The internal local forces are according to @ref linkkin2DFrame2.
        Eigen::VectorXd ComputeLocalForces() const;

        ///The Linear Strain Displacement matrix: 
        ///@param B The linearized strain-displacement matrix.
        ///@note The strain-displacement matrix are according to @ref linkkin2DFrame2.
        void ComputeLinearStrainDisplacementMatrix(Eigen::MatrixXd &B) const;

        ///The Linear/Geometric Strain Displacement matrix: 
        ///@param B The linearized strain-displacement matrix.
        ///@param z Geometric Axial Strain-displacement vector.
        ///@param r Geometric Bending Strain-displacement vector.
        ///@note The previous matrix and vectors are according to @ref linkkin2DFrame2.
        void ComputeGeometricStrainDisplacementMatrix(Eigen::MatrixXd &B, Eigen::VectorXd &z, Eigen::VectorXd &r) const;

        ///Compute the initial stiffness matrix of the element.
        ///@return Matrix with the initial Tanget stiffness matrix.
        ///@note This matrix is computed when the strain deformation is zero.
        Eigen::MatrixXd ComputeInitialStiffnessMatrix() const;
};

#endif
