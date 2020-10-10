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
//
// Description:
///This file contains the "lin3DShell4" linearized four-node element declarations, 
///which defines a frame element in a finite element mesh.
//------------------------------------------------------------------------------

#ifndef _LIN3DSHELL4_HPP_
#define _LIN3DSHELL4_HPP_

#include <map>
#include <memory>
#include <string>
#include <Eigen/Dense>

#include "Node.hpp"
#include "Load.hpp"
#include "Section.hpp"
#include "Element.hpp"
#include "Damping.hpp"
#include "QuadratureRule.hpp"

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      November 4, 2019
/// @version   1.0
/// @file      lin3DShell4.hpp
/// @class     lin3DShell4
/// @see       Element.hpp Mesh.hpp
/// @brief     Class for creating a 3D linearized four-node shell element in a mesh
class lin3DShell4 : public Element{

    public:
        ///Creates a lin3DShell4 in a finite element Mesh.
        ///@param nodes The Node connectivity array of this Element.
        ///@param section Pointer to the Section that this Element is made out of.
        ///@param nGauss Numner of Gauss-integration points.
        ///@param massform The mass formulation to compute the mass matrix.
        ///@note More details can be found at @ref linklin3DShell4.
        ///@see lin3DShell4::theNodes, lin3DShell4::theMaterial, lin3DShell4::QuadraturePoints.
        lin3DShell4(const std::vector<unsigned int> nodes, std::unique_ptr<Section> &section, const unsigned int nGauss=9, bool massform=true);

        ///Destroys this lin3DShell4 object.
        ~lin3DShell4();

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
        ///@note The mass matrix can be revisited in @ref linklin3DShell4.
        ///@see Assembler::ComputeMassMatrix(), Integrator::ComputeEffectiveStiffness().
        Eigen::MatrixXd ComputeMassMatrix();

        ///Compute the stiffness matrix of the element using gauss-integration.
        ///@return Matrix with the Element stiffness matrix.
        ///@note The stiffness matrix can be revisited in @ref linklin3DShell4.
        ///@see Assembler::ComputeStiffnessMatrix(), Integrator::ComputeEffectiveStiffness().
        Eigen::MatrixXd ComputeStiffnessMatrix();

        ///Compute the damping matrix of the element using gauss-integration.
        ///@return Matrix with the Element damping matrix.
        ///@note The damping matrix can be revisited in @ref linklin3DShell4.
        ///@see Assembler::ComputeDampingMatrix(), Integrator::ComputeEffectiveStiffness().
        Eigen::MatrixXd ComputeDampingMatrix();

        ///Compute the PML history matrix using gauss-integration.
        ///@return Matrix with the PML Element matrix.
        ///@note The PML matrix is none existent for this element.
        ///@see Assembler::ComputePMLHistoryMatrix(), Integrator::ComputeEffectiveStiffness().
        Eigen::MatrixXd ComputePMLMatrix();

        ///Compute the internal (elastic) forces acting on the element.
        ///@return Vector with the Element internal force.
        ///@note The internal force vector can be revisited in @ref linklin3DShell4.
        ///@see Assembler::ComputeInternalForceVector(), Integrator::ComputeEffectiveForce().
        Eigen::VectorXd ComputeInternalForces();

        ///Compute the elastic, inertial, and vicous forces acting on the element.
        ///@return Vector with the Element dynamic internal force.
        ///@note The internal force vector can be revisited in @ref linkElement.
        ///@see Assembler::ComputeDynamicInternalForceVector().
        Eigen::VectorXd ComputeInternalDynamicForces();

        ///Compute the PML history vector using gauss-integration.
        ///@return Vector with the PML Element history values.
        ///@note The PML vector is none existent for this element.
        ///@see Assembler::ComputePMLHistoryMatrix(), Integrator::ComputeEffectiveStiffness().
        Eigen::VectorXd ComputePMLVector();

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
        ///Mass Formulation.
        bool MassForm;

        ///The Damping model.
        std::shared_ptr<Damping> theDamping;

        ///The Element's Nodes.
        std::vector<std::shared_ptr<Node> > theNodes;

        ///The Element's Sections.
        std::vector<std::unique_ptr<Section> > theSection;

        ///Coordinate of Gauss points.
        std::unique_ptr<QuadratureRule> QuadraturePoints;

        ///Computes the total element mass.
        double ComputeTotalMass();

        ///Compute/update the local axis-1-2-3 of the element.
        ///@return Matrix with the local axis transformation.
        ///@note Axes transformation are according to @ref linklin3DShell4.
        Eigen::MatrixXd ComputeLocalAxes() const;

        ///Compute/update the projection axis of the element.
        ///@return Matrix with the local axes x1-x2-x3 transformation.
        ///@note Axes transformation are according to @ref linklin3DShell4.
        Eigen::MatrixXd ComputeRotation() const;

        ///Compute/update the node coordinates in axis-1-2-3 of the element.
        ///return The xyloc Node coordinates in 3D projected into r-s plane.
        Eigen::MatrixXd ComputeLocalCoordinates() const;

        ///Update strain in the element.
        ///@param xyloc Node coordinates in 3D projected into r-s plane.
        ///@param BM12 Matrix to enforce constant Tension (wilson) inside element.
        ///@param ri The i-th Gauss coordinate in the r-axis.
        ///@param si The i-th Gauss coordinate in the s-axis.
        Eigen::VectorXd ComputeStrain(Eigen::MatrixXd &BM12, Eigen::MatrixXd &xyloc, double ri, double si);

        ///Update strain rate in the element.
        ///@param xyloc Node coordinates in 3D projected into r-s plane.
        ///@param BM12 Matrix to enforce constant Tension (wilson) inside element.
        ///@param ri The i-th Gauss coordinate in the r-axis.
        ///@param si The i-th Gauss coordinate in the s-axis.
        Eigen::VectorXd ComputeStrainRate(Eigen::MatrixXd &BM12, Eigen::MatrixXd &xyloc, double ri, double si);

        ///Correction matrix to Enforce Constant Tension for membrane effect (wilson).
        Eigen::MatrixXd ComputeConstantTensionMatrix();

        ///Assembles the membrane and plate effects into a single matrix.
        ///@param A The matrix with the shell behavior.
        ///@param Am The matrix with the membrane contribution.
        ///@param Ap The matrix with the plate contribution.
        void AssemblePlateMembraneEffects(Eigen::MatrixXd &A, const Eigen::MatrixXd &Am, const Eigen::MatrixXd &Ap);

        ///Evaluates the Constant Tension matrix and Jacobian for membrane effect at a given Gauss point.
        ///@param ri The i-th Gauss coordinate in the r-axis.
        ///@param si The i-th Gauss coordinate in the s-axis.
        ///@param xyloc Node coordinates in 3D projected into r-s plane.
        ///@param BM12 Matrix to enforce constant Tension (wilson) inside element.
        ///@param Jij The Jacobian matrix evaluated at (ri,si).
        void ConstantTensionMatrix(const double ri, const double si, const Eigen::MatrixXd &xyloc, Eigen::MatrixXd &BM12, Eigen::MatrixXd &Jij);

        ///Evaluates the shape function matrix for plate at a given Gauss point.
        ///@param ri The i-th Gauss coordinate in the r-axis.
        ///@param si The i-th Gauss coordinate in the s-axis.
        ///@param xyloc Node coordinates in 3D projected into r-s plane.
        ///@param Hij The shape function matrix evaluated at (ri,si).
        ///@param Jij The Jacobian matrix evaluated at (ri,si).
        void ComputeLoadShapeFunctionMatrix(const double ri, const double si, const Eigen::MatrixXd &xyloc, Eigen::MatrixXd &Hij, Eigen::MatrixXd &Jij);

        ///Evaluates the shape function matrix for plate at a given Gauss point.
        ///@param ri The i-th Gauss coordinate in the r-axis.
        ///@param si The i-th Gauss coordinate in the s-axis.
        ///@param xyloc Node coordinates in 3D projected into r-s plane.
        ///@param Hij The plate shape function matrix evaluated at (ri,si).
        ///@param Jij The Jacobian matrix evaluated at (ri,si).
        void ComputePlateShapeFunctionMatrix(const double ri, const double si, const Eigen::MatrixXd &xyloc, Eigen::MatrixXd &Hij, Eigen::MatrixXd &Jij);

        ///Evaluates the shape function matrix for membrane at a given Gauss point.
        ///@param ri The i-th Gauss coordinate in the r-axis.
        ///@param si The i-th Gauss coordinate in the s-axis.
        ///@param xyloc Node coordinates in 3D projected into r-s plane.
        ///@param Hij The membrane shape function matrix evaluated at (ri,si).
        ///@param Jij The Jacobian matrix evaluated at (ri,si).
        void ComputeMembraneShapeFunctionMatrix(const double ri, const double si, const Eigen::MatrixXd &xyloc, Eigen::MatrixXd &Hij, Eigen::MatrixXd &Jij);

        ///Evaluates the strain-displacement matrix for a plate effect at a given Gauss point.
        ///@param ri The i-th Gauss coordinate in the r-axis.
        ///@param si The i-th Gauss coordinate in the s-axis.
        ///@param xyloc Node coordinates in 3D projected into r-s plane.
        ///@param Bij The plate strain-displacement matrix evaluated at (ri,si).
        ///@param Jij The Jacobian matrix evaluated at (ri,si).
        void ComputePlateStrainDisplacementMatrix(const double ri, const double si, const Eigen::MatrixXd &xyloc, Eigen::MatrixXd &Bij, Eigen::MatrixXd &Jij); 

        ///Evaluates the strain-displacement matrix for membrane effect at a given Gauss point.
        ///@param ri The i-th Gauss coordinate in the r-axis.
        ///@param si The i-th Gauss coordinate in the s-axis.
        ///@param xyloc Node coordinates in 3D projected into r-s plane.
        ///@param Bij The membrane strain-displacement matrix evaluated at (ri,si).
        ///@param BM12 Matrix to enforce constant Tension (wilson) inside element.
        ///@param Pij The membrane penalty matrix evaluated at (ri,si).
        ///@param Jij The Jacobian matrix evaluated at (ri,si).
        void ComputeMembraneStrainDisplacementMatrix(const double ri, const double si, const Eigen::MatrixXd &xyloc, const Eigen::MatrixXd &Bij, Eigen::MatrixXd &BM12, Eigen::MatrixXd &Pij, Eigen::MatrixXd &Jij);

        ///Compute the initial stiffness matrix of the element.
        ///@return Matrix with the initial Tanget stiffness matrix.
        ///@note This matrix is computed when the strain deformation is zero.
        Eigen::MatrixXd ComputeInitialStiffnessMatrix();
};

#endif
