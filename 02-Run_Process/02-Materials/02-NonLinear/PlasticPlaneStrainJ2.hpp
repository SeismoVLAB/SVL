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
//  [1] J.C. Simo, T.J.R. Hughes, Computational Inelasticity, Springer, 1997, pp.124
//  [2] R.I. Borja, Plasticity Modeling & Computation, Springer, 2013, pp. 47
//
// Description:
///This file contains the "PlasticPlaneStrainJ2" material declarations, which 
///defines an biaxial plastic material for two-dimensional elements.
//------------------------------------------------------------------------------

#ifndef _PLASTICPLANESTRAINJ2_HPP_
#define _PLASTICPLANESTRAINJ2_HPP_

#include <string>
#include <memory>
#include <Eigen/Dense>

#include "Material.hpp"

/// @author    Elnaz E. Seylabi (elnaze@unr.edu)
/// @date      August 20, 2020
/// @version   1.0
/// @file      PlasticPlaneStrainJ2.hpp
/// @class     PlasticPlaneStrainJ2
/// @see       Material.hpp
/// @brief     Class for creating a biaxial plastic J2 material for two-dimensional elements
class PlasticPlaneStrainJ2 : public Material{

    public:
        ///Creates a PlasticPlaneStrainJ2 material to be specified at a Gauss-point in an Element.
        ///@param K The material bulk modulus.
        ///@param G The material shear modulus.
        ///@param rho The material density.
        ///@param Hbar The kinematic hardening modulus.
        ///@param beta The linear hardening parameter.
        ///@param SigmaY The yield stress.
        ///@see PlasticPlaneStrainJ2::K, PlasticPlaneStrainJ2::G, PlasticPlaneStrainJ2::Rho, PlasticPlaneStrainJ2::H, PlasticPlaneStrainJ2::beta, PlasticPlaneStrainJ2::SigmaY.
        PlasticPlaneStrainJ2(const double K, const double G, const double rho, const double Hbar, const double beta, const double SigmaY);

        ///Destroys this PlasticPlaneStrainJ2 material.
        ~PlasticPlaneStrainJ2();

        //Clone the 'PlasticPlaneStrainJ2' material.
        ///@return A Material pointer to the PlasticPlaneStrainJ2 derived class.
        ///@see Material, PlasticPlaneStrainJ2.
        std::unique_ptr<Material> CopyMaterial();

        ///Access material density.
        ///@return The material density value.
        ///@see PlasticPlaneStrainJ2::Rho.
        double GetDensity() const;

        ///Returns the Poisson's ratio.
        ///@return The material poisson's ratio.
        double GetPoissonRatio() const;

        ///Access bulk modulus.
        ///@return The material bulk's modulus.
        ///@see PlasticPlaneStrainJ2::K.
        double GetBulkModulus() const;

        ///Access shear modulus.
        ///@return The material shear's modulus.
        ///@see PlasticPlaneStrainJ2::G.
        double GetShearModulus() const;

        ///Access modulus of elasticity.
        ///@return The material Elasticity's modulus.
        double GetElasticityModulus() const;

        ///Access the material's energy at current strain.
        ///@return Scalar with the material energy value.
        double GetEnergy() const;

        ///Returns the material viscous damping.
        ///@return Vector with the material damping components.
        Eigen::MatrixXd GetDamping() const;

        ///Returns the material strain.
        ///@return Vector with the material strain components.
        ///@note The strain can be revisited in @ref linkPlasticPlaneStrainJ2.
        ///@see PlasticPlaneStrainJ2::Strain.
        Eigen::VectorXd GetStrain() const;

        ///Returns the material stress. 
        ///@return Vector with the material stress components.
        ///@note The stress can be revisited in @ref linkPlasticPlaneStrainJ2.
        ///@see PlasticPlaneStrainJ2::Stress.
        Eigen::VectorXd GetStress() const;

        ///Returns the material strain rate.
        ///@return Vector with the material strain-rate components.
        ///@note The strain-rate can be revisited in @ref linkPlasticPlaneStrainJ2.
        Eigen::VectorXd GetStrainRate() const;

        ///Computes the material total stress.
        ///@return Vector with the material total stress components.
        ///@see PlasticPlaneStrainJ2::Stress.
        Eigen::VectorXd GetTotalStress() const;

        ///Returns the material stiffness.
        ///@return Matrix with the material consistent tangent stiffness matrix.
        ///@note The stiffness matrix can be revisited in @ref linkPlasticPlaneStrainJ2.
        ///@see PlasticPlaneStrainJ2::TangentStiffness.
        Eigen::MatrixXd GetTangentStiffness() const;

        ///Returns the initial material stiffness.
        ///@return Matrix with the initial material tangent stiffness matrix.
        ///@note The initial tangent stiffness matrix is computed when the strain vector is zero.
        Eigen::MatrixXd GetInitialTangentStiffness() const;

        ///Perform converged material state update.
        ///@note This function sets the trail stress and strain as converged.
        void CommitState();

        ///Reverse the material states to previous converged state.
        ///@note This funtion returns the material states to previous converged states.
        void ReverseState();

        ///Brings the material states to its initial state in the element.
        ///@note This funtion returns the material states to the beginning.
        void InitialState();

        ///Update the material state for this iteration.
        ///@param strain Vector with the strain components at this Gauss-point.
        ///@param cond If the the elastic/plastic material components will be updated.
        ///@note This function computes the strain and tanget stiffness matrix once the trial strain converged.
        void UpdateState(const Eigen::VectorXd strain, const unsigned int cond);

    protected:
        ///Constructs a Second Order Identity Tensor.
        ///@return The identity second-rank tensor.
        Eigen::VectorXd ComputeIdentityVector() const;

        ///Constructs deviatoric tensor.
        ///@return The deviatoric tensor.
        Eigen::MatrixXd ComputeDeviatoricTensor() const;

        ///Constructs fourth-rank identity tensor.
        ///@return A rank-four tensor.
        Eigen::MatrixXd ComputeIdentityTensor() const;

        ///Constructs fourth-rank symmetric identity operator.
        ///@return A rank-four tensor.
        Eigen::MatrixXd ComputeIdentityOperator() const;

        ///Computes the trace of a tensor of second-rank.
        ///@param T A second-rank tensor.
        ///@return An scalar with the trace.
        double ComputeTensorTrace(const Eigen::VectorXd &T);

        ///Computes the norm of a tensor of second-rank.
        ///@param T A second-rank tensor.
        ///@return An scalar with the norm.
        double ComputeTensorNorm(const Eigen::VectorXd &T);

    private:
        ///Bulk modulus.
        double K;

        ///Shear modulus.
        double G;

        ///Material density.
        double Rho;
    
        ///Kinematic hardening modulus.
        double H;

        ///The linear hardening parameter.
        double beta;
    
        ///Yield stress
        double SigmaY;
    
        ///Internal hardening variable
        double alpha;

        ///Internal hardening variable at time n.
        double alpha_n;
    
        ///Trial Strain vector
        Eigen::VectorXd Strain;

        ///Trial Stress vector.
        Eigen::VectorXd Stress;

        ///Trial Back stress.
        Eigen::VectorXd BackStress;

        ///Trial Plastic strain.
        Eigen::VectorXd PlasticStrain;

        ///Trial Consistent tangent stiffness.
        Eigen::MatrixXd TangentStiffness;

        ///Back stress at time n.
        Eigen::VectorXd BackStress_n;

        ///Plastic strain at time n.
        Eigen::VectorXd PlasticStrain_n;

        ///Strain at time n
        Eigen::VectorXd Strain_n;

        ///Stress vector at time n
        Eigen::VectorXd Stress_n;

        ///Consistent tangent stiffness at time n.
        Eigen::MatrixXd TangentStiffness_n;
};

#endif
