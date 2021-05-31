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
//  [1] R.J. Borja, A.P. Amies (1994), Multiaxial cyclic plasticity model for 
//      clays, ASCE.
//
// Description:
///This file contains the "PlasticPlaneStrainBA" material declarations, which 
///defines an biaxial plastic bounding surface material for two-dimensional 
///elements. 
//------------------------------------------------------------------------------

#ifndef _PLASTICPLANESTRAINBA_HPP_
#define _PLASTICPLANESTRAINBA_HPP_

#include <string>
#include <Eigen/Dense>
#include <iostream>

#include "Material.hpp"

/// @author    Elnaz E. Seylabi (elnaze@unr.edu)
/// @date      August 18, 2020
/// @version   1.0
/// @file      PlasticPlaneStrainBA.hpp
/// @class     PlasticPlaneStrainBA
/// @see       Material.hpp
/// @brief     Class for creating a biaxial plastic bounding surface material for two-dimensional elements 
class PlasticPlaneStrainBA : public Material{

    public:
        ///Creates a PlasticPlaneStrainBA material to be specified at a Gauss-point in an Element.
        ///@param K The material bulk modulus.
        ///@param G The material shear modulus.
        ///@param rho The material density.
        ///@param H0 The bounding surface hardening.
        ///@param h The hardening parameter.
        ///@param m The exponential hardening parameter.
        ///@param Su Undrained soil strength.
        ///@param beta The integration factor.
        ///@see PlasticPlaneStrainBA::K, PlasticPlaneStrainBA::G, PlasticPlaneStrainBA::Rho, PlasticPlaneStrainBA::H0, PlasticPlaneStrainBA::h, PlasticPlaneStrainBA::m, PlasticPlaneStrainBA::Su, PlasticPlaneStrainBA::R, PlasticPlaneStrainBA::beta.
        PlasticPlaneStrainBA(const double K, const double G, const double rho, const double H0, const double h, const double m, const double Su, const double beta);

        ///Destroys this PlasticPlaneStrainBA material.
        ~PlasticPlaneStrainBA();

        //Clone the 'PlasticPlaneStrainBA' material.
        ///@return A Material pointer to the PlasticPlaneStrainBA derived class.
        ///@see Material, PlasticPlaneStrainBA.
        std::unique_ptr<Material> CopyMaterial();

        ///Access material density.
        ///@return The material density value.
        ///@see PlasticPlaneStrainBA::Rho.
        double GetDensity() const;

        ///Returns the Poisson's ratio.
        ///@return The material poisson's ratio.
        double GetPoissonRatio() const;

        ///Access bulk modulus.
        ///@return The material bulk's modulus.
        ///@see PlasticPlaneStrainBA::K.
        double GetBulkModulus() const;

        ///Access shear modulus.
        ///@return The material shear's modulus.
        ///@see PlasticPlaneStrainBA::G.
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
        ///@note The damping can be revisited in @ref linkPlasticPlaneStrainBA.
        ///@see PlasticPlaneStrainBA::Strain.
        Eigen::VectorXd GetStrain() const;

        ///Returns the material stress. 
        ///@return Vector with the material stress components.
        ///@note The stress can be revisited in @ref linkPlasticPlaneStrainBA.
        ///@see PlasticPlaneStrainBA::Stress.
        Eigen::VectorXd GetStress() const;

        ///Returns the material strain rate.
        ///@return Vector with the material strain-rate components.
        ///@note The strain-rate can be revisited in @ref linkPlasticPlaneStrainBA.
        Eigen::VectorXd GetStrainRate() const;

        ///Computes the material total stress.
        ///@return Vector with the material total stress components.
        ///@see PlasticPlaneStrainBA::Stress.
        Eigen::VectorXd GetTotalStress() const;

        ///Returns the material stiffness.
        ///@return Matrix with the material consistent tangent stiffness matrix.
        ///@note The stiffness matrix can be revisited in @ref linkPlasticPlaneStrainBA.
        ///@see PlasticPlaneStrainBA::TangentStiffness.
        Eigen::MatrixXd GetTangentStiffness() const;

        ///Returns the initial material stiffness.
        ///@return Matrix with the initial material tangent stiffness matrix.
        ///@note The initial tangent stiffness matrix is computed when the strain vector is zero.
        Eigen::MatrixXd GetInitialTangentStiffness() const;

        ///Perform converged material state update.
        ///@note This function sets the trail stress and strain as converged.
        ///@see PlasticPlaneStrainBA::Strain_n, PlasticPlaneStrainBA::Stress_n.
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

        ///Compute inner product of two tensor of first-rank.
        ///@param V1 A first-rank tensor.
        ///@param V2 A first-rank tensor.
        ///@return An scalar with the inner product.
        double ComputeInnerProduct(const Eigen::VectorXd &V1, const Eigen::VectorXd &V2);

        ///Compute hardening using the Newton-Rhapson method.
        ///@param DeviatoricIncrStrain The deviatoric incremental strain.
        ///@param DeviatoricStress The deviatoric stress.
        ///@return Vector that contains kappa and psi using Newton-Rhapson method.
        Eigen::VectorXd ComputeHardening(const Eigen::VectorXd &_DeviatoricIncrStrain, const Eigen::VectorXd &_DeviatoricStress);

        ///Compute hardening using the bi-section method
        ///@param DeviatoricIncrStrain The deviatoric incremental strain.
        ///@param DeviatoricStress The deviatoric stress.
        ///@return Vector that contains kappa and psi using bisection method
        Eigen::VectorXd ComputeHardeningBisection(const Eigen::VectorXd &_DeviatoricIncrStrain, const Eigen::VectorXd &_DeviatoricStress);

    private:
        ///Bulk modulus.
        double K;

        ///Shear modulus.
        double G;

        ///Material density.
        double Rho;
    
        ///Bounding surface hardening.
        double H0;

        ///exponential hardening parameter.
        double h;
    
        ///Exponential hardening parameter.
        double m;

        ///Undrained soil strength.
        double Su;

        ///Bounding surface radius.
        double R;

        ///integration factor.
        double beta;

        ///Internal hardening parameter.
        double psi;

        ///Internal hardening parameter.
        double kappa;
    
        ///Strain vector.
        Eigen::VectorXd Strain;

        ///Strain at time n.
        Eigen::VectorXd Strain_n;

        ///Trial stress vector.
        Eigen::VectorXd Stress;

        ///Stress vector at time n.
        Eigen::VectorXd Stress_n;

        ///Stress tensor at F0.
        Eigen::VectorXd DeviatoricStress0;

        ///Consistent tangent stiffness.
        Eigen::MatrixXd TangentStiffness;

        ///First loading flag
        int FirstLoadFlag;
        
        ///Flag to control the linear system solution.
        int rootFlag;
        
        ///The hardening function.
        double Hn;
};

#endif