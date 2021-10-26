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
//  [1] J.C. Simo, T.J.R. Hughes, Computational Inelasticity, Springer, 1997, 
//      pp.45
//
// Description:
///This file contains the "Plastic1DJ2" material declarations, which defines 
///an uniaxial isotropic plastic material for one-dimensional elements.
//------------------------------------------------------------------------------

#ifndef _PLASTIC1DJ2_HPP_
#define _PLASTIC1DJ2_HPP_

#include <string>
#include <memory>
#include <Eigen/Dense>

#include "Material.hpp"
 
/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      January 5, 2019
/// @version   1.0
/// @file      Plastic1DJ2.hpp
/// @class     Plastic1DJ2
/// @see       Material.hpp
/// @brief     Class for creating an uniaxial isotropic plastic material for one-dimensional elements
class Plastic1DJ2 : public Material{

    public:
        ///Creates a Plastic1DJ2 material to be specified at a Gauss-point in an Element.
        ///@param E The material elasticity modulus.
        ///@param nu The material Poisson's ratio.
        ///@param rho The material density.
        ///@param K The plastic modulus.
        ///@param H The hardening parameter.
        ///@param SigmaY The yield stress.
        ///@see Plastic1DJ2::E, Plastic1DJ2::nu, Plastic1DJ2::Rho, Plastic1DJ2::H, Plastic1DJ2::K, Plastic1DJ2::SigmaY.
        Plastic1DJ2(const double E, const double nu, const double rho, const double K, const double H, const double SigmaY);

        ///Destroys this Plastic1DJ2 material.
        ~Plastic1DJ2();

        //Clone the 'Plastic1DJ2' material.
        ///@return A Material pointer to the Plastic1DJ2 derived class.
        ///@see Material, Plastic1DJ2.
        std::unique_ptr<Material> CopyMaterial();

        ///Access material density.
        ///@return The material density value.
        ///@see Plastic1DJ2::Rho.
        double GetDensity() const;

        ///Returns the Poisson's ratio.
        ///@return The material poisson's ratio.
        ///@see Plastic1DJ2::nu.
        double GetPoissonRatio() const;

        ///Access bulk modulus.
        ///@return The material bulk's modulus.
        double GetBulkModulus() const;

        ///Access shear modulus.
        ///@return The material shear's modulus.
        double GetShearModulus() const;

        ///Access modulus of elasticity.
        ///@return The material Elasticity's (E) modulus.
        ///@see Plastic1DJ2::E.
        double GetElasticityModulus() const;

        ///Access the material's energy at current strain.
        ///@return Scalar with the material energy value.
        double GetEnergy() const;

        ///Returns the material viscous damping.
        ///@return Vector with the material damping components.
        ///@see Plastic1DJ2::Damping.
        Eigen::MatrixXd GetDamping() const;

        ///Returns the material strain.
        ///@return Vector with the material strain components.
        ///@note The strain can be revisited in @ref linkPlastic1DJ2.
        ///@see Plastic1DJ2::Strain.
        Eigen::VectorXd GetStrain() const;

        ///Returns the material stress. 
        ///@return Vector with the material stress components.
        ///@note The stress can be revisited in @ref linkPlastic1DJ2.
        ///@see Plastic1DJ2::Stress.
        Eigen::VectorXd GetStress() const;

        ///Returns the material strain rate.
        ///@return Vector with the material strain-rate components.
        ///@note The strain-rate can be revisited in @ref linkPlastic1DJ2.
        ///@see Plastic1DJ2::StrainRate.
        Eigen::VectorXd GetStrainRate() const;

        ///Computes the material total stress.
        ///@return Vector with the material total stress components.
        ///@see Plastic1DJ2::Stress.
        Eigen::VectorXd GetTotalStress() const;

        ///Returns the material stiffness.
        ///@return Matrix with the material consistent tangent stiffness matrix.
        ///@note The stiffness matrix can be revisited in @ref linkPlastic1DJ2.
        ///@see Plastic1DJ2::TangentStiffness.
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

    private:
        ///Modulus of elasticity.
        double E;

        ///Poisson's ratio.
        double nu;

        ///Material density.
        double Rho;
    
        ///Plastic modulus.
        double K;
    
        ///Kinematic hardening modulus.
        double H;
    
        ///Yield stress.
        double SigmaY;
    
        ///Internal hardening variable.
        double alpha;

        ///Internal hardening variable at time n.
        double alpha_n;
    
        ///Strain vector.
        Eigen::VectorXd Strain;

        ///Stress vector.
        Eigen::VectorXd Stress;

        ///Back stress.
        Eigen::VectorXd BackStress;

        ///Plastic strain.
        Eigen::VectorXd PlasticStrain;

        ///Consistent tangent stiffness.
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
