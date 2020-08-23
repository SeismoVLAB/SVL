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
//  [1] Tuning of Surface-Acoustic-Wave Dispersion via Magnetically Modulated 
//      Contact Resonances Antonio Palermo, Yifan Wang, Paolo Celli, and Chiara 
//      Daraio Phys. Rev. Applied 11, 044057 – Published 18 April 2019
//
// Description:
///This file contains the "Hertzian1DLinear" material declarations, which 
///defines an uniaxial isotropic linear elastic material for one-dimensional 
///elements.
//------------------------------------------------------------------------------

#ifndef _HERTZIAN1DLINEAR_HPP_
#define _HERTZIAN1DLINEAR_HPP_

#include <string>
#include <memory>
#include <Eigen/Dense>

#include "Material.hpp"

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      January 5, 2020
/// @version   1.0
/// @file      Hertzian1DLinear.hpp
/// @class     Hertzian1DLinear
/// @see       Material.hpp
/// @brief     Class for creating an uniaxial hertzian contact material for one-dimensional elements
class Hertzian1DLinear : public Material{

    public:
        ///Creates a Hertzian1DLinear material to be especified at a Gauss-point in an Element.
        ///@param k1 The linear spring value.
        ///@param k2 The quadratic spring value.
        ///@param k3 The cubic spring value.
        ///@param rho The material density.
        ///@see Hertzian1DLinear::k1, Hertzian1DLinear::k2, Hertzian1DLinear::k3, Hertzian1DLinear::Rho.
        Hertzian1DLinear(const double k1, const double k2, const double k3, const double rho = 0.0);

        ///Destroys this Hertzian1DLinear material.
        ~Hertzian1DLinear();

        //Clone the 'Hertzian1DLinear' material.
        std::unique_ptr<Material> CopyMaterial();

        ///Access material density.
        ///@return The material density value.
        ///@see Hertzian1DLinear::Rho.
        double GetDensity() const;

        ///Returns the Poisson's ratio.
        ///@return The material poisson's ratio.
        double GetPoissonRatio() const;

        ///Access bulk modulus.
        ///@return The material bulk's (K) modulus.
        double GetBulkModulus() const;

        ///Access shear modulus.
        ///@return The material shear's (G) modulus.
        double GetShearModulus() const;

        ///Access modulus of elasticity.
        ///@return The material Elasticity's (E) modulus.
        double GetElasticityModulus() const;

        ///Returns the material viscous damping.
        ///@return Vector with the material damping components.
        Eigen::MatrixXd GetDamping() const;

        ///Returns the material strain.
        ///@return Vector with the material strain components.
        ///@note The strain can be revisited in @ref linkHertzian1DLinear.
        Eigen::VectorXd GetStrain() const;

        ///Returns the material stress. 
        ///@return Vector with the material stress components.
        ///@note The stress can be revisited in @ref linkHertzian1DLinear.
        Eigen::VectorXd GetStress() const;

        ///Returns the material strain rate.
        ///@return Vector with the material strain-rate components.
        ///@note The strain-rate can be revisited in @ref linkHertzian1DLinear.
        Eigen::VectorXd GetStrainRate() const;

        ///Computes the material total stress.
        ///@return Vector with the material total stress components.
        Eigen::VectorXd GetTotalStress() const;

        ///Returns the material stiffness.
        ///@return Matrix with the material tangent stiffness matrix.
        ///@note The stiffness matrix can be revisited in @ref linkHertzian1DLinear.
        Eigen::MatrixXd GetTangentStiffness() const;

        ///Returns the initial material stiffness.
        ///@return Matrix with the initial material tangent stiffness matrix.
        ///@note The initial tangent stiifness matrix is computed when the strain vector is zero.
        Eigen::MatrixXd GetInitialTangentStiffness() const;

        ///Perform converged material state update.
        ///@note This function sets the trail stress and strain as converged.
        void CommitState();

        ///Update the material state for this iteration.
        ///@param strain Vector with the strain components at this Gauss-point.
        ///@param cond If the the elatic/platic material components will be updated.
        ///@note This function computes the strain and tanget stiffness matrix once the trial strain converged.
        void UpdateState(const Eigen::VectorXd strain, const unsigned int cond);

    private:
        ///The linear Spring Constant.
        double k1;

        ///The quadratic Spring Constant.
        double k2;

        ///The cubic Spring Constant.
        double k3;

        ///Material density.
        double Rho;

        ///Strain vector.
        Eigen::VectorXd Strain;

        ///Stress vector.
        Eigen::VectorXd Stress;

        ///Tangent stiffness matrix.
        Eigen::MatrixXd TangentStiffness;
};

#endif
