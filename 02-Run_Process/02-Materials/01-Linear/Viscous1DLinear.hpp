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
//  [1] 
//
// Description:
///This file contains the "Viscous1DLinear" material declarations, which 
///defines an uniaxial viscous material in for one-dimensional elements.
//------------------------------------------------------------------------------

#ifndef _VISCOUS1DLINEAR_HPP_
#define _VISCOUS1DLINEAR_HPP_

#include <string>
#include <memory>
#include <Eigen/Dense>

#include "Material.hpp"

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      February 14, 2019
/// @version   1.0
/// @file      Viscous1DLinear.hpp
/// @class     Viscous1DLinear
/// @see       Material.hpp
/// @brief     Class for creating an uniaxial viscous material in for one-dimensional elements
class Viscous1DLinear : public Material{

    public:
        ///Creates a Viscous1DLinear material to be especified at a Gauss-point in an Element.
        ///@param eta The material viscous coefficient.
        ///@see Viscous1DLinear::eta.
        Viscous1DLinear(const double eta);

        ///Destroys this Viscous1DLinear material.
        ~Viscous1DLinear();

        //Clone the 'Viscous1DLinear' material.
        ///@return A Material pointer to the Viscous1DLinear derived class.
        ///@see Material, Viscous1DLinear.
        std::unique_ptr<Material> CopyMaterial();

        ///Access material density.
        ///@return The material density value.
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
        ///@see Viscous1DLinear::Damping.
        Eigen::MatrixXd GetDamping() const;

        ///Returns the material strain.
        ///@return Vector with the material strain components.
        ///@note The strain can be revisited in @ref linkViscous1DLinear.
        Eigen::VectorXd GetStrain() const;

        ///Returns the material stress. 
        ///@return Vector with the material stress components.
        ///@note The stress can be revisited in @ref linkViscous1DLinear.
        Eigen::VectorXd GetStress() const;

        ///Returns the material strain rate.
        ///@return Vector with the material strain-rate components.
        ///@note The strain-rate can be revisited in @ref linkViscous1DLinear.
        ///@see Viscous1DLinear::StrainRate.
        Eigen::VectorXd GetStrainRate() const;

        ///Computes the material total stress.
        ///@return Vector with the material total stress components.
        ///@see Viscous1DLinear::Stress.
        Eigen::VectorXd GetTotalStress() const;

        ///Returns the material stiffness.
        ///@return Matrix with the material tangent stiffness matrix.
        ///@note The stiffness matrix can be revisited in @ref linkViscous1DLinear.
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
        ///The material viscosity.
        double eta;

        ///The material strain-rate vector.
        Eigen::VectorXd StrainRate;

        ///The material stress vector.
        Eigen::VectorXd Stress;

        ///The material damping matrix.
        Eigen::MatrixXd Damping;
};

#endif
