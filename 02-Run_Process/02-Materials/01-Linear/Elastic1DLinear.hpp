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
//  [1] Finite Element Procedures, Bathe, K.J., Chapter 4: pages 194, Table 4.3,  
//      Prentice-Hall, 1996. 
//
// Description:
///This file contains the "Elastic1DLinear" material declarations, which defines an 
///uniaxial isotropic linear elastic material for one-dimensional elements. 
//------------------------------------------------------------------------------

#ifndef _ELASTIC1DLINEAR_HPP_
#define _ELASTIC1DLINEAR_HPP_

#include <string>
#include <memory>
#include <Eigen/Dense>

#include "Material.hpp"

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      January 17, 2018
/// @version   1.0
/// @file      Elastic1DLinear.hpp
/// @class     Elastic1DLinear
/// @see       Material.hpp
/// @brief     Class for creating an uniaxial isotropic linear elastic material for one-dimensional elements
class Elastic1DLinear : public Material{

    public:
        ///Creates a Elastic1DLinear material to be specified at a Gauss-point in an Element.
        ///@param E The material elasticity modulus.
        ///@param nu The material Poisson's ratio.
        ///@param rho The material density.
        ///@see Elastic1DLinear::E, Elastic1DLinear::nu, Elastic1DLinear::Rho.
        Elastic1DLinear(const double E, const double nu = 0, const double rho = 0.0);

        ///Destroys this Elastic1DLinear material.
        ~Elastic1DLinear();

        //Clone the 'Elastic1DLinear' material.
        ///@return A Material pointer to the Elastic1DLinear derived class.
        ///@see Material, Elastic1DLinear.
        std::unique_ptr<Material> CopyMaterial();

        ///Access material density.
        ///@return The material density value.
        ///@see Elastic1DLinear::Rho.
        double GetDensity() const;

        ///Returns the Poisson's ratio.
        ///@return The material poisson's ratio.
        ///@see Elastic1DLinear::nu.
        double GetPoissonRatio() const;

        ///Access bulk modulus.
        ///@return The material bulk's (K) modulus.
        double GetBulkModulus() const;

        ///Access shear modulus.
        ///@return The material shear's (G) modulus.
        double GetShearModulus() const;

        ///Access modulus of elasticity.
        ///@return The material Elasticity's (E) modulus.
        ///@see Elastic1DLinear::E.
        double GetElasticityModulus() const;

        ///Access the material's energy at current strain.
        ///@return Scalar with the material energy value.
        double GetEnergy() const;

        ///Returns the material viscous damping.
        ///@return Vector with the material damping components.
        Eigen::MatrixXd GetDamping() const;

        ///Returns the material strain.
        ///@return Vector with the material strain components.
        ///@note The strain can be revisited in @ref linkElastic1DLinear.
        ///@see Elastic1DLinear::Strain.
        Eigen::VectorXd GetStrain() const;

        ///Returns the material stress. 
        ///@return Vector with the material stress components.
        ///@note The stress can be revisited in @ref linkElastic1DLinear.
        ///@see Elastic1DLinear::Stress.
        Eigen::VectorXd GetStress() const;

        ///Returns the material strain rate.
        ///@return Vector with the material strain-rate components.
        ///@note The stress-rate can be revisited in @ref linkElastic1DLinear.
        Eigen::VectorXd GetStrainRate() const;

        ///Computes the material total stress.
        ///@return Vector with the material total stress components.
        ///@note The damping can be revisited in @ref linkElastic1DLinear.
        ///@see Elastic1DLinear::Stress.
        Eigen::VectorXd GetTotalStress() const;

        ///Returns the material stiffness.
        ///@return Matrix with the material tangent stiffness matrix.
        ///@note The stiffness matrix can be revisited in @ref linkElastic1DLinear.
        ///@see Elastic1DLinear::TangentStiffness.
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

        ///Strain vector.
        Eigen::VectorXd Strain;

        ///Stress vector.
        Eigen::VectorXd Stress;

        ///Tangent stiffness matrix.
        Eigen::MatrixXd TangentStiffness;
};

#endif