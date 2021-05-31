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
//   [1] Bathe, K. Jurgen, "Finite Element Procedures", Chapter 4: pages 194, 
//       Table 4.3, Prentice-Hall, 1996.
//
// Description:
///This file contains the "Steel1DFiber" material declarations, which defines 
///a linear material for 1D, 2D, or 3D analysis to be used along with
///zero-length element or fiber section.
//------------------------------------------------------------------------------

#ifndef _ELASTIC1DFIBER_HPP_
#define _ELASTIC1DFIBER_HPP_

#include "Material.hpp"

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      April 17, 2021
/// @version   1.0
/// @file      Elastic1DFiber.hpp
/// @class     Elastic1DFiber
/// @see       Material.hpp
/// @brief     Class for creating an uniaxial linear elastic material compatible with fiber section
class Elastic1DFiber : public Material {

    public:
        ///Creates a Elastic1DFiber material to be specified at a Gauss-point in an Element.
        ///@param E The material elasticity modulus.
        ///@param rho The material density.
        ///@see Elastic1DFiber::E, Elastic1DFiber::nu, Elastic1DFiber::Rho.
        Elastic1DFiber(const double E, const double nu, const double rho = 0.0);

        ///Destroys this Elastic1DFiber material.
        ~Elastic1DFiber();

        //Clone the 'Elastic1DFiber' material.
        ///@return A Material pointer to the Elastic1DFiber derived class.
        ///@see Material, Elastic1DFiber.
        std::unique_ptr<Material> CopyMaterial();

        ///Access material density.
        ///@return The material density value.
        ///@see Elastic1DFiber::Rho.
        double GetDensity() const;

        ///Returns the Poisson's ratio.
        ///@return The material poisson's ratio.
        ///@see Elastic1DFiber::nu.
        double GetPoissonRatio() const;

        ///Access bulk modulus.
        ///@return The material bulk's (K) modulus.
        double GetBulkModulus() const;

        ///Returns linear shear modulus.
        ///@return The material density value.
        double GetShearModulus() const;

        ///Access modulus of elasticity.
        ///@return The material Elasticity's (E) modulus.
        ///@see Elastic1DFiber::E.
        double GetElasticityModulus() const;

        ///Access the material's energy at current strain.
        ///@return Scalar with the material energy value.
        double GetEnergy() const;

        ///Returns the material viscous damping.
        ///@return Vector with the material damping components.
        Eigen::MatrixXd GetDamping() const;

        ///Returns the material strain.
        ///@return Vector with the material strain components.
        ///@note The strain can be revisited in @ref linkElastic1DFiber.
        ///@see Elastic1DFiber::Strain.
        Eigen::VectorXd GetStrain() const;

        ///Returns the material stress. 
        ///@return Vector with the material stress components.
        ///@note The stress can be revisited in @ref linkElastic1DFiber.
        ///@see Elastic1DFiber::Stress.
        Eigen::VectorXd GetStress() const;

        ///Returns the material strain rate.
        ///@return Vector with the material strain-rate components.
        ///@note The stress-rate can be revisited in @ref linkElastic1DFiber.
        Eigen::VectorXd GetStrainRate() const;

        ///Computes the material total stress.
        ///@return Vector with the material total stress components.
        ///@note The damping can be revisited in @ref linkElastic1DFiber.
        ///@see Elastic1DFiber::Stress.
        Eigen::VectorXd GetTotalStress() const;

        ///Returns the material stiffness.
        ///@return Matrix with the material tangent stiffness matrix.
        ///@note The stiffness matrix can be revisited in @ref linkElastic1DFiber.
        ///@see Elastic1DFiber::TangentStiffness.
        Eigen::MatrixXd GetTangentStiffness() const;

        ///Returns the initial material stiffness.
        ///@return Matrix with the initial material tangent stiffness matrix.
        ///@note The initial tangent stiffness matrix is computed when the strain vector is zero.
        Eigen::MatrixXd GetInitialTangentStiffness() const;

        ///Perform converged material state update.
        ///@note This function sets the trail stress and strain as converged.
        void CommitState();

        ///Reverse the material states to previous converged state.
        ///@note This function returns the material states to previous converged states.
        void ReverseState();

        ///Brings the material states to its initial state in the element.
        ///@note This function returns the material states to the beginning.
        void InitialState();

        ///Update the material state for this iteration.
        ///@param strain Vector with the strain components at this Gauss-point.
        ///@param cond If the the elastic/plastic material components will be updated.
        ///@note This function computes the strain and tanget stiffness matrix once the trial strain converged.
        void UpdateState(const Eigen::VectorXd strain, const unsigned int cond=0);

    private:
        ///Modulus of elasticity.
        double E;

        ///Modulus of elasticity.
        double nu;

        ///Material density.
        double Rho;

        ///Strain vector.
        Eigen::VectorXd Strain;
};

#endif