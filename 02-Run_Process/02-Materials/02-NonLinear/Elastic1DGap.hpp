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
//   [1]
//
// Description:
///This file contains the "Elastic1DGap" material declarations, which defines 
///an elastic gap material for 1D, 2D, or 3D analysis to be used along with
///zero-length element.
//------------------------------------------------------------------------------

#ifndef _ELASTIC1DGAP_HPP_
#define _ELASTIC1DGAP_HPP_

#include "Material.hpp"

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      April 26, 2021
/// @version   1.0
/// @file      Elastic1DGap.hpp
/// @class     Elastic1DGap
/// @see       Material.hpp
/// @brief     Class for creating an uniaxial elastic gap material compatible with zero-length element
class Elastic1DGap : public Material {

    public:
        ///Creates a Elastic1DGap material to be specified at a Gauss-point in an Element.
        ///@param E The elastic moduli of this material
        ///@param gap The separation from which the material start reacting
        ///@param behavior The material behavior (true: Tension, false: Compression)
        ///@see Elastic1DGap::E, Elastic1DGap::gap.
        Elastic1DGap(double E, double gap, bool behavior=false);

        ///Destroys this Elastic1DGap material.
        ~Elastic1DGap();

        //Clone the 'Elastic1DGap' material.
        ///@return A Material pointer to the Elastic1DGap derived class.
        ///@see Material, Elastic1DGap.
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

        ///Returns linear shear modulus.
        ///@return The material density value.
        double GetShearModulus() const;

        ///Access modulus of elasticity.
        ///@return The material Elasticity's (E) modulus.
        ///@see Elastic1DGap::E.
        double GetElasticityModulus() const;

        ///Access the material's energy at current strain.
        ///@return Scalar with the material energy value.
        double GetEnergy() const;

        ///Returns the material viscous damping.
        ///@return Vector with the material damping components.
        Eigen::MatrixXd GetDamping() const;

        ///Returns the material strain.
        ///@return Vector with the material strain components.
        ///@note The strain can be revisited in @ref linkElastic1DGap.
        ///@see Elastic1DGap::oldStrain.
        Eigen::VectorXd GetStrain() const;

        ///Returns the material stress. 
        ///@return Vector with the material stress components.
        ///@note The stress can be revisited in @ref linkElastic1DGap.
        ///@see Elastic1DGap::oldStress.
        Eigen::VectorXd GetStress() const;

        ///Returns the material strain rate.
        ///@return Vector with the material strain-rate components.
        ///@note The stress-rate can be revisited in @ref linkElastic1DGap.
        Eigen::VectorXd GetStrainRate() const;

        ///Computes the material total stress.
        ///@return Vector with the material total stress components.
        ///@note The damping can be revisited in @ref linkElastic1DGap.
        ///@see Elastic1DGap::oldStress.
        Eigen::VectorXd GetTotalStress() const;

        ///Returns the material stiffness.
        ///@return Matrix with the material tangent stiffness matrix.
        ///@note The stiffness matrix can be revisited in @ref linkElastic1DGap.
        ///@see Elastic1DGap::oldTangentStiffness.
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
		//Elastic Gap properties: 
		double E;
		double Gap;
		bool Behavior;

		//Elastic Gap history variables:
		double oldStrain;
};

#endif
