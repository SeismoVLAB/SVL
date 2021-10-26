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
//   [1] Filippou, F. C., Popov, E. P., Bertero, V. V. (1983). "Effects of Bond 
//       Deterioration on Hysteretic Behavior of Reinforced Concrete Joints". 
//       Report EERC 83-19, Earthquake Engineering Research Center, University 
//       of California, Berkeley. 
//
// Description:
///This file contains the "Steel1DFiber" material declarations, which defines 
///a steel material for 1D, 2D, or 3D analysis to be used along with
///zero-length element or fiber section.
//------------------------------------------------------------------------------

#ifndef _STEEL1DFIBER_HPP_
#define _STEEL1DFIBER_HPP_

#include "Material.hpp"

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      April 17, 2021
/// @version   1.0
/// @file      Steel1DFiber.hpp
/// @class     Steel1DFiber
/// @see       Material.hpp
/// @brief     Class for creating an uniaxial steel (Menegotto and Pinto 1973) material compatible with fiber section
class Steel1DFiber : public Material {

    public:
        ///Creates a Steel1DFiber material to be specified at a Gauss-point in an Element.
        ///@param fy Represents the steel yield strength.
        ///@param E Represents the elasticity modulus.
        ///@param b Represents the strain-hardening ratio.
        ///@param R0 Represents the parameters to control the transition from elastic to plastic branches.
        ///@param cR1 Represents the parameters to control the transition from elastic to plastic branches.
        ///@param cR2 Represents the parameters to control the transition from elastic to plastic branches.
        ///@param a1 Represents the isotropic hardening parameter.
        ///@param a2 Represents the isotropic hardening parameter.
        ///@param a3 Represents the isotropic hardening parameter.
        ///@param a4 Represents the isotropic hardening parameter.
        ///@param nu Represents the Poisson's ratio.
        ///@param rho Represents the material density. 
        ///@see Steel1DFiber::E, Steel1DFiber::fy, Steel1DFiber::nu, Steel1DFiber::Rho.
        Steel1DFiber(double fy, double E, double b, double R0=15.00, double cR1=0.925, double cR2=0.150, double a1=0.0, double a2=1.0, double a3=0.0, double a4=1.0, double nu=0.33, double rho=0.0);

        ///Destroys this Steel1DFiber material.
        ~Steel1DFiber();

        //Clone the 'Steel1DFiber' material.
        ///@return A Material pointer to the Steel1DFiber derived class.
        ///@see Material, Steel1DFiber.
        std::unique_ptr<Material> CopyMaterial();

        ///Access material density.
        ///@return The material density value.
        ///@see Steel1DFiber::Rho.
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
        ///@see Steel1DFiber::E.
        double GetElasticityModulus() const;

        ///Access the material's energy at current strain.
        ///@return Scalar with the material energy value.
        double GetEnergy() const;

        ///Returns the material viscous damping.
        ///@return Vector with the material damping components.
        Eigen::MatrixXd GetDamping() const;

        ///Returns the material strain.
        ///@return Vector with the material strain components.
        ///@note The strain can be revisited in @ref linkSteel1DFiber.
        ///@see Steel1DFiber::oldStrain.
        Eigen::VectorXd GetStrain() const;

        ///Returns the material stress. 
        ///@return Vector with the material stress components.
        ///@note The stress can be revisited in @ref linkSteel1DFiber.
        ///@see Steel1DFiber::oldStress.
        Eigen::VectorXd GetStress() const;

        ///Returns the material strain rate.
        ///@return Vector with the material strain-rate components.
        ///@note The stress-rate can be revisited in @ref linkSteel1DFiber.
        Eigen::VectorXd GetStrainRate() const;

        ///Computes the material total stress.
        ///@return Vector with the material total stress components.
        ///@note The damping can be revisited in @ref linkSteel1DFiber.
        ///@see Steel1DFiber::oldStress.
        Eigen::VectorXd GetTotalStress() const;

        ///Returns the material stiffness.
        ///@return Matrix with the material tangent stiffness matrix.
        ///@note The stiffness matrix can be revisited in @ref linkSteel1DFiber.
        ///@see Steel1DFiber::oldTangentStiffness.
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
		//Steel properties:
        double nu;
		double fy; 
		double E;
		double b;
		double R0;
		double cR1;
		double cR2;
		double a1;
		double a2;
		double a3;
		double a4;
        double Rho;

		//Steel Commited history variables:
		double newStrain;
		double newStress;
		double newTangentStiffness;
		double newMinStrain;
		double newMaxStrain;
		double newPlasticStrain;
		double newStrainInterception;
		double newStressInterception;
		double newStrainLastInversion;
		double newStressLastInversion;
		int newLoadUnloadIndex;

		//Steel Updated history variables:
		double oldStrain;
		double oldStress;
		double oldTangentStiffness;
		double oldMinStrain;
		double oldMaxStrain;
		double oldPlasticStrain;
		double oldStrainInterception;
		double oldStressInterception;
		double oldStrainLastInversion;
		double oldStressLastInversion;
		int oldLoadUnloadIndex;
};

#endif
