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
//   [1] Mohd Hisham Mohd Yassin, "Nonlinear Analysis of Prestressed Concrete 
//       Structures under Monotonic and Cycling Loads", PhD dissertation, 
//       University of California, Berkeley, 1994.
//
// Description:
///This file contains the "Concrete1DFiber" material declarations, which defines 
///a concrete material for 1D, 2D, or 3D analysis to be used along with
///zero-length element or fiber section.
//------------------------------------------------------------------------------

#ifndef _CONCRETE1DFIBER_HPP_
#define _CONCRETE1DFIBER_HPP_

#include "Material.hpp"

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      April 17, 2021
/// @version   1.0
/// @file      Concrete1DFiber.hpp
/// @class     Concrete1DFiber
/// @see       Material.hpp
/// @brief     Class for creating an uniaxial concrete material compatible with fiber section
class Concrete1DFiber : public Material {

    public:
        ///Creates a Concrete1DFiber material to be specified at a Gauss-point in an Element.
        ///@param fc The concrete compressive strength at 28 days.
        ///@param ecc The concrete strain at maximum strength.
        ///@param fcu The concrete crushing strength.
        ///@param ecu The concrete strain at crushing strength.
        ///@param ratio Ratio between unloading slope at and initial slope.
        ///@param ft The concrete tensile strength.
        ///@param Et The tension softening stiffness (absolute value).
        ///@param nu Represents the Poisson's ratio.
        ///@param Rho Represents the material density.
        ///@see Concrete1DFiber::E, Concrete1DFiber::nu, Concrete1DFiber::Rho.
        Concrete1DFiber(double fc, double ecc, double fcu, double ecu, double ratio, double ft=0.0, double Et=0.0, double nu= 0.25, double rho=0.0);

        ///Destroys this Concrete1DFiber material.
        ~Concrete1DFiber();

        //Clone the 'Concrete1DFiber' material.
        ///@return A Material pointer to the Concrete1DFiber derived class.
        ///@see Material, Concrete1DFiber.
        std::unique_ptr<Material> CopyMaterial();

        ///Access material density.
        ///@return The material density value.
        ///@see Concrete1DFiber::Rho.
        double GetDensity() const;

        ///Returns the Poisson's ratio.
        ///@return The material poisson's ratio.
        ///@see Concrete1DFiber::nu.
        double GetPoissonRatio() const;

        ///Access bulk modulus.
        ///@return The material bulk's (K) modulus.
        double GetBulkModulus() const;

        ///Returns linear shear modulus.
        ///@return The material density value.
        double GetShearModulus() const;

        ///Access modulus of elasticity.
        ///@return The material Elasticity's (E) modulus.
        ///@see Concrete1DFiber::fc.
        double GetElasticityModulus() const;

        ///Access the material's energy at current strain.
        ///@return Scalar with the material energy value.
        double GetEnergy() const;

        ///Returns the material viscous damping.
        ///@return Vector with the material damping components.
        Eigen::MatrixXd GetDamping() const;

        ///Returns the material strain.
        ///@return Vector with the material strain components.
        ///@note The strain can be revisited in @ref linkConcrete1DFiber.
        ///@see Concrete1DFiber::Strain.
        Eigen::VectorXd GetStrain() const;

        ///Returns the material stress. 
        ///@return Vector with the material stress components.
        ///@note The stress can be revisited in @ref linkConcrete1DFiber.
        ///@see Concrete1DFiber::Stress.
        Eigen::VectorXd GetStress() const;

        ///Returns the material strain rate.
        ///@return Vector with the material strain-rate components.
        ///@note The stress-rate can be revisited in @ref linkConcrete1DFiber.
        Eigen::VectorXd GetStrainRate() const;

        ///Computes the material total stress.
        ///@return Vector with the material total stress components.
        ///@note The damping can be revisited in @ref linkConcrete1DFiber.
        ///@see Concrete1DFiber::Stress.
        Eigen::VectorXd GetTotalStress() const;

        ///Returns the material stiffness.
        ///@return Matrix with the material tangent stiffness matrix.
        ///@note The stiffness matrix can be revisited in @ref linkConcrete1DFiber.
        ///@see Concrete1DFiber::TangentStiffness.
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
		//Concrete properties:
        double nu;
		double Et;
		double fc; 
		double ft;
		double fcu;
		double ecu;
		double ecc;
        double Rho;
		double lambda;

		//Concrete history variables:
		double newStrain;
		double newStress;
		double newTangentStiffness;
		double newMinStrain;
		double newMaxStrain;

		double oldStrain;
		double oldStress;
		double oldTangentStiffness;
		double oldMinStrain;
		double oldMaxStrain;

		//Envelope stress functions:
		void Tensile_Envelope(double eccn, double &fccn, double &Etn);
		void Compressive_Envelope(double eccn, double &fccn, double &Ecn);
};

#endif