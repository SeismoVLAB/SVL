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
///This file contains the abstract "Material object" declarations, which computes 
///the strain, stress, and material tangent stiffness matrix.
//------------------------------------------------------------------------------

#ifndef _MATERIAL_HPP_
#define _MATERIAL_HPP_

#include <string>
#include <memory>
#include <Eigen/Dense>

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      January 5, 2018
/// @version   1.0
/// @file      Material.hpp
/// @class     Material
/// @see       Mesh.hpp
/// @brief     Virtual class for creating a material object
class Material{

    public:
        ///Creates a Material to be specified at a Gauss-point in an Element.
        ///@param name Name of the derived class.
        ///@param model If the material is viscous.
        ///@see Material::Name Material::isViscous.
        Material(std::string name, bool model);

        ///Destroys this Material object.
        virtual ~Material() = 0;

        ///Clone the selected material.
        ///@return A Material pointer to the derived class.
        ///@see Material.
        virtual std::unique_ptr<Material> CopyMaterial() = 0;

        ///Access material density.
        ///@return The material density value.
        virtual double GetDensity() const = 0;

        ///Returns the Poisson's ratio.
        ///@return The material Poisson's ratio (nu) value.
        virtual double GetPoissonRatio() const = 0;

        ///Access bulk modulus.
        ///@return The material bulk's (K) moduli value.
        virtual double GetBulkModulus() const = 0;

        ///Access shear modulus.
        ///@return The material shear's (G) moduli value.
        virtual double GetShearModulus() const = 0;

        ///Access modulus of elasticity.
        ///@return The material elasticity's (E) moduli value.
        virtual double GetElasticityModulus() const = 0;

        ///Computes the material energy for a given strain.
        ///@return Scalar with the material energy value.
        virtual double GetEnergy() const = 0;

        ///Returns the material viscous damping.
        ///@return Vector with the material damping components.
        virtual Eigen::MatrixXd GetDamping() const = 0;

        ///Computes the material strain in Voight notation.
        ///@return Vector with the material strain components.
        virtual Eigen::VectorXd GetStrain() const = 0;

        ///Computes the material stress in Voight notation.
        ///@return Vector with the material stress components.
        virtual Eigen::VectorXd GetStress() const = 0;

        ///Computes the material strain rate in Voight notation.
        ///@return Vector with the material strain-rate components.
        virtual Eigen::VectorXd GetStrainRate() const = 0;

        ///Computes the material total stress in Voight notation.
        ///@return Vector with the material total stress components.
        virtual Eigen::VectorXd GetTotalStress() const = 0;

        ///Computes the material stiffness in Voight notation.
        ///@return Matrix with the material tangent stiffness matrix.
        virtual Eigen::MatrixXd GetTangentStiffness() const = 0;

        ///Computes the initial material stiffness in Voight notation.
        ///@return Matrix with the initial material tangent stiffness matrix.
        ///@note The initial tangent stiffness matrix is computed when the strain vector is zero.
        virtual Eigen::MatrixXd GetInitialTangentStiffness() const = 0;

        ///Perform converged material state update.
        ///@note This function sets the trail stress and strain as converged.
        virtual void CommitState() = 0;

        ///Reverse the material states to previous converged state.
        ///@note This function returns the material states to previous converged states.
        virtual void ReverseState() = 0;

        ///Brings the material states to its initial state in the element.
        ///@note This function returns the material states to the beginning.
        virtual void InitialState() = 0;

        ///Update the material state for this iteration.
        ///@param strain Vector with the strain components at this Gauss-point.
        ///@param cond If the the elastic/plastic material components will be updated.
        ///@note This function computes the strain and tanget stiffness matrix once the trial strain converged.
        virtual void UpdateState(const Eigen::VectorXd strain, const unsigned int cond=0) = 0;

        ///Gets material information.
        ///@return The name of the derived class.
        ///@see Material::Name.
        std::string GetName();

        ///Gets material stress model.
        ///@return If the material is viscous.
        ///@see Material::isViscous.
        bool IsViscous();

    private:
        ///Material Name.
        std::string Name;

        ///Material Stress Model.
        bool isViscous;
};

#endif
