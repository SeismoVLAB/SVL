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
///This file contains the abstract "Quadrature Rule" declarations to integrate 
///quantities in the iso-parametric elements. 
//------------------------------------------------------------------------------

#ifndef _QUADRATURERULE_HPP_
#define _QUADRATURERULE_HPP_

#include <string>
#include <Eigen/Dense>

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      March 18, 2018
/// @version   1.0
/// @file      QuadratureRule.hpp
/// @class     QuadratureRule
/// @see       Element.hpp GaussQuadrature.hpp LobattoQuadrature.hpp
/// @brief     Virtual class for defining a Quadrature Rule to integrate quantities in the iso-parametric elements
class QuadratureRule {

    public:
        ///Creates a QuadratureRule for Element integration.
        ///@param name Generic name of the Element for this quadrature rule.
        ///@param nQp Number of quadrature points.
        ///@note More details can be found at @ref linkQuadratureRule.
        ///@see QuadratureRule::Name.
        QuadratureRule(std::string name);

        ///Destroys this QuadratureRule object.
        virtual ~QuadratureRule() = 0;

        ///Gets Number of Integration Points.
        ///@return The number of Quadrature Points.
        virtual unsigned int GetNumberOfQuadraturePoints() = 0;

        ///Gets Gauss Integration Points.
        ///@param name Generic name of the Element for this quadrature rule.
        ///@param Weights The Gauss quadratue weights.
        ///@param Points The Gauss quadratue coordinates.
        virtual void GetQuadraturePoints(std::string name, Eigen::VectorXd &Weights, Eigen::MatrixXd &Points) = 0;

        ///Gets the Element's Quadrature Name. 
        std::string GetQuadratureName();

    private:
        ///Element's Quadrature Name.
        std::string Name;
};

#endif
