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
///This file contains the "Lobatto Quadrature" declarations to integrate quantities 
///in the iso-parametric elements.
//------------------------------------------------------------------------------

#ifndef _LOBATTOQUADRATURE_HPP_
#define _LOBATTOQUADRATURE_HPP_

#include <string>
#include <Eigen/Dense>

#include "QuadratureRule.hpp"

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      November 7, 2020
/// @version   1.0
/// @file      LobattoQuadrature.hpp
/// @class     LobattoQuadrature
/// @see       Element.hpp QuadratureRule.hpp
/// @brief     Class for creating a Gauss-Lobatto Quadrature declarations to integrate quantities in the iso-parametric elements
class LobattoQuadrature : public QuadratureRule {

    public:
        ///Creates a LobattoQuadrature for Element integration.
        ///@param name Generic name of the Element for this quadrature rule.
        ///@param nQp Number of quadrature points.
        ///@note More details can be found at @ref linkLobattoQuadrature.
        ///@see LobattoQuadrature::nQuadraturePoints, LobattoQuadrature::nOrderQuadrature.
        LobattoQuadrature(std::string name, unsigned int nQp);

        ///Destroys this LobattoQuadrature object.
        ~LobattoQuadrature();

        ///Gets Number of Gauss Integration Points.
        ///@return The number of Gauss Integration Points.
        unsigned int GetNumberOfQuadraturePoints();

        ///Gets Gauss Integration Points.
        ///@param name Generic name of the Element for this quadrature rule.
        ///@param wi The Gauss quadratue weights.
        ///@param xi The Gauss quadratue coordinates.
        void GetQuadraturePoints(std::string name, Eigen::VectorXd &wi, Eigen::MatrixXd &xi);

    private:
        ///Number of Quadrature Points.
        unsigned int nQuadraturePoints;

        ///Number of Surface Quadrature Points.
        unsigned int nOrderQuadrature;

        ///Sets Gauss Quadrature for Line Elements.
        ///@param wi The Gauss quadratue weights.
        ///@param xi The Gauss quadratue coordinates.
        ///@note This function returns the Hexahedron quadrature information, see linkLobattoQuadrature.
        void SetLineQuadraturePoints(Eigen::VectorXd &wi, Eigen::MatrixXd &xi);

        ///Sets Gauss Quadrature for Triangular Elements.
        ///@param wi The Gauss quadratue weights.
        ///@param xi The Gauss quadratue coordinates.
        ///@note This function returns the Hexahedron quadrature information, see linkLobattoQuadrature.
        void SetTriaQuadraturePoints(Eigen::VectorXd &wi, Eigen::MatrixXd &xi);

        ///Sets Gauss Quadrature for Quadrilateral Elements.
        ///@param wi The Gauss quadratue weights.
        ///@param xi The Gauss quadratue coordinates.
        ///@note This function returns the Hexahedron quadrature information, see linkLobattoQuadrature.
        void SetQuadQuadraturePoints(Eigen::VectorXd &wi, Eigen::MatrixXd &xi);

        ///Sets Gauss Quadrature for Tetrahedron Elements.
        ///@param wi The Gauss quadratue weights.
        ///@param xi The Gauss quadratue coordinates.
        ///@note This function returns the Hexahedron quadrature information, see linkLobattoQuadrature.
        void SetTetraQuadraturePoints(Eigen::VectorXd &wi, Eigen::MatrixXd &xi);

        ///Sets Gauss Quadrature for Hexahedron Elements.
        ///@param wi The Gauss quadratue weights.
        ///@param xi The Gauss quadratue coordinates.
        ///@note This function returns the Hexahedron quadrature information, see linkLobattoQuadrature.
        void SetHexaQuadraturePoints(Eigen::VectorXd &wi, Eigen::MatrixXd &xi);
};

#endif
