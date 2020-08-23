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
//   [1]
//
// Description:
///This file contains the "Load Combination" class declarations, which defines 
///how the loads are going to be combined for an specific analysis.
//------------------------------------------------------------------------------

#ifndef _LOADCOMBO_HPP_
#define _LOADCOMBO_HPP_

#include <vector>
#include <string>

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      June 25, 2018
/// @version   1.0
/// @file      LoadCombo.hpp
/// @class     LoadCombo
/// @see       Load.hpp StaticAnalysis.hpp DynamicAnalysis.hpp
/// @brief     Class for defining how the loads are going to be combined for an specific analysis
class LoadCombo{

    public:
        ///Creates a LoadCombo object.
        ///@param name Name of the load combination object.
        ///@param loads The load identifiers to be combined.
        ///@param factors Numerical combinational values.
        ///@note More details can be found at @ref linkLoadCombination
        ///@see Load.
        LoadCombo(std::string name, std::vector<unsigned int> loads, std::vector<double> factors);

        ///Destroys this LoadCombo object.
        ~LoadCombo();

        ///Returns the loads to combine. 
        ///@return The name of the combination.
        ///@see LoadCombo::Name.
        std::string GetCombinationName() const;

        ///Returns the loads factors.
        ///@return List with the combinational factors.
        ///@note More details can be found at @ref linkLoadCombination
        ///@see LoadCombo::Factors.
        std::vector<double> GetLoadFactors() const;

        ///Returns the Load identifiers to be combined,
        ///@return List with the Load identifiers. 
        ///@note More details can be found at @ref linkLoadCombination
        ///@see LoadCombo::Loads.
        std::vector<unsigned int> GetLoadCombination() const;

    private:
        ///The combination's name.
        std::string Name;

        ///The loads to be combined.
        std::vector<unsigned int> Loads;

        ///The factors of combination.
        std::vector<double> Factors;
};

#endif
