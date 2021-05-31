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
//   [1] Carlos Felippa. "Introduction to Finite Element Methods (ASEN 5007)": 
//       Multifreedom Constraints I, Fall 2005.
//
// Description:
///This file contains the "Constraint object" declarations, which impose 
///linear kinematic relations between degree-of-freedom in the model.
//------------------------------------------------------------------------------

#ifndef _CONSTRAINT_HPP_
#define _CONSTRAINT_HPP_

#include <vector>

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      March 15, 2018
/// @version   1.0
/// @file      Constraint.hpp
/// @class     Constraint
/// @see       Node.hpp Mesh.hpp
/// @brief     Class for creating a linear kinematic relations between degree-of-freedom in the model
class Constraint{

    public:
        ///Creates a empty Constraint to a degree-of-freedom.
        Constraint();

        ///Creates a Constraint to be applied to a degree-of-freedom.
        ///@param slave 'Total' degree-of-freedom number to which the constraint will be applied.
        ///@param master List of 'free' degree-of-freedom a.k.a master (or combinational) numbering.
        ///@param factors List of combinational values.
        ///@note More details can be found at @ref linkConstraint.
        ///@see Constraint::Slave, Constraint::Master, Constraint::Coefficients.
        Constraint(unsigned int slave, std::vector<unsigned int> master, std::vector<double> factors);

        ///Destroys this Constraint.
        ~Constraint();

        ///Sets the slave's node degree-of-freedom.
        ///@param slave 'Total' degree-of-freedom number to which the constraint will be applied.
        ///@note More details can be found at @ref linkConstraint.
        ///@see Constraint::Slave.
        void SetSlaveInformation(unsigned int slave);

        ///Sets the master's free degree-of-freedom list numbering.
        ///@param master List of 'free' degree-of-freedom a.k.a master (or combinational) numbering.
        ///@note More details can be found at @ref linkConstraint.
        ///@see Constraint::Master.
        void SetMasterInformation(std::vector<unsigned int> master);

        ///Sets the master's combinational factors applied to degree-of-freedom list numbering.
        ///@param factors List of combinational values.
        ///@note More details can be found at @ref linkConstraint.
        ///@see Constraint::Coefficients.
        void SetCombinationFactors(std::vector<double> factors);

        ///Gets the number of combinations applied to this slave degree-of-freedom.
        ///@return Number of combinational factors.
        unsigned int GetNumberOfConstraints() const;

        ///Gets the slave free degree-of-freedom of this constraint.
        ///@return The slave total degree-of-freedom number.
        ///@see Constraint::Slave.
        unsigned int GetSlaveInformation() const;

        ///Gets the master list of total degree-of-freedom to be combined.
        ///@return List of 'free' degree-of-freedom a.k.a master (or combinational) numbering.
        ///@note More details can be found at @ref linkConstraint.
        ///@see Constraint::Master.
        const std::vector<unsigned int> GetMasterInformation() const;

        ///Gets the combinational factor list for each degree-of-freedom.
        ///@return The coefficient values applied to the master nodes.
        ///@note More details can be found at @ref linkConstraint.
        ///@see Constraint::Coefficients.
        const std::vector<double> GetCombinationFactors() const;

    private:
        ///The slave total degree-of-freedom.
        unsigned int Slave;

        ///The list of master free degree-of-freedom.
        std::vector<unsigned int> Master;

        ///The combination coefficient vector.
        std::vector<double> Coefficients;
};

#endif