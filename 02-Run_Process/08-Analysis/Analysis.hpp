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
///This file contains the pure virtual "Analysis object" declarations, so far, 
///static and dynamic analysis are available. 
//------------------------------------------------------------------------------

#ifndef _ANALYSIS_HPP_
#define _ANALYSIS_HPP_

#include <mpi.h>
#include <ctime>
#include <memory>
#include <iostream>
#include <algorithm>

#include "Mesh.hpp"
#include "Recorder.hpp"
#include "LoadCombo.hpp"
#include "Integrator.hpp"

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      November 22, 2018
/// @version   1.0
/// @file      Analysis.hpp
/// @class     Analysis
/// @see       StaticAnalysis.hpp DynamicAnalysis.hpp
/// @brief     Virtual class for defining the analysis type to be performed
class Analysis{

    public:
        ///Creates a Analysis object.
        ///@param loadcombo Pointer to the LoadCombo object.
        ///@param nteps The number of time steps to evolve solution.
        Analysis(std::shared_ptr<LoadCombo> &loadcombo, unsigned int nteps);

        ///Destroys this Analysis object.
        virtual ~Analysis() = 0;

        ///Performs the required analysis on the domain.
        ///@return Whether or not the analysis was successful.
        virtual bool Analyze() = 0;

        ///Update internal variables in Mesh according to form of simulation.
        ///@param mesh The finite element mesh object where changes are performed
        ///@param integrator The integrator object where changed are taken from
        ///@see Mesh::Nodes Mesh::Elements.
        void UpdateMesh(std::shared_ptr<Mesh> &mesh, std::shared_ptr<Integrator> &integrator);

        ///Sets the recorder for the analysis.
        ///@param recorder Pointer to the recorder where solution is stored.
        ///@see Analysis::theRecorders.
        void SetRecorder(std::shared_ptr<Recorder> &recorder);

        ///Returns the combination name
        ///@return The analysis combination name.
        std::string GetCombinationName();

        ///Construct the reaction vector force from each processor.
        ///@param Reaction The reaction vector in this partition.
        ///@param numberOfTotalDofs The number of total degree-of-freedom.
        ///@note In parallel execution Reaction is gathered across processors.
        void ReducedParallelReaction(Eigen::VectorXd &Reaction);

    protected:
        ///Initialize recorder.
        ///@param mesh Pointer to the Mesh object.
        ///@param nsteps The number of time steps to be recorded.
        void StartRecorders(std::shared_ptr<Mesh> &mesh, unsigned int nsteps);

        ///Writes information on the recorders.
        ///@param mesh Pointer to the Mesh object.
        ///@param step The time step to be recorded.
        void WriteRecorders(std::shared_ptr<Mesh> &mesh, unsigned int step);

        ///Finalize recorder.
        void EndRecorders();

        ///Prints out solving bar for the analysis.
        ///@param percent The progress bar percentage in analysis.
        void PrintProgress(unsigned int percent);

    private:
        ///Total number of time increments.
        unsigned int NumberOfSteps;

        ///The load combination to be used.
        std::shared_ptr<LoadCombo> theLoadCombo;

        ///The recorder associated to analysis.
        std::vector<std::unique_ptr<Recorder> > theRecorders;
};

#endif
