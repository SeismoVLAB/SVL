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
///This file contains the "Recorder object" declarations, for which Node, Element,
///Section responses are supported.
//------------------------------------------------------------------------------

#ifndef _RECORDER_HPP_
#define _RECORDER_HPP_

#include <map> 
#include <vector>
#include <memory>
#include <string>
#include <fstream> 
#include <Eigen/Dense>

#include "Mesh.hpp"

/// @author    Danilo S. Kusanovic (dkusanov@caltech.edu)
/// @date      March 11, 2019
/// @version   1.0
/// @file      Recorder.hpp
/// @class     Recorder
/// @see       Node.hpp Element.hpp Section.hpp
/// @brief     Class that stores the node, element, and section responses for a given analysis
class Recorder{

    public:
        ///Creates a Recorder object to store Paraview solutions.
        ///@param file The file name where the results are going to be stored.
        ///@param name The object name where the response is computed.
        ///@param nparaview The number of features (for elements) to be written in Paraview.
        ///@param nsample The number of sampling points to record the solution.
        ///@param precision The precision number to store the results.
        ///@note More details can be found at @ref linkRecorder.
        ///@see Recorder::Path, Recorder::Name, Recorder::Response, Recorder::IDs.
        Recorder(std::string file, std::string name, unsigned int nparaview, unsigned int nsample = 1, unsigned int precision = 10);

        ///Creates a Recorder object to store Node / Element solutions.
        ///@param file The file name where the results are going to be stored.
        ///@param name The object name where the response is computed.
        ///@param type The response to be stored.
        ///@param index The object identifiers that are going to be stored.
        ///@param nsample The number of sampling points to record the solution.
        ///@param precision The precision number to store the results.
        ///@note More details can be found at @ref linkRecorder.
        ///@see Recorder::Path, Recorder::Name, Recorder::Response, Recorder::IDs.
        Recorder(std::string file, std::string name, std::string type, std::vector<unsigned int> index, unsigned int nsample = 1, unsigned int precision = 10);

        ///Creates a Recorder object to store Section solutions.
        ///@param file The file name where the results are going to be stored.
        ///@param name The object name where the response is computed.
        ///@param type The response to be stored.
        ///@param coordinates The position where the section response will be computed.
        ///@param index The object identifiers that are going to be stored.
        ///@param nsample The number of sampling points to record the solution.
        ///@param precision The precision number to store the results.
        ///@note More details can be found at @ref linkRecorder.
        ///@see Recorder::Path, Recorder::Name, Recorder::Response, Recorder::IDs.
        Recorder(std::string file, std::string name, std::string type, std::vector<double> coordinates, std::vector<unsigned int> index, unsigned int nsample = 1, unsigned int precision = 10);

        ///Destroys this Recorder object.
        ~Recorder();

        ///Clone the 'Recorder' object.
        ///@return A Recorder pointer to store solution.
        std::unique_ptr<Recorder> CopyRecorder();

        ///Initialize the recorder.
        ///@param mesh The finite element Mesh object.  
        ///@param nsteps The number of time step of this simulation.  
        void Initialize(std::shared_ptr<Mesh> &mesh, unsigned int nsteps);

        ///Write information in the recorder.
        ///@param mesh The finite element Mesh object.
        ///@note More details can be found at @ref linkRecorder.
        void WriteResponse(std::shared_ptr<Mesh> &mesh, unsigned int step);

        ///Finalize the recorder.
        void Finalize();

        ///Sets the combination's name.
        ///@param name The LoadCombo name.  
        void SetComboName(std::string name);
        
    private:
        ///The name of the file to record.    
        std::string File;

        ///The object name to record solution.    
        std::string Name;

        ///The name of the LoadCombo.    
        std::string Combo;

        ///The object response to be stored.    
        std::string Response;

        ///The number of sampling points to record the solution.  
        unsigned int nSample;

        ///The sampling rate counter.    
        unsigned int Counter;

        ///The number of features to be written in paraview.
        unsigned int nParaview;

        ///The precision for results.    
        unsigned int Precision;

        ///Section position at strain/stress is computed.
        std::vector<double> Position; 

        ///Nodal/Element indexes to be recorded.
        std::vector<unsigned int> IDs;

        ///Nodal Local indexes for this partition.
        std::map<unsigned int, unsigned int> Tag;

        ///Name of the recorder handler.
        std::ofstream OutputFile;

        ///Writes state variable data to the file.
        ///@param mesh The finite element Mesh object.  
        void WriteNodalResponse(std::shared_ptr<Mesh> &mesh);

        ///Writes element data to the file.
        ///@param mesh The finite element Mesh object.  
        void WriteElementResponse(std::shared_ptr<Mesh> &mesh);

        ///Writes section element data to the file.
        ///@param mesh The finite element Mesh object.  
        void WriteSectionResponse(std::shared_ptr<Mesh> &mesh);

        ///Writes output data in VTK format to the file.
        ///@param mesh The finite element Mesh object.  
        ///@param step The time step to be stored.
        void WriteVTKFiles(std::shared_ptr<Mesh> &mesh, unsigned int step);

        ///Fix blank spaces provided by user in path.
        ///@param theFile String with full path.
        ///@param toReplace The string pattern to be replaced with. 
        ///@return String with the replaced pattern.
        std::string GetSpacedName(std::string theFile, std::string toReplace);
};

#endif