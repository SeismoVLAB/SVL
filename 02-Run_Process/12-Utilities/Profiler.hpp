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
//   [1] https://gist.github.com/TheCherno/Instrumentor.h
//
// Description:
// This file contains the "Profiler object" declarations to compute how long it 
// takes for a program to be executed, and sequentially measures the time of 
// Execution of each function.
//------------------------------------------------------------------------------

#ifndef _PROFILER_HHP_
#define _PROFILER_HHP_

#include <omp.h>
#include <string>
#include <thread>
#include <chrono>
#include <algorithm>
#include <fstream>
#include <sstream>

#include "Definitions.hpp"

///Structure that holds profiling parameters.
struct ProfileResult{
    ///Name of the function to be profiled.
    std::string Name;

    ///Time when Profiler starts.
    long long Start;

    ///Time when Profiler ends.
    long long End;

    ///The identifier of this thread.
    int ThreadID;
};

///Structure that holds function name to be profile.
struct InstrumentationSession{
    ///Name of the function to be profiled.
    std::string Name;
};

///This file contains the "Profiler object" declarations to compute how long it 
///takes for a program to be executed, and sequentially measures the time of 
///Execution of each function.

/// @author    The Cherno (https://www.patreon.com/thecherno)
/// @date      July 2, 2018
/// @version   1.0
/// @file      Profiler.hpp
/// @class     Profiler
/// @see       main.cpp
/// @brief     Class that creates a JSON file with the time of execution of each function, files can be open at chrome://tracing/
class Profiler{

    public:
        ///Creates the Profiler object.
        Profiler() : m_CurrentSession(nullptr), m_ProfileCount(0){
        }

        ///Starts the Profiler output file.
        ///@param name Name o the function to be timed.
        void BeginSession(const std::string& name){
            std::stringstream filename;
            filename << filePath << "/Profiler." << rank << ".json";
            m_OutputStream.open(filename.str().c_str());
            WriteHeader();
            m_CurrentSession = new InstrumentationSession{ name };
        }

        ///Ends the Profiler output file.
        void EndSession(){
            WriteFooter();
            m_OutputStream.close();
            delete m_CurrentSession;
            m_CurrentSession = nullptr;
            m_ProfileCount = 0;
        }

        ///Starts the Profiler output file.
        ///@param result Structure that contains the profiler results.
        void WriteProfile(const ProfileResult& result){
            if (m_ProfileCount++ > 0)
                m_OutputStream << ",";

            std::string name = result.Name;
            std::replace(name.begin(), name.end(), '"', '\'');

            //The JSON file format to be opened with Google-Chrome 
            m_OutputStream << "{";
            m_OutputStream << "\"cat\":\"function\",";
            m_OutputStream << "\"dur\":" << (result.End - result.Start) << ',';
            m_OutputStream << "\"name\":\"" << name << "\",";
            m_OutputStream << "\"ph\":\"X\",";
            m_OutputStream << "\"pid\":\"" << result.ThreadID << "\",";
            m_OutputStream << "\"tid\":"   << rank << ",";
            m_OutputStream << "\"ts\":"    << result.Start;
            m_OutputStream << "}";

            //Flush in case the program abruptly stops.
            m_OutputStream.flush();
        }

        ///Writes the header file for the JSON-file.
        void WriteHeader(){
            m_OutputStream << "{\"otherData\": {},\"traceEvents\":[";
            m_OutputStream.flush();
        }

        ///Writes the footer file for the JSON-file.
        void WriteFooter(){
            m_OutputStream << "]}";
            m_OutputStream.flush();
        }

        ///Gets the Profiler instance.
        static Profiler& Get(){
            static Profiler instance;
            return instance;
        }

    private:
        ///The name of the function to be Profile.
        InstrumentationSession* m_CurrentSession;

        ///The output JSON file where data is written.
        std::ofstream m_OutputStream;

        ///Counter for array separator.
        int m_ProfileCount;
};

/// @author    The Cherno (https://www.patreon.com/thecherno)
/// @date      July 2, 2018
/// @version   1.0
/// @file      Profiler.hpp
/// @class     Timer
/// @see       
/// @brief     Class that compute how long it takes for a program to be executed, and sequentially measures the time of execution of each function.
class Timer{

    public:
        ///Creates the Scope Timer.
        ///@param name Name o the function to be timed.
        Timer(const char* name) : m_Name(name), m_Stopped(false){
            std::this_thread::sleep_for (std::chrono::microseconds(3));
            m_StartTimepoint = std::chrono::high_resolution_clock::now();
        }

        ///Destroys this timer.
        ~Timer(){
            if (!m_Stopped)
                Stop();
        }

        ///Stops this timer.
        void Stop(){
            auto endTimepoint = std::chrono::high_resolution_clock::now();

            long long start = std::chrono::time_point_cast<std::chrono::microseconds>(m_StartTimepoint).time_since_epoch().count();
            long long end = std::chrono::time_point_cast<std::chrono::microseconds>(endTimepoint).time_since_epoch().count();
            int threadID = omp_get_num_threads();
            Profiler::Get().WriteProfile({ m_Name, start, end, threadID });

            m_Stopped = true;
        }

    private:
        ///Name of the function to be timed.
        const char* m_Name;

        ///Time when the timer starts.
        std::chrono::time_point<std::chrono::high_resolution_clock> m_StartTimepoint;

        ///Whether or not the timer is stoped.
        bool m_Stopped;
};

#if PROFILING
///Macro to compute time for profiler
///@param name Name of the function to be timed.
#define PROFILE_SCOPE(name) Timer timer(name)
///Macro to profile a function
///@param __FUNCTION__ Name macro of the function to be timed.
#define PROFILE_FUNCTION() PROFILE_SCOPE( __PRETTY_FUNCTION__ )
#else
#define PROFILE_FUNCTION()
#endif

#endif
