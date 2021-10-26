#include <iostream>
#include "Load.hpp"
#include "Analysis.hpp"
#include "Definitions.hpp"
#include "Profiler.hpp"

//Default constructor.
Analysis::Analysis(std::shared_ptr<LoadCombo> &loadcombo, unsigned int nteps) : NumberOfSteps(nteps){
    theLoadCombo = loadcombo;
}

//Virtual destructor.
Analysis::~Analysis(){
    //Does nothing.
}

//Update internal variables according to simulation specifications.
void 
Analysis::UpdateMesh(std::shared_ptr<Mesh> &mesh, std::shared_ptr<Integrator> &integrator){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    if(strcasecmp(UpdateOption.c_str(),"RESTARTABLE") == 0){
        //Gets element information from the mesh.
        std::map<unsigned int, std::shared_ptr<Node> > Nodes = mesh->GetNodes();

        //Gets element information from the mesh.
        std::map<unsigned int, std::shared_ptr<Element> > Elements = mesh->GetElements();

        //RESTART ANALYSIS TO ITS INITIAL STATE
        //Nodes internal variables are set to initial state
        for(auto it : Nodes){
            auto &Tag = it.first;
            Nodes[Tag]->InitialState();
        }

        //Elements internal variables are set to initial state
        for(auto it : Elements){
            auto &Tag = it.first;
            Elements[Tag]->InitialState();
        }
    }
    else if(strcasecmp(UpdateOption.c_str(),"PROGRESSIVE") == 0){
        //Gets element information from the mesh.
        std::map<unsigned int, std::shared_ptr<Node> > Nodes = mesh->GetNodes();

        Eigen::VectorXd F = integrator->ComputeProgressiveForce(mesh, NumberOfSteps-1);

        //SEQUENTIAL ANALYSIS WITH HOMOGENEOUS BOUNDARIES
        for(auto it : Nodes){
            auto &Tag = it.first;

            //Gets the associated nodal degree-of-freedom.
            std::vector<int> TotalDofs = Nodes[Tag]->GetTotalDegreeOfFreedom();
            unsigned int nTotalDofs = TotalDofs.size();

            //Creates the nodal vector state.
            Eigen::VectorXd Fj(nTotalDofs); 
            Fj.fill(0.0);

            for(unsigned int j = 0; j < nTotalDofs; j++){
                Fj(j) = F(TotalDofs[j]);
            }
                
            //Sets the nodal progressive force vector.
            Nodes[Tag]->SetProgressiveForces(Fj);   
        }
    }
    else if(strcasecmp(UpdateOption.c_str(),"TRANSMISSIVE") == 0){
        //Gets element information from the mesh.
        std::map<unsigned int, std::shared_ptr<Node> > Nodes = mesh->GetNodes();

        Eigen::VectorXd R = integrator->ComputeReactionForce(mesh, NumberOfSteps-1);
        Eigen::VectorXd F = integrator->ComputeProgressiveForce(mesh, NumberOfSteps-1);

        //SEQUENTIAL ANALYSIS WITH INHOMOGENEOUS BOUNDARIES
        for(auto it : Nodes){
            auto &Tag = it.first;

            //Gets the associated nodal degree-of-freedom.
            std::vector<int> FreeDofs = Nodes[Tag]->GetFreeDegreeOfFreedom();
            std::vector<int> TotalDofs = Nodes[Tag]->GetTotalDegreeOfFreedom();
            unsigned int nTotalDofs = TotalDofs.size();

            //Creates the nodal vector state.
            Eigen::VectorXd Fj(nTotalDofs); 
            Fj.fill(0.0);

            for(unsigned int j = 0; j < nTotalDofs; j++){
                Fj(j) += F(TotalDofs[j]);

                //Only restrained degree-of-freedom is added
                if(FreeDofs[j] == -1){
                    Fj(j) += R(TotalDofs[j]);
                }
            }
                
            //Sets the nodal progressive force vector.
            Nodes[Tag]->SetProgressiveForces(Fj);   
        }
    }
}

//Sets the recorders for the analysis
void 
Analysis::SetRecorder(std::shared_ptr<Recorder> &recorder){
    theRecorders.push_back(recorder->CopyRecorder());
}

//Returns the combination name
std::string 
Analysis::GetCombinationName(){
    return theLoadCombo->GetCombinationName();
}

//Construct the residual vector force from each processor.
void
Analysis::ReducedParallelReaction(Eigen::VectorXd &vector){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Pointer to add residual contribution for each processor.
    double *reduced = new double[numberOfTotalDofs];

    //Adds the residual force contribution from each processor.
    MPI_Allreduce(vector.data(), reduced, numberOfTotalDofs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    
    //Transform the pointer into Eigen-vector.
    for(unsigned int k = 0; k < numberOfTotalDofs; k++)
        vector(k) = reduced[k];

    //Erase auxiliary variables. 
    delete[] reduced;
}

//Prints progress bar in analysis.
void 
Analysis::PrintProgress(unsigned int percent){
    if(rank == 0){
        //Progress bar frame.
        std::cout << "\r RUNNING (" << GetCombinationName() << ") : " << percent << "%" << std::flush;
    }
}

//Initialize all recorders.
void 
Analysis::StartRecorders(std::shared_ptr<Mesh> &mesh, unsigned int nsteps){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Gets name of combination.
    std::string name = theLoadCombo->GetCombinationName();

    //Initialize the recorders.
    for(unsigned int k = 0; k < theRecorders.size(); k++){
        theRecorders[k]->SetComboName(name);
        theRecorders[k]->Initialize(mesh, nsteps);

        //PATCH: This section is to remove the high values of stress/strain generated in the DRM element interface.
        if(strcasecmp(theRecorders[k]->GetName().c_str(), "PARAVIEW") == 0){
            std::vector<unsigned int> lTags = theLoadCombo->GetLoadCombination();
            std::map<unsigned int, std::shared_ptr<Load> > Loads = mesh->GetLoads();

            //There is a DRM Load Pattern Applied
            for(unsigned int m = 0; m < lTags.size(); m++){
                if( Loads[lTags[m]]->GetClassification() == 7 ){
                    //The DRM Elements applied to this load pattern
                    std::vector<unsigned int> DRMElems = Loads[lTags[m]]->GetElements();

                    //Gets the element information from the mesh.
                    std::map<unsigned int, std::shared_ptr<Element> > Elements = mesh->GetElements();
                    std::map<unsigned int, bool> DRMConditions;
                    
                    for(auto it : Elements){
                        auto &eTag = it.first;
                        if( std::find(DRMElems.begin(), DRMElems.end(), eTag) != DRMElems.end() ){
                            DRMConditions[eTag] = true;
                        }
                        else{
                            DRMConditions[eTag] = false;
                        }
                    }
                    theRecorders[k]->SetDRMParaviewInterface(DRMConditions);
                }
            }           
        }
    }
}

//Writes information in plain text format.
void 
Analysis::WriteRecorders(std::shared_ptr<Mesh> &mesh, unsigned int step){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Goes over all recorders and write solution.
    for(unsigned int k = 0; k < theRecorders.size(); k++){
        theRecorders[k]->WriteResponse(mesh, step);
    }
}

//Finalize all recorders.
void 
Analysis::EndRecorders(){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Skips one line for next simulation.
    if(rank == 0)
        std::cout << "\n";

    //Finalize the recorders.
    for(unsigned int k = 0; k < theRecorders.size(); k++){
        theRecorders[k]->Finalize();
    }
}
