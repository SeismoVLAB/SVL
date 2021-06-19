#include "Load.hpp"
#include "Analysis.hpp"
#include "Definitions.hpp"
#include "Profiler.hpp"

//Default constructor.
Analysis::Analysis(std::shared_ptr<LoadCombo> &loadcombo){
    theLoadCombo = loadcombo;
}

//Virtual destructor.
Analysis::~Analysis(){
    //Does nothing.
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
Analysis::ReducedParallelReaction(Eigen::VectorXd &Reaction){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Pointer to Eigen residual vector.
    double *vector = new double[numberOfTotalDofs];

    //Pointer to add residual contribution for each processor.
    double *reduced = new double[numberOfTotalDofs];

    //Transform Eigen-vector into pointer.
    Eigen::Map<Eigen::VectorXd>(vector, numberOfTotalDofs) = Reaction;

    //Adds the residual force contribution from each processor.
    MPI_Allreduce(vector, reduced, numberOfTotalDofs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    
    //Transform the pointer into Eigen-vector.
    for(unsigned int k = 0; k < numberOfTotalDofs; k++)
        Reaction(k) = reduced[k];

    //Erase auxiliary variables. 
    delete[] reduced;
    delete[] vector;
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