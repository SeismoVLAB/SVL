#include "StaticAnalysis.hpp"
#include "Definitions.hpp"
#include "Profiler.hpp"

//Default constructor.
StaticAnalysis::StaticAnalysis(std::shared_ptr<Mesh> &mesh, std::shared_ptr<Algorithm> &algorithm, std::shared_ptr<Integrator> &integrator, std::shared_ptr<LoadCombo> &loadcombo, unsigned int nSteps) :
Analysis(loadcombo), theMesh(mesh), NumberOfSteps(nSteps){
    //Moves the linear system algorithm to this analysis.
    theAlgorithm = std::move(algorithm); 

    //Moves the static integrator to this analysis.
    theIntegrator = std::move(integrator); 

    //Sets the load combination to be used.
    theIntegrator->SetLoadCombination(loadcombo);
}

//Virtual destructor.
StaticAnalysis::~StaticAnalysis(){
    //Does nothing.
}

//Perform the analysis in incremental steps.
bool 
StaticAnalysis::Analyze(){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Initialize all recorders of this analysis.
    StartRecorders(theMesh, NumberOfSteps);

    //Perform the incremental anlaysis in nStep increments.
    for(unsigned int k = 0; k < NumberOfSteps; k++){
        //Solves until converge at k-th time steps.
        bool stop = theIntegrator->ComputeNewStep(theMesh, k); 

        //Checks the solution has issues.
        if(stop){
            EndRecorders();
            return stop;
        }

        //Updates the domain state.
        UpdateDomain(k); 

        //Store information in recorders.
        WriteRecorders(theMesh, k);

        //Prints the analysis progress.
        PrintProgressBar((k+1)*100/NumberOfSteps);
    }

    //Finilize all recorders of this analysis.
    EndRecorders();

    //Return the analysis status.
    return false;
}

//Update state variables and reset incremental state in the mesh.
void 
StaticAnalysis::UpdateDomain(unsigned int k){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Obtains the total current response state.  
    Eigen::VectorXd U = theIntegrator->GetDisplacements();

    //Gets all nodes from the mesh.
    std::map<unsigned int, std::shared_ptr<Node> >  Nodes = theMesh->GetNodes();    

    //TODO: Use OpenMP to accelerate here.
    for(auto it : Nodes){
        auto &Tag = it.first;

        //Gets the associated nodal degree-of-freedom.
        std::vector<int> TotalDofs = Nodes[Tag]->GetTotalDegreeOfFreedom();
        unsigned int nTotalDofs = TotalDofs.size();

        //Creates the nodal/incremental vector state.
        Eigen::VectorXd Uij(nTotalDofs); Uij.fill(0.0);
        Eigen::VectorXd dUij(nTotalDofs); dUij.fill(0.0);

        for(unsigned int j = 0; j < nTotalDofs; j++){
            Uij(j) = U(TotalDofs[j]);
        }
    
        //Sets the nodal states.
        Nodes[Tag]->SetDisplacements(Uij);
        Nodes[Tag]->SetIncrementalDisplacements(dUij);
    }

    //Gets all elements from the mesh.
    std::map<unsigned int, std::shared_ptr<Element> > Elements = theMesh->GetElements();

    for(auto it : Elements){
        auto &Tag = it.first;

        //Saves the material states for each element.
        Elements[Tag]->CommitState();
    }

    //Transform to total degree of freedom the displacement vector.
    Eigen::VectorXd R = theIntegrator->ComputeReactionForce(theMesh, k);
    ReducedParallelReaction(R);

    for(auto it : Nodes){
        auto &Tag = it.first;

        //Only if the node is fixed.
        if(Nodes[Tag]->IsFixed()){
            //Gets the associated nodal degree-of-freedom.
            std::vector<int> TotalDofs = Nodes[Tag]->GetTotalDegreeOfFreedom();
            unsigned int nTotalDofs = TotalDofs.size();

            //Creates the nodal/incremental vector state.
            Eigen::VectorXd Rij(nTotalDofs); Rij.fill(0.0);

            for(unsigned int j = 0; j < nTotalDofs; j++)
                Rij(j) = R(TotalDofs[j]);
                
            //Sets the nodal reaction.
            Nodes[Tag]->SetReaction(Rij);  
        }
    }
}
