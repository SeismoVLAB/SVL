#include "Algorithm.hpp"
#include "Definitions.hpp"
#include "Profiler.hpp"

//Defaul constructor.
Algorithm::Algorithm(const std::shared_ptr<Mesh> &mesh, unsigned int flag, double NormFactor) : flag(flag), NormFactor(NormFactor){
    //Operator that enforced restrain/constraint. 
    Total2FreeMatrix = mesh->GetTotalToFreeMatrix();
}

//Virtual destructor.
Algorithm::~Algorithm(){
    //Does nothing.
}

//Update the incremental state variables in the mesh.
void 
Algorithm::UpdateStatesIncrements(std::shared_ptr<Mesh> &mesh, const Eigen::VectorXd &DeltaU){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Update the incremental state variables in the mesh.
    Eigen::VectorXd dU = Total2FreeMatrix*DeltaU;

    //Gets all nodes from the mesh.
    std::map<unsigned int, std::shared_ptr<Node> >  Nodes = mesh->GetNodes();    

    //TODO: Use OpenMP to accelerate here.
    for(auto it : Nodes){
        auto &Tag = it.first;

        //Gets the associated nodal degree-of-freedom.
        std::vector<int> TotalDofs = Nodes[Tag]->GetTotalDegreeOfFreedom();
        unsigned int nTotalDofs = TotalDofs.size();

        //Creates the nodal/incremental vector state.
        Eigen::VectorXd dUij(nTotalDofs); dUij.fill(0.0);

        for(unsigned int j = 0; j < nTotalDofs; j++){
            dUij(j) = dU(TotalDofs[j]);
        }    
        
        //Sets the nodal incremental states.
        Nodes[Tag]->SetIncrementalDisplacements(dUij);
    }

    //Gets all elements from the mesh.
    std::map<unsigned int, std::shared_ptr<Element> > Elements = mesh->GetElements();

    for(auto it : Elements){
        auto &Tag = it.first;

        //Updates the material states for each element.
        Elements[Tag]->UpdateState();
    }
}

//Revert to the previous converged state in the mesh.
void 
Algorithm::ReverseStatesIncrements(std::shared_ptr<Mesh> &mesh){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Gets all nodes from the mesh.
    std::map<unsigned int, std::shared_ptr<Node> >  Nodes = mesh->GetNodes();    

    //TODO: Use OpenMP to accelerate here.
    for(auto it : Nodes){
        auto &Tag = it.first;

        //Gets the associated nodal degree-of-freedom.
        unsigned int nTotalDofs = Nodes[Tag]->GetNumberOfDegreeOfFreedom();

        //Creates the nodal/incremental vector state.
        Eigen::VectorXd dUij(nTotalDofs); 
        dUij.fill(0.0);

        //Sets the nodal incremental states.
        Nodes[Tag]->SetIncrementalDisplacements(dUij);
    }

    //Gets all elements from the mesh.
    std::map<unsigned int, std::shared_ptr<Element> > Elements = mesh->GetElements();

    for(auto it : Elements){
        auto &Tag = it.first;

        //Updates the material states for each element.
        Elements[Tag]->ReverseState();
    }
}

//Construct the residual vector force from each processor.
void
Algorithm::ReducedParallelResidual(const Eigen::VectorXd &vector, double &Residual){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Pointer to add residual contribution for each processor.
    double *reduced = new double[numberOfFreeDofs];

    //Adds the residual force contribution from each processor.
    MPI_Reduce(vector.data(), reduced, numberOfFreeDofs, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    //Computes the norm of the residual force verctor.
    double residual = 0.0;
    for(unsigned int k = 0; k < numberOfFreeDofs; k++)
        residual += reduced[k]*reduced[k];
    Residual = sqrt(residual);

    //Pass the residual error to each partition.
    MPI_Bcast(&Residual, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    //Erase auxiliary variables. 
    delete[] reduced;
}

//Computes convergence tests for this algorithm.
double 
Algorithm::ComputeConvergence(const Eigen::VectorXd &Force, const Eigen::VectorXd &Delta, const Eigen::VectorXd &delta, unsigned int k){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    double Residual = 0.0;

    //Convergence tests possibilities for this algorithm.
    if(flag == 1){
        //Unbalanced Force Norm.
        ReducedParallelResidual(Force, Residual);  
    }
    else if(flag == 2){
        //Increment Displacement Norm.
        Residual = delta.norm();
    }
    else if(flag == 3){
        //Energy Increment Norm
        Eigen::VectorXd Energy = delta.cwiseProduct(Force);
        ReducedParallelResidual(Energy, Residual);
    }
    else if(flag == 4){
        //Relative Unbalanced Force Norm. 
        if(k != 0){
            ReducedParallelResidual(Force, Residual); 
            Residual = Residual/NormFactor;
        }
        else{
            ReducedParallelResidual(Force, NormFactor);
            Residual = NormFactor;
        }
    }
    else if(flag == 5){
        //Relative Increment Displacement Norm
        if(k != 0){
            Residual = delta.norm()/NormFactor;
        }
        else{
            NormFactor = Delta.norm();
            Residual = delta.norm();
        }
    }
    else if(flag == 6){
        //Relative Energy Increment Norm
        Eigen::VectorXd Energy = delta.cwiseProduct(Force);
        if(k != 0){
            ReducedParallelResidual(Energy, Residual);
            Residual = Residual/NormFactor;
        }
        else{
            ReducedParallelResidual(Energy, NormFactor);
            Residual = NormFactor;
        }
    }
    else if(flag == 7){
        //Total Relative Increment Displacement Norm
        Residual = delta.norm()/Delta.norm();
    }
    else if(flag == 8){
        //Maximum Number of Iterations
        Residual = -1.00;
    }
    else{
        //TODO: Set a warning that Convergence Test is not defined.
    }

    return Residual;
}
