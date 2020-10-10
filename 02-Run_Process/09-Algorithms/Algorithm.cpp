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

//Construct the residual vector force from each processor.
void
Algorithm::ReducedParallelResidual(const Eigen::VectorXd &Feff, double &Residual){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Pointer to Eigen residual vector.
    double *vec = new double[numberOfFreeDofs];

    //Pointer to add residual contribution for each processor.
    double *red = new double[numberOfFreeDofs];

    //Transform Eigen-vector into pointer.
    Eigen::Map<Eigen::VectorXd>(vec, numberOfFreeDofs) = Feff;

    //Adds the residual force contribution from each processor.
    MPI_Reduce(vec, red, numberOfFreeDofs, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    
    //Transform the pointer into Eigen-vector.
    Eigen::Map<Eigen::VectorXd> feff(red,numberOfFreeDofs);

    //Computes the norm of the residual force verctor.
    Residual = feff.norm();

    //Pass the residual error to each partition.
    MPI_Bcast(&Residual, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    //Erase auxiiliar variables. 
    delete[] red;
    delete[] vec;
}

//Computes convergence tests for this algorithm.
double 
Algorithm::ComputeConvergence(const Eigen::VectorXd &Force, double Delta, bool isFirstIteration){
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

    double Residual = 0.0;

    //Convergence tests possibilities for this algorithm.
    if(flag == 1){
        //Unbalanced Force Norm.
        if(!isFirstIteration)
            ReducedParallelResidual(Force, Residual);
    }
    else if(flag == 2){
        //Increment Displacement Norm.
        if(isFirstIteration){
            NormFactor = Delta;
        }
        else{
            Residual = Delta - NormFactor;
            NormFactor = Delta;
        }
    }
    else if(flag == 3){
        //Relative Unbalanced Force Norm. 
        if(isFirstIteration){
            ReducedParallelResidual(Force, NormFactor);
        }
        else{
            ReducedParallelResidual(Force, Residual); 
            Residual = Residual/NormFactor;
        }
    }
    else if(flag == 4){
        //Relative Increment Displacement Norm.
        if(isFirstIteration){
            NormFactor = Delta;
        }
        else{
            Residual = (Delta - NormFactor)/Delta;
            NormFactor = Delta;
        }
    }
    else{
        //TODO: Set a warning that Convergence Test is not defined.
    }

    return Residual;
}
