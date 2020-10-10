#include "Mesh.hpp"
#include "Definitions.hpp"

typedef Eigen::Triplet<double> T;

//Overload constructor.
Mesh::Mesh() { 
    //Initialize global variables.
    LumpedStorage = 0;
    ConsistentStorage = 0;
    numberOfConstrainedDofs = 0;
}
 
//Destructor.
Mesh::~Mesh(){
    //Does nothing.
}

//Compute mapping and storage.
void 
Mesh::Initialize(){
    //Initialize the mesh objects.
    for(auto it : Elements){
        auto &Tag = it.first;

        //Sets the finite element dependance among objects.
        Elements[Tag]->SetDomain(Nodes);

        //Computes maximum memory allocation for sparse matrix.
        unsigned int numDofs = Elements[Tag]->GetNumberOfDegreeOfFreedom();
        LumpedStorage += numDofs;
        ConsistentStorage += numDofs*numDofs;
    }

    //Computes the number of constraints components.
    for(std::map<int, std::unique_ptr<Constraint> >::iterator it = Constraints.begin(); it != Constraints.end(); ++it){
        //Gets the node identifier.
        int Tag = it->first;
        numberOfConstrainedDofs += Constraints[Tag]->GetNumberOfConstraints();
    }
}

//Add node to the mesh.
void 
Mesh::AddNode(unsigned int tag, std::shared_ptr<Node> &node){
    Nodes[tag] = std::move(node);
}

//Add constraint to the mesh.
void 
Mesh::AddConstraint(unsigned int tag, std::unique_ptr<Constraint> &constraint){
    Constraints[tag] = std::move(constraint);
}

//Add a material to the mesh.
void 
Mesh::AddMaterial(unsigned int tag, std::unique_ptr<Material> &material){
    Materials[tag] = std::move(material);
}

//Add a material to the mesh.
void 
Mesh::AddSection(unsigned int tag, std::unique_ptr<Section> &section){
    Sections[tag] = std::move(section);
}

//Add classic damping to the mesh.
void
Mesh::AddDamping(unsigned int tag, std::shared_ptr<Damping> &damping){
    Dampings[tag] = std::move(damping);
}

//Add an element to the mesh.
void 
Mesh::AddElement(unsigned int tag, std::shared_ptr<Element>& element){
    Elements[tag] = std::move(element);
}

//Add a load to the mesh.
void 
Mesh::AddLoad(unsigned int tag, std::shared_ptr<Load>& load){
    Loads[tag] = std::move(load);
}

//Add point mass to node.
void 
Mesh::AddMass(unsigned int tag, Eigen::VectorXd& mass){
    Nodes[tag]->SetMass(mass);
}

//Sets damping to elements.
void 
Mesh::SetDamping(unsigned int tag, std::vector<unsigned int> &group){
    //List of elements to apply damping.
    for(unsigned int k = 0; k < group.size(); k++){
        Elements[group[k]]->SetDamping(Dampings[tag]);
    }
}

//Specifies the initial condition.
void 
Mesh::SetInitialCondition(unsigned int Tag, int Condition, Eigen::VectorXd& Xo){
   //Sets the initial state.
    if(Condition == 1)
        Nodes[Tag]->SetDisplacements(Xo);
    else if(Condition == 2)
        Nodes[Tag]->SetVelocities(Xo);
    else if(Condition == 3)
        Nodes[Tag]->SetAccelerations(Xo);
}

//Specifies the support motion for a certain node.
void 
Mesh::SetSupportMotion(unsigned int Tag, unsigned int dof, std::vector<double> &Xo){
    Nodes[Tag]->SetSupportMotion(dof, Xo);
}

//Gets a material from the mesh.
std::unique_ptr<Material>&
Mesh::GetMaterial(unsigned int tag){
    return Materials[tag];
}

//Gets a section from the mesh.
std::unique_ptr<Section>&
Mesh::GetSection(unsigned int tag){
    return Sections[tag];
}

//Gets classic damping from the mesh.
std::shared_ptr<Damping>&
Mesh::GetDamping(unsigned int tag){
    return Dampings[tag];
}

//Gets the nodes from the mesh.
std::map<unsigned int, std::shared_ptr<Node> >&
Mesh::GetNodes(){
    return Nodes;
}

//Gets the elements from the mesh.
std::map<unsigned int, std::shared_ptr<Element> >&
Mesh::GetElements(){
    return Elements;
}

//Gets the loads from the mesh.
std::map<unsigned int, std::shared_ptr<Load> >& 
Mesh::GetLoads(){
    return Loads;
}

//Gets operator that impose restrain/constrain on the model. 
Eigen::SparseMatrix<double>
Mesh::GetTotalToFreeMatrix(){
    //Total memory storage required for sparse matrix.
    unsigned int TotalStorage = numberOfConstrainedDofs + numberOfTotalDofs;

    //Construct the assembler matrix.
    Eigen::SparseMatrix<double> Free2TotalMatrix(numberOfTotalDofs, numberOfFreeDofs);

    //Sparse matrix format.
    std::vector<T> tripletList;
    tripletList.reserve(TotalStorage); 

    //Assembly transformation matrix process.
    unsigned int sum = 0;
    for(auto it : Nodes){
        auto &Tag = it.first;

        //Free degree-of fredom index list for this node.
        std::vector<int> free = Nodes[Tag]->GetFreeDegreeOfFreedom();

        //Total degree-of fredom index list for this node.
        std::vector<int> total = Nodes[Tag]->GetTotalDegreeOfFreedom();

        //Construct the total-to-free degree-of-freedom list matrix.
        for(unsigned int j = 0; j < free.size(); j++){

            //The selected degree-of-freedom is free.
            if (free[j] > -1){
                tripletList[sum] = T(total[j], free[j], 1.0);
                sum++;
            }

            //The selected degree-of-freedom is constrained.
            if (free[j] < -1){
                //Gets the slave degree-of-freedom to be constrained.
                unsigned int Slave = Constraints[free[j]]->GetSlaveInformation();

                //Gets the combinational factors for the constrained degree-of-freedom.
                std::vector<double> Factors = Constraints[free[j]]->GetCombinationFactors();

                //Gets the master degree-of-freedoms that the slave node must satisfy.
                std::vector<unsigned int> Master = Constraints[free[j]]->GetMasterInformation();

                //Assign combination factors to constrained degree-of-freedom.
                for(unsigned int i = 0; i < Factors.size(); i++){
                    tripletList[sum] = T(Slave, Master[i], Factors[i]);
                    sum++;
                }
            }
        }
    }

    //Builds the transformation sparse matrix.
    Free2TotalMatrix.setFromTriplets(tripletList.begin(), tripletList.begin() + sum);

    return Free2TotalMatrix;
}
