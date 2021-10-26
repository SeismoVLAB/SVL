#include "Node.hpp"

//Overload constructor.
Node::Node(unsigned int numberDofs, Eigen::VectorXd coordinates, bool isFixed) : 
Fixed(isFixed), NumberOfDegreeOfFreedom(numberDofs), Coordinates(coordinates){
    //Allocate memory.
    FreeDegreeOfFreedom.resize(numberDofs);
    TotalDegreeOfFreedom.resize(numberDofs);

    Displacements.resize(numberDofs);
    Velocities.resize(numberDofs);
    Accelerations.resize(numberDofs);
    IncrementalDisplacements.resize(numberDofs);

    CumulatedForce.resize(numberDofs);
    PMLIntegratedVector.resize(numberDofs);

    //Assign zero initial conditions.
    Displacements.fill(0.0);
    Velocities.fill(0.0);
    Accelerations.fill(0.0);
    IncrementalDisplacements.fill(0.0);
    CumulatedForce.fill(0.0);

    //The node's history vector for PML 3D.
    PMLIntegratedVector.fill(0.0);
}

//Default Destructor.
Node::~Node(){
    //Does nothing.
}

//Whether the node is fixed/free.
bool 
Node::IsFixed() const{
    return Fixed;
}

//Set all node variables to be zero
void 
Node::InitialState(){
    Displacements.fill(0.0);
    Velocities.fill(0.0);
    Accelerations.fill(0.0);
    IncrementalDisplacements.fill(0.0);
    PMLIntegratedVector.fill(0.0);
}

//Set the node's mass.
void 
Node::SetMass(Eigen::VectorXd &mass){
    Mass = mass;
}

//Set whether the node is fixed/free.
void 
Node::SetAsFixed(bool isFixed){
    Fixed = isFixed;
}

//Set the node's reaction.
void 
Node::SetReaction(Eigen::VectorXd &reaction){
    Reaction = reaction;
}

//Set the node's external force from previous analysis.
void 
Node::SetProgressiveForces(Eigen::VectorXd &force){
    CumulatedForce = force;
}

//Set the reduced degree of freedom list.
void 
Node::SetCoordinates(Eigen::VectorXd &coordinates){
    //Assign New Coordinates.
    Coordinates = coordinates;
}

//Set the free degree of freedom number list.
void Node::SetFreeDegreeOfFreedom(std::vector<int> &freeDOF){
    //Assign Degree of freedom numbering to Free Degree of Freedom list.
    FreeDegreeOfFreedom = freeDOF;
}

//Set the global degree of freedom number list.
void 
Node::SetTotalDegreeOfFreedom(std::vector<int> &totalDOF){
    //Assign Degree of freedom numbering to Total Degree of Freedom list.
    TotalDegreeOfFreedom = totalDOF;
}

//Sets the current displacement state of this node.
void 
Node::SetDisplacements(Eigen::VectorXd &Uo){
    Displacements = Uo;
}

//Sets the current velocity state of this node.
void 
Node::SetVelocities(Eigen::VectorXd &Vo){
    Velocities = Vo;
}

//Sets the current acceleration state of this node.
void 
Node::SetAccelerations(Eigen::VectorXd &Ao){
    Accelerations = Ao;
}

//Sets the current incremental displacement of this node.
void 
Node::SetIncrementalDisplacements(Eigen::VectorXd &dU){
    IncrementalDisplacements = dU;
}

//Sets the domain-reduction nodal displacements.
void 
Node::SetDomainReductionMotion(Eigen::MatrixXd &Uo){
    DomainReductionMotion = Uo;
}

//Sets the nodal absorbent integrated history values.
void 
Node::SetPMLVector(Eigen::VectorXd &Ub){
    PMLIntegratedVector = Ub;
}

//Sets the domain-reduction nodal displacements.
void 
Node::SetSupportMotion(unsigned int k, std::vector<double> &Uo){
    SupportMotion[k] = Uo;
}

//Remove the displacements support motion associated with this Node.
void 
Node::DelSupportMotion(){
    SupportMotion.clear();
}

//Return the node's mass as a vector.
Eigen::VectorXd 
Node::GetMass() const{
    return Mass;
}

//Gets the node reaction force
Eigen::VectorXd 
Node::GetReaction() const{
    //if(!Fixed){
    if(Reaction.size() == 0){
        Eigen::VectorXd ZeroReaction(NumberOfDegreeOfFreedom);
        ZeroReaction.fill(0.0);

        return ZeroReaction;
    }
    return Reaction;
}

//Return the node's coordinate as a vector.
const Eigen::VectorXd&
Node::GetCoordinates() const{
    return Coordinates;
}

//Returns the current displacement state of this node.
const Eigen::VectorXd&
Node::GetDisplacements() const{
    return Displacements;
}

//Returns the current velocity state of this node.
const Eigen::VectorXd&
Node::GetVelocities() const{
    return Velocities;
}

//Returns the current acceleration state of this node.
const Eigen::VectorXd&
Node::GetAccelerations() const{
    return Accelerations;
}

//Returns the inertial forces associated to this node.
Eigen::VectorXd 
Node::GetInertialForces() const{
    if(Mass.size() != 0){
        //Initialize the node inertial force vector.
        Eigen::VectorXd Forces(NumberOfDegreeOfFreedom);

        //Computes the nodal inertial force if there is mass.
        for(unsigned int k = 0; k < NumberOfDegreeOfFreedom; k++)
            Forces(k) = Mass(k)*Accelerations(k);
        
        return Forces;
    }

    return Mass;
}

//Returns the external force from previous analysis.
const Eigen::VectorXd&
Node::GetProgressiveForces() const{
    return CumulatedForce;
}

//Gets the current PML history vector of this node.
const Eigen::VectorXd&
Node::GetPMLVector() const{
    return PMLIntegratedVector;
}

//Returns the current incremental displacement of this node.
const Eigen::VectorXd&
Node::GetIncrementalDisplacements() const{
    return IncrementalDisplacements;
}

//Returns the current incremental displacement of this node.
Eigen::VectorXd
Node::GetDomainReductionMotion(unsigned int k) const{
    return DomainReductionMotion.row(k);
}

//Returns the current incremental displacement of this node.
Eigen::VectorXd 
Node::GetSupportMotion(unsigned int k){
    //Construct the support motion vector.
    Eigen::VectorXd Lg(NumberOfDegreeOfFreedom);
    Lg.fill(0.0); 

    for(auto it : SupportMotion){
        auto &dof = it.first;

        if(k < SupportMotion[dof].size()){
            //The time step is within (assumes dynamic load)
            Lg(dof) = SupportMotion[dof][k];
        }
        else{
            //The time step is out-of-bound (assumes static load)
            Lg(dof) = SupportMotion[dof][0];
        }
    }

    return Lg;
}

//Gets the number of degree-of-freedom with support motion in this node.
unsigned int 
Node::GetNumberOfSupportMotion(){
    return SupportMotion.size();
}

//Return the index of this node in the mesh.
unsigned int 
Node::GetNumberOfDegreeOfFreedom() const{
    return NumberOfDegreeOfFreedom;
}

//Return the node's restrained degree of freedom as a vector.
const std::vector<int>& 
Node::GetFreeDegreeOfFreedom() const{
    return FreeDegreeOfFreedom;
}

//Return the node's total degree of freedom as a vector.
const std::vector<int>&
Node::GetTotalDegreeOfFreedom() const{
    return TotalDegreeOfFreedom;
}
