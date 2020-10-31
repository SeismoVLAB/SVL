#include <cmath>
#include <iostream>
#include "null3DFrame2.hpp"
#include "Definitions.hpp"

//Define constant tolerance value:
const double TOL = 0.9999995;

//Define VTK cell value for Paraview:
const unsigned int VTKCELL = 3;

//Overload constructor.
null3DFrame2::null3DFrame2(const std::vector<unsigned int> nodes) :
Element("null3DFrame2", nodes, 12, VTKCELL){
    //The element nodes.
    theNodes.resize(2);
}

//Destructor.
null3DFrame2::~null3DFrame2(){
    //Does nothing.
}

//Save the material states in the element.
void 
null3DFrame2::CommitState(){
    //Does nothing.
}

//Update the material states in the element.
void 
null3DFrame2::UpdateState(){
    //Does nothing.
}

//Sets the finite element dependance among objects.
void 
null3DFrame2::SetDomain(std::map<unsigned int, std::shared_ptr<Node> > &nodes){
    //Gets the global element connectivity.
    std::vector<unsigned int> conn = GetNodes();

    //Assign the element to mesh node pointer.  
    for(unsigned int i = 0; i < GetNumberOfNodes(); i++){
        theNodes[i] = nodes[conn[i]];
    }
}

//Sets the damping model.
void 
null3DFrame2::SetDamping(const std::shared_ptr<Damping> &damping){
    //The damping model
    theDamping = damping;
}

//Gets the list of total-degree of freedom of this element.
std::vector<unsigned int> 
null3DFrame2::GetTotalDegreeOfFreedom() const{
    //Total number of degree-of-freedom.
    unsigned int nDofs = GetNumberOfDegreeOfFreedom();

    //Reserve memory for the element list of degree-of-freedom.
    std::vector<unsigned int> dofs(nDofs);

    //Construct the element list of degree-of-freedom for assembly.
    for(unsigned int j = 0; j < 2; j++){    
        unsigned int LengthDofs = theNodes[j]->GetNumberOfDegreeOfFreedom();
        std::vector<int> totalDofs = theNodes[j]->GetTotalDegreeOfFreedom();

        for(unsigned int i = 0; i < LengthDofs; i++)
            dofs[i + LengthDofs*j] = totalDofs[i];    
    }

    return dofs;
}

//Returns the section generalised strain at integration point.
Eigen::MatrixXd 
null3DFrame2::GetStrain() const{
    //Considers one integration point.
    Eigen::MatrixXd theStrain(1,6);
    theStrain.fill(0.0);

    return theStrain;
}

//Returns the section generalised stress at integration point.
Eigen::MatrixXd 
null3DFrame2::GetStress() const{
    //Considers one integration point.
    Eigen::MatrixXd theStress(1,6);
    theStress.fill(0.0);

    return theStress;
}

//Returns the section generalised strain rate at integration point.
Eigen::MatrixXd 
null3DFrame2::GetStrainRate() const{
    //Considers one integration point.
    Eigen::MatrixXd theStrainRate(1,6);
    theStrainRate.fill(0.0);

    return theStrainRate;
}

//Gets the material strain in section at  coordinate (x3,x2).
Eigen::MatrixXd 
null3DFrame2::GetStrainAt(double x3, double x2) const{
    //Stress at coordinate is define within section.
    Eigen::MatrixXd theStrain(1, 6);
    theStrain.fill(0.0);

    return theStrain;
}

//Gets the material stress in section at  coordinate (x3,x2).
Eigen::MatrixXd 
null3DFrame2::GetStressAt(double x3, double x2) const{
    //Stress at coordinate is define within section.
    Eigen::MatrixXd theStress(1, 6);
    theStress.fill(0.0);

    return theStress;
}

//Gets the element internal response in VTK format.
Eigen::VectorXd 
null3DFrame2::GetVTKResponse(std::string response) const{   
    //The VTK response vector.
    Eigen::VectorXd theResponse(6);
    theResponse.fill(0.0);

    return theResponse;
}

//Compute the mass matrix of the element.
Eigen::MatrixXd 
null3DFrame2::ComputeMassMatrix(){
    //Consistent mass definition.
    Eigen::MatrixXd MassMatrix(12,12);
    MassMatrix.fill(0.0);

    return MassMatrix;
}

//Compute the stiffness matrix of the element using gauss-integration.
Eigen::MatrixXd 
null3DFrame2::ComputeStiffnessMatrix(){
    //Stiffness matrix definition:
    Eigen::MatrixXd StiffnessMatrix(12, 12);
    StiffnessMatrix.fill(0.0);

    return StiffnessMatrix;
}

//Compute the damping matrix of the element.
Eigen::MatrixXd 
null3DFrame2::ComputeDampingMatrix(){
    //Damping matrix definition.
    Eigen::MatrixXd DampingMatrix(12, 12);
    DampingMatrix.fill(0.0);
    
    return DampingMatrix;
}

//Compute the PML history matrix for Perfectly-Matched Layer (PML).
Eigen::MatrixXd 
null3DFrame2::ComputePMLMatrix(){
    Eigen::MatrixXd Kpml;
    return Kpml;
}

//Compute the element the internal forces acting on the element.
Eigen::VectorXd 
null3DFrame2::ComputeInternalForces(){
    //Stiffness matrix definition:
    Eigen::VectorXd InternalForces(12);
    InternalForces.fill(0.0);

    return InternalForces;
}

//Compute the elastic, inertial, and vicous forces acting on the element.
Eigen::VectorXd 
null3DFrame2::ComputeInternalDynamicForces(){
    //The Internal dynamic force vector
    Eigen::VectorXd InternalForces(12);
    InternalForces.fill(0.0);

    return InternalForces;
}

//Compute the PML history vector using gauss-integration.
Eigen::VectorXd 
null3DFrame2::ComputePMLVector(){
    //Empty PML vector.
    Eigen::VectorXd Fpml;
    return Fpml;
}

//Compute the surface forces acting on the element.
Eigen::VectorXd 
null3DFrame2::ComputeSurfaceForces(const std::shared_ptr<Load> &surfaceLoad, unsigned int face){
    //Local surface load vector:
    Eigen::VectorXd surfaceForces(12);
    surfaceForces.fill(0.0);

    return surfaceForces;
}

//Compute the body forces acting on the element.
Eigen::VectorXd 
null3DFrame2::ComputeBodyForces(const std::shared_ptr<Load> &bodyLoad, unsigned int k){
    //Local body load vector:
    Eigen::VectorXd bodyForces(12);
    bodyForces.fill(0.0);

    return bodyForces;
}

//Compute the domain reduction forces acting on the element.
Eigen::VectorXd 
null3DFrame2::ComputeDomainReductionForces(const std::shared_ptr<Load> &drm, unsigned int k){
    //TODO: Domain reduction forces not implemented for frame.
    Eigen::VectorXd DRMForces(12);
    DRMForces.fill(0.0);

    return DRMForces;
}
