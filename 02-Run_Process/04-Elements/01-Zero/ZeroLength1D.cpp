#include <iostream>
#include "ZeroLength1D.hpp"
#include "Definitions.hpp"
#include "Profiler.hpp"

//Define VTK cell value for Paraview:
const unsigned int VTKCELL = 3; 

//Overload constructor.
ZeroLength1D::ZeroLength1D(const std::vector<unsigned int> nodes, std::unique_ptr<Material> &material, const unsigned int dir) :
Element("ZeroLength1D", nodes, 2*nDimensions, VTKCELL), theDirection(dir){
    //The element nodes.
    theNodes.resize(2);

    //The element material.
    theMaterial = material->CopyMaterial();
}

//Destructor:
ZeroLength1D::~ZeroLength1D(){
    //Does nothing.
}

//Save the material states in the element.
void 
ZeroLength1D::CommitState(){
    if(theMaterial->IsViscous()){
        //Computes strain-rate vector.
        Eigen::VectorXd strainrate = ComputeStrainRate();

        //Update material states.
        theMaterial->UpdateState(strainrate, 2);
    }

    theMaterial->CommitState();
}

//Update the material states in the element.
void 
ZeroLength1D::UpdateState(){
    //Computes strain vector.
    Eigen::VectorXd strain = ComputeStrain();

    //Update material states.
    theMaterial->UpdateState(strain, 1);
}

//Sets the finite element dependance among objects.
void 
ZeroLength1D::SetDomain(std::map<unsigned int, std::shared_ptr<Node> > &nodes){
    //Gets the global element connectivity.
    std::vector<unsigned int> conn = GetNodes();

    //Assign the element to mesh node pointer.  
    for(unsigned int i = 0; i < GetNumberOfNodes(); i++){
        theNodes[i] = nodes[conn[i]];
    }
}

//Sets the damping model.
void 
ZeroLength1D::SetDamping(const std::shared_ptr<Damping>& UNUSED(damping)){
    //does nothing.
}

//Gets the list of total-degree of freedom of this element.
std::vector<unsigned int> 
ZeroLength1D::GetTotalDegreeOfFreedom() const{
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

//Returns the material strain.
Eigen::MatrixXd 
ZeroLength1D::GetStrain() const{
    Eigen::VectorXd Strain = theMaterial->GetStrain();
    Eigen::MatrixXd theStrain(1,1);
    theStrain << Strain(0);

    return theStrain;
}

//Returns the material stress.
Eigen::MatrixXd 
ZeroLength1D::GetStress() const{
    Eigen::VectorXd Stress = theMaterial->GetTotalStress();
    Eigen::MatrixXd theStress(1,1);
    theStress << Stress(0);

    return theStress;
}

//Gets the material strain rate.
Eigen::MatrixXd 
ZeroLength1D::GetStrainRate() const{
    Eigen::VectorXd strainrate = theMaterial->GetStrainRate();
    Eigen::MatrixXd theStrainRate(1,1);
    theStrainRate << strainrate(0);

    return theStrainRate;
}

//Gets the material strain in section at  coordinate (x3,x2).
Eigen::MatrixXd 
ZeroLength1D::GetStrainAt(double UNUSED(x3), double UNUSED(x2)) const{
    //Stress at coordinate is define within section.
    unsigned int ndim = nDimensions*(nDimensions + 1)/2;
    Eigen::MatrixXd theStrain(1, ndim);
    theStrain.fill(0.0);

    return theStrain;
}

//Gets the material stress in section at  coordinate (x3,x2).
Eigen::MatrixXd 
ZeroLength1D::GetStressAt(double UNUSED(x3), double UNUSED(x2)) const{
    //Stress at coordinate is define within section.
    unsigned int ndim = nDimensions*(nDimensions + 1)/2;
    Eigen::MatrixXd theStress(1, ndim);
    theStress.fill(0.0);

    return theStress;
}

//Gets the element internal response in VTK format.
Eigen::VectorXd 
ZeroLength1D::GetVTKResponse(std::string response) const{
    //The VTK response vector.
    Eigen::VectorXd theResponse(6);

    if (strcasecmp(response.c_str(),"Strain") == 0){
        Eigen::VectorXd Strain = theMaterial->GetStrain();
        theResponse << Strain(0), 0.0, 0.0, 0.0, 0.0, 0.0;
    }
    else if(strcasecmp(response.c_str(),"Stress") == 0){
        Eigen::VectorXd Stress = theMaterial->GetTotalStress();
        theResponse << Stress(0), 0.0, 0.0, 0.0, 0.0, 0.0;
    }

    return theResponse;
}

//Computes the element energy for a given deformation.
double 
ZeroLength1D::ComputeEnergy(){
    //TODO: Integrate over element volume to compute the energy
    return 0.0;
}

//Compute the mass matrix of the element.
Eigen::MatrixXd 
ZeroLength1D::ComputeMassMatrix(){
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

    //The matrix dimension.
    unsigned int nDim = 2*nDimensions;

    //Consistent mass definition.
    Eigen::MatrixXd MassMatrix(nDim, nDim);
    MassMatrix.fill(0.0);

    return MassMatrix;
}

//Compute the stiffness matrix of the element.
Eigen::MatrixXd 
ZeroLength1D::ComputeStiffnessMatrix(){
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

    //Gets material tangent matrix.
    Eigen::MatrixXd E = theMaterial->GetTangentStiffness();

    //The global axes transformation: 
    Eigen::VectorXd localAxes = ComputeLocalAxes();

    //zerolength element stiffness.
    double k = E(0,0);

    //Stiffness matrix definition.
    Eigen::MatrixXd StiffnessMatrix = k*localAxes*localAxes.transpose();

    return StiffnessMatrix;
}

//Compute damping matrix of the element.
Eigen::MatrixXd 
ZeroLength1D::ComputeDampingMatrix(){
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

    //Gets material damping matrix.
    Eigen::MatrixXd eta = theMaterial->GetDamping();
    
    //The global axes transformation:
    Eigen::VectorXd localAxes = ComputeLocalAxes();
    
    //zerolength element stiffness.
    double c = eta(0,0);
    
    //Define damping matrix
    Eigen::MatrixXd DampingMatrix = c*localAxes*localAxes.transpose();

    return DampingMatrix;
}

//Compute the PML history matrix for Perfectly-Matched Layer (PML).
Eigen::MatrixXd 
ZeroLength1D::ComputePMLMatrix(){
    Eigen::MatrixXd Kpml;
    return Kpml;
}

//Compute the element internal forces acting on the element.
Eigen::VectorXd 
ZeroLength1D::ComputeInternalForces(){
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

    //The global axes transformation.
    Eigen::VectorXd localAxes = ComputeLocalAxes(); 

    //The nodal internal force vector in global coordinates.
    Eigen::VectorXd stress = theMaterial->GetStress();
 
    //Compute the internal elastic force contribution.
    Eigen::VectorXd InternalForces = stress(0)*localAxes;

    return InternalForces;
}

//Compute the elastic, inertial, and vicous forces acting on the element.
Eigen::VectorXd 
ZeroLength1D::ComputeInternalDynamicForces(){
    //The Internal dynamic force vector
    Eigen::VectorXd InternalForces;

    if( HasFixedNode(theNodes) ){
        unsigned int ndims = 2*nDimensions;

        //Allocate memory for velocity/acceleraton. 
        Eigen::VectorXd V(ndims); 
        Eigen::VectorXd A(ndims);

        //Fills the response vectors with velocity/acceleraton values.
        V << theNodes[0]->GetVelocities(), theNodes[1]->GetVelocities();
        A << theNodes[0]->GetAccelerations(), theNodes[1]->GetAccelerations();

        //Compute the inertial/viscous/elastic dynamic force contribution.
        InternalForces = ComputeInternalForces() + ComputeDampingMatrix()*V + ComputeMassMatrix()*A;
    }

    return InternalForces;
}

//Compute the surface forces acting on the element.
Eigen::VectorXd 
ZeroLength1D::ComputeSurfaceForces(const std::shared_ptr<Load>& UNUSED(surface), unsigned int UNUSED(face)){
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

    //Local surface load vector.
    Eigen::VectorXd surfaceForces(2*nDimensions);
    surfaceForces.fill(0.0);

    return surfaceForces;
}

//Compute the body forces acting on the element.
Eigen::VectorXd 
ZeroLength1D::ComputeBodyForces(const std::shared_ptr<Load>& UNUSED(body), unsigned int UNUSED(k)){
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

    //Local body load vector.
    Eigen::VectorXd bodyForces(2*nDimensions);
    bodyForces.fill(0.0);

    return bodyForces;
}

//Compute the domain reduction forces acting on the element.
Eigen::VectorXd 
ZeroLength1D::ComputeDomainReductionForces(const std::shared_ptr<Load>& UNUSED(drm), unsigned int UNUSED(k)){
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

    //Domain reduction force vector.
    unsigned int nDofs = GetNumberOfDegreeOfFreedom();
    Eigen::VectorXd DRMForces(nDofs);
    DRMForces.fill(0.0);

    return DRMForces;
}

//Compute the local axis of the element.
Eigen::VectorXd
ZeroLength1D::ComputeLocalAxes() const{
    //The global axes transformation:
    Eigen::VectorXd localAxes;

    if(nDimensions == 1){
        localAxes.resize(2);
        localAxes << -1.0, 1.0;
    }
    else if(nDimensions == 2){
        localAxes.resize(4);
        if (theDirection == 0){
            localAxes << -1.0, 0.0, 1.0, 0.0;
        }
        else if (theDirection == 1){
            localAxes << 0.0, -1.0, 0.0, 1.0;
        }
    }
    else if(nDimensions == 3){
        localAxes.resize(6);
        if (theDirection == 0){
            localAxes << -1.0, 0.0, 0.0, 1.0, 0.0, 0.0;
        }
        else if (theDirection == 1){
            localAxes << 0.0, -1.0, 0.0, 0.0, 1.0, 0.0;
        }
        else if (theDirection == 2){
            localAxes << 0.0, 0.0, -1.0, 0.0, 0.0, 1.0;
        }
    }
    
    return localAxes;
}

//Update strain in the element.
Eigen::VectorXd 
ZeroLength1D::ComputeStrain() const{

    //Strain vector definition:
    Eigen::VectorXd strain(1);
    
    //Gets the element displacements.
    Eigen::VectorXd Ui = theNodes[0]->GetDisplacements() + theNodes[0]->GetIncrementalDisplacements();
    Eigen::VectorXd Uj = theNodes[1]->GetDisplacements() + theNodes[1]->GetIncrementalDisplacements();
    
    strain << Uj(theDirection) - Ui(theDirection);

    return strain;
}

//Update strain rate in the element.
Eigen::VectorXd 
ZeroLength1D::ComputeStrainRate() const{
    //Gets the element velocities in undeformed configuration.  
    Eigen::VectorXd Vi = theNodes[0]->GetVelocities();
    Eigen::VectorXd Vj = theNodes[1]->GetVelocities();

    //Strain vector definition:
    Eigen::VectorXd strainrate(1);
    strainrate << Vj(theDirection) - Vi(theDirection);

    return strainrate;
}
