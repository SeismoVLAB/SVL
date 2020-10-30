#include <cmath>
#include <cfloat>
#include <iostream>
#include "HDRBYamamoto3DLink.hpp"
#include "Definitions.hpp"

//Define constant tolerance value:
const double TOL = 0.9999995;

//Define VTK cell value for Paraview:
const unsigned int VTKCELL = 3; 

//Overload constructor.
HDRBYamamoto3DLink::HDRBYamamoto3DLink(const std::vector<unsigned int> nodes, double De, double Di, double hr, unsigned int dim, bool massform) :
Element("HDRBYamamoto3DLink", nodes, 2*dim, VTKCELL), Hr(hr), Dimension(dim), MassForm(massform){
    //The element nodes.
    theNodes.resize(2);

    //Default Bridgestone rubber values.
    n = 0.7;
    cr = 1.0;
    cs = 1.0;

    //HDRB Parameters
    alpha = 0.7*Hr;
    Ar    = (De*De - Di*Di)*3.141592653589793238/4.0;
    Krb   = 1.0E6*(0.22*cr + 1.0*cs)*Ar/Hr; 

    //Assign internal variables.
    Pn.resize(2); Pn.fill(0.0);
    Qn.resize(2); Qn.fill(0.0);
    Fn.resize(2); Fn.fill(0.0);

    Paux.resize(2); Paux.fill(0.0);
    Qaux.resize(2); Qaux.fill(0.0);

    //Number of Newton-Raphson iterations
    nMax  = 50;
}

//Destructor:
HDRBYamamoto3DLink::~HDRBYamamoto3DLink(){
    //Does nothing.
}

//Save the internal variable states in the element.
void 
HDRBYamamoto3DLink::CommitState(){
    Qn = Qaux;
    Pn = Paux;
}

//Update the material states in the element.
void 
HDRBYamamoto3DLink::UpdateState(){
    //Transformation matrix from global to local.
    Eigen::MatrixXd Tr = ComputeRotationMatrix();

    //Computes deformation vector in local coordinates.
    Eigen::VectorXd U = Tr*ComputeRelativeDeformation();

    //Total (Relative) deformation at n-th time step
    Eigen::VectorXd deformation(2); 
    deformation << U(1), U(2);       
        
    //Integration for internal variables.
    Paux = deformation;
    
    //Incremental trajectory vector.      
    Eigen::VectorXd dP = Paux - Pn;
    Eigen::VectorXd strain = deformation/Hr;

    double Pnorm = dP.norm();
    double Qnorm = Qn.norm();

    if (Pnorm < DBL_EPSILON){
        Qaux = Qn;
    }
    else if ((Pnorm > DBL_EPSILON) &&  (Qnorm < DBL_EPSILON)){
        Qaux = Qn + dP/alpha;
    }
    else{
        Qaux = Qn + Pnorm/alpha*(dP/Pnorm - pow(Qnorm, n)*Qn/Qnorm);
    }
        
    double gamma = strain.norm();
        
    //Radial Restoring Force (taur: [MPa], Fr: [N])
    double taur = 0.22*gamma + 0.2*pow(gamma - 1.8, 2.0)*(gamma > 1.8);
    Eigen::VectorXd Fr = 1.0E6*cr*taur*Ar*strain/gamma;

    //Non-Linear Restoring Force (taus: [MPa], Fs: [N])
    double taus = 0.25 + 0.02*gamma + 0.016*pow(gamma, 3.0);
    Eigen::VectorXd Fs = 1.0E6*cs*taus*Ar*Qaux;

    //Total Non-Linear Restoring Force
    Fn = Fr + Fs;
}

//Sets the finite element dependance among objects.
void 
HDRBYamamoto3DLink::SetDomain(std::map<unsigned int, std::shared_ptr<Node> > &nodes){
    //Gets the global element connectivity.
    std::vector<unsigned int> conn = GetNodes();

    //Assign the element to mesh node pointer.  
    for(unsigned int i = 0; i < GetNumberOfNodes(); i++){
        theNodes[i] = nodes[conn[i]];
    }
}

//Sets the damping model.
void 
HDRBYamamoto3DLink::SetDamping(const std::shared_ptr<Damping> &damping){
    //does nothing.
}

//Gets the list of total-degree of freedom of this element.
std::vector<unsigned int> 
HDRBYamamoto3DLink::GetTotalDegreeOfFreedom() const{
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
HDRBYamamoto3DLink::GetStrain() const{
    //TODO: Compute the correct strain.
    Eigen::MatrixXd theStrain(1,2);
    theStrain << Pn(0), Pn(1);

    return theStrain;
}

//Returns the material stress.
Eigen::MatrixXd 
HDRBYamamoto3DLink::GetStress() const{
    //TODO: Compute the correct stress.
    Eigen::MatrixXd theStress(1,2);
    theStress << Fn(0), Fn(1);

    return theStress;
}

//Gets the material strain rate.
Eigen::MatrixXd 
HDRBYamamoto3DLink::GetStrainRate() const{
    Eigen::MatrixXd theStrainRate(1,1);
    theStrainRate << 0.0;

    return theStrainRate;
}

//Gets the material strain in section at  coordinate (x3,x2).
Eigen::MatrixXd 
HDRBYamamoto3DLink::GetStrainAt(double x3, double x2) const{
    //Stress at coordinate is define within section.
    Eigen::MatrixXd theStrain(1, 6); 
    theStrain.fill(0.0);

    return theStrain;
}

//Gets the material stress in section at  coordinate (x3,x2).
Eigen::MatrixXd 
HDRBYamamoto3DLink::GetStressAt(double x3, double x2) const{
    //Stress at coordinate is define within section.
    Eigen::MatrixXd theStress(1, 6); 
    theStress.fill(0.0);

    return theStress;
}

//Gets the element internal response in VTK format.
Eigen::VectorXd 
HDRBYamamoto3DLink::GetVTKResponse(std::string response) const{
    //TODO: Stress/Strain responses
    //The VTK response vector.
    Eigen::VectorXd theResponse(6);
    theResponse.fill(0.0);

    return theResponse;
}

//Compute the mass matrix of the element.
Eigen::MatrixXd 
HDRBYamamoto3DLink::ComputeMassMatrix(){
    //The matrix dimension.
    unsigned int nDim = 2*Dimension;

    //Consistent mass definition.
    Eigen::MatrixXd MassMatrix(nDim, nDim);
    MassMatrix.fill(0.0);

    return MassMatrix;
}

//Compute the stiffness matrix of the element.
Eigen::MatrixXd 
HDRBYamamoto3DLink::ComputeStiffnessMatrix(){
    //The vector dimension.
    unsigned int nDim = 2*Dimension;

    //The global axes transformation.
    Eigen::MatrixXd localStiffness(nDim, nDim);
    localStiffness.fill(0.0);

    localStiffness(1, 1) = Krb;
    localStiffness(2, 2) = Krb;
    localStiffness(Dimension + 1, Dimension + 1) = Krb;
    localStiffness(Dimension + 2, Dimension + 2) = Krb;

    //The global axes transformation: 
    Eigen::MatrixXd localAxes = ComputeLocalAxes();

    //Stiffness matrix definition.
    Eigen::MatrixXd StiffnessMatrix = localAxes.transpose()*localStiffness*localAxes;

    return StiffnessMatrix;
}

//Compute damping matrix of the element.
Eigen::MatrixXd 
HDRBYamamoto3DLink::ComputeDampingMatrix(){
    //The matrix dimension.
    unsigned int nDim = 2*Dimension;

    //Consistent mass definition.
    Eigen::MatrixXd DampingMatrix(nDim, nDim);
    DampingMatrix.fill(0.0);

    return DampingMatrix;
}

//Compute the PML history matrix for Perfectly-Matched Layer (PML).
Eigen::MatrixXd 
HDRBYamamoto3DLink::ComputePMLMatrix(){
    Eigen::MatrixXd Kpml;
    return Kpml;
}

//Compute the element internal forces acting on the element.
Eigen::VectorXd 
HDRBYamamoto3DLink::ComputeInternalForces(){
    //The vector dimension.
    unsigned int nDim = 2*Dimension;

    //The global axes transformation.
    Eigen::VectorXd localForces(nDim);
    localForces.fill(0.0);

    localForces(1) = -Fn(0);
    localForces(2) = -Fn(1);
    localForces(Dimension + 1) = Fn(0);
    localForces(Dimension + 2) = Fn(1);

    //The global axes transformation.
    Eigen::MatrixXd localAxes = ComputeLocalAxes(); 
 
    //Compute the internal elastic force contribution.
    Eigen::VectorXd InternalForces = localAxes.transpose()*localForces;

    return InternalForces;
}

//Compute the elastic, inertial, and vicous forces acting on the element.
Eigen::VectorXd 
HDRBYamamoto3DLink::ComputeInternalDynamicForces(){
    //The Internal dynamic force vector
    Eigen::VectorXd InternalForces;

    if( HasFixedNode(theNodes) ){
        unsigned int ndims = 2*Dimension;

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

//Compute the PML history vector using gauss-integration.
Eigen::VectorXd 
HDRBYamamoto3DLink::ComputePMLVector(){
    //Empty PML vector.
    Eigen::VectorXd Fpml;
    return Fpml;
}

//Compute the surface forces acting on the element.
Eigen::VectorXd 
HDRBYamamoto3DLink::ComputeSurfaceForces(const std::shared_ptr<Load> &surface, unsigned int face){
    //Local surface load vector.
    Eigen::VectorXd surfaceForces(2*Dimension);
    surfaceForces.fill(0.0);

    return surfaceForces;
}

//Compute the body forces acting on the element.
Eigen::VectorXd 
HDRBYamamoto3DLink::ComputeBodyForces(const std::shared_ptr<Load> &body, unsigned int k){
    //Local body load vector.
    Eigen::VectorXd bodyForces(2*Dimension);
    bodyForces.fill(0.0);

    return bodyForces;
}

//Compute the domain reduction forces acting on the element.
Eigen::VectorXd 
HDRBYamamoto3DLink::ComputeDomainReductionForces(const std::shared_ptr<Load> &drm, unsigned int k){
    //Domain reduction force vector.
    unsigned int nDofs = GetNumberOfDegreeOfFreedom();
    Eigen::VectorXd DRMForces(nDofs);
    DRMForces.fill(0.0);

    return DRMForces;
}

//Compute the local axis of the element.
Eigen::MatrixXd
HDRBYamamoto3DLink::ComputeLocalAxes() const{
    //Gets the element coordinates in undeformed configuration.  
    Eigen::VectorXd Xi = theNodes[0]->GetCoordinates();
    Eigen::VectorXd Xj = theNodes[1]->GetCoordinates();

    //Local axis definition.
    Eigen::Vector3d v1;
    Eigen::Vector3d v2;
    Eigen::Vector3d v3;

    //Local axis 1.
    v1 = Xj - Xi;
    v1 = v1/v1.norm();

    //Local Axis 3.
    if(abs(v1(2)) > TOL){
        v3 << 0.0, v1(2), -v1(1);
        v3 = v3/v3.norm();
    }
    else{
        v3 << v1(1), -v1(0), 0.0;
        v3 = v3/v3.norm();
    }

    //Local Axis 2.
    v2 = v3.cross(v1);
    v2 = v2/v2.norm();

    //The matrix dimension.
    unsigned int nDim = 2*Dimension;

    //The global axes transformation:
    Eigen::MatrixXd localAxes(nDim, nDim);

    if(Dimension == 3){
        //For Solid element nodes
        localAxes << v1(0), v1(1), v1(2),  0.0,   0.0,   0.0,
                     v2(0), v2(1), v2(2),  0.0,   0.0,   0.0,
                     v3(0), v3(1), v3(2),  0.0,   0.0,   0.0,
                       0.0,   0.0,   0.0, v1(0), v1(1), v1(2),
                       0.0,   0.0,   0.0, v2(0), v2(1), v2(2),
                       0.0,   0.0,   0.0, v3(0), v3(1), v3(2);
    }
    else if(Dimension == 6){
        //For Structural element nodes
        localAxes << v1(0), v1(1), v1(2),  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
                     v2(0), v2(1), v2(2),  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
                     v3(0), v3(1), v3(2),  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
                       0.0,   0.0,   0.0, v1(0), v1(1), v1(2), 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
                       0.0,   0.0,   0.0, v2(0), v2(1), v2(2), 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
                       0.0,   0.0,   0.0, v3(0), v3(1), v3(2), 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
                       0.0,   0.0,   0.0,  0.0,   0.0,   0.0, v1(0), v1(1), v1(2), 0.0,   0.0,   0.0,
                       0.0,   0.0,   0.0,  0.0,   0.0,   0.0, v2(0), v2(1), v2(2), 0.0,   0.0,   0.0,
                       0.0,   0.0,   0.0,  0.0,   0.0,   0.0, v3(0), v3(1), v3(2), 0.0,   0.0,   0.0,
                       0.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0, v1(0), v1(1), v1(2),
                       0.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0, v2(0), v2(1), v2(2),
                       0.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0, v3(0), v3(1), v3(2);
    }
    else{
        //TODO: For non-standard nodes
        localAxes.setIdentity();
    }
    
    return localAxes;
}

//Compute the local axis of the element.
Eigen::MatrixXd
HDRBYamamoto3DLink::ComputeRotationMatrix() const{
    //Gets the element coordinates in undeformed configuration.  
    Eigen::VectorXd Xi = theNodes[0]->GetCoordinates();
    Eigen::VectorXd Xj = theNodes[1]->GetCoordinates();

    //Local axis definition.
    Eigen::Vector3d v1;
    Eigen::Vector3d v2;
    Eigen::Vector3d v3;

    //Local axis 1.
    v1 = Xj - Xi;
    v1 = v1/v1.norm();

    //Local Axis 3.
    if(abs(v1(2)) > TOL){
        v3 << 0.0, v1(2), -v1(1);
        v3 = v3/v3.norm();
    }
    else{
        v3 << v1(1), -v1(0), 0.0;
        v3 = v3/v3.norm();
    }

    //Local Axis 2.
    v2 = v3.cross(v1);
    v2 = v2/v2.norm();

    //The global axes transformation:
    Eigen::MatrixXd localAxes(Dimension, Dimension);

    if(Dimension == 3){
        //For Solid element nodes
        localAxes << v1(0), v1(1), v1(2),
                     v2(0), v2(1), v2(2),
                     v3(0), v3(1), v3(2);
    }
    else if(Dimension == 6){
        //For Structural element nodes
        localAxes << v1(0), v1(1), v1(2),  0.0,   0.0,   0.0,
                     v2(0), v2(1), v2(2),  0.0,   0.0,   0.0,
                     v3(0), v3(1), v3(2),  0.0,   0.0,   0.0,
                       0.0,   0.0,   0.0, v1(0), v1(1), v1(2),
                       0.0,   0.0,   0.0, v2(0), v2(1), v2(2),
                       0.0,   0.0,   0.0, v3(0), v3(1), v3(2);
    }
    else{
        //TODO: For non-standard nodes
        localAxes.setIdentity();
    }
    
    return localAxes;
}

//Update strain in the element.
Eigen::VectorXd 
HDRBYamamoto3DLink::ComputeRelativeDeformation() const{  
    //Gets the element displacements.
    Eigen::VectorXd Ui = theNodes[0]->GetDisplacements() + theNodes[0]->GetIncrementalDisplacements();
    Eigen::VectorXd Uj = theNodes[1]->GetDisplacements() + theNodes[1]->GetIncrementalDisplacements();
    
    Eigen::VectorXd deformation = Uj - Ui;

    return deformation;
}
