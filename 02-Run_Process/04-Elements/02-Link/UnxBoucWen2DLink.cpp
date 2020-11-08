#include <cmath>
#include <iostream>
#include "UnxBoucWen2DLink.hpp"
#include "Definitions.hpp"

//Define VTK cell value for Paraview:
const unsigned int VTKCELL = 3; 

//Overload constructor.
UnxBoucWen2DLink::UnxBoucWen2DLink(const std::vector<unsigned int> nodes, std::vector<double> params, std::vector<double> vars, const unsigned int dim, const unsigned int dir, bool massform) :
Element("UnxBoucWen2DLink", nodes, 2*dim, VTKCELL), MassForm(massform), Dimension(dim), Direction(dir){
    //The element nodes.
    theNodes.resize(2);

    //Assign internal variables.
    A      = params[0];
    mu     = params[1];
    eta    = params[2];
    beta   = params[3];
    gamma  = params[4];
    Tol    = params[5];

    //Auxiliar variables.
    double fy = vars[0];
    double Kinit  = vars[1];
    double alpha1 = vars[2];
    double alpha2 = vars[3];

    //Compute Bouc-Wen backbone curve parameters
    z    = 0.0;
    nMax = 50;
    qY   = fy*(1.0 - alpha1 - alpha2*pow(fy/Kinit, mu - 1.0));
    uY   = qY/Kinit;
    k0   = (1.0 - alpha1)*Kinit;
    k2   = alpha1*Kinit;
    k3   = alpha2*Kinit;

    //Initialize stiffnes and internal force.
    qbw = 0.0;
    kbw = Kinit;
}

//Destructor:
UnxBoucWen2DLink::~UnxBoucWen2DLink(){
    //Does nothing.
}

//Save the internal variable states in the element.
void 
UnxBoucWen2DLink::CommitState(){
    zn = z;
    Un = U;
}

//Update the material states in the element.
void 
UnxBoucWen2DLink::UpdateState(){
    //Transformation matrix from global to local.
    Eigen::MatrixXd Tr = ComputeRotationMatrix();

    //Computes deformation vector in local coordinates.
    Eigen::VectorXd Ubw = Tr*ComputeRelativeDeformation();
    U = Ubw(Direction);

    //Relative link deformation.
    double dUn = U - Un;

    //Newton-Raphson iteration.
    double f, df, dz;
    unsigned int k = 0;

    do{
        f  = z - zn - dUn/uY*(A - pow(abs(z), eta)*(gamma + beta*sign(z*dUn)));
        df = 1.0 + dUn/uY*eta*pow(abs(z), eta - 1.0)*sign(z)*(gamma + beta*sign(z*dUn));
        dz = f/df;
        z  = z - dz;
        k++;
    } while( (fabs(dz) > Tol) & (k < nMax) );

    //Derivative of internal variable w.r.t displacement.
    double dzdu = A - pow(fabs(z), eta)*(gamma + beta*sign(z*dUn));

    //Compute consistent stiggnes matrix and internal force vector.
    qbw = qY*z + k2*U + k3*sign(U)*pow(fabs(U), mu - 1.0);
    kbw = k0*dzdu + k2 + k3*mu*pow(abs(U), mu - 1.0);
}

//Sets the finite element dependance among objects.
void 
UnxBoucWen2DLink::SetDomain(std::map<unsigned int, std::shared_ptr<Node> > &nodes){
    //Gets the global element connectivity.
    std::vector<unsigned int> conn = GetNodes();

    //Assign the element to mesh node pointer.  
    for(unsigned int i = 0; i < GetNumberOfNodes(); i++){
        theNodes[i] = nodes[conn[i]];
    }
}

//Sets the damping model.
void 
UnxBoucWen2DLink::SetDamping(const std::shared_ptr<Damping> &damping){
    //does nothing.
    UNUNSED_PARAMETER(damping);
}

//Gets the list of total-degree of freedom of this element.
std::vector<unsigned int> 
UnxBoucWen2DLink::GetTotalDegreeOfFreedom() const{
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
UnxBoucWen2DLink::GetStrain() const{
    //The relative deformation.
    Eigen::MatrixXd theStrain(1,1);
    theStrain << Un;

    return theStrain;
}

//Returns the material stress.
Eigen::MatrixXd 
UnxBoucWen2DLink::GetStress() const{
    //The non-linear internal force.
    Eigen::MatrixXd theStress(1,1);
    theStress << qbw;

    return theStress;
}

//Gets the material strain rate.
Eigen::MatrixXd 
UnxBoucWen2DLink::GetStrainRate() const{
    Eigen::MatrixXd theStrainRate(1,1);
    theStrainRate << 0.0;

    return theStrainRate;
}

//Gets the material strain in section at  coordinate (x3,x2).
Eigen::MatrixXd 
UnxBoucWen2DLink::GetStrainAt(double x3, double x2) const{
    UNUNSED_PARAMETER(x3);
    UNUNSED_PARAMETER(x2);

    //Stress at coordinate is define within section.
    Eigen::MatrixXd theStrain(1, 3); 
    theStrain.fill(0.0);

    return theStrain;
}

//Gets the material stress in section at  coordinate (x3,x2).
Eigen::MatrixXd 
UnxBoucWen2DLink::GetStressAt(double x3, double x2) const{
    UNUNSED_PARAMETER(x3);
    UNUNSED_PARAMETER(x2);

    //Stress at coordinate is define within section.
    Eigen::MatrixXd theStress(1, 3); 
    theStress.fill(0.0);

    return theStress;
}

//Gets the element internal response in VTK format.
Eigen::VectorXd 
UnxBoucWen2DLink::GetVTKResponse(std::string response) const{
    UNUNSED_PARAMETER(response);

    //TODO: Stress/Strain responses
    //The VTK response vector.
    Eigen::VectorXd theResponse(6);
    theResponse.fill(0.0);

    return theResponse;
}

//Compute the mass matrix of the element.
Eigen::MatrixXd 
UnxBoucWen2DLink::ComputeMassMatrix(){
    //The matrix dimension.
    unsigned int nDim = 2*Dimension;

    //Consistent mass definition.
    Eigen::MatrixXd MassMatrix(nDim, nDim);
    MassMatrix.fill(0.0);

    return MassMatrix;
}

//Compute the stiffness matrix of the element.
Eigen::MatrixXd 
UnxBoucWen2DLink::ComputeStiffnessMatrix(){
    //The vector dimension.
    unsigned int nDim = 2*Dimension;

    //The global axes transformation.
    Eigen::MatrixXd localStiffness(nDim, nDim);
    localStiffness.fill(0.0);

    localStiffness(Direction , Direction) = kbw;
    localStiffness(Dimension + Direction, Dimension + Direction) = kbw;

    //The global axes transformation: 
    Eigen::MatrixXd localAxes = ComputeLocalAxes();

    //Stiffness matrix definition.
    Eigen::MatrixXd StiffnessMatrix = localAxes.transpose()*localStiffness*localAxes;

    return StiffnessMatrix;
}

//Compute damping matrix of the element.
Eigen::MatrixXd 
UnxBoucWen2DLink::ComputeDampingMatrix(){
    //The matrix dimension.
    unsigned int nDim = 2*Dimension;

    //Consistent mass definition.
    Eigen::MatrixXd DampingMatrix(nDim, nDim);
    DampingMatrix.fill(0.0);

    return DampingMatrix;
}

//Compute the PML history matrix for Perfectly-Matched Layer (PML).
Eigen::MatrixXd 
UnxBoucWen2DLink::ComputePMLMatrix(){
    Eigen::MatrixXd Kpml;
    return Kpml;
}

//Compute the element internal forces acting on the element.
Eigen::VectorXd 
UnxBoucWen2DLink::ComputeInternalForces(){
    //The vector dimension.
    unsigned int nDim = 2*Dimension;

    //The global axes transformation.
    Eigen::VectorXd localForces(nDim);
    localForces.fill(0.0);

    localForces(Direction) = -qbw;
    localForces(Dimension + Direction) = qbw;

    //The global axes transformation.
    Eigen::MatrixXd localAxes = ComputeLocalAxes(); 
 
    //Compute the internal elastic force contribution.
    Eigen::VectorXd InternalForces = localAxes.transpose()*localForces;

    return InternalForces;
}

//Compute the elastic, inertial, and vicous forces acting on the element.
Eigen::VectorXd 
UnxBoucWen2DLink::ComputeInternalDynamicForces(){
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

//Compute the surface forces acting on the element.
Eigen::VectorXd 
UnxBoucWen2DLink::ComputeSurfaceForces(const std::shared_ptr<Load> &surface, unsigned int face){
    UNUNSED_PARAMETER(face);
    UNUNSED_PARAMETER(surface);

    //Local surface load vector.
    Eigen::VectorXd surfaceForces(2*Dimension);
    surfaceForces.fill(0.0);

    return surfaceForces;
}

//Compute the body forces acting on the element.
Eigen::VectorXd 
UnxBoucWen2DLink::ComputeBodyForces(const std::shared_ptr<Load> &body, unsigned int k){
    UNUNSED_PARAMETER(k);
    UNUNSED_PARAMETER(body);

    //Local body load vector.
    Eigen::VectorXd bodyForces(2*Dimension);
    bodyForces.fill(0.0);

    return bodyForces;
}

//Compute the domain reduction forces acting on the element.
Eigen::VectorXd 
UnxBoucWen2DLink::ComputeDomainReductionForces(const std::shared_ptr<Load> &drm, unsigned int k){
    UNUNSED_PARAMETER(k);
    UNUNSED_PARAMETER(drm);

    //Domain reduction force vector.
    unsigned int nDofs = GetNumberOfDegreeOfFreedom();
    Eigen::VectorXd DRMForces(nDofs);
    DRMForces.fill(0.0);

    return DRMForces;
}

//Compute the local axis of the element.
Eigen::MatrixXd
UnxBoucWen2DLink::ComputeLocalAxes() const{
    //Gets the element coordinates in undeformed configuration.  
    Eigen::VectorXd Xi = theNodes[0]->GetCoordinates();
    Eigen::VectorXd Xj = theNodes[1]->GetCoordinates();

    //Local axis definition.
    Eigen::Vector2d v1 = Xj - Xi;
    v1 = v1/v1.norm();

    //The matrix dimension.
    unsigned int nDim = 2*Dimension;

    //The global axes transformation:
    Eigen::MatrixXd localAxes(nDim, nDim);

    if(Dimension == 2){
        //For Solid element nodes
        localAxes <<  v1(0), v1(1),   0.0,    0.0,
                     -v1(1), v1(0),   0.0,    0.0,
                        0.0,   0.0,  v1(0), v1(1),
                        0.0,   0.0, -v1(1), v1(0);
    }
    else if(Dimension == 3){
        //For Structural element nodes
        localAxes <<  v1(0), v1(1), 0.0,   0.0,    0.0, 0.0,
                     -v1(1), v1(0), 0.0,   0.0,    0.0, 0.0,
                        0.0,   0.0, 1.0,   0.0,    0.0, 0.0,
                        0.0,   0.0, 0.0,  v1(0), v1(1), 0.0,
                        0.0,   0.0, 0.0, -v1(1), v1(0), 0.0,
                        0.0,   0.0, 0.0,   0.0,    0.0, 1.0;
    }
    else{
        //TODO: For non-standard nodes
        localAxes.setIdentity();
    }
    
    return localAxes;
}

//Compute the local axis of the element.
Eigen::MatrixXd
UnxBoucWen2DLink::ComputeRotationMatrix() const{
    //Gets the element coordinates in undeformed configuration.  
    Eigen::VectorXd Xi = theNodes[0]->GetCoordinates();
    Eigen::VectorXd Xj = theNodes[1]->GetCoordinates();

    //Local axis definition.
    Eigen::Vector2d v1 = Xj - Xi;
    v1 = v1/v1.norm();

    //The global axes transformation:
    Eigen::MatrixXd localAxes(Dimension, Dimension);

    if(Dimension == 2){
        //For Solid element nodes
        localAxes <<  v1(0), v1(1),
                     -v1(1), v1(0);
    }
    else if(Dimension == 3){
        //For Structural element nodes
        localAxes <<  v1(0), v1(1), 0.0,
                     -v1(1), v1(0), 0.0,
                        0.0,   0.0, 1.0;
    }
    else{
        //TODO: For non-standard nodes
        localAxes.setIdentity();
    }
    
    return localAxes;
}

//Sign function implementation.
double
UnxBoucWen2DLink::sign(double x) const{
    return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
}

//Update strain in the element.
Eigen::VectorXd 
UnxBoucWen2DLink::ComputeRelativeDeformation() const{  
    //Gets the element displacements.
    Eigen::VectorXd Ui = theNodes[0]->GetDisplacements() + theNodes[0]->GetIncrementalDisplacements();
    Eigen::VectorXd Uj = theNodes[1]->GetDisplacements() + theNodes[1]->GetIncrementalDisplacements();
    
    Eigen::VectorXd deformation = Uj - Ui;

    return deformation;
}
