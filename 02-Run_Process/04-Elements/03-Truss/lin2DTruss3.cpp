#include <cmath>
#include <iostream>
#include "Material.hpp"
#include "lin2DTruss3.hpp"
#include "GaussQuadrature.hpp"
#include "Definitions.hpp"

//Define constant tolerance value:
const double TOL = 0.9999995;

//Define VTK cell value for Paraview:
const unsigned int VTKCELL = 21;

//Overload constructor.
lin2DTruss3::lin2DTruss3(const std::vector<unsigned int> nodes, std::unique_ptr<Material> &material, const double area, const unsigned int nGauss, bool massform) :
Element("lin2DTruss3", nodes, 6, VTKCELL), A(area), MassForm(massform){
    //The element nodes.
    theNodes.resize(3);

    //Numerical integration rule.
    QuadraturePoints = std::make_unique<GaussQuadrature>("Line", nGauss);

    //The element material. 
    theMaterial.resize(nGauss);
    for(unsigned int i = 0; i < nGauss; i++)
        theMaterial[i] = material->CopyMaterial();
}

//Destructor.
lin2DTruss3::~lin2DTruss3(){
    //Does nothing.
}

//Save the material states in the element.
void 
lin2DTruss3::CommitState(){
    //Updates the viscous material components.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    if(theMaterial[0]->IsViscous()){
        //Gets the quadrature information.    
        Eigen::VectorXd wi;
        Eigen::MatrixXd xi;
        QuadraturePoints->GetQuadraturePoints("Line", wi, xi);

        //Update material states.
        for(unsigned int k = 0; k < wi.size(); k++){
            //Compute Strain-Displacement Matrix at Gauss Point.
            Eigen::MatrixXd Bij = ComputeStrainDisplacementMatrix(xi(k,0));

            //Computes Strain Rate vector.
            Eigen::VectorXd strainrate = ComputeStrainRate(Bij);

            //Update the material state.
            theMaterial[k]->UpdateState(strainrate, 2);
        }
    }

    for(unsigned int k = 0; k < nPoints; k++)
        theMaterial[k]->CommitState();
}

//Update the material states in the element.
void 
lin2DTruss3::UpdateState(){
    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Line", wi, xi);

    //Update material states.
    for(unsigned int k = 0; k < wi.size(); k++){
        //Compute Strain-Displacement Matrix at Gauss Point.
        Eigen::MatrixXd Bij = ComputeStrainDisplacementMatrix(xi(k,0));

        //Computes strain vector.
        Eigen::VectorXd strain = ComputeStrain(Bij);

        //Update the material state.
        theMaterial[k]->UpdateState(strain, 1);
    }
}

//Sets the finite element dependance among objects.
void 
lin2DTruss3::SetDomain(std::map<unsigned int, std::shared_ptr<Node> > &nodes){
    //Gets the global element connectivity.
    std::vector<unsigned int> conn = GetNodes();

    //Assign the element to mesh node pointer.  
    for(unsigned int i = 0; i < GetNumberOfNodes(); i++){
        theNodes[i] = nodes[conn[i]];
    }

    //Computes length of element. 
    Eigen::VectorXd Xi = theNodes[0]->GetCoordinates();
    Eigen::VectorXd Xj = theNodes[1]->GetCoordinates();

    Lo = (Xj - Xi).norm();
}

//Sets the damping model.
void 
lin2DTruss3::SetDamping(const std::shared_ptr<Damping> &damping){
    //The damping model
    theDamping = damping;
}

//Gets the list of total-degree of freedom of this element.
std::vector<unsigned int> 
lin2DTruss3::GetTotalDegreeOfFreedom() const{
    //Total number of degree-of-freedom.
    unsigned int nDofs = GetNumberOfDegreeOfFreedom();

    //Reserve memory for the element list of degree-of-freedom.
    std::vector<unsigned int> dofs(nDofs);

    //Construct the element list of degree-of-freedom for assembly.
    for(unsigned int j = 0; j < 3; j++){    
        unsigned int LengthDofs = theNodes[j]->GetNumberOfDegreeOfFreedom();
        std::vector<int> totalDofs = theNodes[j]->GetTotalDegreeOfFreedom();

        for(unsigned int i = 0; i < LengthDofs; i++)
            dofs[i + LengthDofs*j] = totalDofs[i];    
    }

    return dofs;
}

//Returns the material strain at integration points.
Eigen::MatrixXd 
lin2DTruss3::GetStrain() const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    Eigen::MatrixXd theStrain(nPoints,1);
    for(unsigned int k = 0; k < nPoints; k++)
        theStrain.row(k) = theMaterial[k]->GetStrain();

    return theStrain;
}

//Returns the material stress at integration points.
Eigen::MatrixXd 
lin2DTruss3::GetStress() const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    Eigen::MatrixXd theStress(nPoints,1);
    for(unsigned int k = 0; k < nPoints; k++)
        theStress.row(k) = theMaterial[k]->GetTotalStress();

    return theStress;
}

//Returns the material strain-rate at integration points.
Eigen::MatrixXd 
lin2DTruss3::GetStrainRate() const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    Eigen::MatrixXd theStrainRate(nPoints,1);
    for(unsigned int k = 0; k < nPoints; k++)
        theStrainRate.row(k) = theMaterial[k]->GetStrainRate();

    return theStrainRate;
}

//Gets the material strain in section at  coordinate (x3,x2).
Eigen::MatrixXd 
lin2DTruss3::GetStrainAt(double x3, double x2) const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    //Stress at coordinate is define within section.
    Eigen::MatrixXd theStrain(nPoints, 3);
    theStrain.fill(0.0);

    return theStrain;
}

//Gets the material stress in section at  coordinate (x3,x2).
Eigen::MatrixXd 
lin2DTruss3::GetStressAt(double x3, double x2) const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    //Stress at coordinate is define within section.
    Eigen::MatrixXd theStress(nPoints, 3);
    theStress.fill(0.0);

    return theStress;
}

//Gets the element internal response in VTK format.
Eigen::VectorXd 
lin2DTruss3::GetVTKResponse(std::string response) const{
    //The VTK response vector.
    Eigen::VectorXd theResponse(6);

    if (strcasecmp(response.c_str(),"Strain") == 0){
        double nu = theMaterial[0]->GetPoissonRatio();
        Eigen::MatrixXd strain = GetStrain();
        Eigen::VectorXd Strain = strain.colwise().mean();
        theResponse << Strain(0), -nu*Strain(0), -nu*Strain(0), 0.0, 0.0, 0.0;
    }
    else if(strcasecmp(response.c_str(),"Stress") == 0){
        Eigen::MatrixXd stress = GetStress();
        Eigen::VectorXd Stress = stress.colwise().mean();
        theResponse << Stress(0), 0.0, 0.0, 0.0, 0.0, 0.0;
    }

    return theResponse;
}

//Compute the mass matrix of the element using a consistent definition.
Eigen::MatrixXd 
lin2DTruss3::ComputeMassMatrix(){
    //Use consistent mass definition:
    Eigen::MatrixXd MassMatrix(6,6);
    MassMatrix.fill(0.0);

    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Line", wi, xi);

    //Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
        //Gets material properties:
        double rho = theMaterial[i]->GetDensity();

        //Compute Interpolation Function Matrix at Gauss Point:
        Eigen::MatrixXd Hij = ComputeShapeFunctionMatrix(xi(i,0));

        //Numerical integration:
        MassMatrix += wi(i)*rho*A*Lo/2.0*Hij.transpose()*Hij;
    }

    //Lumped Mass Formulation
    if(MassForm){
        //Lumped Mass in diagonal terms.
        double m11 = MassMatrix(0,0) + MassMatrix(0,2) + MassMatrix(0,4);
        double m22 = MassMatrix(2,0) + MassMatrix(2,2) + MassMatrix(2,4);
        double m33 = MassMatrix(4,0) + MassMatrix(4,2) + MassMatrix(4,4);

        MassMatrix << m11, 0.0, 0.0, 0.0, 0.0, 0.0,
                      0.0, m11, 0.0, 0.0, 0.0, 0.0,
                      0.0, 0.0, m22, 0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0, m22, 0.0, 0.0,
                      0.0, 0.0, 0.0, 0.0, m33, 0.0,
                      0.0, 0.0, 0.0, 0.0, 0.0, m33;
    }

    //Computes the transformation matrix.
    Eigen::MatrixXd localAxes = ComputeTransformationAxes();

    //Transform Mass matrix into Global Coordinates.
    MassMatrix = localAxes.transpose()*MassMatrix*localAxes;



    return MassMatrix;
}

//Compute the stiffness matrix of the element.
Eigen::MatrixXd 
lin2DTruss3::ComputeStiffnessMatrix(){
    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Line", wi, xi);

    //Stiffness matrix definition:
    Eigen::MatrixXd ke(3,3);
    ke.fill(0.0);

    //Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
        //Compute Strain-Displacement Matrix at Gauss Point.
        Eigen::MatrixXd Bij = ComputeStrainDisplacementMatrix(xi(i,0));

        //Gets material tangent matrix at Gauss point.
        Eigen::MatrixXd Cij = theMaterial[i]->GetTangentStiffness();

        //Numerical integration.
        ke += wi(i)*A*Lo/2.0*Bij.transpose()*Cij*Bij;
    }

    //The global axes transformation: 
    Eigen::MatrixXd localAxes = ComputeLocalAxes();

    //Stiffness matrix definition.
    Eigen::MatrixXd StiffnessMatrix = localAxes.transpose()*ke*localAxes;

    return StiffnessMatrix;
}

//Compute the damping matrix of the element.
Eigen::MatrixXd 
lin2DTruss3::ComputeDampingMatrix(){
    //Damping matrix definition.
    Eigen::MatrixXd DampingMatrix(6,6);
    DampingMatrix.fill(0.0);
    
    //TODO: check if this is initial stiffness
    std::string dampName = theDamping->GetName();
    std::vector<double> dampParam = theDamping->GetParameters();

    if(strcasecmp(dampName.c_str(),"Free") == 0){
        //Does nothing.
    }
    else if(strcasecmp(dampName.c_str(),"Rayleigh") == 0){    
        //Compute stiffness and mass matrix.
        Eigen::MatrixXd MassMatrix = ComputeMassMatrix();
        Eigen::MatrixXd StiffnessMatrix = ComputeInitialStiffnessMatrix();

        DampingMatrix += dampParam[0]*MassMatrix + dampParam[1]*StiffnessMatrix;
    }
    else if(strcasecmp(dampName.c_str(),"Caughey") == 0){
        //TODO: implement Caughey damping
    }

    //Adds material damping contribution.
    if(theMaterial[0]->IsViscous()){ 
        //Gets the quadrature information.
        Eigen::VectorXd wi;
        Eigen::MatrixXd xi;
        QuadraturePoints->GetQuadraturePoints("Line", wi, xi);

        //Stiffness matrix definition:
        Eigen::MatrixXd ce(3,3);
        ce.fill(0.0);

        //Numerical integration.
        for(unsigned int i = 0; i < wi.size(); i++){
            //Compute Strain-Displacement Matrix at Gauss Point.
            Eigen::MatrixXd Bij = ComputeStrainDisplacementMatrix(xi(i,0));

            //Gets material tangent matrix at Gauss point.
            Eigen::MatrixXd eta = theMaterial[i]->GetDamping();

            //Numerical integration.
            ce += wi(i)*A*Lo/2.0*Bij.transpose()*eta*Bij;
        }

        //The global axes transformation:
        Eigen::MatrixXd localAxes = ComputeLocalAxes();
        
        //update damping matrix for material damping
        DampingMatrix += localAxes.transpose()*ce*localAxes;
    }  

    return DampingMatrix;
}

//Compute the element the internal forces acting on the element.
Eigen::VectorXd 
lin2DTruss3::ComputeInternalForces(){
    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Line", wi, xi);

    //Stiffness matrix definition:
    Eigen::VectorXd Finternal(3);
    Finternal.fill(0.0);

    //Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
        //Compute Strain-Displacement Matrix at Gauss Point.
        Eigen::MatrixXd Bij = ComputeStrainDisplacementMatrix(xi(i,0));

        //Gets material strain at Gauss point.
        Eigen::VectorXd Stress = theMaterial[i]->GetStress();

        //Numerical integration.
        Finternal += wi(i)*A*Lo/2.0*Bij.transpose()*Stress;
    }

    //The global axes transformation.
    Eigen::MatrixXd localAxes = ComputeLocalAxes(); 

    //Internal force vector definition.
    Eigen::VectorXd InternalForces = localAxes.transpose()*Finternal;

    return InternalForces;
}

//Compute the elastic, inertial, and vicous forces acting on the element.
Eigen::VectorXd 
lin2DTruss3::ComputeInternalDynamicForces(){
    //The Internal dynamic force vector
    Eigen::VectorXd InternalForces;

    if( HasFixedNode(theNodes) ){
        //Allocate memory for velocity/acceleraton. 
        Eigen::VectorXd V(6); 
        Eigen::VectorXd A(6);

        //Fills the response vectors with velocity/acceleraton values.
        V << theNodes[0]->GetVelocities(), theNodes[1]->GetVelocities(), theNodes[2]->GetVelocities();
        A << theNodes[0]->GetAccelerations(), theNodes[1]->GetAccelerations(), theNodes[2]->GetAccelerations();

        //Compute the inertial/viscous/elastic dynamic force contribution.
        InternalForces = ComputeInternalForces() + ComputeDampingMatrix()*V + ComputeMassMatrix()*A;
    }

    return InternalForces;
}

//Compute the surface forces acting on the element.
Eigen::VectorXd 
lin2DTruss3::ComputeSurfaceForces(const std::shared_ptr<Load> &surfaceLoad, unsigned int face){
    //Local surface load vector:
    Eigen::VectorXd surfaceForces(6);

    //Gets the surface force:
    Eigen::VectorXd qs = surfaceLoad->GetLoadVector();

    //Transformation matrix to local coordinates.
    Eigen::MatrixXd localAxes = ComputeTransformationAxes();
    Eigen::MatrixXd RotationMatrix = localAxes.block(0,0,2,2);

    //Transform load from global to local coordinates.
    qs = RotationMatrix*qs;

    surfaceForces << 1.0/6.0*qs(0)*Lo, 
                     1.0/6.0*qs(1)*Lo, 
                     1.0/6.0*qs(0)*Lo, 
                     1.0/6.0*qs(1)*Lo, 
                     2.0/3.0*qs(0)*Lo, 
                     2.0/3.0*qs(1)*Lo;

    //Node load vector in global coordinates.
    surfaceForces = localAxes.transpose()*surfaceForces;

    return surfaceForces;
}

//Compute the body forces acting on the element.
Eigen::VectorXd 
lin2DTruss3::ComputeBodyForces(const std::shared_ptr<Load> &bodyLoad, unsigned int k){
    //Local body load vector:
    Eigen::VectorXd bodyForces(6);

    //Gets material properties:
    double rho = theMaterial[0]->GetDensity();

    //Gets the body force:
    Eigen::VectorXd qb = A*rho*bodyLoad->GetLoadVector(k);

    //Transformation matrix to local coordinates.
    Eigen::MatrixXd localAxes = ComputeTransformationAxes();
    Eigen::MatrixXd RotationMatrix = localAxes.block(0,0,2,2);

    //Transform load into local coordinates.
    qb = RotationMatrix*qb;

    bodyForces << 1.0/6.0*qb(0)*Lo, 
                  1.0/6.0*qb(1)*Lo, 
                  1.0/6.0*qb(0)*Lo, 
                  1.0/6.0*qb(1)*Lo, 
                  2.0/3.0*qb(0)*Lo, 
                  2.0/3.0*qb(1)*Lo;

    //Node load vector in global coordinates.
    bodyForces = localAxes.transpose()*bodyForces;

    return bodyForces;
}

//Compute the domain reduction forces acting on the element.
Eigen::VectorXd 
lin2DTruss3::ComputeDomainReductionForces(const std::shared_ptr<Load> &drm, unsigned int k){
    //TODO: Domain reduction forces are not implemented for Truss.
    Eigen::VectorXd DRMForces(6);
    DRMForces.fill(0.0);

    return DRMForces;
}

//Compute the current length of the element.
double
lin2DTruss3::ComputeLength() const{
    //Gets the element coordinates in undeformed configuration.  
    Eigen::VectorXd Xi = theNodes[0]->GetCoordinates();
    Eigen::VectorXd Xj = theNodes[1]->GetCoordinates();

    //Gets the element coordinates in deformed configuration. 
    Eigen::VectorXd xi = Xi + theNodes[0]->GetDisplacements() + theNodes[0]->GetIncrementalDisplacements();
    Eigen::VectorXd xj = Xj + theNodes[1]->GetDisplacements() + theNodes[1]->GetIncrementalDisplacements();

    //Local axis-1.
    Eigen::Vector2d v1;
    v1 = Xj - Xi;
    v1 = v1/v1.norm();

    //Current length of element:
    double L = v1.dot(xj - xi); 

    return L;
}

//Compute/update the local axis of the element.
Eigen::MatrixXd
lin2DTruss3::ComputeLocalAxes() const{
    //Gets the element node coordinates. 
    Eigen::VectorXd Xi = theNodes[0]->GetCoordinates();
    Eigen::VectorXd Xj = theNodes[1]->GetCoordinates();

    //The global axes transformation: 
    Eigen::MatrixXd localAxes(3,6);
    Eigen::Vector2d v1;

    //Local axis 1:
    v1 = Xj - Xi;
    v1 = v1/v1.norm();

    localAxes << v1(0), v1(1),  0.0 ,  0.0 ,  0.0 ,  0.0 ,
                  0.0 ,  0.0 , v1(0), v1(1),  0.0 ,  0.0 ,
                  0.0 ,  0.0 ,  0.0 ,  0.0 , v1(0), v1(1);

    return localAxes;
}

//Compute/update the local axis of the element.
Eigen::MatrixXd
lin2DTruss3::ComputeTransformationAxes() const{
    //Gets the element node coordinates. 
    Eigen::VectorXd Xi = theNodes[0]->GetCoordinates();
    Eigen::VectorXd Xj = theNodes[1]->GetCoordinates();

    //Local axis 1:
    Eigen::Vector2d v1;
    v1 = Xj - Xi;
    v1 = v1/v1.norm();

    //The global axes transformation: 
    Eigen::MatrixXd transformationAxes(6,6);

    transformationAxes <<  v1(0), v1(1),    0.0,   0.0,    0.0,   0.0,
                          -v1(1), v1(0),    0.0,   0.0,    0.0,   0.0,
                             0.0,   0.0,  v1(0), v1(1),    0.0,   0.0,
                             0.0,   0.0, -v1(1), v1(0),    0.0,   0.0,
                             0.0,   0.0,    0.0,   0.0,  v1(0), v1(1),
                             0.0,   0.0,    0.0,   0.0, -v1(1), v1(0);
     
    return transformationAxes;
}

//Update strain in the element.
Eigen::VectorXd 
lin2DTruss3::ComputeStrain(const Eigen::MatrixXd &Bij) const{
    //Gets the element coordinates in deformed configuration. 
    Eigen::VectorXd U1 = theNodes[0]->GetDisplacements() + theNodes[0]->GetIncrementalDisplacements();
    Eigen::VectorXd U2 = theNodes[1]->GetDisplacements() + theNodes[1]->GetIncrementalDisplacements();
    Eigen::VectorXd U3 = theNodes[2]->GetDisplacements() + theNodes[2]->GetIncrementalDisplacements();

    Eigen::VectorXd nodalDisplacement(6);
    nodalDisplacement << U1, U2, U3;

    //Compute the local axis transformation.
    Eigen::MatrixXd localAxes = ComputeLocalAxes();

    //Compute strain rate:
    Eigen::MatrixXd strain = Bij*localAxes*nodalDisplacement;

    //Strain vector definition:
    Eigen::VectorXd Strain(1);
    Strain << strain(0,0);

    return Strain;
}

//Update strain rate in the element.
Eigen::VectorXd 
lin2DTruss3::ComputeStrainRate(const Eigen::MatrixXd &Bij) const{
    //Gets the element velocities in undeformed configuration.  
    Eigen::VectorXd V1 = theNodes[0]->GetVelocities();
    Eigen::VectorXd V2 = theNodes[1]->GetVelocities();
    Eigen::VectorXd V3 = theNodes[2]->GetVelocities();

    Eigen::VectorXd nodalVelocities(6);
    nodalVelocities << V1, V2, V3;

    //Compute the local axis transformation.
    Eigen::MatrixXd localAxes = ComputeLocalAxes();

    //Compute strain rate:
    Eigen::MatrixXd rate = Bij*localAxes*nodalVelocities;

    //Strain vector definition:
    Eigen::VectorXd strainrate(1);
    strainrate << rate(0,0);

    return strainrate;
}

//Evaluates the shape function matrix at a given Gauss point.
Eigen::MatrixXd 
lin2DTruss3::ComputeShapeFunctionMatrix(const double ri) const{
    //Shape function coefficients:
    double H11 = ri/2.0*(ri - 1.0);
    double H22 = ri/2.0*(ri + 1.0);
    double H33 = 1.0 - ri*ri;

    //Shape function matrix:
    Eigen::MatrixXd Hij(2,6);
    Hij << H11, 0.0, H22, 0.0, H33, 0.0,
           0.0, H11, 0.0, H22, 0.0, H33;

    return Hij;
}

//Evaluates the strain-displacement matrix at a given Gauss point.
Eigen::MatrixXd 
lin2DTruss3::ComputeStrainDisplacementMatrix(const double ri) const{
    //Strain-displacement matrix coefficients:
    double B11 =  ri - 1.0/2.0;
    double B21 =  ri + 1.0/2.0;
    double B31 = -2.0*ri;

    //Shape function matrix:
    Eigen::MatrixXd Bij(1,3);
    Bij << 2.0/Lo*B11, 2.0/Lo*B21, 2.0/Lo*B31;

    return Bij;
}

//Compute the initial stiffness matrix of the element
Eigen::MatrixXd 
lin2DTruss3::ComputeInitialStiffnessMatrix() const{
    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Line", wi, xi);

    //Stiffness matrix definition:
    Eigen::MatrixXd ke(3,3);
    ke.fill(0.0);

    //Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
        //Compute Strain-Displacement Matrix at Gauss Point.
        Eigen::MatrixXd Bij = ComputeStrainDisplacementMatrix(xi(i,0));

        //Gets material tangent matrix at Gauss point.
        Eigen::MatrixXd Cij = theMaterial[i]->GetInitialTangentStiffness();

        //Numerical integration.
        ke += wi(i)*A*Lo/2.0*Bij.transpose()*Cij*Bij;
    }

    //The global axes transformation: 
    Eigen::MatrixXd localAxes = ComputeLocalAxes();

    //Stiffness matrix definition.
    Eigen::MatrixXd StiffnessMatrix = localAxes.transpose()*ke*localAxes;

    return StiffnessMatrix;
}
