#include <cmath>
#include "lin2DFrame2.hpp"
#include "GaussQuadrature.hpp"
#include "LobattoQuadrature.hpp"
#include "Definitions.hpp"
#include "Profiler.hpp"

//Define constant tolerance value:
const double TOL = 0.9999995;

//Define VTK cell value for Paraview:
const unsigned int VTKCELL = 3;

//Overload constructor.
lin2DFrame2::lin2DFrame2(const std::vector<unsigned int> nodes, std::unique_ptr<Section> &section, bool formulation, const std::string quadrature, unsigned int nGauss) :
Element("lin2DFrame2", nodes, 6, VTKCELL, GROUPFRAME), Formulation(formulation), Phi(0.0){
    //The element nodes.
    theNodes.resize(2);

    //Numerical integration rule.
    if(strcasecmp(quadrature.c_str(),"GAUSS") == 0)
        QuadraturePoints = std::make_unique<GaussQuadrature>("Line", nGauss);
    else if(strcasecmp(quadrature.c_str(),"LOBATTO") == 0)
        QuadraturePoints = std::make_unique<LobattoQuadrature>("Line", nGauss);

    //The element material. 
    theSection.resize(nGauss);
    for(unsigned int i = 0; i < nGauss; i++)
        theSection[i] = section->CopySection();

    //Section Shear Factor.
    if(Formulation){
        Eigen::MatrixXd Cs = section->GetInitialTangentStiffness();
        Phi = Cs(1,1)/Cs(2,2);
    }
}

//Destructor.
lin2DFrame2::~lin2DFrame2(){
    //Does nothing.
}

//Save the section states in the element.
void 
lin2DFrame2::CommitState(){
    //It considers only linear elastic material, thus viscous material is not allowed.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    for(unsigned int k = 0; k < nPoints; k++){
        theSection[k]->CommitState();
    }
}

//Reverse the section states to previous converged state in this element.
void 
lin2DFrame2::ReverseState(){
    //Reverse the material components.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    for(unsigned int k = 0; k < nPoints; k++)
        theSection[k]->ReverseState();
}

//Brings the section state to its initial state in this element.
void 
lin2DFrame2::InitialState(){
    //Brings the material components to initial state.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    for(unsigned int k = 0; k < nPoints; k++)
        theSection[k]->InitialState();
}

//Update the section states in the element.
void 
lin2DFrame2::UpdateState(){
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
        theSection[k]->UpdateState(strain, 1);
    }
}

//Sets the finite element dependance among objects.
void 
lin2DFrame2::SetDomain(std::map<unsigned int, std::shared_ptr<Node> > &nodes){
    //Gets the global element connectivity.
    std::vector<unsigned int> conn = GetNodes();

    //Assign the element to mesh node pointer.  
    for(unsigned int i = 0; i < GetNumberOfNodes(); i++){
        theNodes[i] = nodes[conn[i]];
    }

    //Computes length of element. 
    Eigen::VectorXd Xi = theNodes[0]->GetCoordinates();
    Eigen::VectorXd Xj = theNodes[1]->GetCoordinates();

    L = (Xj - Xi).norm();

    //Section Shear Factor.
    if(Formulation)
        Phi = 12.0*Phi/L/L;
}

//Sets the damping model.
void 
lin2DFrame2::SetDamping(const std::shared_ptr<Damping> &damping){
    //The damping model
    theDamping = damping;
}

//Gets the list of total-degree of freedom of this element.
std::vector<unsigned int> 
lin2DFrame2::GetTotalDegreeOfFreedom() const{
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
lin2DFrame2::GetStrain() const{
    //Number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    Eigen::MatrixXd theStrain(nPoints,3);
    for(unsigned int k = 0; k < nPoints; k++)
        theStrain.row(k) = theSection[k]->GetStrain();

    return theStrain;
}

//Returns the section generalised stress at integration point.
Eigen::MatrixXd 
lin2DFrame2::GetStress() const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    Eigen::MatrixXd theStress(nPoints,3);
    for(unsigned int k = 0; k < nPoints; k++)
        theStress.row(k) = theSection[k]->GetStress();

    return theStress;
}

//Returns the section generalised strain-rate at integration point.
Eigen::MatrixXd 
lin2DFrame2::GetStrainRate() const{
    //TODO: No strain-rate is employed for frame section.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    Eigen::MatrixXd theStrainRate(nPoints,3);
    theStrainRate.fill(0.0);

    return theStrainRate;
}

//Gets the material strain in section at  coordinate (x3,x2).
Eigen::MatrixXd 
lin2DFrame2::GetStrainAt(double x3, double x2) const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    //Strain at coordinate is define within section.
    Eigen::MatrixXd theStrain(nPoints, 3);
    for(unsigned int k = 0; k < nPoints; k++)
        theStrain.row(k) = theSection[k]->GetStrainAt(x3, x2);

    return theStrain;
}

//Gets the material stress in section at  coordinate (x3,x2).
Eigen::MatrixXd 
lin2DFrame2::GetStressAt(double x3, double x2) const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    //Stress at coordinate is define within section.
    Eigen::MatrixXd theStress(nPoints, 3);
    for(unsigned int k = 0; k < nPoints; k++)
        theStress.row(k) = theSection[k]->GetStressAt(x3, x2);

    return theStress;
}

//Gets the element internal response in VTK format.
Eigen::VectorXd 
lin2DFrame2::GetVTKResponse(std::string response) const{   
    //The VTK response vector.
    Eigen::VectorXd theResponse(18);

    if (strcasecmp(response.c_str(),"Strain") == 0){
        Eigen::MatrixXd strain = GetStrain();
        Eigen::VectorXd Strain = strain.colwise().mean();

        //[exx, kzz, txy] = [u,x, v,xx, u,y+v,x]
        theResponse << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Strain(0), 0.0, Strain(2), 0.0, 0.0, Strain(1);
    }
    else if(strcasecmp(response.c_str(),"Stress") == 0){
        Eigen::MatrixXd stress = GetStress();
        Eigen::VectorXd Stress = stress.colwise().mean();

        //[Fx, Mzz, Qy] = [EA u,x, EI v,xx, GA (u,y+v,x)]
        theResponse << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Stress(0), 0.0, Stress(2), 0.0, 0.0, Stress(1);
    }

    return theResponse;
}

//Computes the element energy for a given deformation.
double 
lin2DFrame2::ComputeEnergy(){
    //TODO: Integrate over element volume to compute the energy
    return 0.0;
}

//Compute the mass matrix of the element.
Eigen::MatrixXd 
lin2DFrame2::ComputeMassMatrix(){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Line", wi, xi);

    //Consistent mass definition.
    Eigen::MatrixXd MassMatrix(6,6);
    MassMatrix.fill(0.0);

    //Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
        //Gets material properties:
        Eigen::MatrixXd Rho = theSection[i]->GetDensity();

        //Compute Interpolation Function Matrix at Gauss Point:
        Eigen::MatrixXd Hij = ComputeShapeFunctionMatrix(xi(i,0));

        //Numerical integration:
        MassMatrix += wi(i)*L/2.0*Hij.transpose()*Rho*Hij;
    }

    //Lumped Mass Formulation
    if(MassFormulation){
        //Lumped Mass in diagonal terms.
        double m11 = MassMatrix(0,0) + MassMatrix(0,3);
        double m22 = MassMatrix(1,1) + MassMatrix(1,4);
        double m33 = MassMatrix(2,2);

        MassMatrix << m11,  0.0,  0.0,  0.0,  0.0,  0.0,
                      0.0,  m22,  0.0,  0.0,  0.0,  0.0,
                      0.0,  0.0,  m33,  0.0,  0.0,  0.0,
                      0.0,  0.0,  0.0,  m11,  0.0,  0.0,
                      0.0,  0.0,  0.0,  0.0,  m22,  0.0,
                      0.0,  0.0,  0.0,  0.0,  0.0,  m33;
    }

    //Computes the transformation matrix.
    Eigen::MatrixXd localAxes = ComputeLocalAxes();

    //Transform Mass matrix into Global Coordinates.
    MassMatrix = localAxes.transpose()*MassMatrix*localAxes;

    return MassMatrix;
}

//Compute the stiffness matrix of the element.
Eigen::MatrixXd 
lin2DFrame2::ComputeStiffnessMatrix(){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Line", wi, xi);

    //Stiffness matrix definition:
    Eigen::MatrixXd StiffnessMatrix(6,6);
    StiffnessMatrix.fill(0.0);

    //Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
        //Compute Strain-Displacement Matrix at Gauss Point.
        Eigen::MatrixXd Bij = ComputeStrainDisplacementMatrix(xi(i,0));

        //Gets material tangent matrix at Gauss point.
        Eigen::MatrixXd Cij = theSection[i]->GetTangentStiffness();

        //Numerical integration.
        StiffnessMatrix += wi(i)*L/2.0*Bij.transpose()*Cij*Bij;
    }

    //The global axes transformation: 
    Eigen::MatrixXd localAxes = ComputeLocalAxes();

    //Stiffness matrix definition.
    StiffnessMatrix = localAxes.transpose()*StiffnessMatrix*localAxes;

    return StiffnessMatrix;
}

//Compute the damping matrix of the element.
Eigen::MatrixXd 
lin2DFrame2::ComputeDampingMatrix(){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Damping matrix definition.
    Eigen::MatrixXd DampingMatrix;
    DampingMatrix.resize(6,6);
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

    //No material damping contribution is allowed.
    return DampingMatrix;
}

//Compute the PML history matrix for Perfectly-Matched Layer (PML).
Eigen::MatrixXd 
lin2DFrame2::ComputePMLMatrix(){
    Eigen::MatrixXd Kpml;
    return Kpml;
}

//Compute the element the internal forces acting on the element.
Eigen::VectorXd 
lin2DFrame2::ComputeInternalForces(){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Line", wi, xi);

    //Stiffness matrix definition:
    Eigen::VectorXd InternalForces(6);
    InternalForces.fill(0.0);

    //Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
        //Compute Strain-Displacement Matrix at Gauss Point.
        Eigen::MatrixXd Bij = ComputeStrainDisplacementMatrix(xi(i,0));

        //Gets material strain at Gauss point.
        Eigen::VectorXd Stress = theSection[i]->GetStress();

        //Numerical integration.
        InternalForces += wi(i)*L/2.0*Bij.transpose()*Stress;
    }

    //The global axes transformation.
    Eigen::MatrixXd localAxes = ComputeLocalAxes(); 

    //Internal load vector in global coordinates.
    InternalForces = localAxes.transpose()*InternalForces;

    return InternalForces;
}

//Compute the elastic, inertial, and viscous forces acting on the element.
Eigen::VectorXd 
lin2DFrame2::ComputeInternalDynamicForces(){
    //The Internal dynamic force vector
    Eigen::VectorXd InternalForces;

    if( HasFixedNode(theNodes) ){
        //Allocate memory for velocity/acceleraton. 
        Eigen::VectorXd V(6); 
        Eigen::VectorXd A(6);

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
lin2DFrame2::ComputeSurfaceForces(const std::shared_ptr<Load> &surfaceLoad, unsigned int UNUSED(face)){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Local surface load vector:
    Eigen::VectorXd surfaceForces(6);

    //Gets the surface force:
    Eigen::VectorXd qs = surfaceLoad->GetLoadVector();

    //Transformation matrix to local coordinates.
    Eigen::MatrixXd localAxes = ComputeLocalAxes();
    Eigen::MatrixXd rotationMatrix = localAxes.block(0,0,2,2);

    //Transform load from global to local coordinates.
    qs = rotationMatrix*qs;

    surfaceForces << qs(0)*L/2.0, 
                     qs(1)*L/2.0,
                     qs(1)*L*L/12.0, 
                     qs(0)*L/2.0, 
                     qs(1)*L/2.0,
                    -qs(1)*L*L/12.0;

    //Node load vector in global coordinates.
    surfaceForces = localAxes.transpose()*surfaceForces;

    return surfaceForces;
}

//Compute the body forces acting on the element.
Eigen::VectorXd 
lin2DFrame2::ComputeBodyForces(const std::shared_ptr<Load> &bodyLoad, unsigned int k){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Local body load vector:
    Eigen::VectorXd bodyForces(6);

    //Gets section density .
    Eigen::MatrixXd rho = theSection[0]->GetDensity();

    //Gets the body force:
    Eigen::VectorXd qb = rho(0,0)*bodyLoad->GetLoadVector(k);

    //Transformation matrix to local coordinates.
    Eigen::MatrixXd localAxes = ComputeLocalAxes();
    Eigen::MatrixXd rotationMatrix = localAxes.block(0,0,2,2);

    //Transform load into local coordinates.
    qb = rotationMatrix*qb;

    bodyForces << qb(0)*L/2.0, 
                  qb(1)*L/2.0,
                  qb(1)*L*L/12.0, 
                  qb(0)*L/2.0, 
                  qb(1)*L/2.0,
                 -qb(1)*L*L/12.0;

    //Node load vector in global coordinates.
    bodyForces = localAxes.transpose()*bodyForces;

    return bodyForces;
}

//Compute the domain reduction forces acting on the element.
Eigen::VectorXd 
lin2DFrame2::ComputeDomainReductionForces(const std::shared_ptr<Load>& UNUSED(drm), unsigned int UNUSED(k)){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //TODO: Domain reduction forces not implemented for frame.
    Eigen::VectorXd DRMForces(6);
    DRMForces.fill(0.0);

    return DRMForces;
}

//Compute the current length of the element.
double
lin2DFrame2::ComputeLength() const{
    //Gets the element coordinates in undeformed configuration.  
    Eigen::VectorXd Xi = theNodes[0]->GetCoordinates();
    Eigen::VectorXd Xj = theNodes[1]->GetCoordinates();

    //Gets the element coordinates in deformed configuration. 
    Eigen::VectorXd xi = Xi + (theNodes[0]->GetDisplacements() + theNodes[0]->GetIncrementalDisplacements()).head(2);
    Eigen::VectorXd xj = Xj + (theNodes[1]->GetDisplacements() + theNodes[1]->GetIncrementalDisplacements()).head(2);

    //Local axis-1.
    Eigen::Vector2d v1;
    v1 = Xj - Xi;
    v1 = v1/v1.norm();

    //Current length of element:
    double l = v1.dot(xj - xi); 

    return l;
}

//Compute/update the local axis of the element.
Eigen::MatrixXd
lin2DFrame2::ComputeLocalAxes() const{
    //Gets the element coordinates in undeformed configuration.  
    Eigen::VectorXd Xi = theNodes[0]->GetCoordinates();
    Eigen::VectorXd Xj = theNodes[1]->GetCoordinates();

    //Local axis definition.
    Eigen::Vector2d v1 = Xj - Xi;
    v1 = v1/v1.norm();

    //The global axes transformation. 
    Eigen::MatrixXd localAxes(6, 6);

    localAxes <<  v1(0), v1(1), 0.0,   0.0,    0.0, 0.0,
                 -v1(1), v1(0), 0.0,   0.0,    0.0, 0.0,
                    0.0,   0.0, 1.0,   0.0,    0.0, 0.0,
                    0.0,   0.0, 0.0,  v1(0), v1(1), 0.0,
                    0.0,   0.0, 0.0, -v1(1), v1(0), 0.0,
                    0.0,   0.0, 0.0,   0.0,    0.0, 1.0;

    return localAxes;
}

//Update strain in the element.
Eigen::VectorXd 
lin2DFrame2::ComputeStrain(const Eigen::MatrixXd &Bij) const{
    //The global axes transformation.
    Eigen::MatrixXd localAxes    = ComputeLocalAxes(); 
    Eigen::MatrixXd rotationAxes = localAxes.block(0,0,3,3);
    
    //Gets the element displacements in local coordinates. 
    Eigen::VectorXd U1 = rotationAxes*(theNodes[0]->GetDisplacements() + theNodes[0]->GetIncrementalDisplacements());
    Eigen::VectorXd U2 = rotationAxes*(theNodes[1]->GetDisplacements() + theNodes[1]->GetIncrementalDisplacements());

    Eigen::VectorXd nodalDisplacement(6);
    nodalDisplacement << U1, U2;

    //Strain vector:
    //[exx, kzz, txy] = [u,x, v,xx, u,y+v,x]
    Eigen::VectorXd Strain(3); 
    Strain = Bij*nodalDisplacement;

    return Strain;
}

//Update strain rate in the element.
Eigen::VectorXd 
lin2DFrame2::ComputeStrainRate(const Eigen::MatrixXd &Bij) const{
    //The global axes transformation.
    Eigen::MatrixXd localAxes    = ComputeLocalAxes(); 
    Eigen::MatrixXd rotationAxes = localAxes.block(0,0,3,3);

    //Gets the element velocities in local coordinates. 
    Eigen::VectorXd V1 = rotationAxes*theNodes[0]->GetVelocities();
    Eigen::VectorXd V2 = rotationAxes*theNodes[1]->GetVelocities();

    Eigen::VectorXd nodalVelocities(6);
    nodalVelocities << V1, V2;

    //Strain rate vector: 
    Eigen::VectorXd StrainRate(3); 
    StrainRate = Bij*nodalVelocities;

    return StrainRate;
}

//Evaluates the shape function matrix at a given Gauss point.
Eigen::MatrixXd 
lin2DFrame2::ComputeShapeFunctionMatrix(const double ri) const{
    //Shape function coefficients:
    double N11 = 1.0/2.0*(1.0 - ri);
    double N21 = 1.0/2.0*(1.0 + ri);
    double N31 = 1.0/(4.0*(1.0 + Phi))*(1.0 - ri)*(1.0 - ri)*(2.0 + ri) + Phi/(2.0*(1.0 + Phi))*(1.0 - ri);
    double N41 = L/(8.0*(1.0 + Phi))*(1.0 - ri)*(1.0 - ri)*(ri + 1.0) + L*Phi/(8.0*(1.0 + Phi))*(1.0 - ri*ri);
    double N51 = 1.0/(4.0*(1.0 + Phi))*(1.0 + ri)*(1.0 + ri)*(2.0 - ri) + Phi/(2.0*(1.0 + Phi))*(1.0 + ri);
    double N61 = L/(8.0*(1.0 + Phi))*(1.0 + ri)*(1.0 + ri)*(ri - 1.0) - L*Phi/(8.0*(1.0 + Phi))*(1.0 - ri*ri);

    //Shape function matrix:
    Eigen::MatrixXd Nij(2,6);
    Nij << N11, 0.0, 0.0, N21, 0.0, 0.0,
           0.0, N31, N41, 0.0, N51, N61;

    return Nij;
}

//Evaluates the strain-displacement matrix at a given Gauss point.
Eigen::MatrixXd 
lin2DFrame2::ComputeStrainDisplacementMatrix(const double ri) const{
    //Strain-displacement matrix coefficients for bending.
    double B11 = -1.0/L;
    double B21 =  1.0/L;
    double B31 =  6.0*ri/(L*L*(1.0 + Phi));
    double B41 =  (3.0*ri - 1.0)/(L*(1.0 + Phi)) - Phi/(L*(1.0 + Phi));
    double B51 = -6.0*ri/(L*L*(1.0 + Phi));
    double B61 =  (3.0*ri + 1.0)/(L*(1.0 + Phi)) + Phi/(L*(1.0 + Phi));

    //Strain-displacement matrix coefficients for shearing.
    double S11 = -Phi/(L*(1.0 + Phi));
    double S21 = -Phi/(2.0*(1.0 + Phi));
    double S31 =  Phi/(L*(1.0 + Phi));
    double S41 = -Phi/(2.0*(1.0 + Phi));

    //Shape function matrix:
    Eigen::MatrixXd Bij(3,6);
    Bij << B11,  0.0,  0.0, B21,  0.0,  0.0,
           0.0, -B31, -B41, 0.0, -B51, -B61,
           0.0,  S11,  S21, 0.0,  S31,  S41;

    return Bij;
}

//Compute the initial stiffness matrix of the element
Eigen::MatrixXd 
lin2DFrame2::ComputeInitialStiffnessMatrix() const{
    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Line", wi, xi);

    //Stiffness matrix definition:
    Eigen::MatrixXd StiffnessMatrix(6,6);
    StiffnessMatrix.fill(0.0);

    //Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
        //Compute Strain-Displacement Matrix at Gauss Point.
        Eigen::MatrixXd Bij = ComputeStrainDisplacementMatrix(xi(i,0));

        //Gets material tangent matrix at Gauss point.
        Eigen::MatrixXd Cij = theSection[i]->GetInitialTangentStiffness();

        //Numerical integration.
        StiffnessMatrix += wi(i)*L/2.0*Bij.transpose()*Cij*Bij;
    }

    //The global axes transformation: 
    Eigen::MatrixXd localAxes = ComputeLocalAxes();

    //Stiffness matrix definition.
    StiffnessMatrix = localAxes.transpose()*StiffnessMatrix*localAxes;

    return StiffnessMatrix;
}
