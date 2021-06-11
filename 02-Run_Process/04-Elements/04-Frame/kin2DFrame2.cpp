#include <cmath>
#include "kin2DFrame2.hpp"
#include "Definitions.hpp"
#include "Profiler.hpp"

//Define constant tolerance value:
const double TOL = 0.9999995;

//Define VTK cell value for Paraview:
const unsigned int VTKCELL = 3;

//Overload constructor.
kin2DFrame2::kin2DFrame2(const std::vector<unsigned int> nodes, std::unique_ptr<Section> &section) :
Element("kin2DFrame2", nodes, 6, VTKCELL, GROUPFRAME){
    //The element nodes.
    theNodes.resize(2);

    //The element section.
    theSection = section->CopySection();
}

//Destructor.
kin2DFrame2::~kin2DFrame2(){
    //Does nothing.
}

//Save the section states in the element.
void 
kin2DFrame2::CommitState(){
    theSection->CommitState();
}

//Reverse the section states to previous converged state in this element.
void 
kin2DFrame2::ReverseState(){
    theSection->ReverseState();
}

//Brings the section state to its initial state in this element.
void 
kin2DFrame2::InitialState(){
    theSection->InitialState();
}

//Update the section states in the element.
void 
kin2DFrame2::UpdateState(){
    //Computes strain vector.
    Eigen::VectorXd strain = ComputeStrain();

    //Update material states.
    theSection->UpdateState(strain, 1);
}

//Sets the finite element dependance among objects.
void 
kin2DFrame2::SetDomain(std::map<unsigned int, std::shared_ptr<Node> > &nodes){
    //Gets the global element connectivity.
    std::vector<unsigned int> conn = GetNodes();

    //Assign the element to mesh node pointer.  
    for(unsigned int i = 0; i < GetNumberOfNodes(); i++){
        theNodes[i] = nodes[conn[i]];
    }

    //Computes length of element. 
    Eigen::VectorXd Xi = theNodes[0]->GetCoordinates();
    Eigen::VectorXd Xj = theNodes[1]->GetCoordinates();

    //Compute the initial angle (radians) of element.
    double dx = Xj(0) - Xi(0);
    double dz = Xj(1) - Xi(1);
    Alpha = atan2(dz, dx);

    //Undeformed (initial) Beam Length.
    Lo = (Xj - Xi).norm();
}

//Sets the damping model.
void 
kin2DFrame2::SetDamping(const std::shared_ptr<Damping> &damping){
    //The damping model
    theDamping = damping;
}

//Gets the list of total-degree of freedom of this element.
std::vector<unsigned int> 
kin2DFrame2::GetTotalDegreeOfFreedom() const{
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
kin2DFrame2::GetStrain() const{
    //Gets the section srains from section.
    Eigen::MatrixXd Strain = theSection->GetStrain().transpose();

    return Strain;
}

//Returns the section generalised stress at integration point.
Eigen::MatrixXd 
kin2DFrame2::GetStress() const{
    //Gets the section stress from section.
    Eigen::MatrixXd Stress = theSection->GetStress().transpose();

    return Stress;
}

//Returns the section generalised strain-rate at integration point.
Eigen::MatrixXd 
kin2DFrame2::GetStrainRate() const{
    //TODO: No strain-rate is employed for frame section.
    Eigen::MatrixXd StrainRate(1,3);
    StrainRate << 0.0, 0.0, 0.0;
    return StrainRate;
}

//Gets the material strain in section at  coordinate (x3,x2).
Eigen::MatrixXd 
kin2DFrame2::GetStrainAt(double x3, double x2) const{
    //Stress at coordinate is define within section.
    Eigen::MatrixXd theStrain = theSection->GetStrainAt(x3, x2);

    return theStrain;
}

//Gets the material stress in section at  coordinate (x3,x2).
Eigen::MatrixXd 
kin2DFrame2::GetStressAt(double x3, double x2) const{
    //Stress at coordinate is define within section.
    Eigen::MatrixXd theStress = theSection->GetStressAt(x3, x2);

    return theStress;
}

//Gets the element internal response in VTK format.
Eigen::VectorXd 
kin2DFrame2::GetVTKResponse(std::string response) const{   
    //The VTK response vector.
    Eigen::VectorXd theResponse(18);

    if (strcasecmp(response.c_str(),"Strain") == 0){
        Eigen::MatrixXd Strain = GetStrain();

        //[exx, kzz, txy] = [u,x, v,xx, u,y+v,x]
        theResponse << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Strain(0,0), 0.0, Strain(0,2), 0.0, 0.0, Strain(0,1);
    }
    else if(strcasecmp(response.c_str(),"Stress") == 0){
        Eigen::MatrixXd Stress = GetStress();

        //[Fx, Mzz, Qy] = [EA u,x, EI v,xx, GA (u,y+v,x)]
        theResponse << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Stress(0,0), 0.0, Stress(0,2), 0.0, 0.0, Stress(0,1);
    }

    return theResponse;
}

//Computes the element energy for a given deformation.
double 
kin2DFrame2::ComputeEnergy(){
    //TODO: Integrate over element volume to compute the energy
    return 0.0;
}

//Compute the mass matrix of the element.
Eigen::MatrixXd 
kin2DFrame2::ComputeMassMatrix(){
    PROFILE_FUNCTION();

    //Gets section density.
    Eigen::MatrixXd rho = theSection->GetDensity();

    //Gets the total mass:
    double mass = rho(0,0)*Lo;

    //Consistent mass definition.
    Eigen::MatrixXd MassMatrix(6,6);

    //Construct mass matrix according to formulation.
    if(MassFormulation){
        //Gets the total mass.
        double m11 = mass/2.0;
        double m33 = mass*Lo*Lo/105.0;

        MassMatrix << m11,  0.0,  0.0,  0.0,  0.0,  0.0,
                      0.0,  m11,  0.0,  0.0,  0.0,  0.0,
                      0.0,  0.0,  m33,  0.0,  0.0,  0.0,
                      0.0,  0.0,  0.0,  m11,  0.0,  0.0,
                      0.0,  0.0,  0.0,  0.0,  m11,  0.0,
                      0.0,  0.0,  0.0,  0.0,  0.0,  m33;
    }
    else{
        //Gets the total mass.
        double m11 = mass/3.0;
        double m22 = mass*13.0/35.0;
        double m33 = mass*Lo*Lo/105.0;
        double m14 = mass/6.0;
        double m23 = mass*Lo*11.0/210.0;
        double m25 = mass*9.0/70.0;
        double m26 = mass*Lo*13.0/420.0;
        double m36 = mass*Lo*Lo/140.0;

        MassMatrix << m11,  0.0,  0.0,  m14,  0.0,  0.0,
                      0.0,  m22,  m23,  0.0,  m25, -m26,
                      0.0,  m23,  m33,  0.0,  m26, -m36,
                      m14,  0.0,  0.0,  m11,  0.0,  0.0,
                      0.0,  m25,  m26,  0.0,  m22, -m23,
                      0.0, -m26, -m36,  0.0, -m23,  m33;
    }

    //Computes the transformation matrix.
    Eigen::MatrixXd localAxes = ComputeLocalAxes();

    //Transform Mass matrix into Global Coordinates.
    MassMatrix = localAxes.transpose()*MassMatrix*localAxes;

    return MassMatrix;
}

//Compute the stiffness matrix of the element.
Eigen::MatrixXd 
kin2DFrame2::ComputeStiffnessMatrix(){
    PROFILE_FUNCTION();

    //Computes the current element length.
    double L = ComputeLength();

    //Computes the current axial force.
    double N = ComputeAxialForce();

    //Computes the current bending moments.
    Eigen::VectorXd M = ComputeBendingMoment();

    //Cross-sectional stiffness.
    Eigen::MatrixXd EI = theSection->GetTangentStiffness();

    //Beam Tangent Section Stiffness matrix;
    Eigen::MatrixXd C(3,3);
    C << EI(0,0)/Lo,      0.0      ,        0.0    ,
            0.0    , 4.0*EI(1,1)/Lo, 2.0*EI(1,1)/Lo,
            0.0    , 2.0*EI(1,1)/Lo, 4.0*EI(1,1)/Lo;

    //The Linear/Geometric Strain Displacement matrix: 
    Eigen::VectorXd z(6);
    Eigen::VectorXd r(6);
    Eigen::MatrixXd B(3,6); 
    ComputeGeometricStrainDisplacementMatrix(B, z, r);

    //The material stiffness matrix.
    Eigen::MatrixXd KL = B.transpose()*C*B;

    //The material stiffness matrix.
    Eigen::MatrixXd KNL = N/L*z*z.transpose() + (M(0) + M(1))/L/L*(r*z.transpose() + z*r.transpose());

    //Stiffness matrix definition:
    Eigen::MatrixXd StiffnessMatrix = KL + KNL;

    return StiffnessMatrix;
}

//Compute the damping matrix of the element.
Eigen::MatrixXd 
kin2DFrame2::ComputeDampingMatrix(){
    PROFILE_FUNCTION();

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

    //No material damping contribution is allowed.
    return DampingMatrix;
}

//Compute the PML history matrix for Perfectly-Matched Layer (PML).
Eigen::MatrixXd 
kin2DFrame2::ComputePMLMatrix(){
    Eigen::MatrixXd Kpml;
    return Kpml;
}

//Compute the element the internal forces acting on the element.
Eigen::VectorXd 
kin2DFrame2::ComputeInternalForces(){
    PROFILE_FUNCTION();

    //Computes the current bending moments.
    Eigen::VectorXd Qlocal = ComputeLocalForces();

    //The Linear/Geometric Strain Displacement matrix: 
    Eigen::MatrixXd B(3,6);
    ComputeLinearStrainDisplacementMatrix(B);

    //The nodal internal force vector in global coordinates.
    Eigen::VectorXd InternalForces = B.transpose()*Qlocal;

    return InternalForces;
}

//Compute the elastic, inertial, and viscous forces acting on the element.
Eigen::VectorXd 
kin2DFrame2::ComputeInternalDynamicForces(){
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
kin2DFrame2::ComputeSurfaceForces(const std::shared_ptr<Load>& surfaceLoad, unsigned int UNUSED(face)){
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

    surfaceForces << qs(0)*Lo/2.0, 
                     qs(1)*Lo/2.0,
                     qs(1)*Lo*Lo/12.0, 
                     qs(0)*Lo/2.0, 
                     qs(1)*Lo/2.0,
                    -qs(1)*Lo*Lo/12.0;

    //Node load vector in global coordinates.
    surfaceForces = localAxes.transpose()*surfaceForces;

    return surfaceForces;
}

//Compute the body forces acting on the element.
Eigen::VectorXd 
kin2DFrame2::ComputeBodyForces(const std::shared_ptr<Load> &bodyLoad, unsigned int k){
    PROFILE_FUNCTION();

    //Local body load vector:
    Eigen::VectorXd bodyForces(6);

    //Gets section density .
    Eigen::MatrixXd rho = theSection->GetDensity();

    //Gets the body force:
    Eigen::VectorXd qb = rho(0,0)*bodyLoad->GetLoadVector(k);

    //Transformation matrix to local coordinates.
    Eigen::MatrixXd localAxes = ComputeLocalAxes();
    Eigen::MatrixXd rotationMatrix = localAxes.block(0,0,2,2);

    //Transform load into local coordinates.
    qb = rotationMatrix*qb;

    bodyForces << qb(0)*Lo/2.0, 
                  qb(1)*Lo/2.0,
                  qb(1)*Lo*Lo/12.0, 
                  qb(0)*Lo/2.0, 
                  qb(1)*Lo/2.0,
                 -qb(1)*Lo*Lo/12.0;

    //Node load vector in global coordinates.
    bodyForces = localAxes.transpose()*bodyForces;

    return bodyForces;
}

//Compute the domain reduction forces acting on the element.
Eigen::VectorXd 
kin2DFrame2::ComputeDomainReductionForces(const std::shared_ptr<Load>& UNUSED(drm), unsigned int UNUSED(k)){
    PROFILE_FUNCTION();

    //TODO: Domain reduction forces not implemented for frame.
    Eigen::VectorXd DRMForces(6);
    DRMForces.fill(0.0);

    return DRMForces;
}

//Compute the current length of the element.
double
kin2DFrame2::ComputeLength() const{
    //Gets the element node displacement.
    Eigen::VectorXd ui = theNodes[0]->GetDisplacements() + theNodes[0]->GetIncrementalDisplacements();
    Eigen::VectorXd uj = theNodes[1]->GetDisplacements() + theNodes[1]->GetIncrementalDisplacements();

    //Gets the element node coordinates. 
    Eigen::VectorXd xi = theNodes[0]->GetCoordinates() + ui.head(2);
    Eigen::VectorXd xj = theNodes[1]->GetCoordinates() + uj.head(2);

    //Current length of element:
    double L = (xj - xi).norm(); 

    return L;
}

//Compute/update the local axis of the element.
Eigen::MatrixXd
kin2DFrame2::ComputeLocalAxes() const{
    //Gets the element node displacement.
    Eigen::VectorXd ui = theNodes[0]->GetDisplacements() + theNodes[0]->GetIncrementalDisplacements();
    Eigen::VectorXd uj = theNodes[1]->GetDisplacements() + theNodes[1]->GetIncrementalDisplacements();

    //Gets the element coordinates in deformed configuration.  
    Eigen::VectorXd xi = theNodes[0]->GetCoordinates() + ui.head(2);
    Eigen::VectorXd xj = theNodes[1]->GetCoordinates() + uj.head(2);

    //Local axis definition.
    Eigen::Vector2d v1 = xj - xi;
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

//Update strain in the element compatible with section.
Eigen::VectorXd 
kin2DFrame2::ComputeStrain() const{
    //Computes the curent element length.
    double L = ComputeLength();

    //Cross-sectional stiffness.
    Eigen::MatrixXd EI = theSection->GetTangentStiffness();

    //Computes the current axial and bending moments.
    Eigen::VectorXd Qlocal = ComputeLocalForces();

    //Strain vector:
    //[exx, kzz, txy] = [u,x, v,xx, u,y+v,x]
    Eigen::VectorXd Strain(3); 
    Strain << Qlocal(0)/EI(0,0), 
              (Qlocal(2) - Qlocal(1))/2.0/EI(1,1),
              (Qlocal(1) + Qlocal(2))/L/EI(2,2);

    return Strain;
}

//Update strain rate in the element.
Eigen::VectorXd 
kin2DFrame2::ComputeStrainRate() const{
    //Strain rate vector: 
    Eigen::VectorXd StrainRate(3); 
    StrainRate.fill(0.0);

    return StrainRate;
}

//Computes the current axial strain.
double 
kin2DFrame2::ComputeAxialStrain() const{
    //Current length of element:
    double L = ComputeLength();

    //Strain vector definition:
    double strain = (L - Lo)/Lo;

    return strain;
}

//Computes the current axial force.
double 
kin2DFrame2::ComputeAxialForce() const{
    //Gets the element node displacement.
    Eigen::VectorXd ui = theNodes[0]->GetDisplacements();
    Eigen::VectorXd uj = theNodes[1]->GetDisplacements();

    //Gets the element node coordinates. 
    Eigen::VectorXd xi = theNodes[0]->GetCoordinates() + ui.head(2);
    Eigen::VectorXd xj = theNodes[1]->GetCoordinates() + uj.head(2);

    //Longitudinall cross-sectional stiffness.
    Eigen::MatrixXd EA = theSection->GetTangentStiffness();  

    //Current length of element:
    double L = (xj - xi).norm();

    //Compute the axial strain.     
    Eigen::VectorXd strain(1);
    strain << (L - Lo)/Lo;

    //Computes the axial force.
    double N = EA(0,0)*strain(0);

    return N;
}

//Computes the current bending moments.
Eigen::VectorXd 
kin2DFrame2::ComputeBendingMoment() const{
    //Gets the element coordinates in deformed configuration.  
    Eigen::VectorXd Xi = theNodes[0]->GetCoordinates();
    Eigen::VectorXd Xj = theNodes[1]->GetCoordinates();

    //Gets the element coordinates in deformed configuration.  
    Eigen::VectorXd Ui = theNodes[0]->GetDisplacements();
    Eigen::VectorXd Uj = theNodes[1]->GetDisplacements();

    //Flexural cross-sectional stiffness.
    Eigen::MatrixXd EI = theSection->GetTangentStiffness();

    //Compute the rotation angles at each node.
    double dx = (Xj(0) + Uj(0)) - (Xi(0) + Ui(0));
    double dz = (Xj(1) + Uj(1)) - (Xi(1) + Ui(1));

    double alpha = atan2(dz, dx);
    double beta1 = Ui(2) + Alpha;
    double beta2 = Uj(2) + Alpha;

    double theta1 = atan((cos(alpha)*sin(beta1) - sin(alpha)*cos(beta1))/(cos(alpha)*cos(beta1) + sin(alpha)*sin(beta1)));
    double theta2 = atan((cos(alpha)*sin(beta2) - sin(alpha)*cos(beta2))/(cos(alpha)*cos(beta2) + sin(alpha)*sin(beta2)));

    //Local internal force vector.
    Eigen::VectorXd M(2);
    M << EI(1,1)/Lo*(4.0*theta1 + 2.0*theta2),
         EI(1,1)/Lo*(2.0*theta1 + 4.0*theta2);

    return M;
}

//Computes the current axial and bending moments.
Eigen::VectorXd 
kin2DFrame2::ComputeLocalForces() const{
    //Gets the element coordinates in deformed configuration.  
    Eigen::VectorXd Xi = theNodes[0]->GetCoordinates();
    Eigen::VectorXd Xj = theNodes[1]->GetCoordinates();

    //Gets the element coordinates in deformed configuration.  
    Eigen::VectorXd Ui = theNodes[0]->GetDisplacements() + theNodes[0]->GetIncrementalDisplacements();
    Eigen::VectorXd Uj = theNodes[1]->GetDisplacements() + theNodes[1]->GetIncrementalDisplacements();

    //Cross-sectional stiffness.
    Eigen::MatrixXd EI = theSection->GetTangentStiffness();

    //Compute the current strain.
    double strain = ComputeAxialStrain();

    //Compute the rotation angles at each node.
    double dx = (Xj(0) + Uj(0)) - (Xi(0) + Ui(0));
    double dz = (Xj(1) + Uj(1)) - (Xi(1) + Ui(1));
    double alpha = atan2(dz, dx);

    double beta1  = Ui(2) + Alpha;
    double beta2  = Uj(2) + Alpha;
    double theta1 = atan((cos(alpha)*sin(beta1)-sin(alpha)*cos(beta1))/(cos(alpha)*cos(beta1)+sin(alpha)*sin(beta1)));
    double theta2 = atan((cos(alpha)*sin(beta2)-sin(alpha)*cos(beta2))/(cos(alpha)*cos(beta2)+sin(alpha)*sin(beta2)));

    //Local internal force vector.
    Eigen::VectorXd Qlocal(3);
    Qlocal << EI(0,0)*strain, 
              EI(1,1)/Lo*(4.0*theta1 + 2.0*theta2),
              EI(1,1)/Lo*(2.0*theta1 + 4.0*theta2);

    return Qlocal;
}

//The Linear Strain Displacement matrix: 
void 
kin2DFrame2::ComputeLinearStrainDisplacementMatrix(Eigen::MatrixXd &B) const{
    //Gets the element node displacement.
    Eigen::VectorXd ui = theNodes[0]->GetDisplacements() + theNodes[0]->GetIncrementalDisplacements();
    Eigen::VectorXd uj = theNodes[1]->GetDisplacements() + theNodes[1]->GetIncrementalDisplacements();

    //Gets the element node coordinates. 
    Eigen::VectorXd xi = theNodes[0]->GetCoordinates() + ui.head(2);
    Eigen::VectorXd xj = theNodes[1]->GetCoordinates() + uj.head(2);

    //Current beam length.
    double L = (xj - xi).norm();

    //Local longitudinal axis definition.
    Eigen::Vector2d v1 = (xj - xi)/L;

    //Linear Strain-displacement matrix.
    B << -v1(0)  , -v1(1)  , 0.0, v1(0)  ,  v1(1)  , 0.0,
         -v1(1)/L,  v1(0)/L, 1.0, v1(1)/L, -v1(0)/L, 0.0,        
         -v1(1)/L,  v1(0)/L, 0.0, v1(1)/L, -v1(0)/L, 1.0;
}

//The Linear/Geometric Strain Displacement matrix: 
void 
kin2DFrame2::ComputeGeometricStrainDisplacementMatrix(Eigen::MatrixXd &B, Eigen::VectorXd &z, Eigen::VectorXd &r) const{
    //Gets the element node displacement.
    Eigen::VectorXd ui = theNodes[0]->GetDisplacements() + theNodes[0]->GetIncrementalDisplacements();
    Eigen::VectorXd uj = theNodes[1]->GetDisplacements() + theNodes[1]->GetIncrementalDisplacements();

    //Gets the element node coordinates. 
    Eigen::VectorXd xi = theNodes[0]->GetCoordinates() + ui.head(2);
    Eigen::VectorXd xj = theNodes[1]->GetCoordinates() + uj.head(2);

    //Current beam length.
    double L = (xj - xi).norm();

    //Local longitudinal axis definition.
    Eigen::Vector2d v1 = (xj - xi)/L;

    //Linear Strain-displacement matrix.
    B << -v1(0)  , -v1(1)  , 0.0, v1(0)  ,  v1(1)  , 0.0,
         -v1(1)/L,  v1(0)/L, 1.0, v1(1)/L, -v1(0)/L, 0.0,        
         -v1(1)/L,  v1(0)/L, 0.0, v1(1)/L, -v1(0)/L, 1.0;

    //Geometric Axial Strain-displacement vector.    
    z << v1(1), -v1(0), 0.0, -v1(1),  v1(0), 0.0;

    //Geometric Bending Strain-displacement vector.
    r << -v1(0), -v1(1), 0.0, v1(0),  v1(1), 0.0;
}

//Compute the initial stiffness matrix of the element.
Eigen::MatrixXd 
kin2DFrame2::ComputeInitialStiffnessMatrix() const{
    //Cross-sectional stiffness.
    Eigen::MatrixXd EI = theSection->GetInitialTangentStiffness();

    //Beam Tangent Section Stiffness matrix;
    Eigen::MatrixXd C(3,3);
    C << EI(0,0)/Lo,      0.0      ,        0.0    ,
            0.0    , 4.0*EI(1,1)/Lo, 2.0*EI(1,1)/Lo,
            0.0    , 2.0*EI(1,1)/Lo, 4.0*EI(1,1)/Lo;

    //Gets the element node coordinates. 
    Eigen::VectorXd Xi = theNodes[0]->GetCoordinates();
    Eigen::VectorXd Xj = theNodes[1]->GetCoordinates();

    //Local longitudinal axis definition.
    Eigen::Vector2d v1 = (Xj - Xi)/Lo;

    //Linear Strain-displacement matrix.
    Eigen::MatrixXd B(3,6);
    B << -v1(0)   , -v1(1)   , 0.0, v1(0)   ,  v1(1)   , 0.0,
         -v1(1)/Lo,  v1(0)/Lo, 1.0, v1(1)/Lo, -v1(0)/Lo, 0.0,        
         -v1(1)/Lo,  v1(0)/Lo, 0.0, v1(1)/Lo, -v1(0)/Lo, 1.0;

    //The material stiffness matrix.
    Eigen::MatrixXd StiffnessMatrix = B.transpose()*C*B;

    return StiffnessMatrix;
}
