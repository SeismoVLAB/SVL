#include <cmath>
#include "Section.hpp"
#include "lin3DShell4.hpp"
#include "GaussQuadrature.hpp"
#include "LobattoQuadrature.hpp"
#include "Definitions.hpp"
#include "Profiler.hpp"

//Define constant tolerance value:
const double TOL = 0.9999995;

//Define VTK cell value for Paraview:
const unsigned int VTKCELL = 9;

//Overload constructor.
lin3DShell4::lin3DShell4(const std::vector<unsigned int> nodes, std::unique_ptr<Section> &section, const std::string quadrature, const unsigned int nGauss) :
Element("lin3DShell4", nodes, 24, VTKCELL, GROUPSHELL){
    //The element nodes.
    theNodes.resize(4);

    //Numerical integration rule.
    if(strcasecmp(quadrature.c_str(),"GAUSS") == 0)
        QuadraturePoints = std::make_unique<GaussQuadrature>("Quad", nGauss);
    else if(strcasecmp(quadrature.c_str(),"LOBATTO") == 0)
        QuadraturePoints = std::make_unique<LobattoQuadrature>("Quad", nGauss);

    //The element material. 
    theSection.resize(nGauss);
    for(unsigned int i = 0; i < nGauss; i++)
        theSection[i] = section->CopySection();
}

//Destructor.
lin3DShell4::~lin3DShell4(){
    //Does nothing.
}

//Save the section states in the element.
void 
lin3DShell4::CommitState(){
    //TODO: Viscous material in shell element is not allowed.
    //Updates the viscous material components.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    for(unsigned int k = 0; k < nPoints; k++)
        theSection[k]->CommitState();
}

//Reverse the section states to previous converged state in this element.
void 
lin3DShell4::ReverseState(){
    //Reverse the section components.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    for(unsigned int k = 0; k < nPoints; k++)
        theSection[k]->ReverseState();
}

//Brings the section state to its initial state in this element.
void 
lin3DShell4::InitialState(){
    //Brings the material components to initial state.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    for(unsigned int k = 0; k < nPoints; k++)
        theSection[k]->InitialState();
}

//Update the section states in the element.
void 
lin3DShell4::UpdateState(){
    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Quad", wi, xi);

    //Gets Node Local coordinates.
    Eigen::MatrixXd xyloc = ComputeLocalCoordinates();

    //Computes the Matrix to enforce constant Tension (wilson).
    Eigen::MatrixXd BM12 = ComputeConstantTensionMatrix();

    //Update section states.
    for(unsigned int k = 0; k < wi.size(); k++){
        //Computes strain vector.
        Eigen::VectorXd strain = ComputeStrain(BM12, xyloc, xi(k,0), xi(k,1));

        //Update the material state.
        theSection[k]->UpdateState(strain, 1);
    }
}

//Sets the finite element dependance among objects.
void 
lin3DShell4::SetDomain(std::map<unsigned int, std::shared_ptr<Node> > &nodes){
    //Gets the global element connectivity.
    std::vector<unsigned int> conn = GetNodes();

    //Assign the element to mesh node pointer.  
    for(unsigned int i = 0; i < GetNumberOfNodes(); i++){
        theNodes[i] = nodes[conn[i]];
    }
}

//Sets the damping model.
void 
lin3DShell4::SetDamping(const std::shared_ptr<Damping> &damping){
    //The damping model
    theDamping = damping;
}

//Gets the list of total-degree of freedom of this element.
std::vector<unsigned int> 
lin3DShell4::GetTotalDegreeOfFreedom() const{
    //Total number of degree-of-freedom.
    unsigned int nDofs = GetNumberOfDegreeOfFreedom();

    //Reserve memory for the element list of degree-of-freedom.
    std::vector<unsigned int> dofs(nDofs);

    //Construct the element list of degree-of-freedom for assembly.
    for(unsigned int j = 0; j < 4; j++){    
        unsigned int LengthDofs = theNodes[j]->GetNumberOfDegreeOfFreedom();
        std::vector<int> totalDofs = theNodes[j]->GetTotalDegreeOfFreedom();

        for(unsigned int i = 0; i < LengthDofs; i++)
            dofs[i + LengthDofs*j] = totalDofs[i];    
    }

    return dofs;
}

//Returns the section generalised strain at integration point.
Eigen::MatrixXd 
lin3DShell4::GetStrain() const{
    //Number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    Eigen::MatrixXd theStrain(nPoints,6);
    for(unsigned int k = 0; k < nPoints; k++)
        theStrain.row(k) = theSection[k]->GetStrain();

    return theStrain;
}

//Returns the section generalised stress at integration point.
Eigen::MatrixXd 
lin3DShell4::GetStress() const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    Eigen::MatrixXd theStress(nPoints,6);
    for(unsigned int k = 0; k < nPoints; k++)
        theStress.row(k) = theSection[k]->GetStress();

    return theStress;
}

//Returns the section generalised strain rate at integration point.
Eigen::MatrixXd 
lin3DShell4::GetStrainRate() const{
    //TODO: No strain-rate is employed for shell section.
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    //Material Strain at eah Gauss-Coordinate.
    Eigen::MatrixXd theStrainRate(nPoints,6);
    theStrainRate.fill(0.0);

    return theStrainRate;
}

//Gets the material strain in section at  coordinate (x3,x2).
Eigen::MatrixXd 
lin3DShell4::GetStrainAt(double x3, double UNUSED(x2)) const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    //Strain at coordinate is define within section.
    Eigen::MatrixXd theStrain(nPoints, 6);
    for(unsigned int k = 0; k < nPoints; k++)
        theStrain.row(k) = theSection[k]->GetStrainAt(x3);

    return theStrain;
}

//Gets the material stress in section at  coordinate (x3,x2).
Eigen::MatrixXd 
lin3DShell4::GetStressAt(double x3, double UNUSED(x2)) const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    //Stress at coordinate is define within section.
    Eigen::MatrixXd theStress(nPoints, 6);
    for(unsigned int k = 0; k < nPoints; k++)
        theStress.row(k) = theSection[k]->GetStressAt(x3);

    return theStress;
}

//Gets the element internal response in VTK format.
Eigen::VectorXd 
lin3DShell4::GetVTKResponse(std::string response) const{
    //The VTK response vector.
    Eigen::VectorXd theResponse(18);

    if (strcasecmp(response.c_str(),"Strain") == 0){
        Eigen::MatrixXd strain = GetStrain();
        Eigen::MatrixXd Strain = strain.colwise().mean();

        theResponse << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Strain(0), Strain(1), Strain(2), Strain(3), Strain(4), Strain(5), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    }
    else if(strcasecmp(response.c_str(),"Stress") == 0){
        Eigen::MatrixXd stress = GetStress();
        Eigen::MatrixXd Stress = stress.colwise().mean();

        theResponse << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Stress(0), Stress(1), Stress(2), Stress(3), Stress(4), Stress(5), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    }

    return theResponse;
}

//Computes the element energy for a given deformation.
double 
lin3DShell4::ComputeEnergy(){
    //TODO: Integrate over element volume to compute the energy
    return 0.0;
}

//Compute the mass matrix of the element gauss-integration.
Eigen::MatrixXd 
lin3DShell4::ComputeMassMatrix(){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Computes the transformation matrix.
    Eigen::MatrixXd localAxes = ComputeLocalAxes();

    //Gets Node Local coordinates.
    Eigen::MatrixXd xyloc = ComputeLocalCoordinates();

    //Plate mass matrix component.
    Eigen::MatrixXd PlateMass(12,12);
    Eigen::MatrixXd MembraneMass(12,12);

    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Quad", wi, xi);

    //Numerical integration.
    PlateMass.fill(0.0);
    MembraneMass.fill(0.0);
    for(unsigned int i = 0; i < wi.size(); i++){
        //Compute Shape Function Matrix and Jacobian at Gauss Point.
        Eigen::MatrixXd Jij(2,2);
        Eigen::MatrixXd Hij(3,12);

        //Gets material tangent matrix at Gauss point.
        Eigen::MatrixXd Rho = theSection[i]->GetDensity();

        //Plate Numerical integration.
        ComputePlateShapeFunctionMatrix(xi(i,0), xi(i,1), xyloc, Hij, Jij);
        PlateMass += wi(i)*fabs(Jij.determinant())*Hij.transpose()*Rho.block(3,3,3,3)*Hij;

        //Numerical integration.
        ComputeMembraneShapeFunctionMatrix(xi(i,0), xi(i,1), xyloc, Hij, Jij);
        MembraneMass += wi(i)*fabs(Jij.determinant())*Hij.transpose()*Rho.block(0,0,3,3)*Hij;
    }

    //Assembles the total the total stiffness matrix;
    Eigen::MatrixXd MassMatrix(24,24);
    AssemblePlateMembraneEffects(MassMatrix, MembraneMass, PlateMass);

    //Lumped Mass Formulation
    if(MassFormulation){
        double mtot = ComputeTotalMass();
        double MassX = MassMatrix(0,0) + MassMatrix(6,6) + MassMatrix(12,12) + MassMatrix(18,18);
        double MassY = MassMatrix(1,1) + MassMatrix(7,7) + MassMatrix(13,13) + MassMatrix(19,19);
        double MassZ = MassMatrix(2,2) + MassMatrix(8,8) + MassMatrix(14,14) + MassMatrix(20,20);
        
        //Lumped mass matrix definition.
        for(unsigned int i = 0; i < 24; i++){
            for(unsigned int j = 0; j < 24; j++){
                if(i != j)
                    MassMatrix(i,j) = 0.0;
            }
        }

        MassMatrix( 0, 0) = mtot*MassMatrix( 0, 0)/MassX; MassMatrix( 1, 1) = mtot*MassMatrix( 1, 1)/MassY; MassMatrix( 2, 2) = mtot*MassMatrix( 2, 2)/MassZ; 
        MassMatrix( 6, 6) = mtot*MassMatrix( 6, 6)/MassX; MassMatrix( 7, 7) = mtot*MassMatrix( 7, 7)/MassY; MassMatrix( 8, 8) = mtot*MassMatrix( 8, 8)/MassZ;
        MassMatrix(12,12) = mtot*MassMatrix(12,12)/MassX; MassMatrix(13,13) = mtot*MassMatrix(13,13)/MassY; MassMatrix(14,14) = mtot*MassMatrix(14,14)/MassZ;
        MassMatrix(18,18) = mtot*MassMatrix(18,18)/MassX; MassMatrix(19,19) = mtot*MassMatrix(19,19)/MassY; MassMatrix(20,20) = mtot*MassMatrix(20,20)/MassZ;
    }

    //Transform Mass matrix into Global Coordinates.
    MassMatrix = localAxes.transpose()*MassMatrix*localAxes;

    return MassMatrix;
}

//Compute the stiffness matrix of the element using gauss-integration.
Eigen::MatrixXd 
lin3DShell4::ComputeStiffnessMatrix(){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Computes the transformation matrix.
    Eigen::MatrixXd localAxes = ComputeLocalAxes();

    //Gets Node Local coordinates.
    Eigen::MatrixXd xyloc = ComputeLocalCoordinates();

    //Computes the Matrix to enforce constant Tension (wilson).
    Eigen::MatrixXd BM12 = ComputeConstantTensionMatrix();

    //Membrane stiffness matrix component.
    Eigen::MatrixXd PlateStiffness(12,12); 
    Eigen::MatrixXd MembraneStiffness(12,12); 
    
    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Quad", wi, xi);

    //Numerical integration.
    PlateStiffness.fill(0.0);
    MembraneStiffness.fill(0.0);

    for(unsigned int i = 0; i < wi.size(); i++){
        //Compute Strain-Displacement Matrix and Jacobian at Gauss Point.
        Eigen::MatrixXd Jij(2,2);
        Eigen::MatrixXd Bij(3,12);
        Eigen::MatrixXd Pij(1,12);

        //Gets material tangent matrix at Gauss point.
        Eigen::MatrixXd Cij = theSection[i]->GetTangentStiffness();

        //Plate Numerical integration.
        ComputePlateStrainDisplacementMatrix(xi(i,0), xi(i,1), xyloc, Bij, Jij);
        PlateStiffness += wi(i)*fabs(Jij.determinant())*Bij.transpose()*Cij.block(3,3,3,3)*Bij;

        //Membrane Numerical integration.
        ComputeMembraneStrainDisplacementMatrix(xi(i,0), xi(i,1), xyloc, BM12, Bij, Pij, Jij);
        MembraneStiffness += wi(i)*fabs(Jij.determinant())*Bij.transpose()*Cij.block(0,0,3,3)*Bij;
        MembraneStiffness += wi(i)*fabs(Jij.determinant())*Cij(2,2)*Pij.transpose()*Pij;
    }

    //Assembles the total the total stiffness matrix;
    Eigen::MatrixXd StiffnessMatrix(24,24);
    AssemblePlateMembraneEffects(StiffnessMatrix, MembraneStiffness, PlateStiffness);

    StiffnessMatrix = localAxes.transpose()*StiffnessMatrix*localAxes;

    return StiffnessMatrix;
}

//Compute the damping matrix of the element gauss-integration.
Eigen::MatrixXd 
lin3DShell4::ComputeDampingMatrix(){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Damping matrix definition.
    Eigen::MatrixXd DampingMatrix(24, 24);
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
lin3DShell4::ComputePMLMatrix(){
    Eigen::MatrixXd Kpml;
    return Kpml;
}

//Compute the element the internal forces acting on the element.
Eigen::VectorXd 
lin3DShell4::ComputeInternalForces(){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //TODO: compute internal forces for shell should be done by integration of stresses.
    //Gets the element coordinates in deformed configuration. 
    Eigen::VectorXd U1 = theNodes[0]->GetDisplacements() + theNodes[0]->GetIncrementalDisplacements();
    Eigen::VectorXd U2 = theNodes[1]->GetDisplacements() + theNodes[1]->GetIncrementalDisplacements();
    Eigen::VectorXd U3 = theNodes[2]->GetDisplacements() + theNodes[2]->GetIncrementalDisplacements();
    Eigen::VectorXd U4 = theNodes[3]->GetDisplacements() + theNodes[3]->GetIncrementalDisplacements();

    Eigen::VectorXd nodalDisplacement(24);
    nodalDisplacement << U1, U2, U3, U4;

    //Gets the stiffness matrix.
    Eigen::MatrixXd StiffnessMatrix = ComputeStiffnessMatrix();

    //The nodal internal force vector in local coordinates.
    Eigen::VectorXd InternalForces(24);
    InternalForces = StiffnessMatrix*nodalDisplacement;

    return InternalForces;
}

//Compute the elastic, inertial, and vicous forces acting on the element.
Eigen::VectorXd 
lin3DShell4::ComputeInternalDynamicForces(){
    //The Internal dynamic force vector
    Eigen::VectorXd InternalForces;

    if( HasFixedNode(theNodes) ){
        //Allocate memory for velocity/acceleraton. 
        Eigen::VectorXd V(24); 
        Eigen::VectorXd A(24);

        //Fills the response vectors with velocity/acceleraton values.
        V << theNodes[0]->GetVelocities(), theNodes[1]->GetVelocities(), theNodes[2]->GetVelocities(), theNodes[3]->GetVelocities();
        A << theNodes[0]->GetAccelerations(), theNodes[1]->GetAccelerations(), theNodes[2]->GetAccelerations(), theNodes[3]->GetAccelerations();

        //Compute the inertial/viscous/elastic dynamic force contribution.
        InternalForces = ComputeInternalForces() + ComputeDampingMatrix()*V + ComputeMassMatrix()*A;
    }

    return InternalForces;
}

//Compute the surface forces acting on the element.
Eigen::VectorXd 
lin3DShell4::ComputeSurfaceForces(const std::shared_ptr<Load> &surfaceLoad, unsigned int UNUSED(face)){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Gets Node Local coordinates.
    Eigen::MatrixXd xyloc = ComputeLocalCoordinates();

    //Transformation matrix to local coordinates.
    Eigen::MatrixXd localAxes = ComputeLocalAxes();
    Eigen::MatrixXd rotationMatrix = localAxes.block(0,0,3,3);

    //Gets the surface force.
    Eigen::VectorXd qs = surfaceLoad->GetLoadVector();

    //Transform load into local coordinates.
    qs = rotationMatrix*qs;

    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Quad", wi, xi);

    //Local body load vector.
    Eigen::VectorXd Qs(12);
    Qs.fill(0.0);

    //Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
        //Compute Strain-Displacement Matrix and Jacobian at Gauss Point.
        Eigen::MatrixXd Jij(2,2);
        Eigen::MatrixXd Nij(3,12);
        ComputeLoadShapeFunctionMatrix(xi(i,0), xi(i,1), xyloc, Nij, Jij);

        //Numerical integration.
        Qs += wi(i)*fabs(Jij.determinant())*Nij.transpose()*qs;
    }

    //The global body force vector.
    Eigen::VectorXd surfaceForces(24);
    surfaceForces << Qs(0), Qs(1), Qs(2), 0.0, 0.0, 0.0, Qs(3), Qs(4), Qs(5), 0.0, 0.0, 0.0, Qs(6), Qs(7), Qs(8), 0.0, 0.0, 0.0, Qs(9), Qs(10), Qs(11), 0.0, 0.0, 0.0;

    //Node load vector in global coordinates.
    surfaceForces = localAxes.transpose()*surfaceForces;

    return surfaceForces;
}

//Compute the body forces acting on the element.
Eigen::VectorXd 
lin3DShell4::ComputeBodyForces(const std::shared_ptr<Load> &bodyLoad, unsigned int k){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Gets Node Local coordinates.
    Eigen::MatrixXd xyloc = ComputeLocalCoordinates();

    //Transformation matrix to local coordinates.
    Eigen::MatrixXd localAxes = ComputeLocalAxes();
    Eigen::MatrixXd rotationMatrix = localAxes.block(0,0,3,3);

    //Gets the body force.
    Eigen::VectorXd qb = bodyLoad->GetLoadVector(k);

    //Transform load into local coordinates.
    qb = rotationMatrix*qb;

    //Local body load vector.
    Eigen::VectorXd Qb(12);
    Qb.fill(0.0);

    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Quad", wi, xi);

    //Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
        //Compute Strain-Displacement Matrix and Jacobian at Gauss Point.
        Eigen::MatrixXd Jij(2,2);
        Eigen::MatrixXd Nij(3,12);
        ComputeLoadShapeFunctionMatrix(xi(i,0), xi(i,1), xyloc, Nij, Jij);

        //Gets material tangent matrix at Gauss point.
        Eigen::MatrixXd Rho = theSection[i]->GetDensity();

        //Numerical integration.
        Qb += wi(i)*fabs(Jij.determinant())*Nij.transpose()*Rho.block(0,0,3,3)*qb;
    }

    //The global body force vector.
    Eigen::VectorXd bodyForces(24);
    bodyForces << Qb(0), Qb(1), Qb(2), 0.0, 0.0, 0.0, Qb(3), Qb(4), Qb(5), 0.0, 0.0, 0.0, Qb(6), Qb(7), Qb(8), 0.0, 0.0, 0.0, Qb(9), Qb(10), Qb(11), 0.0, 0.0, 0.0;

    //Node load vector in global coordinates.
    bodyForces = localAxes.transpose()*bodyForces;

    return bodyForces;
}

//Compute the domain reduction forces acting on the element.
Eigen::VectorXd 
lin3DShell4::ComputeDomainReductionForces(const std::shared_ptr<Load>& UNUSED(drm), unsigned int UNUSED(k)){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //TODO: Domain reduction forces are not implemented for shell.
    Eigen::VectorXd DRMForces(24);
    DRMForces.fill(0.0);

    return DRMForces;
}

//Compute/update the Rotation axis of the element.
Eigen::MatrixXd 
lin3DShell4::ComputeRotation() const{
    //Gets the element coordinates in undeformed configuration.  
    Eigen::VectorXd Xi = theNodes[0]->GetCoordinates();
    Eigen::VectorXd Xj = theNodes[1]->GetCoordinates();
    Eigen::VectorXd Xk = theNodes[2]->GetCoordinates();
    Eigen::VectorXd Xl = theNodes[3]->GetCoordinates();

    //Gets the side mid-point coordinates of the element.
    Eigen::VectorXd Xmij = 0.5*(Xi + Xj);
    Eigen::VectorXd Xmjk = 0.5*(Xj + Xk);
    Eigen::VectorXd Xmkl = 0.5*(Xk + Xl);
    Eigen::VectorXd Xmli = 0.5*(Xl + Xi);

    //In-plane vectors.
    Eigen::Vector3d vn1 = Xmjk - Xmli;
    Eigen::Vector3d vn2 = Xmkl - Xmij;

    //Local axis definition.
    Eigen::Vector3d v1;
    Eigen::Vector3d v2;
    Eigen::Vector3d v3;

    //Local axis 3.
    v3 = vn1.cross(vn2);
    v3 = v3/v3.norm();

    //Local Axis 1.
    if(fabs(v3(2)) > TOL){
        v1 << v3(2), 0.0, -v3(0);
        v1 = v1/v1.norm();
    }
    else{
        v1 << -v3(1), v3(0), 0.0;
        v1 = v1/v1.norm();
    }

    //Local Axis 2.
    v2 = v3.cross(v1);
    v2 = v2/v2.norm();

    //Element projection axis.
    Eigen::MatrixXd rotationAxes(3,3);

    rotationAxes << v1(0), v1(1), v1(2),
                    v2(0), v2(1), v2(2),
                    v3(0), v3(1), v3(2); 

    return rotationAxes;
}

//Compute/update the node coordinates in axis-1-2-3 of the element.
Eigen::MatrixXd 
lin3DShell4::ComputeLocalCoordinates() const{
    //Rotated coordinates.
    Eigen::MatrixXd xyz(3,4);
    xyz << theNodes[0]->GetCoordinates(), 
           theNodes[1]->GetCoordinates(), 
           theNodes[2]->GetCoordinates(), 
           theNodes[3]->GetCoordinates(); 

    Eigen::MatrixXd Tr = ComputeRotation();
    Eigen::MatrixXd xylocal = Tr*xyz;

    return xylocal;
}

//Compute/update the local axis-1-2-3 of the element.
Eigen::MatrixXd 
lin3DShell4::ComputeLocalAxes() const{
    //Compute the projection axis.
    Eigen::MatrixXd Tr = ComputeRotation();

    //The global axes transformation. 
    Eigen::MatrixXd localAxes(24,24);

    localAxes << Tr(0,0), Tr(0,1), Tr(0,2),   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                 Tr(1,0), Tr(1,1), Tr(1,2),   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                 Tr(2,0), Tr(2,1), Tr(2,2),   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                   0.0,   0.0,   0.0, Tr(0,0), Tr(0,1), Tr(0,2),   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                   0.0,   0.0,   0.0, Tr(1,0), Tr(1,1), Tr(1,2),   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                   0.0,   0.0,   0.0, Tr(2,0), Tr(2,1), Tr(2,2),   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, Tr(0,0), Tr(0,1), Tr(0,2),   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, Tr(1,0), Tr(1,1), Tr(1,2),   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, Tr(2,0), Tr(2,1), Tr(2,2),   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, Tr(0,0), Tr(0,1), Tr(0,2),   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, Tr(1,0), Tr(1,1), Tr(1,2),   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, Tr(2,0), Tr(2,1), Tr(2,2),   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, Tr(0,0), Tr(0,1), Tr(0,2),   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, Tr(1,0), Tr(1,1), Tr(1,2),   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, Tr(2,0), Tr(2,1), Tr(2,2),   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, Tr(0,0), Tr(0,1), Tr(0,2),   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, Tr(1,0), Tr(1,1), Tr(1,2),   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, Tr(2,0), Tr(2,1), Tr(2,2),   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, Tr(0,0), Tr(0,1), Tr(0,2),   0.0,   0.0,   0.0,
                   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, Tr(1,0), Tr(1,1), Tr(1,2),   0.0,   0.0,   0.0,
                   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, Tr(2,0), Tr(2,1), Tr(2,2),   0.0,   0.0,   0.0, 
                   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, Tr(0,0), Tr(0,1), Tr(0,2),
                   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, Tr(1,0), Tr(1,1), Tr(1,2),
                   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, Tr(2,0), Tr(2,1), Tr(2,2);

    return localAxes;
}

//Update strain in the element.
Eigen::VectorXd 
lin3DShell4::ComputeStrain(Eigen::MatrixXd &BM12, Eigen::MatrixXd &xyloc, double ri, double si){
    //Compute the projection axis.
    Eigen::MatrixXd Tr = ComputeRotation();

    //The global to local axes transformation. 
    Eigen::MatrixXd localAxes(6,6); 
    localAxes << Tr(0,0), Tr(0,1), Tr(0,2),   0.0,     0.0,     0.0,
                 Tr(1,0), Tr(1,1), Tr(1,2),   0.0,     0.0,     0.0, 
                 Tr(2,0), Tr(2,1), Tr(2,2),   0.0,     0.0,     0.0, 
                     0.0,   0.0,     0.0,   Tr(0,0), Tr(0,1), Tr(0,2),
                     0.0,   0.0,     0.0,   Tr(1,0), Tr(1,1), Tr(1,2),
                     0.0,   0.0,     0.0,   Tr(2,0), Tr(2,1), Tr(2,2);

    //Gets the element coordinates in deformed configuration. 
    Eigen::VectorXd U1 = localAxes*(theNodes[0]->GetDisplacements() + theNodes[0]->GetIncrementalDisplacements());
    Eigen::VectorXd U2 = localAxes*(theNodes[1]->GetDisplacements() + theNodes[1]->GetIncrementalDisplacements());
    Eigen::VectorXd U3 = localAxes*(theNodes[2]->GetDisplacements() + theNodes[2]->GetIncrementalDisplacements());
    Eigen::VectorXd U4 = localAxes*(theNodes[3]->GetDisplacements() + theNodes[3]->GetIncrementalDisplacements());

    //Membrane local displacements.
    Eigen::VectorXd Um(12);
    Um << U1(0), U1(1), U1(5), U2(0), U2(1), U2(5), U3(0), U3(1), U3(5), U4(0), U4(1), U4(5); 

    //Plate local displacements.
    Eigen::VectorXd Up(12);
    Up << U1(2), U1(3), U1(4), U2(2), U2(3), U2(4), U3(2), U3(3), U3(4), U4(2), U4(3), U4(4); 

    //Strain values for membrane and plate effects.
    Eigen::MatrixXd Jij(2, 2);
    Eigen::MatrixXd Bij(3,12);
    Eigen::MatrixXd Pij(1,12);

    //Membrane generalised stresses at Gauss point.
    ComputeMembraneStrainDisplacementMatrix(ri, si, xyloc, BM12, Bij, Pij, Jij);
    Eigen::VectorXd Sm = Bij*Um; 

    //Plate generalised stresses at each Gauss point.
    ComputePlateStrainDisplacementMatrix(ri, si, xyloc, Bij, Jij);
    Eigen::VectorXd Sp = Bij*Up; 

    //The strain vector Gauss point.
    Eigen::VectorXd Strain(6);
    Strain << Sm, Sp;

    return Strain;
}

//Update strain rate in the element.
Eigen::VectorXd 
lin3DShell4::ComputeStrainRate(Eigen::MatrixXd &BM12, Eigen::MatrixXd &xyloc, double ri, double si){
    //Compute the projection axis.
    Eigen::MatrixXd Tr = ComputeRotation();

    //The global to local axes transformation. 
    Eigen::MatrixXd localAxes(6,6); 
    localAxes << Tr(0,0), Tr(0,1), Tr(0,2),   0.0,     0.0,     0.0,
                 Tr(1,0), Tr(1,1), Tr(1,2),   0.0,     0.0,     0.0, 
                 Tr(2,0), Tr(2,1), Tr(2,2),   0.0,     0.0,     0.0, 
                     0.0,   0.0,     0.0,   Tr(0,0), Tr(0,1), Tr(0,2),
                     0.0,   0.0,     0.0,   Tr(1,0), Tr(1,1), Tr(1,2),
                     0.0,   0.0,     0.0,   Tr(2,0), Tr(2,1), Tr(2,2);

    //Gets the element coordinates in deformed configuration. 
    Eigen::VectorXd V1 = localAxes*theNodes[0]->GetVelocities();
    Eigen::VectorXd V2 = localAxes*theNodes[1]->GetVelocities();
    Eigen::VectorXd V3 = localAxes*theNodes[2]->GetVelocities();
    Eigen::VectorXd V4 = localAxes*theNodes[3]->GetVelocities();

    //Membrane local displacements.
    Eigen::VectorXd Vm(12);
    Vm << V1(0), V1(1), V1(5), V2(0), V2(1), V2(5), V3(0), V3(1), V3(5), V4(0), V4(1), V4(5); 

    //Plate local displacements.
    Eigen::VectorXd Vp(12);
    Vp << V1(2), V1(3), V1(4), V2(2), V2(3), V2(4), V3(2), V3(3), V3(4), V4(2), V4(3), V4(4); 

    //Strain values for membrane and plate effects.
    Eigen::MatrixXd Jij(2, 2);
    Eigen::MatrixXd Bij(3,12);
    Eigen::MatrixXd Pij(1,12);

    //Membrane generalised stresses at Gauss point.
    ComputeMembraneStrainDisplacementMatrix(ri, si, xyloc, BM12, Bij, Pij, Jij);
    Eigen::VectorXd Sm = Bij*Vm; 

    //Plate generalised stresses at each Gauss point.
    ComputePlateStrainDisplacementMatrix(ri, si, xyloc, Bij, Jij);
    Eigen::VectorXd Sp = Bij*Vp; 

    //The strain vector Gauss point.
    Eigen::VectorXd StrainRate(6);
    StrainRate << Sm, Sp;

    return StrainRate;
}

//Correction matrix to Enforce Constant Tension for membrane effect (wilson).
double 
lin3DShell4::ComputeTotalMass(){
    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Quad", wi, xi);

    //Gets Node Local coordinates.
    Eigen::MatrixXd xyloc = ComputeLocalCoordinates();

    //Constant Tension matrix component.
    double Mass = 0.0;
    double x1 = xyloc(0,0); double y1 = xyloc(1,0);
    double x2 = xyloc(0,1); double y2 = xyloc(1,1);
    double x3 = xyloc(0,2); double y3 = xyloc(1,2);
    double x4 = xyloc(0,3); double y4 = xyloc(1,3);

    //Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
        //Components Node Coordinates.
        double p = xi(i,0); 
        double q = xi(i,1);

        //Jacobian Matrix.
        Eigen::MatrixXd Jij(2,2);
        Jij(0,0) = -1.0/4.0*x1 + 1.0/4.0*x1*q + 1.0/4.0*x2 - 1.0/4.0*x2*q + 1.0/4.0*x3 + 1.0/4.0*x3*q - 1.0/4.0*x4*q - 1.0/4.0*x4;
        Jij(0,1) = -1.0/4.0*y1 + 1.0/4.0*y1*q + 1.0/4.0*y2 - 1.0/4.0*q*y2 + 1.0/4.0*y3 + 1.0/4.0*q*y3 - 1.0/4.0*y4*q - 1.0/4.0*y4;
        Jij(1,0) = -1.0/4.0*x1 + 1.0/4.0*x1*p - 1.0/4.0*x2 - 1.0/4.0*x2*p + 1.0/4.0*x3 + 1.0/4.0*x3*p - 1.0/4.0*x4*p + 1.0/4.0*x4;
        Jij(1,1) = -1.0/4.0*y1 + 1.0/4.0*y1*p - 1.0/4.0*y2 - 1.0/4.0*y2*p + 1.0/4.0*y3 + 1.0/4.0*y3*p - 1.0/4.0*y4*p + 1.0/4.0*y4;

        //Gets material tangent matrix at Gauss point.
        Eigen::MatrixXd Rho = theSection[i]->GetDensity();

        //Numerical integration.
        Mass += wi(i)*Rho(0,0)*fabs(Jij.determinant());
    }

    return Mass;
}

//Correction matrix to Enforce Constant Tension for membrane effect (wilson).
Eigen::MatrixXd 
lin3DShell4::ComputeConstantTensionMatrix() {
    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Quad", wi, xi);

    //Gets Node Local coordinates.
    Eigen::MatrixXd xyloc = ComputeLocalCoordinates();

    //Constant Tension matrix component.
    double Area = 0.0;
    Eigen::MatrixXd BM12(3,4); 
    BM12.fill(0.0);

    //Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
        //Compute Constant Tension Matrix and Jacobian at Gauss Point.
        Eigen::MatrixXd Jij(2,2);
        Eigen::MatrixXd Bij(3,4);
        ConstantTensionMatrix(xi(i,0), xi(i,1), xyloc, Bij, Jij);

        //Numerical integration.
        Area += wi(i)*fabs(Jij.determinant());
        BM12 += wi(i)*fabs(Jij.determinant())*Bij;
    }

    //Correction by Area.
    BM12 = 1.0/Area*BM12;

    return BM12;
}

//Assembles the membrane and plate effects into a single matrix.
//Plate (dofs) : [3 4 5 9 10 11 15 16 17 21 22 23], and Membrane (dofs) : [1 2 6 7 8 12 13 14 18 19 20 24]
void
lin3DShell4::AssemblePlateMembraneEffects(Eigen::MatrixXd &A, const Eigen::MatrixXd &Am, const Eigen::MatrixXd &Ap){
    //The assembled shell matrix.
    A <<  Am( 0,0), Am( 0,1),      0.0,      0.0,      0.0, Am( 0,2), Am( 0,3), Am( 0,4),      0.0,      0.0,      0.0, Am( 0,5), Am( 0,6), Am( 0,7),      0.0,      0.0,      0.0, Am( 0,8), Am( 0,9),Am( 0,10),      0.0,      0.0,      0.0, Am( 0,11), 
          Am( 1,0), Am( 1,1),      0.0,      0.0,      0.0, Am( 1,2), Am( 1,3), Am( 1,4),      0.0,      0.0,      0.0, Am( 1,5), Am( 1,6), Am( 1,7),      0.0,      0.0,      0.0, Am( 1,8), Am( 1,9),Am( 1,10),      0.0,      0.0,      0.0, Am( 1,11), 
               0.0,      0.0, Ap( 0,0), Ap( 0,1), Ap( 0,2),      0.0,      0.0,      0.0, Ap( 0,3), Ap( 0,4), Ap( 0,5),      0.0,      0.0,      0.0, Ap( 0,6), Ap( 0,7), Ap( 0,8),      0.0,      0.0,      0.0, Ap( 0,9), Ap(0,10), Ap(0,11),       0.0, 
               0.0,      0.0, Ap( 1,0), Ap( 1,1), Ap( 1,2),      0.0,      0.0,      0.0, Ap( 1,3), Ap( 1,4), Ap( 1,5),      0.0,      0.0,      0.0, Ap( 1,6), Ap( 1,7), Ap( 1,8),      0.0,      0.0,      0.0, Ap( 1,9), Ap(1,10), Ap(1,11),       0.0, 
               0.0,      0.0, Ap( 2,0), Ap( 2,1), Ap( 2,2),      0.0,      0.0,      0.0, Ap( 2,3), Ap( 2,4), Ap( 2,5),      0.0,      0.0,      0.0, Ap( 2,6), Ap( 2,7), Ap( 2,8),      0.0,      0.0,      0.0, Ap( 2,9), Ap(2,10), Ap(2,11),       0.0, 
          Am( 2,0), Am( 2,1),      0.0,      0.0,      0.0, Am( 2,2), Am( 2,3), Am( 2,4),      0.0,      0.0,      0.0, Am( 2,5), Am( 2,6), Am( 2,7),      0.0,      0.0,      0.0, Am( 2,8), Am( 2,9),Am( 2,10),      0.0,      0.0,      0.0, Am( 2,11), 
          Am( 3,0), Am( 3,1),      0.0,      0.0,      0.0, Am( 3,2), Am( 3,3), Am( 3,4),      0.0,      0.0,      0.0, Am( 3,5), Am( 3,6), Am( 3,7),      0.0,      0.0,      0.0, Am( 3,8), Am( 3,9),Am( 3,10),      0.0,      0.0,      0.0, Am( 3,11), 
          Am( 4,0), Am( 4,1),      0.0,      0.0,      0.0, Am( 4,2), Am( 4,3), Am( 4,4),      0.0,      0.0,      0.0, Am( 4,5), Am( 4,6), Am( 4,7),      0.0,      0.0,      0.0, Am( 4,8), Am( 4,9),Am( 4,10),      0.0,      0.0,      0.0, Am( 4,11), 
               0.0,      0.0, Ap( 3,0), Ap( 3,1), Ap( 3,2),      0.0,      0.0,      0.0, Ap( 3,3), Ap( 3,4), Ap( 3,5),      0.0,      0.0,      0.0, Ap( 3,6), Ap( 3,7), Ap( 3,8),      0.0,      0.0,      0.0, Ap( 3,9), Ap(3,10), Ap(3,11),       0.0, 
               0.0,      0.0, Ap( 4,0), Ap( 4,1), Ap( 4,2),      0.0,      0.0,      0.0, Ap( 4,3), Ap( 4,4), Ap( 4,5),      0.0,      0.0,      0.0, Ap( 4,6), Ap( 4,7), Ap( 4,8),      0.0,      0.0,      0.0, Ap( 4,9), Ap(4,10), Ap(4,11),       0.0, 
               0.0,      0.0, Ap( 5,0), Ap( 5,1), Ap( 5,2),      0.0,      0.0,      0.0, Ap( 5,3), Ap( 5,4), Ap( 5,5),      0.0,      0.0,      0.0, Ap( 5,6), Ap( 5,7), Ap( 5,8),      0.0,      0.0,      0.0, Ap( 5,9), Ap(5,10), Ap(5,11),       0.0, 
          Am( 5,0), Am( 5,1),      0.0,      0.0,      0.0, Am( 5,2), Am( 5,3), Am( 5,4),      0.0,      0.0,      0.0, Am( 5,5), Am( 5,6), Am( 5,7),      0.0,      0.0,      0.0, Am( 5,8), Am( 5,9),Am( 5,10),      0.0,      0.0,      0.0, Am( 5,11), 
          Am( 6,0), Am( 6,1),      0.0,      0.0,      0.0, Am( 6,2), Am( 6,3), Am( 6,4),      0.0,      0.0,      0.0, Am( 6,5), Am( 6,6), Am( 6,7),      0.0,      0.0,      0.0, Am( 6,8), Am( 6,9),Am( 6,10),      0.0,      0.0,      0.0, Am( 6,11), 
          Am( 7,0), Am( 7,1),      0.0,      0.0,      0.0, Am( 7,2), Am( 7,3), Am( 7,4),      0.0,      0.0,      0.0, Am( 7,5), Am( 7,6), Am( 7,7),      0.0,      0.0,      0.0, Am( 7,8), Am( 7,9),Am( 7,10),      0.0,      0.0,      0.0, Am( 7,11), 
               0.0,      0.0, Ap( 6,0), Ap( 6,1), Ap( 6,2),      0.0,      0.0,      0.0, Ap( 6,3), Ap( 6,4), Ap( 6,5),      0.0,      0.0,      0.0, Ap( 6,6), Ap( 6,7), Ap( 6,8),      0.0,      0.0,      0.0, Ap( 6,9), Ap(6,10), Ap(6,11),       0.0, 
               0.0,      0.0, Ap( 7,0), Ap( 7,1), Ap( 7,2),      0.0,      0.0,      0.0, Ap( 7,3), Ap( 7,4), Ap( 7,5),      0.0,      0.0,      0.0, Ap( 7,6), Ap( 7,7), Ap( 7,8),      0.0,      0.0,      0.0, Ap( 7,9), Ap(7,10), Ap(7,11),       0.0, 
               0.0,      0.0, Ap( 8,0), Ap( 8,1), Ap( 8,2),      0.0,      0.0,      0.0, Ap( 8,3), Ap( 8,4), Ap( 8,5),      0.0,      0.0,      0.0, Ap( 8,6), Ap( 8,7), Ap( 8,8),      0.0,      0.0,      0.0, Ap( 8,9), Ap(8,10), Ap(8,11),       0.0, 
          Am( 8,0), Am( 8,1),      0.0,      0.0,      0.0, Am( 8,2), Am( 8,3), Am( 8,4),      0.0,      0.0,      0.0, Am( 8,5), Am( 8,6), Am( 8,7),      0.0,      0.0,      0.0, Am( 8,8), Am( 8,9),Am( 8,10),      0.0,      0.0,      0.0, Am( 8,11), 
          Am( 9,0), Am( 9,1),      0.0,      0.0,      0.0, Am( 9,2), Am( 9,3), Am( 9,4),      0.0,      0.0,      0.0, Am( 9,5), Am( 9,6), Am( 9,7),      0.0,      0.0,      0.0, Am( 9,8), Am( 9,9),Am( 9,10),      0.0,      0.0,      0.0, Am( 9,11), 
          Am(10,0), Am(10,1),      0.0,      0.0,      0.0, Am(10,2), Am(10,3), Am(10,4),      0.0,      0.0,      0.0, Am(10,5), Am(10,6), Am(10,7),      0.0,      0.0,      0.0, Am(10,8), Am(10,9),Am(10,10),      0.0,      0.0,      0.0, Am(10,11), 
               0.0,      0.0, Ap( 9,0), Ap( 9,1), Ap( 9,2),      0.0,      0.0,      0.0, Ap( 9,3), Ap( 9,4), Ap( 9,5),      0.0,      0.0,      0.0, Ap( 9,6), Ap( 9,7), Ap( 9,8),      0.0,      0.0,      0.0, Ap( 9,9), Ap(9,10), Ap(9,11),       0.0, 
               0.0,      0.0, Ap(10,0), Ap(10,1), Ap(10,2),      0.0,      0.0,      0.0, Ap(10,3), Ap(10,4), Ap(10,5),      0.0,      0.0,      0.0, Ap(10,6), Ap(10,7), Ap(10,8),      0.0,      0.0,      0.0, Ap(10,9),Ap(10,10),Ap(10,11),       0.0, 
               0.0,      0.0, Ap(11,0), Ap(11,1), Ap(11,2),      0.0,      0.0,      0.0, Ap(11,3), Ap(11,4), Ap(11,5),      0.0,      0.0,      0.0, Ap(11,6), Ap(11,7), Ap(11,8),      0.0,      0.0,      0.0, Ap(11,9),Ap(11,10),Ap(11,11),       0.0, 
          Am(11,0), Am(11,1),      0.0,      0.0,      0.0, Am(11,2), Am(11,3), Am(11,4),      0.0,      0.0,      0.0, Am(11,5), Am(11,6), Am(11,7),      0.0,      0.0,      0.0, Am(11,8), Am(11,9),Am(11,10),      0.0,      0.0,      0.0, Am(11,11);
}

//Evaluates the shape function matrix for load at a given Gauss point.
void
lin3DShell4::ComputePlateShapeFunctionMatrix(const double ri, const double si, const Eigen::MatrixXd &xyloc, Eigen::MatrixXd &Hij, Eigen::MatrixXd &Jij) {
    //Shape function 'N' coefficients.
    double N1 = 1.0/4.0*(1.0 - ri)*(1.0 - si);
    double N2 = 1.0/4.0*(1.0 + ri)*(1.0 - si);
    double N3 = 1.0/4.0*(1.0 + ri)*(1.0 + si);
    double N4 = 1.0/4.0*(1.0 - ri)*(1.0 + si);

    //Jacobian Matrix definition:
    double J11 = -1.0/4.0*(1.0 - si)*xyloc(0,0) + 1.0/4.0*(1.0 - si)*xyloc(0,1) + 1.0/4.0*(1.0 + si)*xyloc(0,2) - 1.0/4.0*(1.0 + si)*xyloc(0,3);
    double J12 = -1.0/4.0*(1.0 - si)*xyloc(1,0) + 1.0/4.0*(1.0 - si)*xyloc(1,1) + 1.0/4.0*(1.0 + si)*xyloc(1,2) - 1.0/4.0*(1.0 + si)*xyloc(1,3); 
    double J21 = -1.0/4.0*(1.0 - ri)*xyloc(0,0) - 1.0/4.0*(1.0 + ri)*xyloc(0,1) + 1.0/4.0*(1.0 + ri)*xyloc(0,2) + 1.0/4.0*(1.0 - ri)*xyloc(0,3); 
    double J22 = -1.0/4.0*(1.0 - ri)*xyloc(1,0) - 1.0/4.0*(1.0 + ri)*xyloc(1,1) + 1.0/4.0*(1.0 + ri)*xyloc(1,2) + 1.0/4.0*(1.0 - ri)*xyloc(1,3); 

    //Jacobian matrix:
    Jij << J11, J12,
           J21, J22;

    //Shape function matrix.
    Hij <<  N1, 0.0, 0.0,  N2, 0.0, 0.0,  N3, 0.0, 0.0,  N4, 0.0, 0.0,
           0.0, 0.0,  N1, 0.0, 0.0,  N2, 0.0, 0.0,  N3, 0.0, 0.0,  N4,
           0.0, -N1, 0.0, 0.0, -N2, 0.0, 0.0, -N3, 0.0, 0.0, -N4, 0.0;
}

//Evaluates the shape function matrix for plate at a given Gauss point.
void
lin3DShell4::ComputeLoadShapeFunctionMatrix(const double ri, const double si, const Eigen::MatrixXd &xyloc, Eigen::MatrixXd &Nij, Eigen::MatrixXd &Jij) {
    //Shape function 'N' coefficients.
    double N1 = 1.0/4.0*(1.0 - ri)*(1.0 - si);
    double N2 = 1.0/4.0*(1.0 + ri)*(1.0 - si);
    double N3 = 1.0/4.0*(1.0 + ri)*(1.0 + si);
    double N4 = 1.0/4.0*(1.0 - ri)*(1.0 + si);

    //Jacobian Matrix definition:
    double J11 = -1.0/4.0*(1.0 - si)*xyloc(0,0) + 1.0/4.0*(1.0 - si)*xyloc(0,1) + 1.0/4.0*(1.0 + si)*xyloc(0,2) - 1.0/4.0*(1.0 + si)*xyloc(0,3);
    double J12 = -1.0/4.0*(1.0 - si)*xyloc(1,0) + 1.0/4.0*(1.0 - si)*xyloc(1,1) + 1.0/4.0*(1.0 + si)*xyloc(1,2) - 1.0/4.0*(1.0 + si)*xyloc(1,3); 
    double J21 = -1.0/4.0*(1.0 - ri)*xyloc(0,0) - 1.0/4.0*(1.0 + ri)*xyloc(0,1) + 1.0/4.0*(1.0 + ri)*xyloc(0,2) + 1.0/4.0*(1.0 - ri)*xyloc(0,3); 
    double J22 = -1.0/4.0*(1.0 - ri)*xyloc(1,0) - 1.0/4.0*(1.0 + ri)*xyloc(1,1) + 1.0/4.0*(1.0 + ri)*xyloc(1,2) + 1.0/4.0*(1.0 - ri)*xyloc(1,3); 

    //Jacobian matrix:
    Jij << J11, J12,
           J21, J22;

    //Shape function matrix.
    Nij <<  N1, 0.0, 0.0,  N2, 0.0, 0.0,  N3, 0.0, 0.0,  N4, 0.0, 0.0,
           0.0,  N1, 0.0, 0.0,  N2, 0.0, 0.0,  N3, 0.0, 0.0,  N4, 0.0,
           0.0, 0.0,  N1, 0.0, 0.0,  N2, 0.0, 0.0,  N3, 0.0, 0.0,  N4;
}

//Evaluates the shape function matrix for membrane at a given Gauss point.
void 
lin3DShell4::ComputeMembraneShapeFunctionMatrix(const double ri, const double si, const Eigen::MatrixXd &xyloc, Eigen::MatrixXd &Hij, Eigen::MatrixXd &Jij) {
    //Shape function 'N' coefficients.
    double N1 = 1.0/4.00*(1.0 - ri)*(1.0 - si);
    double N2 = 1.0/4.00*(1.0 + ri)*(1.0 - si);
    double N3 = 1.0/4.00*(1.0 + ri)*(1.0 + si);
    double N4 = 1.0/4.00*(1.0 - ri)*(1.0 + si);
    double N5 = 1.0/16.0*(1.0 - ri*ri)*(1.0 - si);
    double N6 = 1.0/4.00*(1.0 + ri)*(1.0 - si*si);
    double N7 = 1.0/16.0*(1.0 - ri*ri)*(1.0 + si);
    double N8 = 1.0/4.00*(1.0 - ri)*(1.0 - si*si);

    //Shape function 'H' coefficients.
    double dx1 = xyloc(0,1) - xyloc(0,0); 
    double dy1 = xyloc(1,1) - xyloc(1,0);
    double dx2 = xyloc(0,2) - xyloc(0,1); 
    double dy2 = xyloc(1,2) - xyloc(1,1);
    double dx3 = xyloc(0,3) - xyloc(0,2); 
    double dy3 = xyloc(1,3) - xyloc(1,2);
    double dx4 = xyloc(0,0) - xyloc(0,3); 
    double dy4 = xyloc(1,0) - xyloc(1,3);

    double H1 = N8*dy4 - N5*dy1;  
    double H3 = N5*dy1 - N6*dy2;  
    double H5 = N6*dy2 - N7*dy3;  
    double H7 = N7*dy3 - N8*dy4;
    double H2 = N5*dx1 - N8*dx4;  
    double H4 = N6*dx2 - N5*dx1;  
    double H6 = N7*dx3 - N6*dx2;  
    double H8 = N8*dx4 - N7*dx3;

    //Jacobian Matrix definition:
    double J11 = -1.0/4.0*(1.0 - si)*xyloc(0,0) + 1.0/4.0*(1.0 - si)*xyloc(0,1) + 1.0/4.0*(1.0 + si)*xyloc(0,2) - 1.0/4.0*(1.0 + si)*xyloc(0,3);
    double J12 = -1.0/4.0*(1.0 - si)*xyloc(1,0) + 1.0/4.0*(1.0 - si)*xyloc(1,1) + 1.0/4.0*(1.0 + si)*xyloc(1,2) - 1.0/4.0*(1.0 + si)*xyloc(1,3); 
    double J21 = -1.0/4.0*(1.0 - ri)*xyloc(0,0) - 1.0/4.0*(1.0 + ri)*xyloc(0,1) + 1.0/4.0*(1.0 + ri)*xyloc(0,2) + 1.0/4.0*(1.0 - ri)*xyloc(0,3); 
    double J22 = -1.0/4.0*(1.0 - ri)*xyloc(1,0) - 1.0/4.0*(1.0 + ri)*xyloc(1,1) + 1.0/4.0*(1.0 + ri)*xyloc(1,2) + 1.0/4.0*(1.0 - ri)*xyloc(1,3); 

    //Jacobian matrix:
    Jij << J11, J12,
           J21, J22;

    //Shape function matrix.
    Hij <<  N1, 0.0,  H1,  N2, 0.0,  H3,  N3, 0.0,  H5,  N4, 0.0,  H7,
           0.0,  N1,  H2, 0.0,  N2,  H4, 0.0,  N3,  H6, 0.0,  N4,  H8,
           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
}

//Evaluates the strain-displacement matrix for a plate effect at a given Gauss point.
void 
lin3DShell4::ComputePlateStrainDisplacementMatrix(const double ri, const double si, const Eigen::MatrixXd &xyloc, Eigen::MatrixXd &Bij, Eigen::MatrixXd &Jij) {
    //Quadrature point coordinates.
    double p = ri; 
    double q = si;

    //Components Node Coordinates.
    double x1 = xyloc(0,0); double y1 = xyloc(1,0);
    double x2 = xyloc(0,1); double y2 = xyloc(1,1);
    double x3 = xyloc(0,2); double y3 = xyloc(1,2);
    double x4 = xyloc(0,3); double y4 = xyloc(1,3);

    //Shape functions in local element coordinates.
    Eigen::MatrixXd HXpq(2,12);
    Eigen::MatrixXd HYpq(2,12);

    HXpq(0, 0) = -3.0/2.0*(-x2+x1)/((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))*p*(1.0-q)+3.0/4.0*(-x1+x4)/((x1-x4)*(x1-x4)+(y1-y4)*(y1-y4))*(1.0-q*q);
    HXpq(0, 1) = -(-3.0/4.0*x2+3.0/4.0*x1)*(y2-y1)/((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))*p*(1.0-q)-1.0/2.0*(-3.0/4.0*x1+3.0/4.0*x4)*(y1-y4)/((x1-x4)*(x1-x4)+(y1-y4)*(y1-y4))*(1.0-q*q);
    HXpq(0, 2) = -1.0/4.0*q+1.0/4.0*q*q-1.0/2.0*p*(1-q)+(1.0/2.0*(y2-y1)*(y2-y1)-1.0/4.0*(x2-x1)*(x2-x1))/((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))*p*(1.0-q)+1.0/2.0*(1.0/2.0*(y1-y4)*(y1-y4)-1.0/4.0*(x1-x4)*(x1-x4))/((x1-x4)*(x1-x4)+(y1-y4)*(y1-y4))*(1.0-q*q);
    HXpq(0, 3) =  3.0/4.0*(-x3+x2)/((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2))*(1.0-q*q)+3.0/2.0*(-x2+x1)/((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))*p*(1.0-q);
    HXpq(0, 4) =  1.0/2.0*(-3.0/4.0*x3+3.0/4.0*x2)*(y3-y2)/((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2))*(1.0-q*q)-(-3.0/4.0*x2+3.0/4.0*x1)*(y2-y1)/((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))*p*(1.0-q);
    HXpq(0, 5) =  1.0/4.0*q-1.0/2.0*p*(1.0-q)-1.0/4.0*q*q-1.0/2.0*(1.0/2.0*(y3-y2)*(y3-y2)-1.0/4.0*(x3-x2)*(x3-x2))/((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2))*(1.0-q*q)+(1.0/2.0*(y2-y1)*(y2-y1)-1.0/4.0*(x2-x1)*(x2-x1))/((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))*p*(1.0-q);
    HXpq(0, 6) = -3.0/2.0*(-x4+x3)/((x4-x3)*(x4-x3)+(y4-y3)*(y4-y3))*p*(1.0+q)-3.0/4.0*(-x3+x2)/((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2))*(1.0-q*q);
    HXpq(0, 7) = -(-3.0/4.0*x4+3.0/4.0*x3)*(y4-y3)/((x4-x3)*(x4-x3)+(y4-y3)*(y4-y3))*p*(1.0+q)+1.0/2.0*(-3.0/4.0*x3+3.0/4.0*x2)*(y3-y2)/((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2))*(1.0-q*q);
    HXpq(0, 8) = -1.0/4.0*q-1.0/4.0*q*q-1.0/2.0*p*(1.0+q)+(1.0/2.0*(y4-y3)*(y4-y3)-1.0/4.0*(x4-x3)*(x4-x3))/((x4-x3)*(x4-x3)+(y4-y3)*(y4-y3))*p*(1.0+q)-1.0/2.0*(1.0/2.0*(y3-y2)*(y3-y2)-1.0/4.0*(x3-x2)*(x3-x2))/((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2))*(1.0-q*q);
    HXpq(0, 9) = -3.0/4.0*(-x1+x4)/((x1-x4)*(x1-x4)+(y1-y4)*(y1-y4))*(1.0-q*q)+3.0/2.0*(-x4+x3)/((x4-x3)*(x4-x3)+(y4-y3)*(y4-y3))*p*(1.0+q);
    HXpq(0,10) = -1.0/2.0*(-3.0/4.0*x1+3.0/4.0*x4)*(y1-y4)/((x1-x4)*(x1-x4)+(y1-y4)*(y1-y4))*(1.0-q*q)-(-3.0/4.0*x4+3.0/4.0*x3)*(y4-y3)/((x4-x3)*(x4-x3)+(y4-y3)*(y4-y3))*p*(1.0+q);
    HXpq(0,11) =  1.0/4.0*q-1.0/2.0*p*(1.0+q)+1.0/4.0*q*q+1.0/2.0*(1.0/2.0*(y1-y4)*(y1-y4)-1.0/4.0*(x1-x4)*(x1-x4))/((x1-x4)*(x1-x4)+(y1-y4)*(y1-y4))*(1.0-q*q)+(1.0/2.0*(y4-y3)*(y4-y3)-1.0/4.0*(x4-x3)*(x4-x3))/((x4-x3)*(x4-x3)+(y4-y3)*(y4-y3))*p*(1.0+q);
    HXpq(1, 0) = -3.0/2.0*(-x2+x1)/((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))*(1.0/2.0-1.0/2.0*p*p)+3*(-x1+x4)/((x1-x4)*(x1-x4)+(y1-y4)*(y1-y4))*(1.0/2.0-1.0/2.0*p)*q;
    HXpq(1, 1) = -(-3.0/4.0*x2+3.0/4.0*x1)*(y2-y1)/((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))*(1.0/2.0-1.0/2.0*p*p)-2*(-3.0/4.0*x1+3.0/4.0*x4)*(y1-y4)/((x1-x4)*(x1-x4)+(y1-y4)*(y1-y4))*(1.0/2.0-1.0/2.0*p)*q;
    HXpq(1, 2) = -1.0/4.0*p-(1.0/2.0-1.0/2.0*p)*q+1.0/4.0*p*p+(1.0/2.0*(y2-y1)*(y2-y1)-1.0/4.0*(x2-x1)*(x2-x1))/((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))*(1.0/2.0-1.0/2.0*p*p)+2*(1.0/2.0*(y1-y4)*(y1-y4)-1.0/4.0*(x1-x4)*(x1-x4))/((x1-x4)*(x1-x4)+(y1-y4)*(y1-y4))*(1.0/2.0-1.0/2.0*p)*q;
    HXpq(1, 3) = -3.0*(-x3+x2)/((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2))*(1.0/2.0+1.0/2.0*p)*q+3.0/2.0*(-x2+x1)/((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))*(1.0/2.0-1.0/2.0*p*p);
    HXpq(1, 4) = -2.0*(-3.0/4.0*x3+3.0/4.0*x2)*(y3-y2)/((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2))*(1.0/2.0+1.0/2.0*p)*q-(-3.0/4.0*x2+3.0/4.0*x1)*(y2-y1)/((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))*(1.0/2.0-1.0/2.0*p*p);
    HXpq(1, 5) =  1.0/4.0*p+1.0/4.0*p*p-(1.0/2.0+1.0/2.0*p)*q+2*(1.0/2.0*(y3-y2)*(y3-y2)-1.0/4.0*(x3-x2)*(x3-x2))/((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2))*(1.0/2.0+1.0/2.0*p)*q+(1.0/2.0*(y2-y1)*(y2-y1)-1.0/4.0*(x2-x1)*(x2-x1))/((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))*(1.0/2.0-1.0/2.0*p*p);
    HXpq(1, 6) =  3.0/2.0*(-x4+x3)/((x4-x3)*(x4-x3)+(y4-y3)*(y4-y3))*(1.0/2.0-1.0/2.0*p*p)+3*(-x3+x2)/((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2))*(1.0/2.0+1.0/2.0*p)*q;
    HXpq(1, 7) = (-3.0/4.0*x4+3.0/4.0*x3)*(y4-y3)/((x4-x3)*(x4-x3)+(y4-y3)*(y4-y3))*(1.0/2.0-1.0/2.0*p*p)-2*(-3.0/4.0*x3+3.0/4.0*x2)*(y3-y2)/((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2))*(1.0/2.0+1.0/2.0*p)*q;
    HXpq(1, 8) = -1.0/4.0*p-(1.0/2.0+1.0/2.0*p)*q-1.0/4.0*p*p-(1.0/2.0*(y4-y3)*(y4-y3)-1.0/4.0*(x4-x3)*(x4-x3))/((x4-x3)*(x4-x3)+(y4-y3)*(y4-y3))*(1.0/2.0-1.0/2.0*p*p)+2*(1.0/2.0*(y3-y2)*(y3-y2)-1.0/4.0*(x3-x2)*(x3-x2))/((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2))*(1.0/2.0+1.0/2.0*p)*q;
    HXpq(1, 9) = -3.0*(-x1+x4)/((x1-x4)*(x1-x4)+(y1-y4)*(y1-y4))*(1.0/2.0-1.0/2.0*p)*q-3.0/2.0*(-x4+x3)/((x4-x3)*(x4-x3)+(y4-y3)*(y4-y3))*(1.0/2.0-1.0/2.0*p*p);
    HXpq(1,10) = -2.0*(-3.0/4.0*x1+3.0/4.0*x4)*(y1-y4)/((x1-x4)*(x1-x4)+(y1-y4)*(y1-y4))*(1.0/2.0-1.0/2.0*p)*q+(-3.0/4.0*x4+3.0/4.0*x3)*(y4-y3)/((x4-x3)*(x4-x3)+(y4-y3)*(y4-y3))*(1.0/2.0-1.0/2.0*p*p);
    HXpq(1,11) =  1.0/4.0*p-1.0/4.0*p*p-(1.0/2.0-1.0/2.0*p)*q+2*(1.0/2.0*(y1-y4)*(y1-y4)-1.0/4.0*(x1-x4)*(x1-x4))/((x1-x4)*(x1-x4)+(y1-y4)*(y1-y4))*(1.0/2.0-1.0/2.0*p)*q-(1.0/2.0*(y4-y3)*(y4-y3)-1.0/4.0*(x4-x3)*(x4-x3))/((x4-x3)*(x4-x3)+(y4-y3)*(y4-y3))*(1.0/2.0-1.0/2.0*p*p);
    
    HYpq(0, 0) = -3.0/2.0*(-y2+y1)/((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))*p*(1.0-q)+3.0/4.0*(-y1+y4)/((x1-x4)*(x1-x4)+(y1-y4)*(y1-y4))*(1.0-q*q);
    HYpq(0, 1) =  1.0/4.0*q-1.0/4.0*q*q+1.0/2.0*p*(1.0-q)-(1.0/2.0*(x2-x1)*(x2-x1)-1.0/4.0*(y2-y1)*(y2-y1))/((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))*p*(1.0-q)-1.0/2.0*(1.0/2.0*(x1-x4)*(x1-x4)-1.0/4.0*(y1-y4)*(y1-y4))/((x1-x4)*(x1-x4)+(y1-y4)*(y1-y4))*(1.0-q*q);
    HYpq(0, 2) = (-3.0/4.0*x2+3.0/4.0*x1)*(y2-y1)/((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))*p*(1.0-q)+1.0/2.0*(-3.0/4.0*x1+3.0/4.0*x4)*(y1-y4)/((x1-x4)*(x1-x4)+(y1-y4)*(y1-y4))*(1.0-q*q);
    HYpq(0, 3) =  3.0/4.0*(-y3+y2)/((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2))*(1.0-q*q)+3.0/2.0*(-y2+y1)/((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))*p*(1.0-q);
    HYpq(0, 4) = -1.0/4.0*q+1.0/2.0*p*(1.0-q)+1.0/4.0*q*q+1.0/2.0*(1.0/2.0*(x3-x2)*(x3-x2)-1.0/4.0*(y3-y2)*(y3-y2))/((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2))*(1.0-q*q)-(1.0/2.0*(x2-x1)*(x2-x1)-1.0/4.0*(y2-y1)*(y2-y1))/((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))*p*(1.0-q);
    HYpq(0, 5) = -1.0/2.0*(-3.0/4.0*x3+3.0/4.0*x2)*(y3-y2)/((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2))*(1.0-q*q)+(-3.0/4.0*x2+3.0/4.0*x1)*(y2-y1)/((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))*p*(1.0-q);
    HYpq(0, 6) = -3.0/2.0*(-y4+y3)/((x4-x3)*(x4-x3)+(y4-y3)*(y4-y3))*p*(1.0+q)-3.0/4.0*(-y3+y2)/((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2))*(1.0-q*q);
    HYpq(0, 7) =  1.0/4.0*q+1.0/4.0*q*q+1.0/2.0*p*(1.0+q)-(1.0/2.0*(x4-x3)*(x4-x3)-1.0/4.0*(y4-y3)*(y4-y3))/((x4-x3)*(x4-x3)+(y4-y3)*(y4-y3))*p*(1.0+q)+1.0/2.0*(1.0/2.0*(x3-x2)*(x3-x2)-1.0/4.0*(y3-y2)*(y3-y2))/((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2))*(1.0-q*q);
    HYpq(0, 8) = (-3.0/4.0*x4+3.0/4.0*x3)*(y4-y3)/((x4-x3)*(x4-x3)+(y4-y3)*(y4-y3))*p*(1.0+q)-1.0/2.0*(-3.0/4.0*x3+3.0/4.0*x2)*(y3-y2)/((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2))*(1.0-q*q);
    HYpq(0, 9) = -3.0/4.0*(-y1+y4)/((x1-x4)*(x1-x4)+(y1-y4)*(y1-y4))*(1.0-q*q)+3.0/2.0*(-y4+y3)/((x4-x3)*(x4-x3)+(y4-y3)*(y4-y3))*p*(1.0+q);
    HYpq(0,10) = -1.0/4.0*q+1.0/2.0*p*(1.0+q)-1.0/4.0*q*q-1.0/2.0*(1.0/2.0*(x1-x4)*(x1-x4)-1.0/4.0*(y1-y4)*(y1-y4))/((x1-x4)*(x1-x4)+(y1-y4)*(y1-y4))*(1.0-q*q)-(1.0/2.0*(x4-x3)*(x4-x3)-1.0/4.0*(y4-y3)*(y4-y3))/((x4-x3)*(x4-x3)+(y4-y3)*(y4-y3))*p*(1.0+q);
    HYpq(0,11) =  1.0/2.0*(-3.0/4.0*x1+3.0/4.0*x4)*(y1-y4)/((x1-x4)*(x1-x4)+(y1-y4)*(y1-y4))*(1.0-q*q)+(-3.0/4.0*x4+3.0/4.0*x3)*(y4-y3)/((x4-x3)*(x4-x3)+(y4-y3)*(y4-y3))*p*(1.0+q);
    HYpq(1, 0) = -3.0/2.0*(-y2+y1)/((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))*(1.0/2.0-1.0/2.0*p*p)+3*(-y1+y4)/((x1-x4)*(x1-x4)+(y1-y4)*(y1-y4))*(1.0/2.0-1.0/2.0*p)*q;
    HYpq(1, 1) =  1.0/4.0*p+(1.0/2.0-1.0/2.0*p)*q-1.0/4.0*p*p-(1.0/2.0*(x2-x1)*(x2-x1)-1.0/4.0*(y2-y1)*(y2-y1))/((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))*(1.0/2.0-1.0/2.0*p*p)-2*(1.0/2.0*(x1-x4)*(x1-x4)-1.0/4.0*(y1-y4)*(y1-y4))/((x1-x4)*(x1-x4)+(y1-y4)*(y1-y4))*(1.0/2.0-1.0/2.0*p)*q;
    HYpq(1, 2) = (-3.0/4.0*x2+3.0/4.0*x1)*(y2-y1)/((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))*(1.0/2.0-1.0/2.0*p*p)+2*(-3.0/4.0*x1+3.0/4.0*x4)*(y1-y4)/((x1-x4)*(x1-x4)+(y1-y4)*(y1-y4))*(1.0/2.0-1.0/2.0*p)*q;
    HYpq(1, 3) = -3.0*(-y3+y2)/((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2))*(1.0/2.0+1.0/2.0*p)*q+3.0/2.0*(-y2+y1)/((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))*(1.0/2.0-1.0/2.0*p*p);
    HYpq(1, 4) = -1.0/4.0*p-1.0/4.0*p*p+(1.0/2.0+1.0/2.0*p)*q-2*(1.0/2.0*(x3-x2)*(x3-x2)-1.0/4.0*(y3-y2)*(y3-y2))/((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2))*(1.0/2.0+1.0/2.0*p)*q-(1.0/2.0*(x2-x1)*(x2-x1)-1.0/4.0*(y2-y1)*(y2-y1))/((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))*(1.0/2.0-1.0/2.0*p*p);
    HYpq(1, 5) =  2.0*(-3.0/4.0*x3+3.0/4.0*x2)*(y3-y2)/((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2))*(1.0/2.0+1.0/2.0*p)*q+(-3.0/4.0*x2+3.0/4.0*x1)*(y2-y1)/((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))*(1.0/2.0-1.0/2.0*p*p);
    HYpq(1, 6) =  3.0/2.0*(-y4+y3)/((x4-x3)*(x4-x3)+(y4-y3)*(y4-y3))*(1.0/2.0-1.0/2.0*p*p)+3*(-y3+y2)/((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2))*(1.0/2.0+1.0/2.0*p)*q;
    HYpq(1, 7) =  1.0/4.0*p+(1.0/2.0+1.0/2.0*p)*q+1.0/4.0*p*p+(1.0/2.0*(x4-x3)*(x4-x3)-1.0/4.0*(y4-y3)*(y4-y3))/((x4-x3)*(x4-x3)+(y4-y3)*(y4-y3))*(1.0/2.0-1.0/2.0*p*p)-2*(1.0/2.0*(x3-x2)*(x3-x2)-1.0/4.0*(y3-y2)*(y3-y2))/((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2))*(1.0/2.0+1.0/2.0*p)*q;
    HYpq(1, 8) = -1.0*(-3.0/4.0*x4+3.0/4.0*x3)*(y4-y3)/((x4-x3)*(x4-x3)+(y4-y3)*(y4-y3))*(1.0/2.0-1.0/2.0*p*p)+2*(-3.0/4.0*x3+3.0/4.0*x2)*(y3-y2)/((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2))*(1.0/2.0+1.0/2.0*p)*q;
    HYpq(1, 9) = -3.0*(-y1+y4)/((x1-x4)*(x1-x4)+(y1-y4)*(y1-y4))*(1.0/2.0-1.0/2.0*p)*q-3.0/2.0*(-y4+y3)/((x4-x3)*(x4-x3)+(y4-y3)*(y4-y3))*(1.0/2.0-1.0/2.0*p*p);
    HYpq(1,10) = -1.0/4.0*p+1.0/4.0*p*p+(1.0/2.0-1.0/2.0*p)*q-2*(1.0/2.0*(x1-x4)*(x1-x4)-1.0/4.0*(y1-y4)*(y1-y4))/((x1-x4)*(x1-x4)+(y1-y4)*(y1-y4))*(1.0/2.0-1.0/2.0*p)*q+(1.0/2.0*(x4-x3)*(x4-x3)-1.0/4.0*(y4-y3)*(y4-y3))/((x4-x3)*(x4-x3)+(y4-y3)*(y4-y3))*(1.0/2.0-1.0/2.0*p*p);
    HYpq(1,11) =  2.0*(-3.0/4.0*x1+3.0/4.0*x4)*(y1-y4)/((x1-x4)*(x1-x4)+(y1-y4)*(y1-y4))*(1.0/2.0-1.0/2.0*p)*q-(-3.0/4.0*x4+3.0/4.0*x3)*(y4-y3)/((x4-x3)*(x4-x3)+(y4-y3)*(y4-y3))*(1.0/2.0-1.0/2.0*p*p);

    //Jacobian Matrix.
    Jij(0,0) = -1.0/4.0*x1 + 1.0/4.0*x1*q + 1.0/4.0*x2 - 1.0/4.0*x2*q + 1.0/4.0*x3 + 1.0/4.0*x3*q - 1.0/4.0*x4*q - 1.0/4.0*x4;
    Jij(0,1) = -1.0/4.0*y1 + 1.0/4.0*y1*q + 1.0/4.0*y2 - 1.0/4.0*y2*q + 1.0/4.0*y3 + 1.0/4.0*y3*q - 1.0/4.0*y4*q - 1.0/4.0*y4;
    Jij(1,0) = -1.0/4.0*x1 + 1.0/4.0*x1*p - 1.0/4.0*x2 - 1.0/4.0*x2*p + 1.0/4.0*x3 + 1.0/4.0*x3*p - 1.0/4.0*x4*p + 1.0/4.0*x4;
    Jij(1,1) = -1.0/4.0*y1 + 1.0/4.0*y1*p - 1.0/4.0*y2 - 1.0/4.0*y2*p + 1.0/4.0*y3 + 1.0/4.0*y3*p - 1.0/4.0*y4*p + 1.0/4.0*y4;

    //Shape functions in global element coordinates.
    Eigen::MatrixXd HXxy = -1.0*Jij.inverse()*HXpq;
    Eigen::MatrixXd HYxy = -1.0*Jij.inverse()*HYpq;
    Eigen::MatrixXd HZxy = HXxy.row(1) + HYxy.row(0);

    Bij << HXxy(0,0), HXxy(0,1), HXxy(0,2), HXxy(0,3), HXxy(0,4), HXxy(0,5), HXxy(0,6), HXxy(0,7), HXxy(0,8), HXxy(0,9), HXxy(0,10), HXxy(0,11),
           HYxy(1,0), HYxy(1,1), HYxy(1,2), HYxy(1,3), HYxy(1,4), HYxy(1,5), HYxy(1,6), HYxy(1,7), HYxy(1,8), HYxy(1,9), HYxy(1,10), HYxy(1,11),
           HZxy(0,0), HZxy(0,1), HZxy(0,2), HZxy(0,3), HZxy(0,4), HZxy(0,5), HZxy(0,6), HZxy(0,7), HZxy(0,8), HZxy(0,9), HZxy(0,10), HZxy(0,11);
}

//Evaluates the strain-displacement matrix for membrane effect at a given Gauss point.
void
lin3DShell4::ComputeMembraneStrainDisplacementMatrix(const double ri, const double si, const Eigen::MatrixXd &xyloc, const Eigen::MatrixXd &BM12, Eigen::MatrixXd &Bij, Eigen::MatrixXd &Pij, Eigen::MatrixXd &Jij) {
    //Quadrature point coordinates.
    double p = ri; 
    double q = si;

    //Components Node Coordinates.
    double x1 = xyloc(0,0); double y1 = xyloc(1,0);
    double x2 = xyloc(0,1); double y2 = xyloc(1,1);
    double x3 = xyloc(0,2); double y3 = xyloc(1,2);
    double x4 = xyloc(0,3); double y4 = xyloc(1,3);

    //Strain-Displacement Deformation Matrix for Membrane.
    double D = x3*q*y4 + x2*q*y1 - x2*y4*p + x1*q*y3 - x1*y3*p - x1*q*y2 + x2*y3*p - x4*q*y3 + x4*y2*p - x4*y1*p - x3*q*y1 - x3*y2*p + x3*y1*p - x2*q*y4 + x4*q*y2 - x1*y4 + x1*y2 + x1*y4*p - x2*y1 + x2*y3 - x4*y3 + x4*y1 - x3*y2 + x3*y4;

    Bij(0, 0) = (-y4+y2-y3*p+y4*p-y2*q+y3*q)/D;
    Bij(0, 1) = 0.0;
    Bij(0, 2) = 1.0/8.0*(-2.0*y1*y2*p*q-y2*y2-y2*y3-2.0*y2*y2*p+y2*y2*q+2.0*y4*y4*q+y4*y4*q*q-2.0*y4*y4*p*q+y4*y4+y4*y3-y4*y4*p+y1*y1*q*q-3.0*y1*y4+3.0*y1*y2-y1*y1*p*q*q-y4*y4*p*q*q+y1*y2*p-3.0*y1*y3*p*p+3.0*y1*y4*p*p+y1*y3*p*p*q-y1*y4*p*p*q-y2*y2*p*p+y2*y2*p*p*q+y1*y1*p-y1*y1*q+y1*y1*p*p*q+2.0*y2*y3*p-y2*y3*q+2.0*y2*y1*p*p+3.0*y2*y3*p*p-3.0*y2*y4*p*p-y1*y1*p*p-y4*y2*p*q*q+y4*y3*p*q*q-2.0*y2*y3*p*q-2.0*y2*y1*p*p*q-y2*y3*p*p*q+y2*y4*p*p*q+2.0*y4*y3*p*q-2.0*y4*y3*q+3.0*y4*y2*q*q-3.0*y4*y3*q*q+2.0*y2*y2*p*q-2.0*y1*y4*q*q+y4*y2*p+y4*y3*p-y4*y2*q+2.0*y1*y4*p*q+y1*y2*p*q*q-y1*y3*p*q*q+2.0*y1*y4*p*q*q+3.0*y1*y3*q*q-3.0*y1*y3*p+3.0*y1*y3*q-y1*y4*q-3.0*y1*y2*q*q)/D - BM12(0,0);
    Bij(0, 3) = -1.0*(y1-y3-y3*p+y4*p-y1*q+y4*q)/D;
    Bij(0, 4) = 0.0;
    Bij(0, 5) = -1.0/8.0*(-y1*y1-y1*y4+2.0*y1*y1*p+y1*y1*q-y1*y1*p*p+y3*y3+y3*y3*p+2.0*y3*y3*q+y3*y3*q*q+3.0*y3*y1*q*q+3.0*y2*y4*q*q-3.0*y3*y4*q*q-2.0*y3*y2*q*q-2.0*y3*y4*q-y3*y4*p+y2*y2*p*q*q+y3*y3*p*q*q+2.0*y3*y3*p*q+3.0*y1*y4*p*p-3.0*y1*y3*p*p-y1*y4*q-y1*y3*q-2.0*y1*y4*p-y1*y3*p-3.0*y2*y4*p*p+3.0*y2*y3*p*p+2.0*y2*y1*p*p+3.0*y2*y4*q-y2*y3*q+3.0*y2*y4*p-y2*y1*p+y1*y1*p*p*q-2.0*y1*y1*p*q+y2*y2*p*p*q-y2*y2*q-y2*y2*p*p+3.0*y2*y1-3.0*y2*y3-y2*y2*p+y2*y2*q*q-2.0*y3*y2*p*q*q+y3*y1*p*q*q-2.0*y3*y4*p*q-y1*y4*p*p*q+y1*y3*p*p*q+2.0*y1*y4*p*q+y3*y4+y2*y4*p*p*q-y2*y3*p*p*q-2.0*y2*y1*p*p*q-2.0*y2*y3*p*q+2.0*y2*y1*p*q+y2*y4*p*q*q-y2*y1*p*q*q-3.0*y2*y1*q*q-y3*y4*p*q*q)/D - BM12(0,1);
    Bij(0, 6) = (y4-y2+y1*p-y2*p-y1*q+y4*q)/D;
    Bij(0, 7) = 0.0;
    Bij(0, 8) = 1.0/8.0*(-y1*y4-y3*y3*p+y3*y3*q+y3*y3*q*q+3.0*y3*y1*q*q+3.0*y2*y4*q*q-3.0*y3*y4*q*q-2.0*y3*y2*q*q-y3*y4*p+y2*y2*p*q*q+y3*y3*p*q*q+3.0*y1*y4*p*p-3.0*y1*y3*p*p+y1*y4*q-3.0*y1*y3*q-2.0*y1*y4*p+3.0*y1*y3*p-3.0*y2*y4*p*p+3.0*y2*y3*p*p+y2*y4*q+y2*y3*q+2.0*y2*y1*q-y2*y4*p-y2*y1*p-2.0*y2*y2*p*q-2.0*y2*y2*q+y2*y1+y2*y2-3.0*y2*y3+y2*y2*p+y2*y2*q*q-2.0*y3*y2*p*q*q+y3*y1*p*q*q-2.0*y3*y4*p*q+y1*y4*p*p*q-y1*y3*p*p*q-2.0*y1*y4*p*q+3.0*y3*y4-y2*y4*p*p*q+y2*y3*p*p*q+2.0*y2*y3*p*q+2.0*y2*y1*p*q+y2*y4*p*q*q-y2*y1*p*q*q-3.0*y2*y1*q*q-y3*y4*p*q*q+2.0*y3*y4*p*p*q+2.0*y4*y4*p*q-y4*y4*p*p*q+2.0*y3*y4*p*p-y4*y4*p*p-y3*y3*p*p-y4*y4+2.0*y4*y4*p-y4*y4*q-y3*y3*p*p*q)/D - BM12(0,2);
    Bij(0, 9) = -1.0*(-y1+y3+y1*p-y2*p-y2*q+y3*q)/D;
    Bij(0,10) = 0.0;
    Bij(0,11) = -1.0/8.0*(y3*y2*q+y3*y1*q+y3*y4*p+2.0*y3*y2*p+y3*y1*p-y4*y4*p*p*q-y3*y3*p*p*q-2.0*y3*y3*p*q-2.0*y1*y4*q*q+3.0*y1*y3*q*q-3.0*y1*y2*q*q+2.0*y1*y2*q+y1*y2*p-y4*y4*p*q*q-y1*y1*p*q*q+2.0*y1*y1*p*q-3.0*y4*y2*p*p+3.0*y4*y1*p*p-3.0*y4*y2*q+y4*y1*q-3.0*y4*y2*p+2.0*y3*y4*p*p+3.0*y3*y2*p*p-3.0*y3*y1*p*p-3.0*y4*y3*q*q+3.0*y4*y2*q*q-y3*y3-y3*y2+3.0*y3*y4-y3*y3*p*p-y3*y3*q-2.0*y3*y3*p+y3*y2*p*p*q-y3*y1*p*p*q+2.0*y3*y4*p*q+2.0*y3*y2*p*q-3.0*y4*y1-y1*y3*p*q*q+y1*y2*p*q*q-2.0*y1*y2*p*q-y4*y2*p*p*q+y4*y1*p*p*q-2.0*y4*y1*p*q+2.0*y3*y4*p*p*q+2.0*y1*y4*p*q*q+y4*y3*p*q*q-y4*y2*p*q*q-y4*y4*p*p+y4*y4*q+y4*y4*p-2.0*y1*y1*q-y1*y1*p+y1*y2+y1*y1+y1*y1*q*q+y4*y4*q*q)/D - BM12(0,3);

    Bij(1, 0) = 0.0;
    Bij(1, 1) = (-x2+x4+x2*q-x3*q+x3*p-x4*p)/D;
    Bij(1, 2) = 1.0/8.0*(-x1*x1*p*q*q-x4*x4*p*q*q+x1*x2*p-3.0*x1*x2*q*q+3.0*x1*x3*q*q-2.0*x1*x4*q*q-x4*x2*p*q*q+x4*x3*p*q*q-2.0*x1*x2*p*p*q-2.0*x1*x2*p*q+x1*x2*p*q*q-x1*x3*p*q*q+2.0*x1*x4*p*q*q+2.0*x4*x3*p*q+3.0*x1*x2-x1*x1*q+x1*x1*p+x1*x1*q*q+x4*x4*q*q-x1*x1*p*p-x2*x2+x2*x2*q-2.0*x2*x2*p-x2*x2*p*p+x4*x3*p+3.0*x4*x2*q*q-3.0*x4*x3*q*q+x1*x1*p*p*q+2.0*x2*x2*p*q+x2*x2*p*p*q+2.0*x1*x2*p*p+x4*x4+2.0*x4*x4*q-x4*x4*p-2.0*x4*x4*p*q-2.0*x4*x3*q+x4*x2*p+2.0*x4*x1*p*q-3.0*x4*x1+x4*x3-x4*x1*q-x4*x2*q+3.0*x2*x3*p*p-3.0*x2*x4*p*p+x2*x4*p*p*q+2.0*x2*x3*p+x1*x3*p*p*q-x1*x4*p*p*q-2.0*x2*x3*p*q-x2*x3*p*p*q-3.0*x1*x3*p*p+3.0*x1*x4*p*p-x2*x3*q-x2*x3-3.0*x1*x3*p+3.0*x1*x3*q)/D - BM12(1,0);
    Bij(1, 3) = 0.0;
    Bij(1, 4) = -1.0*(-x1+x3+x1*q-x4*q+x3*p-x4*p)/D;
    Bij(1, 5) = -1.0/8.0*(x1*x3*p*p*q+2.0*x1*x4*p*q-2.0*x1*x2*p*p*q+2.0*x1*x2*p*q+x2*x4*p*p*q+3.0*x2*x3*p*p-x2*x3*p*p*q-2.0*x2*x3*p*q-3.0*x2*x4*p*p-x2*x3*q+3.0*x2*x4*q-x1*x4*p*p*q-x2*x2*p*p-x2*x2*p-x2*x2*q+x2*x2*p*p*q-x1*x1*p*p+2.0*x1*x1*p-2.0*x1*x1*p*q-x1*x4+x1*x1*q-x1*x1+3.0*x1*x2-3.0*x3*x4*q*q-x3*x4*p*q*q+x3*x1*p*q*q-2.0*x3*x4*p*q-x3*x4*p+3.0*x3*x1*q*q-2.0*x3*x4*q+x2*x4*p*q*q+3.0*x2*x4*q*q-x2*x1*p*q*q-2.0*x2*x3*p*q*q-3.0*x2*x1*q*q-2.0*x2*x3*q*q+3.0*x2*x4*p+x3*x3*p*q*q+x3*x3*q*q+x3*x3*p+2.0*x3*x3*q+x3*x3+x2*x2*q*q+x2*x2*p*q*q+2.0*x3*x3*p*q-2.0*x1*x4*p-x1*x3*p-x1*x2*p-x1*x4*q+x4*x3-x1*x3*q+x1*x1*p*p*q+2.0*x1*x2*p*p-3.0*x2*x3+3.0*x1*x4*p*p-3.0*x1*x3*p*p)/D - BM12(1,1);
    Bij(1, 6) = 0.0;
    Bij(1, 7) = -1.0*(-x2+x4-x1*q+x4*q+x1*p-x2*p)/D;
    Bij(1, 8) = -1.0/8.0*(x1*x3*p*p*q+2.0*x1*x4*p*q-2.0*x1*x2*p*q+x2*x4*p*p*q-3.0*x2*x3*p*p-x2*x3*p*p*q-2.0*x2*x3*p*q+3.0*x2*x4*p*p-x2*x3*q-x2*x4*q-x1*x4*p*p*q-x2*x2*p+2.0*x2*x2*p*q+2.0*x2*x2*q+x1*x4-x1*x2+3.0*x3*x4*q*q+x3*x4*p*q*q-x3*x1*p*q*q+2.0*x3*x4*p*q+x3*x4*p-3.0*x3*x1*q*q-x2*x4*p*q*q-3.0*x2*x4*q*q+x2*x1*p*q*q+2.0*x2*x3*p*q*q+3.0*x2*x1*q*q+2.0*x2*x3*q*q+x2*x4*p-x3*x3*p*q*q-x3*x3*q*q+x3*x3*p-x3*x3*q-x2*x2*q*q-x2*x2*p*q*q-x2*x2+x4*x4+2.0*x1*x4*p-3.0*x1*x3*p+x1*x2*p-x1*x4*q-3.0*x4*x3+3.0*x1*x3*q-2.0*x1*x2*q-2.0*x4*x3*p*p+x3*x3*p*p+x4*x4*p*p+x4*x4*p*p*q+x4*x4*q-2.0*x4*x4*p-2.0*x4*x3*p*p*q-2.0*x4*x4*p*q+x3*x3*p*p*q+3.0*x2*x3-3.0*x1*x4*p*p+3.0*x1*x3*p*p)/D - BM12(1,2);
    Bij(1, 9) = 0.0;
    Bij(1,10) = (-x1+x3-x2*q+x3*q+x1*p-x2*p)/D;
    Bij(1,11) = 1.0/8.0*(2.0*x3*x3*p-x4*x4*p+x4*x4*p*p+x3*x2+x3*x3+x3*x3*q+2.0*x1*x2*p*q-3.0*x1*x3*q*q+3.0*x4*x1-3.0*x4*x3-x4*x4*q+x1*x1*p*q*q-2.0*x1*x1*p*q-3.0*x3*x2*p*p+3.0*x3*x1*p*p-2.0*x3*x2*p-x3*x1*p-x3*x2*q-x3*x1*q-2.0*x4*x3*p*p+3.0*x4*x2*p*p-3.0*x4*x1*p*p-x4*x3*p+3.0*x4*x2*p+3.0*x4*x2*q-x4*x1*q+x3*x3*p*p*q+2.0*x3*x3*p*q+x4*x4*p*p*q+2.0*x1*x4*q*q-x1*x2*p-2.0*x1*x2*q+x4*x4*p*q*q+3.0*x4*x3*q*q-3.0*x4*x2*q*q-x3*x2*p*p*q+x3*x1*p*p*q-2.0*x3*x2*p*q-x1*x1-x1*x2-2.0*x4*x3*p*p*q+x4*x2*p*p*q-2.0*x1*x4*p*q*q+x1*x3*p*q*q-x1*x2*p*q*q+3.0*x1*x2*q*q-x4*x3*p*q*q+x4*x2*p*q*q+x1*x1*p+2.0*x1*x1*q+2.0*x4*x1*p*q+x3*x3*p*p-x4*x4*q*q-x1*x1*q*q-2.0*x4*x3*p*q-x4*x1*p*p*q)/D - BM12(1,3);

    Bij(2, 0) = (-x2+x4+x2*q-x3*q+x3*p-x4*p)/D;
    Bij(2, 1) = (-y4+y2-y3*p+y4*p-y2*q+y3*q)/D;
    Bij(2, 2) = -1.0/8.0*(-2.0*x3*q*y4+x2*y4*p+3.0*x1*q*y3-3.0*x1*y3*p+2.0*x2*y3*p-2.0*x4*q*y3+x4*y2*p+3.0*x3*q*y1+2.0*x3*y2*p-3.0*x3*y1*p-x2*q*y4-x4*q*y2-3.0*x1*y4+3.0*x1*y2+x2*y4*p*p*q+3.0*x2*y1-x2*y3-x1*y4*p*p*q-2.0*x2*y3*p*q-x2*y3*p*p*q+x4*y3-3.0*x4*y1-x3*y2+x1*y3*p*p*q+x3*y4+x4*y3*p*q*q+2.0*x4*y3*p*q-x4*y2*p*q*q+2.0*y1*x1*p*p*q-2.0*y1*x2*p*p*q+y1*x3*p*p*q-y1*x4*p*p*q+x1*y2*p*q*q-x1*y3*p*q*q+2.0*y2*x2*p*p*q-y2*x3*p*p*q+y2*x4*p*p*q+2.0*y4*x3*p*q-4.0*y4*x4*p*q+2.0*y4*x1*p*q*q-y4*x2*p*q*q+y4*x3*p*q*q-2.0*y4*x4*p*q*q-2.0*y2*x1*p*q+4.0*y2*x2*p*q-2.0*y2*x3*p*q-2.0*y2*x1*p*p*q-2.0*y1*x2*p*q+2.0*y1*x4*p*q-2.0*y1*x1*p*q*q+y1*x2*p*q*q-y1*x3*p*q*q+2.0*y1*x4*p*q*q+3.0*x2*y3*p*p-3.0*x2*y4*p*p-x2*y3*q-3.0*x1*y3*p*p+3.0*x1*y4*p*p+x4*y3*p+3.0*x4*y2*q*q-3.0*x4*y3*q*q-3.0*x1*y2*q*q+3.0*x1*y3*q*q+2.0*y2*x1*p*p-2.0*y2*x2*p*p+3.0*y2*x3*p*p-3.0*y2*x4*p*p-2.0*y1*x1*p*p+2.0*y1*x2*p*p-3.0*y1*x3*p*p+3.0*y1*x4*p*p+y4*x3*p-2.0*y4*x4*p-2.0*y4*x1*q*q+3.0*y4*x2*q*q-3.0*y4*x3*q*q+2.0*y4*x4*q*q+2.0*y2*x2*q-y2*x3*q+y2*x1*p-4.0*y2*x2*p+2.0*y4*x1*p*q-y4*x1*q+4.0*y4*x4*q+y1*x2*p+2.0*y1*x1*q*q-3.0*y1*x2*q*q+3.0*y1*x3*q*q-2.0*y1*x4*q*q-2.0*y2*x2+2.0*y4*x4-2.0*y1*x1*q-y1*x4*q+2.0*y1*x1*p)/D - BM12(2,0);
    Bij(2, 3) = -1.0*(-x1+x3+x1*q-x4*q+x3*p-x4*p)/D;
    Bij(2, 4) = -1.0*(y1-y3-y3*p+y4*p-y1*q+y4*q)/D;
    Bij(2, 5) = 1.0/8.0*(-2.0*x3*q*y4+3.0*x2*y4*p-x1*q*y3-x1*y3*p-2.0*x4*q*y3+3.0*x4*y2*p-2.0*x4*y1*p-x3*q*y1-x3*y1*p+3.0*x2*q*y4+3.0*x4*q*y2-x1*y4+3.0*x1*y2-2.0*x1*y4*p+x2*y4*p*p*q+3.0*x2*y1-3.0*x2*y3-x1*y4*p*p*q-2.0*x2*y3*p*q-x2*y3*p*p*q+x4*y3-x4*y1-3.0*x3*y2+x1*y3*p*p*q+x3*y4-x4*y3*p*q*q-2.0*x4*y3*p*q+x4*y2*p*q*q+2.0*y1*x1*p*p*q-2.0*y1*x2*p*p*q+y1*x3*p*p*q-y1*x4*p*p*q-x1*y2*p*q*q+x1*y3*p*q*q+2.0*y2*x2*p*p*q-y2*x3*p*p*q+y2*x4*p*p*q-2.0*y4*x3*p*q+y4*x2*p*q*q-y4*x3*p*q*q+2.0*y2*x1*p*q-2.0*y2*x3*p*q-2.0*y2*x1*p*p*q-4.0*y1*x1*p*q+2.0*y1*x2*p*q+2.0*y1*x4*p*q-y1*x2*p*q*q+y1*x3*p*q*q+3.0*x2*y3*p*p-3.0*x2*y4*p*p-x2*y3*q-3.0*x1*y3*p*p+3.0*x1*y4*p*p-x4*y3*p+3.0*x4*y2*q*q-3.0*x4*y3*q*q-3.0*x1*y2*q*q+3.0*x1*y3*q*q+2.0*y2*x1*p*p-2.0*y2*x2*p*p+3.0*y2*x3*p*p-3.0*y2*x4*p*p-2.0*y1*x1*p*p+2.0*y1*x2*p*p-3.0*y1*x3*p*p+3.0*y1*x4*p*p-y4*x3*p+3.0*y4*x2*q*q-3.0*y4*x3*q*q-2.0*y2*x2*q-y2*x3*q-y2*x1*p-2.0*y2*x2*p+2.0*y4*x1*p*q-y4*x1*q-y1*x2*p-3.0*y1*x2*q*q+3.0*y1*x3*q*q+2.0*y1*x1*q-y1*x4*q+4.0*y1*x1*p-2.0*y1*x1+2.0*y2*x2*p*q*q-2.0*y2*x3*p*q*q+4.0*y3*x3*p*q-2.0*y3*x2*p*q*q+2.0*y3*x3*p*q*q+2.0*y3*x3+4.0*y3*x3*q+2.0*y3*x3*p-2.0*y3*x2*q*q+2.0*y3*x3*q*q+2.0*y2*x2*q*q-2.0*y2*x3*q*q)/D - BM12(2,1);
    Bij(2, 6) = -1.0*(-x2+x4-x1*q+x4*q+x1*p-x2*p)/D;
    Bij(2, 7) = (y4-y2+y1*p-y2*p-y1*q+y4*q)/D;
    Bij(2, 8) = -1.0/8.0*(2.0*x2*q*y1-x2*y4*p-3.0*x1*q*y3+3.0*x1*y3*p+2.0*x1*q*y2-x4*y2*p-2.0*x4*y1*p-3.0*x3*q*y1+3.0*x3*y1*p+x2*q*y4+x4*q*y2-x1*y4+x1*y2-2.0*x1*y4*p-x2*y4*p*p*q+x2*y1-3.0*x2*y3+x1*y4*p*p*q+2.0*x2*y3*p*q+x2*y3*p*p*q+3.0*x4*y3-x4*y1-3.0*x3*y2-x1*y3*p*p*q+3.0*x3*y4-x4*y3*p*q*q-2.0*x4*y3*p*q+x4*y2*p*q*q-2.0*y3*x3*p*p*q+2.0*y3*x4*p*p*q-y1*x3*p*p*q+y1*x4*p*p*q+2.0*y4*x3*p*p*q-2.0*y4*x4*p*p*q-x1*y2*p*q*q+x1*y3*p*q*q+y2*x3*p*p*q-y2*x4*p*p*q-2.0*y4*x3*p*q+4.0*y4*x4*p*q+y4*x2*p*q*q-y4*x3*p*q*q+2.0*y2*x1*p*q-4.0*y2*x2*p*q+2.0*y2*x3*p*q+2.0*y1*x2*p*q-2.0*y1*x4*p*q-y1*x2*p*q*q+y1*x3*p*q*q+2.0*y4*x3*p*p-2.0*y4*x4*p*p-2.0*y3*x3*p*p+2.0*y3*x4*p*p+3.0*x2*y3*p*p-3.0*x2*y4*p*p+x2*y3*q-3.0*x1*y3*p*p+3.0*x1*y4*p*p-x4*y3*p+3.0*x4*y2*q*q-3.0*x4*y3*q*q-3.0*x1*y2*q*q+3.0*x1*y3*q*q+3.0*y2*x3*p*p-3.0*y2*x4*p*p-3.0*y1*x3*p*p+3.0*y1*x4*p*p-y4*x3*p+4.0*y4*x4*p+3.0*y4*x2*q*q-3.0*y4*x3*q*q-4.0*y2*x2*q+y2*x3*q-y2*x1*p+2.0*y2*x2*p-2.0*y4*x1*p*q+y4*x1*q-2.0*y4*x4*q-y1*x2*p-3.0*y1*x2*q*q+3.0*y1*x3*q*q+2.0*y2*x2-2.0*y4*x4+y1*x4*q+2.0*y2*x2*p*q*q-2.0*y2*x3*p*q*q-2.0*y3*x2*p*q*q+2.0*y3*x3*p*q*q+2.0*y3*x3*q-2.0*y3*x3*p-2.0*y3*x2*q*q+2.0*y3*x3*q*q+2.0*y2*x2*q*q-2.0*y2*x3*q*q)/D - BM12(2,2);
    Bij(2, 9) = (-x1+x3-x2*q+x3*q+x1*p-x2*p)/D;
    Bij(2,10) = -1.0*(-y1+y3+y1*p-y2*p-y2*q+y3*q)/D;
    Bij(2,11) = 1.0/8.0*(2.0*x2*q*y1-3.0*x2*y4*p+x1*q*y3+x1*y3*p+2.0*x1*q*y2+2.0*x2*y3*p-3.0*x4*y2*p+x3*q*y1+2.0*x3*y2*p+x3*y1*p-3.0*x2*q*y4-3.0*x4*q*y2-3.0*x1*y4+x1*y2-x2*y4*p*p*q+x2*y1-x2*y3+x1*y4*p*p*q+2.0*x2*y3*p*q+x2*y3*p*p*q+3.0*x4*y3-3.0*x4*y1-x3*y2-x1*y3*p*p*q+3.0*x3*y4+x4*y3*p*q*q+2.0*x4*y3*p*q-x4*y2*p*q*q-2.0*y3*x3*p*p*q+2.0*y3*x4*p*p*q-y1*x3*p*p*q+y1*x4*p*p*q+2.0*y4*x3*p*p*q-2.0*y4*x4*p*p*q+x1*y2*p*q*q-x1*y3*p*q*q+y2*x3*p*p*q-y2*x4*p*p*q+2.0*y4*x3*p*q+2.0*y4*x1*p*q*q-y4*x2*p*q*q+y4*x3*p*q*q-2.0*y4*x4*p*q*q-2.0*y2*x1*p*q+2.0*y2*x3*p*q+4.0*y1*x1*p*q-2.0*y1*x2*p*q-2.0*y1*x4*p*q-2.0*y1*x1*p*q*q+y1*x2*p*q*q-y1*x3*p*q*q+2.0*y1*x4*p*q*q+2.0*y4*x3*p*p-2.0*y4*x4*p*p-2.0*y3*x3*p*p+2.0*y3*x4*p*p+3.0*x2*y3*p*p-3.0*x2*y4*p*p+x2*y3*q-3.0*x1*y3*p*p+3.0*x1*y4*p*p+x4*y3*p+3.0*x4*y2*q*q-3.0*x4*y3*q*q-3.0*x1*y2*q*q+3.0*x1*y3*q*q+3.0*y2*x3*p*p-3.0*y2*x4*p*p-3.0*y1*x3*p*p+3.0*y1*x4*p*p+y4*x3*p+2.0*y4*x4*p-2.0*y4*x1*q*q+3.0*y4*x2*q*q-3.0*y4*x3*q*q+2.0*y4*x4*q*q+y2*x3*q+y2*x1*p-2.0*y4*x1*p*q+y4*x1*q+2.0*y4*x4*q+y1*x2*p+2.0*y1*x1*q*q-3.0*y1*x2*q*q+3.0*y1*x3*q*q-2.0*y1*x4*q*q-4.0*y1*x1*q+y1*x4*q-2.0*y1*x1*p+2.0*y1*x1-4.0*y3*x3*p*q-2.0*y3*x3-2.0*y3*x3*q-4.0*y3*x3*p)/D - BM12(2,3);

    //The Contant-Strain Penalty Matrix for Membrane.
    double A = (-x1*y4*p+x2*y4*p-x4*q*y2+x3*q*y1-x2*q*y1-x2*y3*p-x1*q*y3+x1*q*y2+x1*y3*p+x4*q*y3-x4*y2*p+x4*y1*p-x3*q*y4+x3*y2*p-x3*y1*p+x2*q*y4+x1*y4-x1*y2-x4*y1-x2*y3-x3*y4+x3*y2+x2*y1+x4*y3);
    Pij(0, 0) =  1.0/2.0*(-x3*q+x2*q-x4*p+x3*p+x4-x2)/A;
    Pij(0, 1) =  1.0/2.0*(-y4*p-y3*q+y2*q+y3*p+y4-y2)/A;
    Pij(0, 2) = -1.0/16.0*(-6.0*x1*y4*p+x2*y4*p-x4*q*y2+7*x3*q*y1-6*x2*q*y1-2.0*x2*y3*p-7.0*x1*q*y3+6.0*x1*q*y2+7.0*x1*y3*p+3*x2*y4*p*p*q+3.0*y4*x2*p*q*q+2.0*x4*q*y3+3.0*x1*y3*p*p*q-x4*y2*p+6.0*x4*y1*p-2.0*x3*q*y4-3.0*x1*y4*p*p*q+2.0*x3*y2*p-7.0*x3*y1*p+x2*q*y4-3.0*x2*y3*p*p*q+5.0*x1*y4-5.0*x1*y2-4.0*y4*x2*p*q+2.0*y4*x3*p*q-5.0*x4*y1-3.0*x2*y3-3.0*x3*y4+2.0*y1*x2*p*q-2.0*y1*x4*p*q+3.0*x3*y2+3.0*y1*x3*p*q*q-3.0*y1*x2*p*q*q+2.0*x2*y3*p*q+2.0*y4*x1*p*q-2*x2*y4+5.0*x2*y1+3.0*x4*y3-3.0*x4*y2*p*q*q+2.0*x4*y2+3.0*x4*y3*p*q*q-2.0*x4*y3*p*q-x1*y2*q*q+x1*y3*q*q-5.0*x4*y3*p+x4*y2*q*q-x4*y3*q*q+5.0*x2*y3*q+x2*y3*p*p-x2*y4*p*p-x1*y3*p*p+x1*y4*p*p-y1*x4*q+y1*x2*p+y1*x2*q*q-y1*x3*q*q+y4*x1*q+5.0*y4*x3*p-y4*x2*q*q+y4*x3*q*q-5.0*y2*x3*q-y2*x1*p-y2*x3*p*p+y2*x4*p*p+y1*x3*p*p-y1*x4*p*p+3.0*x1*y2*p*q*q-3.0*x1*y3*p*q*q-3.0*y1*x3*p*p*q+3.0*y1*x4*p*p*q-2.0*y2*x1*p*q-2.0*y2*x3*p*q+4.0*y2*x4*p*q+3.0*y2*x3*p*p*q-3.0*y2*x4*p*p*q-3.0*y4*x3*p*q*q)/A;
    Pij(0, 3) = -1.0/2.0*(-x4*q+x1*q-x4*p+x3*p-x1+x3)/A;
    Pij(0, 4) = -1.0/2.0*(-y4*p+y1*q+y3*p-y4*q+y3-y1)/A;
    Pij(0, 5) =  1.0/16.0*(2.0*x1*y4*p-7.0*x2*y4*p+7.0*x4*q*y2-3.0*x2*y3*p*p*q+3.0*x2*y4*p*p*q-x3*q*y1+6.0*x2*q*y1+6.0*x2*y3*p+x1*q*y3-6.0*x1*q*y2-x1*y3*p+3.0*x3*y1*p*q*q-3.0*x3*y4*p*q*q-2.0*x4*q*y3+7.0*x4*y2*p-2.0*x4*y1*p+2.0*x3*q*y4-6.0*x3*y2*p+x3*y1*p-7.0*x2*q*y4+3.0*x2*y4*p*q*q-3.0*x2*y1*p*q*q+2.0*x3*y4*p*q-3.0*x1*y4+5.0*x1*y2-3.0*y2*x4*p*q*q+3.0*x1*y3*p*p*q+3.0*x4*y1+5.0*x2*y3+3.0*x3*y4-5.0*x3*y2+3.0*y2*x1*p*q*q-5.0*x2*y1-3.0*x1*y4*p*p*q-3.0*x4*y3+3.0*y3*x4*p*q*q-2.0*x1*y4*p*q-3.0*y3*x1*p*q*q-2.0*y3*x2*p*q-2.0*y3*x4*p*q+2.0*y1*x4*p*q-3.0*y1*x3*p*p*q+3.0*y1*x4*p*p*q-3.0*y2*x4*p*p*q+2*y1*x2*p*q-4.0*y1*x3*p*q-2.0*y2*x1*p*q+2.0*y2*x3*p*q+3.0*y2*x3*p*p*q-2.0*x1*y3+2.0*x3*y1+y2*x1*q*q-y2*x4*q*q+y3*x4*q*q+4.0*y3*x1*p*q+y1*x3*p*p-y1*x4*p*p+y3*x2*q-5.0*y3*x4*p-y3*x1*q*q+y1*x2*p-y2*x3*p*p+y2*x4*p*p-5.0*y1*x4*q-y2*x3*q-y2*x1*p-x2*y1*q*q+x1*y4*p*p+5.0*x3*y4*p+x3*y1*q*q-x3*y4*q*q+5.0*x1*y4*q-x1*y3*p*p+x2*y3*p*p-x2*y4*p*p+x2*y4*q*q)/A;
    Pij(0, 6) = -1.0/2.0*(-x2+x4-x1*q+x4*q+x1*p-x2*p)/A;
    Pij(0, 7) = -1.0/2.0*(-y1*q+y4*q-y2*p+y1*p+y4-y2)/A;
    Pij(0, 8) = -1.0/16.0*(-2.0*x1*y4*p+x2*y4*p-x4*q*y2-2.0*y3*x2*p*q+7.0*x3*q*y1-2.0*x4*y2-2.0*x2*q*y1-6.0*x2*y3*p+2.0*x2*y4-7.0*x1*q*y3+2.0*x1*q*y2+7.0*x1*y3*p-3.0*y3*x1*p*q*q+2.0*y3*x4*p*q+6.0*x4*q*y3-x4*y2*p+2.0*x4*y1*p-6.0*x3*q*y4+6.0*x3*y2*p-7.0*x3*y1*p+x2*q*y4+3.0*y2*x1*p*q*q-4.0*y2*x4*p*q+2.0*y2*x3*p*q+3.0*y3*x4*p*q*q-3.0*y3*x2*p*p*q+3.0*y3*x1*p*p*q+3*y4*x2*p*p*q-2.0*y4*x3*p*q+4.0*y4*x2*p*q-3.0*y4*x1*p*p*q-3.0*y2*x4*p*q*q-y2*x3*q-5.0*y2*x1*p+y3*x4*q*q+2.0*y2*x1*p*q-y3*x1*q*q-2.0*x2*y1*p*q+y3*x2*q-3.0*x3*y4*p*q*q+3.0*x1*y4+3.0*x3*y1*p*q*q-3.0*x1*y2+5.0*y4*x1*q-y3*x4*p+y2*x1*q*q-y2*x4*q*q-5.0*x4*y1*q+x2*y4*q*q-x2*y1*q*q+5.0*x2*y1*p-x3*y4*q*q+x3*y1*q*q+y3*x1*p*p-y3*x2*p*p+y4*x2*p*p-y4*x1*p*p+y4*x3*p-2.0*y4*x1*p*q-3.0*x4*y1-5.0*x2*y3-5.0*x3*y4+5.0*x3*y2-3.0*x2*y1*p*q*q+3.0*x3*y2*p*p*q-3.0*x3*y1*p*p*q-x3*y1*p*p+x3*y2*p*p+3.0*x4*y1*p*p*q-3.0*x4*y2*p*p*q+3.0*x2*y1-x4*y2*p*p+x4*y1*p*p+2.0*x4*y1*p*q+5.0*x4*y3+3.0*x2*y4*p*q*q)/A;
    Pij(0, 9) =  1.0/2.0*(x1*p-x2*p+x3*q-x2*q-x1+x3)/A;
    Pij(0,10) =  1.0/2.0*(-y2*q+y3*q-y2*p+y1*p-y1+y3)/A;
    Pij(0,11) =  1.0/16.0*(6.0*x1*y4*p-7.0*x2*y4*p+7.0*x4*q*y2+2.0*y3*x2*p*q-x3*q*y1+2.0*x1*y3+2.0*x2*q*y1+2.0*x2*y3*p+x1*q*y3-2.0*x1*q*y2-x1*y3*p-4.0*y3*x1*p*q-3.0*y3*x1*p*q*q+2.0*y3*x4*p*q-6.0*x4*q*y3+7.0*x4*y2*p-6.0*x4*y1*p+6.0*x3*q*y4-2.0*x3*y1-2.0*x3*y2*p+x3*y1*p-7.0*x2*q*y4+3.0*y2*x1*p*q*q-2.0*y2*x3*p*q+3.0*y3*x4*p*q*q-3.0*y3*x2*p*p*q+3.0*y3*x1*p*p*q+3.0*y4*x2*p*p*q-2.0*y4*x3*p*q-3.0*y4*x1*p*p*q-3.0*y2*x4*p*q*q-5.0*y2*x3*q-5.0*y2*x1*p-y3*x4*q*q+2.0*y2*x1*p*q+y3*x1*q*q-2.0*x2*y1*p*q+5.0*y3*x2*q-3.0*x3*y4*p*q*q-5.0*x1*y4+3.0*x3*y1*p*q*q+4.0*x3*y1*p*q+3.0*x1*y2+y4*x1*q-y3*x4*p-y2*x1*q*q+y2*x4*q*q-x4*y1*q-x2*y4*q*q+x2*y1*q*q+5.0*x2*y1*p+x3*y4*q*q-x3*y1*q*q+y3*x1*p*p-y3*x2*p*p+y4*x2*p*p-y4*x1*p*p+y4*x3*p+2.0*y4*x1*p*q+5.0*x4*y1+3.0*x2*y3+5.0*x3*y4-3.0*x3*y2-3.0*x2*y1*p*q*q+3.0*x3*y2*p*p*q-3.0*x3*y1*p*p*q-x3*y1*p*p+x3*y2*p*p+3.0*x4*y1*p*p*q-3.0*x4*y2*p*p*q-3.0*x2*y1-x4*y2*p*p+x4*y1*p*p-2.0*x4*y1*p*q-5.0*x4*y3+3.0*x2*y4*p*q*q)/A;

    //Jacobian Matrix.
    Jij(0,0) = -1.0/4.0*x1 + 1.0/4.0*x1*q + 1.0/4.0*x2 - 1.0/4.0*x2*q + 1.0/4.0*x3 + 1.0/4.0*x3*q - 1.0/4.0*x4*q - 1.0/4.0*x4;
    Jij(0,1) = -1.0/4.0*y1 + 1.0/4.0*y1*q + 1.0/4.0*y2 - 1.0/4.0*q*y2 + 1.0/4.0*y3 + 1.0/4.0*q*y3 - 1.0/4.0*y4*q - 1.0/4.0*y4;
    Jij(1,0) = -1.0/4.0*x1 + 1.0/4.0*x1*p - 1.0/4.0*x2 - 1.0/4.0*x2*p + 1.0/4.0*x3 + 1.0/4.0*x3*p - 1.0/4.0*x4*p + 1.0/4.0*x4;
    Jij(1,1) = -1.0/4.0*y1 + 1.0/4.0*y1*p - 1.0/4.0*y2 - 1.0/4.0*y2*p + 1.0/4.0*y3 + 1.0/4.0*y3*p - 1.0/4.0*y4*p + 1.0/4.0*y4;
}

//Evaluates the Constant Tension matrix and Jacobian for membrane effect at a given Gauss point.
void
lin3DShell4::ConstantTensionMatrix(const double ri, const double si, const Eigen::MatrixXd &xyloc, Eigen::MatrixXd &BM12, Eigen::MatrixXd &Jij) {
    //Components Node Coordinates.
    double p = ri; 
    double q = si;

    double x1 = xyloc(0,0); double y1 = xyloc(1,0);
    double x2 = xyloc(0,1); double y2 = xyloc(1,1);
    double x3 = xyloc(0,2); double y3 = xyloc(1,2);
    double x4 = xyloc(0,3); double y4 = xyloc(1,3);

    //Strain-Displacement Deformation Matrix for Membrane.
    double D = x3*q*y4 + x2*q*y1 - x2*y4*p + x1*q*y3 - x1*y3*p - x1*q*y2 + x2*y3*p - x4*q*y3 + x4*y2*p - x4*y1*p - x3*q*y1 - x3*y2*p + x3*y1*p - x2*q*y4 + x4*q*y2 - x1*y4 + x1*y2 + x1*y4*p - x2*y1 + x2*y3 - x4*y3 + x4*y1 - x3*y2 + x3*y4;

    BM12(0,0) = 1.0/8.0*(-2.0*y1*y2*p*q-y2*y2-y2*y3-2.0*y2*y2*p+y2*y2*q+2.0*y4*y4*q+y4*y4*q*q-2.0*y4*y4*p*q+y4*y4+y4*y3-y4*y4*p+y1*y1*q*q-3.0*y1*y4+3.0*y1*y2-y1*y1*p*q*q-y4*y4*p*q*q+y1*y2*p-3.0*y1*y3*p*p+3.0*y1*y4*p*p+y1*y3*p*p*q-y1*y4*p*p*q-y2*y2*p*p+y2*y2*p*p*q+y1*y1*p-y1*y1*q+y1*y1*p*p*q+2.0*y2*y3*p-y2*y3*q+2.0*y2*y1*p*p+3.0*y2*y3*p*p-3.0*y2*y4*p*p-y1*y1*p*p-y4*y2*p*q*q+y4*y3*p*q*q-2.0*y2*y3*p*q-2.0*y2*y1*p*p*q-y2*y3*p*p*q+y2*y4*p*p*q+2.0*y4*y3*p*q-2.0*y4*y3*q+3.0*y4*y2*q*q-3.0*y4*y3*q*q+2.0*y2*y2*p*q-2.0*y1*y4*q*q+y4*y2*p+y4*y3*p-y4*y2*q+2.0*y1*y4*p*q+y1*y2*p*q*q-y1*y3*p*q*q+2.0*y1*y4*p*q*q+3.0*y1*y3*q*q-3.0*y1*y3*p+3.0*y1*y3*q-y1*y4*q-3.0*y1*y2*q*q)/D;
    BM12(0,1) = -1.0/8.0*(-y1*y1-y1*y4+2.0*y1*y1*p+y1*y1*q-y1*y1*p*p+y3*y3+y3*y3*p+2.0*y3*y3*q+y3*y3*q*q+3.0*y3*y1*q*q+3.0*y2*y4*q*q-3.0*y3*y4*q*q-2.0*y3*y2*q*q-2.0*y3*y4*q-y3*y4*p+y2*y2*p*q*q+y3*y3*p*q*q+2.0*y3*y3*p*q+3.0*y1*y4*p*p-3.0*y1*y3*p*p-y1*y4*q-y1*y3*q-2.0*y1*y4*p-y1*y3*p-3.0*y2*y4*p*p+3.0*y2*y3*p*p+2.0*y2*y1*p*p+3.0*y2*y4*q-y2*y3*q+3.0*y2*y4*p-y2*y1*p+y1*y1*p*p*q-2.0*y1*y1*p*q+y2*y2*p*p*q-y2*y2*q-y2*y2*p*p+3.0*y2*y1-3.0*y2*y3-y2*y2*p+y2*y2*q*q-2.0*y3*y2*p*q*q+y3*y1*p*q*q-2.0*y3*y4*p*q-y1*y4*p*p*q+y1*y3*p*p*q+2.0*y1*y4*p*q+y3*y4+y2*y4*p*p*q-y2*y3*p*p*q-2.0*y2*y1*p*p*q-2.0*y2*y3*p*q+2.0*y2*y1*p*q+y2*y4*p*q*q-y2*y1*p*q*q-3.0*y2*y1*q*q-y3*y4*p*q*q)/D;
    BM12(0,2) = 1.0/8.0*(-y1*y4-y3*y3*p+y3*y3*q+y3*y3*q*q+3.0*y3*y1*q*q+3.0*y2*y4*q*q-3.0*y3*y4*q*q-2.0*y3*y2*q*q-y3*y4*p+y2*y2*p*q*q+y3*y3*p*q*q+3.0*y1*y4*p*p-3.0*y1*y3*p*p+y1*y4*q-3.0*y1*y3*q-2.0*y1*y4*p+3.0*y1*y3*p-3.0*y2*y4*p*p+3.0*y2*y3*p*p+y2*y4*q+y2*y3*q+2.0*y2*y1*q-y2*y4*p-y2*y1*p-2.0*y2*y2*p*q-2.0*y2*y2*q+y2*y1+y2*y2-3.0*y2*y3+y2*y2*p+y2*y2*q*q-2.0*y3*y2*p*q*q+y3*y1*p*q*q-2.0*y3*y4*p*q+y1*y4*p*p*q-y1*y3*p*p*q-2.0*y1*y4*p*q+3.0*y3*y4-y2*y4*p*p*q+y2*y3*p*p*q+2.0*y2*y3*p*q+2.0*y2*y1*p*q+y2*y4*p*q*q-y2*y1*p*q*q-3.0*y2*y1*q*q-y3*y4*p*q*q+2.0*y3*y4*p*p*q+2.0*y4*y4*p*q-y4*y4*p*p*q+2.0*y3*y4*p*p-y4*y4*p*p-y3*y3*p*p-y4*y4+2.0*y4*y4*p-y4*y4*q-y3*y3*p*p*q)/D;
    BM12(0,3) = -1.0/8.0*(y3*y2*q+y3*y1*q+y3*y4*p+2.0*y3*y2*p+y3*y1*p-y4*y4*p*p*q-y3*y3*p*p*q-2.0*y3*y3*p*q-2.0*y1*y4*q*q+3.0*y1*y3*q*q-3.0*y1*y2*q*q+2.0*y1*y2*q+y1*y2*p-y4*y4*p*q*q-y1*y1*p*q*q+2.0*y1*y1*p*q-3.0*y4*y2*p*p+3.0*y4*y1*p*p-3.0*y4*y2*q+y4*y1*q-3.0*y4*y2*p+2.0*y3*y4*p*p+3.0*y3*y2*p*p-3.0*y3*y1*p*p-3.0*y4*y3*q*q+3.0*y4*y2*q*q-y3*y3-y3*y2+3.0*y3*y4-y3*y3*p*p-y3*y3*q-2.0*y3*y3*p+y3*y2*p*p*q-y3*y1*p*p*q+2.0*y3*y4*p*q+2.0*y3*y2*p*q-3.0*y4*y1-y1*y3*p*q*q+y1*y2*p*q*q-2.0*y1*y2*p*q-y4*y2*p*p*q+y4*y1*p*p*q-2.0*y4*y1*p*q+2.0*y3*y4*p*p*q+2.0*y1*y4*p*q*q+y4*y3*p*q*q-y4*y2*p*q*q-y4*y4*p*p+y4*y4*q+y4*y4*p-2.0*y1*y1*q-y1*y1*p+y1*y2+y1*y1+y1*y1*q*q+y4*y4*q*q)/D;

    BM12(1,0) = 1.0/8.0*(-x1*x1*p*q*q-x4*x4*p*q*q+x1*x2*p-3.0*x1*x2*q*q+3.0*x1*x3*q*q-2.0*x1*x4*q*q-x4*x2*p*q*q+x4*x3*p*q*q-2.0*x1*x2*p*p*q-2.0*x1*x2*p*q+x1*x2*p*q*q-x1*x3*p*q*q+2.0*x1*x4*p*q*q+2.0*x4*x3*p*q+3.0*x1*x2-x1*x1*q+x1*x1*p+x1*x1*q*q+x4*x4*q*q-x1*x1*p*p-x2*x2+x2*x2*q-2.0*x2*x2*p-x2*x2*p*p+x4*x3*p+3.0*x4*x2*q*q-3.0*x4*x3*q*q+x1*x1*p*p*q+2.0*x2*x2*p*q+x2*x2*p*p*q+2.0*x1*x2*p*p+x4*x4+2.0*x4*x4*q-x4*x4*p-2.0*x4*x4*p*q-2.0*x4*x3*q+x4*x2*p+2.0*x4*x1*p*q-3.0*x4*x1+x4*x3-x4*x1*q-x4*x2*q+3.0*x2*x3*p*p-3.0*x2*x4*p*p+x2*x4*p*p*q+2.0*x2*x3*p+x1*x3*p*p*q-x1*x4*p*p*q-2.0*x2*x3*p*q-x2*x3*p*p*q-3.0*x1*x3*p*p+3.0*x1*x4*p*p-x2*x3*q-x2*x3-3.0*x1*x3*p+3.0*x1*x3*q)/D;
    BM12(1,1) = -1.0/8.0*(x1*x3*p*p*q+2.0*x1*x4*p*q-2.0*x1*x2*p*p*q+2.0*x1*x2*p*q+x2*x4*p*p*q+3.0*x2*x3*p*p-x2*x3*p*p*q-2.0*x2*x3*p*q-3.0*x2*x4*p*p-x2*x3*q+3.0*x2*x4*q-x1*x4*p*p*q-x2*x2*p*p-x2*x2*p-x2*x2*q+x2*x2*p*p*q-x1*x1*p*p+2.0*x1*x1*p-2.0*x1*x1*p*q-x1*x4+x1*x1*q-x1*x1+3.0*x1*x2-3.0*x3*x4*q*q-x3*x4*p*q*q+x3*x1*p*q*q-2.0*x3*x4*p*q-x3*x4*p+3.0*x3*x1*q*q-2.0*x3*x4*q+x2*x4*p*q*q+3.0*x2*x4*q*q-x2*x1*p*q*q-2.0*x2*x3*p*q*q-3.0*x2*x1*q*q-2.0*x2*x3*q*q+3.0*x2*x4*p+x3*x3*p*q*q+x3*x3*q*q+x3*x3*p+2.0*x3*x3*q+x3*x3+x2*x2*q*q+x2*x2*p*q*q+2.0*x3*x3*p*q-2.0*x1*x4*p-x1*x3*p-x1*x2*p-x1*x4*q+x4*x3-x1*x3*q+x1*x1*p*p*q+2.0*x1*x2*p*p-3.0*x2*x3+3.0*x1*x4*p*p-3.0*x1*x3*p*p)/D;
    BM12(1,2) = -1.0/8.0*(x1*x3*p*p*q+2.0*x1*x4*p*q-2.0*x1*x2*p*q+x2*x4*p*p*q-3.0*x2*x3*p*p-x2*x3*p*p*q-2.0*x2*x3*p*q+3.0*x2*x4*p*p-x2*x3*q-x2*x4*q-x1*x4*p*p*q-x2*x2*p+2.0*x2*x2*p*q+2.0*x2*x2*q+x1*x4-x1*x2+3.0*x3*x4*q*q+x3*x4*p*q*q-x3*x1*p*q*q+2.0*x3*x4*p*q+x3*x4*p-3.0*x3*x1*q*q-x2*x4*p*q*q-3.0*x2*x4*q*q+x2*x1*p*q*q+2.0*x2*x3*p*q*q+3.0*x2*x1*q*q+2.0*x2*x3*q*q+x2*x4*p-x3*x3*p*q*q-x3*x3*q*q+x3*x3*p-x3*x3*q-x2*x2*q*q-x2*x2*p*q*q-x2*x2+x4*x4+2.0*x1*x4*p-3.0*x1*x3*p+x1*x2*p-x1*x4*q-3.0*x4*x3+3.0*x1*x3*q-2.0*x1*x2*q-2.0*x4*x3*p*p+x3*x3*p*p+x4*x4*p*p+x4*x4*p*p*q+x4*x4*q-2.0*x4*x4*p-2.0*x4*x3*p*p*q-2.0*x4*x4*p*q+x3*x3*p*p*q+3.0*x2*x3-3.0*x1*x4*p*p+3.0*x1*x3*p*p)/D;
    BM12(1,3) = 1.0/8.0*(2.0*x3*x3*p-x4*x4*p+x4*x4*p*p+x3*x2+x3*x3+x3*x3*q+2.0*x1*x2*p*q-3.0*x1*x3*q*q+3.0*x4*x1-3.0*x4*x3-x4*x4*q+x1*x1*p*q*q-2.0*x1*x1*p*q-3.0*x3*x2*p*p+3.0*x3*x1*p*p-2.0*x3*x2*p-x3*x1*p-x3*x2*q-x3*x1*q-2.0*x4*x3*p*p+3.0*x4*x2*p*p-3.0*x4*x1*p*p-x4*x3*p+3.0*x4*x2*p+3.0*x4*x2*q-x4*x1*q+x3*x3*p*p*q+2.0*x3*x3*p*q+x4*x4*p*p*q+2.0*x1*x4*q*q-x1*x2*p-2.0*x1*x2*q+x4*x4*p*q*q+3.0*x4*x3*q*q-3.0*x4*x2*q*q-x3*x2*p*p*q+x3*x1*p*p*q-2.0*x3*x2*p*q-x1*x1-x1*x2-2.0*x4*x3*p*p*q+x4*x2*p*p*q-2.0*x1*x4*p*q*q+x1*x3*p*q*q-x1*x2*p*q*q+3.0*x1*x2*q*q-x4*x3*p*q*q+x4*x2*p*q*q+x1*x1*p+2.0*x1*x1*q+2.0*x4*x1*p*q+x3*x3*p*p-x4*x4*q*q-x1*x1*q*q-2.0*x4*x3*p*q-x4*x1*p*p*q)/D;

    BM12(2,0) = -1.0/8.0*(-2.0*x3*q*y4+x2*y4*p+3.0*x1*q*y3-3.0*x1*y3*p+2.0*x2*y3*p-2.0*x4*q*y3+x4*y2*p+3.0*x3*q*y1+2.0*x3*y2*p-3.0*x3*y1*p-x2*q*y4-x4*q*y2-3.0*x1*y4+3.0*x1*y2+x2*y4*p*p*q+3.0*x2*y1-x2*y3-x1*y4*p*p*q-2.0*x2*y3*p*q-x2*y3*p*p*q+x4*y3-3.0*x4*y1-x3*y2+x1*y3*p*p*q+x3*y4+x4*y3*p*q*q+2.0*x4*y3*p*q-x4*y2*p*q*q+2.0*y1*x1*p*p*q-2.0*y1*x2*p*p*q+y1*x3*p*p*q-y1*x4*p*p*q+x1*y2*p*q*q-x1*y3*p*q*q+2.0*y2*x2*p*p*q-y2*x3*p*p*q+y2*x4*p*p*q+2.0*y4*x3*p*q-4.0*y4*x4*p*q+2.0*y4*x1*p*q*q-y4*x2*p*q*q+y4*x3*p*q*q-2.0*y4*x4*p*q*q-2.0*y2*x1*p*q+4.0*y2*x2*p*q-2.0*y2*x3*p*q-2.0*y2*x1*p*p*q-2.0*y1*x2*p*q+2.0*y1*x4*p*q-2.0*y1*x1*p*q*q+y1*x2*p*q*q-y1*x3*p*q*q+2.0*y1*x4*p*q*q+3.0*x2*y3*p*p-3.0*x2*y4*p*p-x2*y3*q-3.0*x1*y3*p*p+3.0*x1*y4*p*p+x4*y3*p+3.0*x4*y2*q*q-3.0*x4*y3*q*q-3.0*x1*y2*q*q+3.0*x1*y3*q*q+2.0*y2*x1*p*p-2.0*y2*x2*p*p+3.0*y2*x3*p*p-3.0*y2*x4*p*p-2.0*y1*x1*p*p+2.0*y1*x2*p*p-3.0*y1*x3*p*p+3.0*y1*x4*p*p+y4*x3*p-2.0*y4*x4*p-2.0*y4*x1*q*q+3.0*y4*x2*q*q-3.0*y4*x3*q*q+2.0*y4*x4*q*q+2.0*y2*x2*q-y2*x3*q+y2*x1*p-4.0*y2*x2*p+2.0*y4*x1*p*q-y4*x1*q+4.0*y4*x4*q+y1*x2*p+2.0*y1*x1*q*q-3.0*y1*x2*q*q+3.0*y1*x3*q*q-2.0*y1*x4*q*q-2.0*y2*x2+2.0*y4*x4-2.0*y1*x1*q-y1*x4*q+2.0*y1*x1*p)/D;
    BM12(2,1) = 1.0/8.0*(-2.0*x3*q*y4+3.0*x2*y4*p-x1*q*y3-x1*y3*p-2.0*x4*q*y3+3.0*x4*y2*p-2.0*x4*y1*p-x3*q*y1-x3*y1*p+3.0*x2*q*y4+3.0*x4*q*y2-x1*y4+3.0*x1*y2-2.0*x1*y4*p+x2*y4*p*p*q+3.0*x2*y1-3.0*x2*y3-x1*y4*p*p*q-2.0*x2*y3*p*q-x2*y3*p*p*q+x4*y3-x4*y1-3.0*x3*y2+x1*y3*p*p*q+x3*y4-x4*y3*p*q*q-2.0*x4*y3*p*q+x4*y2*p*q*q+2.0*y1*x1*p*p*q-2.0*y1*x2*p*p*q+y1*x3*p*p*q-y1*x4*p*p*q-x1*y2*p*q*q+x1*y3*p*q*q+2.0*y2*x2*p*p*q-y2*x3*p*p*q+y2*x4*p*p*q-2.0*y4*x3*p*q+y4*x2*p*q*q-y4*x3*p*q*q+2.0*y2*x1*p*q-2.0*y2*x3*p*q-2.0*y2*x1*p*p*q-4.0*y1*x1*p*q+2.0*y1*x2*p*q+2.0*y1*x4*p*q-y1*x2*p*q*q+y1*x3*p*q*q+3.0*x2*y3*p*p-3.0*x2*y4*p*p-x2*y3*q-3.0*x1*y3*p*p+3.0*x1*y4*p*p-x4*y3*p+3.0*x4*y2*q*q-3.0*x4*y3*q*q-3.0*x1*y2*q*q+3.0*x1*y3*q*q+2.0*y2*x1*p*p-2.0*y2*x2*p*p+3.0*y2*x3*p*p-3.0*y2*x4*p*p-2.0*y1*x1*p*p+2.0*y1*x2*p*p-3.0*y1*x3*p*p+3.0*y1*x4*p*p-y4*x3*p+3.0*y4*x2*q*q-3.0*y4*x3*q*q-2.0*y2*x2*q-y2*x3*q-y2*x1*p-2.0*y2*x2*p+2.0*y4*x1*p*q-y4*x1*q-y1*x2*p-3.0*y1*x2*q*q+3.0*y1*x3*q*q+2.0*y1*x1*q-y1*x4*q+4.0*y1*x1*p-2.0*y1*x1+2.0*y2*x2*p*q*q-2.0*y2*x3*p*q*q+4.0*y3*x3*p*q-2.0*y3*x2*p*q*q+2.0*y3*x3*p*q*q+2.0*y3*x3+4.0*y3*x3*q+2.0*y3*x3*p-2.0*y3*x2*q*q+2.0*y3*x3*q*q+2.0*y2*x2*q*q-2.0*y2*x3*q*q)/D;
    BM12(2,2) = -1.0/8.0*(2.0*x2*q*y1-x2*y4*p-3.0*x1*q*y3+3.0*x1*y3*p+2.0*x1*q*y2-x4*y2*p-2.0*x4*y1*p-3.0*x3*q*y1+3.0*x3*y1*p+x2*q*y4+x4*q*y2-x1*y4+x1*y2-2.0*x1*y4*p-x2*y4*p*p*q+x2*y1-3.0*x2*y3+x1*y4*p*p*q+2.0*x2*y3*p*q+x2*y3*p*p*q+3.0*x4*y3-x4*y1-3.0*x3*y2-x1*y3*p*p*q+3.0*x3*y4-x4*y3*p*q*q-2.0*x4*y3*p*q+x4*y2*p*q*q-2.0*y3*x3*p*p*q+2.0*y3*x4*p*p*q-y1*x3*p*p*q+y1*x4*p*p*q+2.0*y4*x3*p*p*q-2.0*y4*x4*p*p*q-x1*y2*p*q*q+x1*y3*p*q*q+y2*x3*p*p*q-y2*x4*p*p*q-2.0*y4*x3*p*q+4.0*y4*x4*p*q+y4*x2*p*q*q-y4*x3*p*q*q+2.0*y2*x1*p*q-4.0*y2*x2*p*q+2.0*y2*x3*p*q+2.0*y1*x2*p*q-2.0*y1*x4*p*q-y1*x2*p*q*q+y1*x3*p*q*q+2.0*y4*x3*p*p-2.0*y4*x4*p*p-2.0*y3*x3*p*p+2.0*y3*x4*p*p+3.0*x2*y3*p*p-3.0*x2*y4*p*p+x2*y3*q-3.0*x1*y3*p*p+3.0*x1*y4*p*p-x4*y3*p+3.0*x4*y2*q*q-3.0*x4*y3*q*q-3.0*x1*y2*q*q+3.0*x1*y3*q*q+3.0*y2*x3*p*p-3.0*y2*x4*p*p-3.0*y1*x3*p*p+3.0*y1*x4*p*p-y4*x3*p+4.0*y4*x4*p+3.0*y4*x2*q*q-3.0*y4*x3*q*q-4.0*y2*x2*q+y2*x3*q-y2*x1*p+2.0*y2*x2*p-2.0*y4*x1*p*q+y4*x1*q-2.0*y4*x4*q-y1*x2*p-3.0*y1*x2*q*q+3.0*y1*x3*q*q+2.0*y2*x2-2.0*y4*x4+y1*x4*q+2.0*y2*x2*p*q*q-2.0*y2*x3*p*q*q-2.0*y3*x2*p*q*q+2.0*y3*x3*p*q*q+2.0*y3*x3*q-2.0*y3*x3*p-2.0*y3*x2*q*q+2.0*y3*x3*q*q+2.0*y2*x2*q*q-2.0*y2*x3*q*q)/D;
    BM12(2,3) = 1.0/8.0*(2.0*x2*q*y1-3.0*x2*y4*p+x1*q*y3+x1*y3*p+2.0*x1*q*y2+2.0*x2*y3*p-3.0*x4*y2*p+x3*q*y1+2.0*x3*y2*p+x3*y1*p-3.0*x2*q*y4-3.0*x4*q*y2-3.0*x1*y4+x1*y2-x2*y4*p*p*q+x2*y1-x2*y3+x1*y4*p*p*q+2.0*x2*y3*p*q+x2*y3*p*p*q+3.0*x4*y3-3.0*x4*y1-x3*y2-x1*y3*p*p*q+3.0*x3*y4+x4*y3*p*q*q+2.0*x4*y3*p*q-x4*y2*p*q*q-2.0*y3*x3*p*p*q+2.0*y3*x4*p*p*q-y1*x3*p*p*q+y1*x4*p*p*q+2.0*y4*x3*p*p*q-2.0*y4*x4*p*p*q+x1*y2*p*q*q-x1*y3*p*q*q+y2*x3*p*p*q-y2*x4*p*p*q+2.0*y4*x3*p*q+2.0*y4*x1*p*q*q-y4*x2*p*q*q+y4*x3*p*q*q-2.0*y4*x4*p*q*q-2.0*y2*x1*p*q+2.0*y2*x3*p*q+4.0*y1*x1*p*q-2.0*y1*x2*p*q-2.0*y1*x4*p*q-2.0*y1*x1*p*q*q+y1*x2*p*q*q-y1*x3*p*q*q+2.0*y1*x4*p*q*q+2.0*y4*x3*p*p-2.0*y4*x4*p*p-2.0*y3*x3*p*p+2.0*y3*x4*p*p+3.0*x2*y3*p*p-3.0*x2*y4*p*p+x2*y3*q-3.0*x1*y3*p*p+3.0*x1*y4*p*p+x4*y3*p+3.0*x4*y2*q*q-3.0*x4*y3*q*q-3.0*x1*y2*q*q+3.0*x1*y3*q*q+3.0*y2*x3*p*p-3.0*y2*x4*p*p-3.0*y1*x3*p*p+3.0*y1*x4*p*p+y4*x3*p+2.0*y4*x4*p-2.0*y4*x1*q*q+3.0*y4*x2*q*q-3.0*y4*x3*q*q+2.0*y4*x4*q*q+y2*x3*q+y2*x1*p-2.0*y4*x1*p*q+y4*x1*q+2.0*y4*x4*q+y1*x2*p+2.0*y1*x1*q*q-3.0*y1*x2*q*q+3.0*y1*x3*q*q-2.0*y1*x4*q*q-4.0*y1*x1*q+y1*x4*q-2.0*y1*x1*p+2.0*y1*x1-4.0*y3*x3*p*q-2.0*y3*x3-2.0*y3*x3*q-4.0*y3*x3*p)/D;

    //Jacobian Matrix.
    Jij(0,0) = -1.0/4.0*x1 + 1.0/4.0*x1*q + 1.0/4.0*x2 - 1.0/4.0*x2*q + 1.0/4.0*x3 + 1.0/4.0*x3*q - 1.0/4.0*x4*q - 1.0/4.0*x4;
    Jij(0,1) = -1.0/4.0*y1 + 1.0/4.0*y1*q + 1.0/4.0*y2 - 1.0/4.0*q*y2 + 1.0/4.0*y3 + 1.0/4.0*q*y3 - 1.0/4.0*y4*q - 1.0/4.0*y4;
    Jij(1,0) = -1.0/4.0*x1 + 1.0/4.0*x1*p - 1.0/4.0*x2 - 1.0/4.0*x2*p + 1.0/4.0*x3 + 1.0/4.0*x3*p - 1.0/4.0*x4*p + 1.0/4.0*x4;
    Jij(1,1) = -1.0/4.0*y1 + 1.0/4.0*y1*p - 1.0/4.0*y2 - 1.0/4.0*y2*p + 1.0/4.0*y3 + 1.0/4.0*y3*p - 1.0/4.0*y4*p + 1.0/4.0*y4;
}

//Compute the initial stiffness matrix of the element
Eigen::MatrixXd 
lin3DShell4::ComputeInitialStiffnessMatrix(){
    //Computes the transformation matrix.
    Eigen::MatrixXd localAxes = ComputeLocalAxes();

    //Gets Node Local coordinates.
    Eigen::MatrixXd xyloc = ComputeLocalCoordinates();

    //Computes the Matrix to enforce constant Tension (wilson).
    Eigen::MatrixXd BM12 = ComputeConstantTensionMatrix();

    //Membrane stiffness matrix component.
    Eigen::MatrixXd StiffnessMatrix(24,24);
    Eigen::MatrixXd PlateStiffness(12,12); 
    Eigen::MatrixXd MembraneStiffness(12,12); 
    
    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Quad", wi, xi);

    //Numerical integration.
    PlateStiffness.fill(0.0);
    MembraneStiffness.fill(0.0);

    for(unsigned int i = 0; i < wi.size(); i++){
        //Compute Strain-Displacement Matrix and Jacobian at Gauss Point.
        Eigen::MatrixXd Jij(2,2);
        Eigen::MatrixXd Bij(3,12);
        Eigen::MatrixXd Pij(1,12);

        //Gets material tangent matrix at Gauss point.
        Eigen::MatrixXd Cij = theSection[i]->GetInitialTangentStiffness();

        //Plate Numerical integration.
        ComputePlateStrainDisplacementMatrix(xi(i,0), xi(i,1), xyloc, Bij, Jij);
        PlateStiffness += wi(i)*fabs(Jij.determinant())*Bij.transpose()*Cij.block(3,3,3,3)*Bij;

        //Membrane Numerical integration.
        ComputeMembraneStrainDisplacementMatrix(xi(i,0), xi(i,1), xyloc, BM12, Bij, Pij, Jij);
        MembraneStiffness += wi(i)*fabs(Jij.determinant())*Bij.transpose()*Cij.block(0,0,3,3)*Bij;
        MembraneStiffness += wi(i)*fabs(Jij.determinant())*Cij(2,2)*Pij.transpose()*Pij;
    }

    //Assembles the total the total stiffness matrix;
    AssemblePlateMembraneEffects(StiffnessMatrix, MembraneStiffness, PlateStiffness);

    StiffnessMatrix = localAxes.transpose()*StiffnessMatrix*localAxes;

    return StiffnessMatrix;
}
