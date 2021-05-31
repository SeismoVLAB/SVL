#include <cmath>
#include <Eigen/LU> 
#include "Material.hpp"
#include "lin2DTria3.hpp"
#include "GaussQuadrature.hpp"
#include "LobattoQuadrature.hpp"
#include "Definitions.hpp"
#include "Profiler.hpp"

//Define VTK cell value for Paraview:
const unsigned int VTKCELL = 5;

//Overload constructor.
lin2DTria3::lin2DTria3(const std::vector<unsigned int> nodes, std::unique_ptr<Material> &material, const double th, const std::string quadrature, const unsigned int nGauss) :
Element("lin2DTria3", nodes, 6, VTKCELL), t(th){
    //The element nodes.
    theNodes.resize(3);

    //Numerical integration rule.
    if(strcasecmp(quadrature.c_str(),"GAUSS") == 0)
        QuadraturePoints = std::make_unique<GaussQuadrature>("Tria", nGauss);
    else if(strcasecmp(quadrature.c_str(),"LOBATTO") == 0)
        QuadraturePoints = std::make_unique<LobattoQuadrature>("Tria", nGauss);

    //The element material. 
    theMaterial = material->CopyMaterial();
}

//Default Destructor.
lin2DTria3::~lin2DTria3(){
    //Does nothing.
}

//Save the material states in the element.
void 
lin2DTria3::CommitState(){
    if(theMaterial->IsViscous()){
        //Jacobian matrix.
        Eigen::MatrixXd Jij = ComputeJacobianMatrix();

        //Compute Strain-Displacement Matrix at Gauss Point.
        Eigen::MatrixXd Bij = ComputeStrainDisplacementMatrix(Jij);

        //Computes Strain Rate vector.
        Eigen::VectorXd strainrate = ComputeStrainRate(Bij);

        //Update the material state.
        theMaterial->UpdateState(strainrate, 2);
    }

    theMaterial->CommitState();
}

//Reverse the material states to previous converged state in this element.
void 
lin2DTria3::ReverseState(){
    //Reverse the material components.
    theMaterial->ReverseState();
}

//Brings the material state to its initial state in this element.
void 
lin2DTria3::InitialState(){
    theMaterial->InitialState();
}

//Update the material states in the element.
void 
lin2DTria3::UpdateState(){
    //Jacobian matrix.
    Eigen::MatrixXd Jij = ComputeJacobianMatrix();

    //Compute Strain-Displacement Matrix at Gauss Point.
    Eigen::MatrixXd Bij = ComputeStrainDisplacementMatrix(Jij);

    //Computes strain vector.
    Eigen::VectorXd strain = ComputeStrain(Bij);

    //Update the material state.
    theMaterial->UpdateState(strain, 1);
}

//Sets the finite element dependance among objects.
void 
lin2DTria3::SetDomain(std::map<unsigned int, std::shared_ptr<Node> > &nodes){
    //Gets the global element connectivity.
    std::vector<unsigned int> conn = GetNodes();

    //Assign the element to mesh node pointer.  
    for(unsigned int i = 0; i < GetNumberOfNodes(); i++){
        theNodes[i] = nodes[conn[i]];
    }
}

//Sets the damping model.
void 
lin2DTria3::SetDamping(const std::shared_ptr<Damping> &damping){
    //The damping model
    theDamping = damping;
}

//Gets the list of total-degree of freedom of this element.
std::vector<unsigned int> 
lin2DTria3::GetTotalDegreeOfFreedom() const{
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
lin2DTria3::GetStrain() const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    Eigen::MatrixXd theStrain(nPoints,3);
    for(unsigned int k = 0; k < nPoints; k++)
        theStrain.row(k) = theMaterial->GetStrain();

    return theStrain;
}

//Returns the material stress at integration points.
Eigen::MatrixXd 
lin2DTria3::GetStress() const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    Eigen::MatrixXd theStress(nPoints,3);
    for(unsigned int k = 0; k < nPoints; k++)
        theStress.row(k) = theMaterial->GetTotalStress();

    return theStress;
}

//Returns the material strain-rate at integration points.
Eigen::MatrixXd 
lin2DTria3::GetStrainRate() const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    Eigen::MatrixXd theStrainRate(nPoints,3);
    for(unsigned int k = 0; k < nPoints; k++)
        theStrainRate.row(k) = theMaterial->GetStrainRate();

    return theStrainRate;
}

//Gets the material strain in section at  coordinate (x3,x2).
Eigen::MatrixXd 
lin2DTria3::GetStrainAt(double UNUSED(x3), double UNUSED(x2)) const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    //Stress at coordinate is define within section.
    Eigen::MatrixXd theStrain(nPoints, 3);
    theStrain.fill(0.0);

    return theStrain;
}

//Gets the material stress in section at  coordinate (x3,x2).
Eigen::MatrixXd 
lin2DTria3::GetStressAt(double UNUSED(x3), double UNUSED(x2)) const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    //Stress at coordinate is define within section.
    Eigen::MatrixXd theStress(nPoints, 3);
    theStress.fill(0.0);

    return theStress;
}

//Gets the element internal response in VTK format.
Eigen::VectorXd 
lin2DTria3::GetVTKResponse(std::string response) const{
    //The VTK response vector.
    Eigen::VectorXd theResponse(6);

    if (strcasecmp(response.c_str(),"Strain") == 0){
        Eigen::MatrixXd strain = GetStrain();
        Eigen::VectorXd Strain = strain.colwise().mean();
        theResponse << Strain(0), Strain(1), 0.0, Strain(2), 0.0, 0.0;
    }
    else if(strcasecmp(response.c_str(),"Stress") == 0){
        Eigen::MatrixXd stress = GetStress();
        Eigen::VectorXd Stress = stress.colwise().mean();
        theResponse << Stress(0), Stress(1), 0.0, Stress(2), 0.0, 0.0; 
    }

    return theResponse;
}

//Computes the element energy for a given deformation.
double 
lin2DTria3::ComputeEnergy(){
    //TODO: Integrate over element volume to compute the energy
    return 0.0;
}

//Compute the mass matrix of the element using gauss-integration.
Eigen::MatrixXd 
lin2DTria3::ComputeMassMatrix(){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Use consistent mass definition:
    Eigen::MatrixXd MassMatrix(6,6);
    MassMatrix.fill(0.0);

    //Gets material properties:
    double rho = theMaterial->GetDensity();

    //Jacobian matrix:
    Eigen::MatrixXd Jij = ComputeJacobianMatrix();

    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Tria", wi, xi);

    //Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
        //Compute Interpolation Function Matrix at Gauss Point:
        Eigen::MatrixXd Hij = ComputeShapeFunctionMatrix(xi(i,0), xi(i,1));

        //Numerical integration:
        MassMatrix += wi(i)*rho*t*fabs(Jij.determinant())*Hij.transpose()*Hij;
    }

    //Lumped Mass Formulation
    if(MassFormulation){
        //Lumped Mass in diagonal terms.
        for (unsigned int i = 0; i < 6; i++){
            for (unsigned int j = 0; j < 6; j++){
                if(i != j){
                    MassMatrix(i,i) += MassMatrix(i,j);
                    MassMatrix(i,j) = 0.0;
                }
            }
        }
    }

    return MassMatrix;
}

//Compute the stiffness matrix of the element using gauss-integration.
Eigen::MatrixXd 
lin2DTria3::ComputeStiffnessMatrix(){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Stiffness matrix definition:
    Eigen::MatrixXd StiffnessMatrix(6,6);
    StiffnessMatrix.fill(0.0);

    //Jacobian matrix.
    Eigen::MatrixXd Jij = ComputeJacobianMatrix();

    //Compute Strain-Displacement Matrix at Gauss Point.
    Eigen::MatrixXd Bij = ComputeStrainDisplacementMatrix(Jij);

    //Gets material tangent matrix at Gauss point.
    Eigen::MatrixXd Cij = theMaterial->GetTangentStiffness();

    //Numerical integration.
    StiffnessMatrix = 0.5*t*fabs(Jij.determinant())*Bij.transpose()*Cij*Bij;

    return StiffnessMatrix;
}

//Compute the damping matrix of the element using gauss-integration.
Eigen::MatrixXd 
lin2DTria3::ComputeDampingMatrix(){
    //Starts profiling this function.
    PROFILE_FUNCTION();
    
    //Damping matrix definition.
    Eigen::MatrixXd DampingMatrix(6,6);
    DampingMatrix.fill(0.0);

    //Material damping contribution.
    if(theMaterial->IsViscous()){        
        //Jacobian matrix.
        Eigen::MatrixXd Jij = ComputeJacobianMatrix();
            
        //Compute Strain-Displacement Matrix at Gauss Point.
        Eigen::MatrixXd Bij = ComputeStrainDisplacementMatrix(Jij);
            
        //Gets material damping matrix at Gauss point.
        Eigen::MatrixXd Dij = theMaterial->GetDamping();
            
        //Numerical integration.
        DampingMatrix += 0.5*t*fabs(Jij.determinant())*Bij.transpose()*Dij*Bij;
    }  

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

    return DampingMatrix;
}

//Compute the PML history matrix for Perfectly-Matched Layer (PML).
Eigen::MatrixXd 
lin2DTria3::ComputePMLMatrix(){
    Eigen::MatrixXd Kpml;
    return Kpml;
}

//Compute the internal forces acting on the element.
Eigen::VectorXd 
lin2DTria3::ComputeInternalForces(){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Jacobian matrix.
    Eigen::MatrixXd Jij = ComputeJacobianMatrix();

    //Compute Strain-Displacement Matrix at Gauss Point.
    Eigen::MatrixXd Bij = ComputeStrainDisplacementMatrix(Jij);

    //Gets material strain at Gauss point.
    Eigen::VectorXd Stress = theMaterial->GetStress();

    //Internal force vector definition:
    Eigen::VectorXd InternalForces = 0.5*t*fabs(Jij.determinant())*Bij.transpose()*Stress;

    return InternalForces;
}

//Compute the elastic, inertial, and vicious forces acting on the element.
Eigen::VectorXd 
lin2DTria3::ComputeInternalDynamicForces(){
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
lin2DTria3::ComputeSurfaceForces(const std::shared_ptr<Load> &surfaceLoad, unsigned int face){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Local surface load vector:
    Eigen::VectorXd surfaceForces(6);
    surfaceForces.fill(0.0);

    //Gets the surface load:
    Eigen::VectorXd qs = surfaceLoad->GetLoadVector();

    if(face == 1){
        //Gets the face coordinates in undeformed configuration. 
        Eigen::VectorXd x1 = theNodes[0]->GetCoordinates();
        Eigen::VectorXd x2 = theNodes[1]->GetCoordinates();

        //Jacobian matrix:
        double detJij = 0.5*t*(x2 - x1).norm();

        //Explicit integration:
        surfaceForces << qs(0)*detJij, qs(1)*detJij, qs(0)*detJij, qs(1)*detJij, 0.0, 0.0;
    }
    if(face == 2){
        //Gets the face coordinates in undeformed configuration. 
        Eigen::VectorXd x2 = theNodes[1]->GetCoordinates();
        Eigen::VectorXd x3 = theNodes[2]->GetCoordinates();

        //Jacobian matrix:
        double detJij = 0.5*t*(x3 - x2).norm();

        //Explicit integration:
        surfaceForces << 0.0, 0.0, qs(0)*detJij, qs(1)*detJij, qs(0)*detJij, qs(1)*detJij;
    }
    if(face == 3){
        //Gets the face coordinates in undeformed configuration. 
        Eigen::VectorXd x1 = theNodes[0]->GetCoordinates();
        Eigen::VectorXd x3 = theNodes[2]->GetCoordinates();

        //Jacobian matrix:
        double detJij = 0.5*t*(x3 - x1).norm();

        //Explicit integration:
        surfaceForces << qs(0)*detJij, qs(1)*detJij, 0.0, 0.0, qs(0)*detJij, qs(1)*detJij;
    }

    return surfaceForces;
}

//Compute the body forces acting on the element.
Eigen::VectorXd 
lin2DTria3::ComputeBodyForces(const std::shared_ptr<Load> &bodyLoad, unsigned int k){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Local body load vector:
    Eigen::VectorXd bodyForces(6);
    bodyForces.fill(0.0);

    //Gets the body force:
    Eigen::VectorXd qb = bodyLoad->GetLoadVector(k);

    //Gets material properties:
    double rho = theMaterial->GetDensity();

    //Jacobian matrix:
    Eigen::MatrixXd Jij = ComputeJacobianMatrix();
    
    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Tria", wi, xi);

    //Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
        //Compute Interpolation Function Matrix at Gauss Point:
        Eigen::MatrixXd Hij = ComputeShapeFunctionMatrix(xi(i,0), xi(i,1));

        //Numerical integration:
        bodyForces += wi(i)*rho*t*fabs(Jij.determinant())*Hij.transpose()*qb;
    }

    return bodyForces;
}

//Compute the domain reduction forces acting on the element.
Eigen::VectorXd 
lin2DTria3::ComputeDomainReductionForces(const std::shared_ptr<Load> &drm, unsigned int k){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Get the Domain-Reduction field motion.
    Eigen::VectorXd x1 = theNodes[0]->GetDomainReductionMotion(k);
    Eigen::VectorXd x2 = theNodes[1]->GetDomainReductionMotion(k);
    Eigen::VectorXd x3 = theNodes[2]->GetDomainReductionMotion(k);

    //Constructs the domain reduction boundary/exterior connectivity.
    std::vector<bool> DRMcond(6);
    std::vector<unsigned int> conn = GetNodes();

    for(unsigned int i = 0; i < conn.size(); i++){
        bool condition = drm->GetDRMCondition(conn[i]);
        DRMcond[2*i  ] = condition;
        DRMcond[2*i+1] = condition;
    }

    //Constructs the displacement, velocity and acceleration vectors. 
    Eigen::VectorXd Uo(6); 
    Eigen::VectorXd Vo(6);
    Eigen::VectorXd Ao(6);
 
    Uo << x1(0), x1(1), x2(0), x2(1), x3(0), x3(1);
    Vo << x1(2), x1(3), x2(2), x2(3), x3(2), x3(3);
    Ao << x1(4), x1(5), x2(4), x2(5), x3(4), x3(5);

    //Computes the mass, damping and stiffness matrices.
    Eigen::MatrixXd MassMatrix = ComputeMassMatrix();
    Eigen::MatrixXd DampingMatrix = ComputeDampingMatrix();
    Eigen::MatrixXd StiffnessMatrix = ComputeStiffnessMatrix();

    //Modifies the Mass, Damping and stiffness matrix.
    for(unsigned int i = 0; i < DRMcond.size(); i++){
        for(unsigned int j = 0; j < DRMcond.size(); j++){
            if(DRMcond[i] == DRMcond[j]){
                MassMatrix(i,j)      = 0.0;
                DampingMatrix(i,j)   = 0.0;
                StiffnessMatrix(i,j) = 0.0;
            }
        }
    }

    //Domain reduction force vector.
    Eigen::VectorXd DRMForces(6);
    DRMForces = MassMatrix*Ao + DampingMatrix*Vo + StiffnessMatrix*Uo;

    return DRMForces;
}

//Update strain in the element.
Eigen::VectorXd 
lin2DTria3::ComputeStrain(const Eigen::MatrixXd &Bij) const{
    //Gets the element coordinates in deformed configuration. 
    Eigen::VectorXd U1 = theNodes[0]->GetDisplacements() + theNodes[0]->GetIncrementalDisplacements();
    Eigen::VectorXd U2 = theNodes[1]->GetDisplacements() + theNodes[1]->GetIncrementalDisplacements();
    Eigen::VectorXd U3 = theNodes[2]->GetDisplacements() + theNodes[2]->GetIncrementalDisplacements();

    Eigen::VectorXd nodalDisplacement(6);
    nodalDisplacement << U1, U2, U3;

    //Strain vector:
    Eigen::VectorXd Strain(3);
    Strain = Bij*nodalDisplacement;

    return Strain;
}

//Update strain rate in the element.
Eigen::VectorXd 
lin2DTria3::ComputeStrainRate(const Eigen::MatrixXd& UNUSED(Bij)) const{
    //TODO: Compute strain rate.
    //Strain vector definition:
    Eigen::VectorXd strainrate(3);
    strainrate.fill(0.0);

    return strainrate;
}

//Computes the jacobian of the transformation. 
Eigen::MatrixXd 
lin2DTria3::ComputeJacobianMatrix() const{
    //Gets the element coordinates in deformed configuration. 
    Eigen::VectorXd X1 = theNodes[0]->GetCoordinates();
    Eigen::VectorXd X2 = theNodes[1]->GetCoordinates();
    Eigen::VectorXd X3 = theNodes[2]->GetCoordinates();

    //Jacobian matrix:
    Eigen::MatrixXd Jij(2,2);
    Jij << X2(0) - X1(0), X2(1) - X1(1),
           X3(0) - X1(0), X3(1) - X1(1);

    return Jij;
}

//Compute Shape Function at Gauss Point:
Eigen::MatrixXd 
lin2DTria3::ComputeShapeFunctionMatrix(const double ri, const double si) const{
    //Shape function coefficients:
    double H11 = 1.0 - ri - si;
    double H22 = ri;
    double H33 = si;

    //Shape function matrix:
    Eigen::MatrixXd Hij(2,6);
    Hij << H11, 0.0, H22, 0.0, H33, 0.0,
           0.0, H11, 0.0, H22, 0.0, H33;

    return Hij;
}

//Evaluates the lumped-mass matrix matrix at a given Gauss point.
Eigen::MatrixXd 
lin2DTria3::ComputeStrainDisplacementMatrix(const Eigen::MatrixXd &Jij) const{
    //Inverse jacobian matrix:
    Eigen::MatrixXd J = Jij.inverse();

    //Strain-displacement matrix coefficients:
    double B11 = -1.0*(J(0,0) + J(0,1)); 
    double B21 = J(0,0); 
    double B31 = J(0,1); 

    double B12 = -1.0*(J(1,1) + J(1,0)); 
    double B22 = J(1,0);
    double B32 = J(1,1);

    //Shape function matrix:
    Eigen::MatrixXd Bij(3,6);
    Bij << B11, 0.0, B21, 0.0, B31, 0.0,
           0.0, B12, 0.0, B22, 0.0, B32,
           B12, B11, B22, B21, B32, B31;

    return Bij;
}

//Compute the initial stiffness matrix of the element using gauss-integration.
Eigen::MatrixXd 
lin2DTria3::ComputeInitialStiffnessMatrix() const{
    //Jacobian matrix.
    Eigen::MatrixXd Jij = ComputeJacobianMatrix();

    //Compute Strain-Displacement Matrix at Gauss Point.
    Eigen::MatrixXd Bij = ComputeStrainDisplacementMatrix(Jij);

    //Gets material tangent matrix at Gauss point.
    Eigen::MatrixXd Cij = theMaterial->GetInitialTangentStiffness();

    //Stiffness matrix definition:
    Eigen::MatrixXd StiffnessMatrix = 0.5*t*fabs(Jij.determinant())*Bij.transpose()*Cij*Bij;

    return StiffnessMatrix;
}
