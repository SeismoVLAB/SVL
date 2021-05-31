#include <cmath>
#include <Eigen/LU> 
#include "Material.hpp"
#include "lin3DTetra4.hpp"
#include "GaussQuadrature.hpp"
#include "LobattoQuadrature.hpp"
#include "Definitions.hpp"
#include "Profiler.hpp"

//Define VTK cell value for Paraview:
const unsigned int VTKCELL = 10;

//Overload constructor.
lin3DTetra4::lin3DTetra4(const std::vector<unsigned int> nodes, std::unique_ptr<Material> &material, const std::string quadrature, const unsigned int nGauss) :
Element("lin3DTetra4", nodes, 12, VTKCELL){
    //The element nodes.
    theNodes.resize(4);

    //Numerical integration rule.
    if(strcasecmp(quadrature.c_str(),"GAUSS") == 0)
        QuadraturePoints = std::make_unique<GaussQuadrature>("Tetra", nGauss);
    else if(strcasecmp(quadrature.c_str(),"LOBATTO") == 0)
        QuadraturePoints = std::make_unique<LobattoQuadrature>("Tetra", nGauss);

    //The element material. 
    theMaterial = material->CopyMaterial();
}

//Destructor.
lin3DTetra4::~lin3DTetra4(){
    //Does nothing.
}

//Save the material states in the element.
void 
lin3DTetra4::CommitState(){
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
lin3DTetra4::ReverseState(){
    //Reverse the material components.
    theMaterial->ReverseState();
}

//Brings the material state to its initial state in this element.
void 
lin3DTetra4::InitialState(){
    //Brings the material components to initial state.
    theMaterial->InitialState();
}

//Update the material states in the element.
void 
lin3DTetra4::UpdateState(){
    //Update material states.
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
lin3DTetra4::SetDomain(std::map<unsigned int, std::shared_ptr<Node> > &nodes){
    //Gets the global element connectivity.
    std::vector<unsigned int> conn = GetNodes();

    //Assign the element to mesh node pointer.  
    for(unsigned int i = 0; i < GetNumberOfNodes(); i++){
        theNodes[i] = nodes[conn[i]];
    }
}

//Sets the damping model.
void 
lin3DTetra4::SetDamping(const std::shared_ptr<Damping> &damping){
    //The damping model
    theDamping = damping;
}

//Gets the list of total-degree of freedom of this element.
std::vector<unsigned int> 
lin3DTetra4::GetTotalDegreeOfFreedom() const{
    //Total number of degree-of-freedom.
    unsigned int nDofs = GetNumberOfDegreeOfFreedom();

    //Reserve memory for the element list of degree-of-freedom.
    std::vector<unsigned int> dofs(nDofs);

    //Construct the element list of degree-of-freedom for assembly.
    for(unsigned int j = 0; j < GetNumberOfNodes(); j++){    
        unsigned int LengthDofs = theNodes[j]->GetNumberOfDegreeOfFreedom();
        std::vector<int> totalDofs = theNodes[j]->GetTotalDegreeOfFreedom();

        for(unsigned int i = 0; i < LengthDofs; i++)
            dofs[i + LengthDofs*j] = totalDofs[i];    
    }

    return dofs;
}

//Returns the material strain at integration points.
Eigen::MatrixXd 
lin3DTetra4::GetStrain() const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    Eigen::MatrixXd theStrain(nPoints,6);
    for(unsigned int k = 0; k < nPoints; k++)
        theStrain.row(k) = theMaterial->GetStrain();

    return theStrain;
}

//Returns the material stress at integration points.
Eigen::MatrixXd 
lin3DTetra4::GetStress() const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    Eigen::MatrixXd theStress(nPoints,6);
    for(unsigned int k = 0; k < nPoints; k++)
        theStress.row(k) = theMaterial->GetTotalStress();

    return theStress;
}

//Returns the material strain-rate at integration points.
Eigen::MatrixXd 
lin3DTetra4::GetStrainRate() const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    Eigen::MatrixXd theStrainRate(nPoints,6);
    for(unsigned int k = 0; k < nPoints; k++)
        theStrainRate.row(k) = theMaterial->GetStrainRate();

    return theStrainRate;
}

//Gets the material strain in section at  coordinate (x3,x2).
Eigen::MatrixXd 
lin3DTetra4::GetStrainAt(double UNUSED(x3), double UNUSED(x2)) const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    //Stress at coordinate is define within section.
    Eigen::MatrixXd theStrain(nPoints, 6);
    theStrain.fill(0.0);

    return theStrain;
}

//Gets the material stress in section at  coordinate (x3,x2).
Eigen::MatrixXd 
lin3DTetra4::GetStressAt(double UNUSED(x3), double UNUSED(x2)) const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    //Stress at coordinate is define within section.
    Eigen::MatrixXd theStress(nPoints, 6);
    theStress.fill(0.0);

    return theStress;
}

//Gets the element internal response in VTK format.
Eigen::VectorXd 
lin3DTetra4::GetVTKResponse(std::string response) const{
    //The VTK response vector.
    Eigen::VectorXd theResponse(6);

    if (strcasecmp(response.c_str(),"Strain") == 0){
        Eigen::MatrixXd strain = GetStrain();
        Eigen::VectorXd Strain = strain.colwise().mean();
        theResponse << Strain(0), Strain(1), Strain(2), Strain(3), Strain(4), Strain(5);
    }
    else if(strcasecmp(response.c_str(),"Stress") == 0){
        Eigen::MatrixXd stress = GetStress();
        Eigen::VectorXd Stress = stress.colwise().mean();
        theResponse << Stress(0), Stress(1), Stress(2), Stress(3), Stress(4), Stress(5);
    }

    return theResponse;
}

//Computes the element energy for a given deformation.
double 
lin3DTetra4::ComputeEnergy(){
    //TODO: Integrate over element volume to compute the energy
    return 0.0;
}

//Compute the mass matrix of the element using gauss-integration.
Eigen::MatrixXd 
lin3DTetra4::ComputeMassMatrix(){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Use consistent mass definition:
    Eigen::MatrixXd MassMatrix(12,12);
    MassMatrix.fill(0.0);

    //Gets material properties:
    double rho = theMaterial->GetDensity();

    //Jacobian matrix:
    Eigen::MatrixXd Jij = ComputeJacobianMatrix();

    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Tetra", wi, xi);

    //Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
        //Compute Interpolation Function Matrix at Gauss Point:
        Eigen::MatrixXd Hij = ComputeShapeFunctionMatrix(xi(i,0), xi(i,1), xi(i,2));

        //Numerical integration:
        MassMatrix += wi(i)*rho*fabs(Jij.determinant())*Hij.transpose()*Hij;
    }

    //Lumped Mass Formulation
    if(MassFormulation){
        //Lumped Mass in diagonal terms.
        for (unsigned int i = 0; i < 12; i++){
            for (unsigned int j = 0; j < 12; j++){
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
lin3DTetra4::ComputeStiffnessMatrix(){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Jacobian matrix.
    Eigen::MatrixXd Jij = ComputeJacobianMatrix();

    //Compute Strain-Displacement Matrix at Gauss Point.
    Eigen::MatrixXd Bij = ComputeStrainDisplacementMatrix(Jij);

    //Gets material tangent matrix at Gauss point.
    Eigen::MatrixXd Cij = theMaterial->GetTangentStiffness();

    //Stiffness matrix definition:
    Eigen::MatrixXd StiffnessMatrix = 1.0/6.0*fabs(Jij.determinant())*Bij.transpose()*Cij*Bij;

    return StiffnessMatrix;
}

//Compute the damping matrix of the element using gauss-integration.
Eigen::MatrixXd 
lin3DTetra4::ComputeDampingMatrix(){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Damping matrix definition
    Eigen::MatrixXd DampingMatrix(12,12);
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
        DampingMatrix += 1.0/6.0*fabs(Jij.determinant())*Bij.transpose()*Dij*Bij;
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
lin3DTetra4::ComputePMLMatrix(){
    Eigen::MatrixXd Kpml;
    return Kpml;
}

//Compute the internal forces acting on the element.
Eigen::VectorXd 
lin3DTetra4::ComputeInternalForces(){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Jacobian matrix.
    Eigen::MatrixXd Jij = ComputeJacobianMatrix();

    //Compute Strain-Displacement Matrix at Gauss Point.
    Eigen::MatrixXd Bij = ComputeStrainDisplacementMatrix(Jij);

    //Gets material strain at Gauss point.
    Eigen::VectorXd Stress = theMaterial->GetStress();

    //Internal Force definition.
    Eigen::VectorXd InternalForces = 1.0/6.0*fabs(Jij.determinant())*Bij.transpose()*Stress;
    
    return InternalForces;
}

//Compute the elastic, inertial, and viscous forces acting on the element.
Eigen::VectorXd 
lin3DTetra4::ComputeInternalDynamicForces(){
    //The Internal dynamic force vector
    Eigen::VectorXd InternalForces;

    if( HasFixedNode(theNodes) ){
        //Allocate memory for velocity/acceleraton. 
        Eigen::VectorXd V(12); 
        Eigen::VectorXd A(12);

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
lin3DTetra4::ComputeSurfaceForces(const std::shared_ptr<Load> &surfaceLoad, unsigned int face){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Local surface load vector:
    Eigen::VectorXd surfaceForces(12);
    surfaceForces.fill(0.0);

    //Gets the surface load:
    Eigen::VectorXd qs = surfaceLoad->GetLoadVector();

    if(face == 1){
        //Gets the face coordinates in undeformed configuration. 
        Eigen::VectorXd x1 = theNodes[0]->GetCoordinates();
        Eigen::VectorXd x2 = theNodes[1]->GetCoordinates();
        Eigen::VectorXd x3 = theNodes[2]->GetCoordinates();

        //vectors along s and t axes.
        Eigen::Vector3d v1, v2;
        v1 = x2 - x1;
        v2 = x3 - x1;

        //Jacobian matrix:
        double detJij = v1.cross(v2).norm();

        //Explicit integration:
        surfaceForces << detJij*qs(0)/6.0, detJij*qs(1)/6.0, detJij*qs(2)/6.0, detJij*qs(0)/6.0, detJij*qs(1)/6.0, detJij*qs(2)/6.0, detJij*qs(0)/6.0, detJij*qs(1)/6.0, detJij*qs(2)/6.0, 0.0, 0.0, 0.0; 
    }
    if(face == 2){
        //Gets the face coordinates in undeformed configuration. 
        Eigen::VectorXd x1 = theNodes[0]->GetCoordinates();
        Eigen::VectorXd x2 = theNodes[1]->GetCoordinates();
        Eigen::VectorXd x4 = theNodes[3]->GetCoordinates();

        //vectors along s and t axes.
        Eigen::Vector3d v1, v2;
        v1 = x2 - x1;
        v2 = x4 - x1;

        //Jacobian matrix:
        double detJij = v1.cross(v2).norm();

        //Explicit integration:
        surfaceForces << detJij*qs(0)/6.0, detJij*qs(1)/6.0, detJij*qs(2)/6.0, detJij*qs(0)/6.0, detJij*qs(1)/6.0, detJij*qs(2)/6.0, 0.0, 0.0, 0.0, detJij*qs(0)/6.0, detJij*qs(1)/6.0, detJij*qs(2)/6.0; 
    }
    if(face == 3){
        //Gets the face coordinates in undeformed configuration. 
        Eigen::VectorXd x2 = theNodes[1]->GetCoordinates();
        Eigen::VectorXd x3 = theNodes[2]->GetCoordinates();
        Eigen::VectorXd x4 = theNodes[3]->GetCoordinates();

        //vectors along s and t axes.
        Eigen::Vector3d v1, v2;
        v1 = x2 - x4;
        v2 = x3 - x4;

        //Jacobian matrix:
        double detJij = v1.cross(v2).norm();

        //Explicit integration:
        surfaceForces << 0.0, 0.0, 0.0, detJij*qs(0)/6.0, detJij*qs(1)/6.0, detJij*qs(2)/6.0, detJij*qs(0)/6.0, detJij*qs(1)/6.0, detJij*qs(2)/6.0, detJij*qs(0)/6.0, detJij*qs(1)/6.0, detJij*qs(2)/6.0; 

    }
    if(face == 4){
        //Gets the face coordinates in undeformed configuration. 
        Eigen::VectorXd x1 = theNodes[0]->GetCoordinates();
        Eigen::VectorXd x3 = theNodes[2]->GetCoordinates();
        Eigen::VectorXd x4 = theNodes[3]->GetCoordinates();

        //vectors along s and t axes.
        Eigen::Vector3d v1, v2;
        v1 = x3 - x1;
        v2 = x4 - x1;

        //Jacobian matrix:
        double detJij = v1.cross(v2).norm();

        //Explicit integration:
        surfaceForces << detJij*qs(0)/6.0, detJij*qs(1)/6.0, detJij*qs(2)/6.0, 0.0, 0.0, 0.0, detJij*qs(0)/6.0, detJij*qs(1)/6.0, detJij*qs(2)/6.0, detJij*qs(0)/6.0, detJij*qs(1)/6.0, detJij*qs(2)/6.0; 
    }

    return surfaceForces;
}

//Compute the body forces acting on the element.
Eigen::VectorXd 
lin3DTetra4::ComputeBodyForces(const std::shared_ptr<Load> &bodyLoad, unsigned int k){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Local body load vector.
    Eigen::VectorXd bodyForces(12);
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
    QuadraturePoints->GetQuadraturePoints("Tetra", wi, xi);

    //Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
        //Compute Interpolation Function Matrix at Gauss Point:
        Eigen::MatrixXd Hij = ComputeShapeFunctionMatrix(xi(i,0), xi(i,1), xi(i,2));

        //Numerical integration:
        bodyForces += wi(i)*rho*fabs(Jij.determinant())*Hij.transpose()*qb;
    }
    
    return bodyForces;
}

//Compute the domain reduction forces acting on the element.
Eigen::VectorXd 
lin3DTetra4::ComputeDomainReductionForces(const std::shared_ptr<Load> &drm, unsigned int k){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Local domain-reduction load vector.
    Eigen::VectorXd DRMForces(12);

    //Get the Domain-Reduction field motion.
    Eigen::VectorXd x1 = theNodes[0]->GetDomainReductionMotion(k);
    Eigen::VectorXd x2 = theNodes[1]->GetDomainReductionMotion(k);
    Eigen::VectorXd x3 = theNodes[2]->GetDomainReductionMotion(k);
    Eigen::VectorXd x4 = theNodes[3]->GetDomainReductionMotion(k);

    //Constructs the domain reduction boundary/exterior connectivity.
    std::vector<bool> DRMcond(12);
    std::vector<unsigned int> conn = GetNodes();

    for(unsigned int i = 0; i < conn.size(); i ++){
        bool condition = drm->GetDRMCondition(conn[i]);
        DRMcond[3*i  ] = condition;
        DRMcond[3*i+1] = condition;
        DRMcond[3*i+2] = condition;
    }

    //Constructs the displacement, velocity and acceleration vectors. 
    Eigen::VectorXd Uo(12); 
    Eigen::VectorXd Vo(12);
    Eigen::VectorXd Ao(12);
 
    Uo << x1(0), x1(1), x1(2), x2(0), x2(1), x2(2), x3(0), x3(1), x3(2), x4(0), x4(1), x4(2);
    Vo << x1(3), x1(4), x1(5), x2(3), x2(4), x2(5), x3(3), x3(4), x3(5), x4(3), x4(4), x4(5);
    Ao << x1(6), x1(7), x1(8), x2(6), x2(7), x2(8), x3(6), x3(7), x3(8), x4(6), x4(7), x4(8);

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
    DRMForces = MassMatrix*Ao + DampingMatrix*Vo + StiffnessMatrix*Uo;

    return DRMForces;
}

//Update strain in the element.
Eigen::VectorXd 
lin3DTetra4::ComputeStrain(const Eigen::MatrixXd &Bij) const{
    //Gets the element coordinates in deformed configuration. 
    Eigen::VectorXd U1 = theNodes[0]->GetDisplacements() + theNodes[0]->GetIncrementalDisplacements();
    Eigen::VectorXd U2 = theNodes[1]->GetDisplacements() + theNodes[1]->GetIncrementalDisplacements();
    Eigen::VectorXd U3 = theNodes[2]->GetDisplacements() + theNodes[2]->GetIncrementalDisplacements();
    Eigen::VectorXd U4 = theNodes[3]->GetDisplacements() + theNodes[3]->GetIncrementalDisplacements();

    Eigen::VectorXd nodalDisplacement(12);
    nodalDisplacement << U1, U2, U3, U4;

    //Strain vector:
    Eigen::VectorXd Strain = Bij*nodalDisplacement;

    return Strain;
}

//Update strain rate in the element.
Eigen::VectorXd 
lin3DTetra4::ComputeStrainRate(const Eigen::MatrixXd& UNUSED(Bij)) const{
    //TODO: Compute strain rate.
    //Strain vector definition:
    Eigen::VectorXd strainrate(6);
    strainrate.fill(0.0);

    return strainrate;
}

//Computes the jacobian of the transformation. 
Eigen::MatrixXd 
lin3DTetra4::ComputeJacobianMatrix() const{
    //Gets the element coordinates in deformed configuration. 
    Eigen::VectorXd X1 = theNodes[0]->GetCoordinates();
    Eigen::VectorXd X2 = theNodes[1]->GetCoordinates();
    Eigen::VectorXd X3 = theNodes[2]->GetCoordinates();
    Eigen::VectorXd X4 = theNodes[3]->GetCoordinates();

    //Jacobian matrix:
    Eigen::MatrixXd Jij(3,3);
    Jij << X2(0) - X1(0), X2(1) - X1(1), X2(2) - X1(2),
           X3(0) - X1(0), X3(1) - X1(1), X3(2) - X1(2),
           X4(0) - X1(0), X4(1) - X1(1), X4(2) - X1(2);

    return Jij;
}

//Compute Shape Function at Gauss Point:
Eigen::MatrixXd 
lin3DTetra4::ComputeShapeFunctionMatrix(const double ri, const double si, const double ti) const{
    //Shape function coefficients:
    double H11 = 1.0 - ri - si - ti;
    double H22 = ri;
    double H33 = si;
    double H44 = ti;

    //Shape function matrix:
    Eigen::MatrixXd Hij(3,12);
    Hij << H11, 0.0, 0.0, H22, 0.0, 0.0, H33, 0.0, 0.0, H44, 0.0, 0.0,
           0.0, H11, 0.0, 0.0, H22, 0.0, 0.0, H33, 0.0, 0.0, H44, 0.0,
           0.0, 0.0, H11, 0.0, 0.0, H22, 0.0, 0.0, H33, 0.0, 0.0, H44;

    return Hij;
}

//Evaluates the deformation matrix at a given Gauss point.
Eigen::MatrixXd 
lin3DTetra4::ComputeStrainDisplacementMatrix(const Eigen::MatrixXd &Jij) const{
    //Inverse jacobian matrix:
    Eigen::MatrixXd J = Jij.inverse();

    //Strain-displacement matrix coefficients:
    double B11 = -J(0,0) - J(0,1) - J(0,2);
    double B21 = J(0,0);
    double B31 = J(0,1);
    double B41 = J(0,2);

    double B12 = -J(1,0) - J(1,1) - J(1,2);
    double B22 = J(1,0);
    double B32 = J(1,1);
    double B42 = J(1,2);

    double B13 = -J(2,0) - J(2,1) - J(2,2);
    double B23 = J(2,0);
    double B33 = J(2,1);
    double B43 = J(2,2);


    //Deformation matrix definition:
    Eigen::MatrixXd Bij(6,12);
    Bij <<  B11, 0.0, 0.0, B21, 0.0, 0.0, B31, 0.0, 0.0, B41, 0.0, 0.0,
            0.0, B12, 0.0, 0.0, B22, 0.0, 0.0, B32, 0.0, 0.0, B42, 0.0,
            0.0, 0.0, B13, 0.0, 0.0, B23, 0.0, 0.0, B33, 0.0, 0.0, B43,
            B12, B11, 0.0, B22, B21, 0.0, B32, B31, 0.0, B42, B41, 0.0,
            0.0, B13, B12, 0.0, B23, B22, 0.0, B33, B32, 0.0, B43, B42,
            B13, 0.0, B11, B23, 0.0, B21, B33, 0.0, B31, B43, 0.0, B41;

    return Bij;
}

//Compute the initial stiffness matrix of the element using gauss-integration.
Eigen::MatrixXd 
lin3DTetra4::ComputeInitialStiffnessMatrix() const{
    //Jacobian matrix.
    Eigen::MatrixXd Jij = ComputeJacobianMatrix();

    //Compute Strain-Displacement Matrix at Gauss Point.
    Eigen::MatrixXd Bij = ComputeStrainDisplacementMatrix(Jij);

    //Gets material tangent matrix at Gauss point.
    Eigen::MatrixXd Cij = theMaterial->GetInitialTangentStiffness();

    //Stiffness matrix definition:
    Eigen::MatrixXd StiffnessMatrix = 1.0/6.0*fabs(Jij.determinant())*Bij.transpose()*Cij*Bij;

    return StiffnessMatrix;
}
