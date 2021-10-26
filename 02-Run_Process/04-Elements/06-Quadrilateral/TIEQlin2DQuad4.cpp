#include <cmath>
#include <string>
#include <Eigen/LU> 
#include "Material.hpp"
#include "TIEQlin2DQuad4.hpp"
#include "GaussQuadrature.hpp"
#include "LobattoQuadrature.hpp"
#include "Definitions.hpp"
#include "Profiler.hpp"

//Define constant value:
const double PI = 3.1415926535897932;

//Overload constructor.
TIEQlin2DQuad4::TIEQlin2DQuad4(const std::vector<unsigned int> nodes, std::unique_ptr<Material> &material, const double th, const std::string quadrature, const unsigned int nGauss, const std::string Type, const double zref, const double cf1, const double cf2, const double eref) :
Element("TIEQlin2DQuad4", nodes, 8, VTK_LINEAR_QUAD, GROUP_ELEMENT_QUAD), t(th), cf1(cf1), cf2(cf2), zref(zref), eref(eref), Type(Type) {
    //The element nodes.
    theNodes.resize(4);

    //Numerical integration rule.
    if(strcasecmp(quadrature.c_str(),"GAUSS") == 0)
        QuadraturePoints = std::make_unique<GaussQuadrature>("Quad", nGauss);
    else if(strcasecmp(quadrature.c_str(),"LOBATTO") == 0)
        QuadraturePoints = std::make_unique<LobattoQuadrature>("Quad", nGauss);

    //The element material. 
    theMaterial.resize(nGauss);
    for(unsigned int i = 0; i < nGauss; i++)
        theMaterial[i] = material->CopyMaterial();
}

//Destructor.
TIEQlin2DQuad4::~TIEQlin2DQuad4(){
    //Does nothing.
}

//Save the material states in the element.
void 
TIEQlin2DQuad4::CommitState(){
    //Updates the viscous material components.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    for(unsigned int k = 0; k < nPoints; k++)
        theMaterial[k]->CommitState();
}

//Reverse the material states to previous converged state in this element.
void 
TIEQlin2DQuad4::ReverseState(){
    //Reverse the material components.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    for(unsigned int k = 0; k < nPoints; k++)
        theMaterial[k]->ReverseState();
}

//Brings the material state to its initial state in this element.
void 
TIEQlin2DQuad4::InitialState(){
    //Brings the material components to initial state.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    for(unsigned int k = 0; k < nPoints; k++)
        theMaterial[k]->InitialState();
}

//Update the material states in the element.
void 
TIEQlin2DQuad4::UpdateState(){
    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Quad", wi, xi);

    //Gets the element coordinates in deformed configuration. 
    Eigen::VectorXd X1 = theNodes[0]->GetCoordinates();
    Eigen::VectorXd X2 = theNodes[1]->GetCoordinates();
    Eigen::VectorXd X3 = theNodes[2]->GetCoordinates();
    Eigen::VectorXd X4 = theNodes[3]->GetCoordinates();

    Eigen::VectorXd X(8);
    X << X1, X2, X3, X4;


    //Update material states.
    for(unsigned int k = 0; k < wi.size(); k++){
        //Jacobian matrix.
        Eigen::MatrixXd Jij = ComputeJacobianMatrix(xi(k,0), xi(k,1));

        //Compute Strain-Displacement Matrix at Gauss Point.
        Eigen::MatrixXd Bij = ComputeStrainDisplacementMatrix(xi(k,0), xi(k,1), Jij);

        //Computes strain vector.
        Eigen::VectorXd strain = ComputeStrain(Bij);

        double Gmax = theMaterial[k]->GetShearModulus();
        double rho  = theMaterial[k]->GetDensity();
        double vs   = sqrt(Gmax/rho);

        Eigen::MatrixXd Hij = ComputeShapeFunctionMatrix(xi(k,0), xi(k,1));

        Eigen::VectorXd XGauss = Hij*X;

        double z = XGauss(1);

        Eigen::VectorXd EqLinParam = ComputeGGmaxDamping(vs,z,rho);

        //double GGmax = EqLinParam(0);  [-Wunused-variable]

        theMaterial[k]->UpdateState(strain, 1);
    }
}

//Sets the finite element dependance among objects.
void 
TIEQlin2DQuad4::SetDomain(std::map<unsigned int, std::shared_ptr<Node> > &nodes){
    //Gets the global element connectivity.
    std::vector<unsigned int> conn = GetNodes();

    //Assign the element to mesh node pointer.  
    for(unsigned int i = 0; i < GetNumberOfNodes(); i++){
        theNodes[i] = nodes[conn[i]];
    }
}

//Sets the damping model.
void 
TIEQlin2DQuad4::SetDamping(const std::shared_ptr<Damping> &damping){
    //The damping model
    theDamping = damping;
}

//Gets the list of total-degree of freedom of this element.
std::vector<unsigned int> 
TIEQlin2DQuad4::GetTotalDegreeOfFreedom() const{
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

//Returns the material strain at integration points.
Eigen::MatrixXd 
TIEQlin2DQuad4::GetStrain() const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    Eigen::MatrixXd theStrain(nPoints,3);
    for(unsigned int k = 0; k < nPoints; k++)
        theStrain.row(k) = theMaterial[k]->GetStrain();

    return theStrain;
}

//Returns the material stress at integration points.
Eigen::MatrixXd 
TIEQlin2DQuad4::GetStress() const{
    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Quad", wi, xi);
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    //Gets the element coordinates in deformed configuration. 
    Eigen::VectorXd X1 = theNodes[0]->GetCoordinates();
    Eigen::VectorXd X2 = theNodes[1]->GetCoordinates();
    Eigen::VectorXd X3 = theNodes[2]->GetCoordinates();
    Eigen::VectorXd X4 = theNodes[3]->GetCoordinates();

    Eigen::VectorXd X(8);
    X << X1, X2, X3, X4;

    Eigen::MatrixXd theStress(nPoints,3);

    //Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
        //Gets material tangent matrix at Gauss point.
        Eigen::MatrixXd Cij = theMaterial[i]->GetTangentStiffness();

        double Gmax = theMaterial[i]->GetShearModulus();
        double rho  = theMaterial[i]->GetDensity();
        double vs   = sqrt(Gmax/rho);

        Eigen::MatrixXd Hij = ComputeShapeFunctionMatrix(xi(i,0), xi(i,1));

        Eigen::VectorXd XGauss = Hij*X;

        double z = XGauss(1);

        Eigen::VectorXd EqLinParam = ComputeGGmaxDamping(vs,z,rho);

        double GGmax = EqLinParam(0);

        theStress.row(i) = GGmax*(theMaterial[i]->GetTotalStress());
    }

    return theStress;
}

//Returns the material strain-rate at integration points.
Eigen::MatrixXd 
TIEQlin2DQuad4::GetStrainRate() const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    Eigen::MatrixXd theStrainRate(nPoints,3);
    for(unsigned int k = 0; k < nPoints; k++)
        theStrainRate.row(k) = theMaterial[k]->GetStrainRate();

    return theStrainRate;
}

//Gets the material strain in section at  coordinate (x3,x2).
Eigen::MatrixXd 
TIEQlin2DQuad4::GetStrainAt(double UNUSED(x3), double UNUSED(x2)) const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    //Stress at coordinate is define within section.
    Eigen::MatrixXd theStrain(nPoints, 3);
    theStrain.fill(0.0);

    return theStrain;
}

//Gets the material stress in section at  coordinate (x3,x2).
Eigen::MatrixXd 
TIEQlin2DQuad4::GetStressAt(double UNUSED(x3), double UNUSED(x2)) const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    //Stress at coordinate is define within section.
    Eigen::MatrixXd theStress(nPoints, 3);
    theStress.fill(0.0);

    return theStress;
}

//Gets the element internal response in VTK format.
Eigen::VectorXd 
TIEQlin2DQuad4::GetVTKResponse(std::string response) const{
    //The VTK response vector.
    Eigen::VectorXd theResponse(18);

    if (strcasecmp(response.c_str(),"Strain") == 0){
        Eigen::MatrixXd strain = GetStrain();
        Eigen::VectorXd Strain = strain.colwise().mean();
        theResponse << Strain(0), Strain(1), 0.0, Strain(2), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    }
    else if(strcasecmp(response.c_str(),"Stress") == 0){
        Eigen::MatrixXd stress = GetStress();
        Eigen::VectorXd Stress = stress.colwise().mean();
        theResponse << Stress(0), Stress(1), 0.0, Stress(2), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0; 
    }

    return theResponse;
}

//Computes the element energy for a given deformation.
double 
TIEQlin2DQuad4::ComputeEnergy(){
    //TODO: Integrate over element volume to compute the energy
    return 0.0;
}

//Compute the mass matrix of the element using gauss-integration.
Eigen::MatrixXd 
TIEQlin2DQuad4::ComputeMassMatrix(){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Use consistent mass definition:
    Eigen::MatrixXd MassMatrix(8, 8);
    MassMatrix.fill(0.0);

    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Quad", wi, xi);

    //Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
        //Gets material properties:
        double rho = theMaterial[i]->GetDensity();

        //Jacobian matrix:
        Eigen::MatrixXd Jij = ComputeJacobianMatrix(xi(i,0), xi(i,1));

        //Compute Interpolation Function Matrix at Gauss Point:
        Eigen::MatrixXd Hij = ComputeShapeFunctionMatrix(xi(i,0), xi(i,1));

        //Numerical integration:
        MassMatrix += wi(i)*rho*t*fabs(Jij.determinant())*Hij.transpose()*Hij;
    }

    //Lumped Mass Formulation
    if(MassFormulation){
        //Lumped Mass in diagonal terms.
        for (unsigned int i = 0; i < 8; i++){
            for (unsigned int j = 0; j < 8; j++){
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
TIEQlin2DQuad4::ComputeStiffnessMatrix(){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Stiffness matrix definition:
    Eigen::MatrixXd StiffnessMatrix(8, 8);
    StiffnessMatrix.fill(0.0);

    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Quad", wi, xi);

    //Gets the element coordinates in deformed configuration. 
    Eigen::VectorXd X1 = theNodes[0]->GetCoordinates();
    Eigen::VectorXd X2 = theNodes[1]->GetCoordinates();
    Eigen::VectorXd X3 = theNodes[2]->GetCoordinates();
    Eigen::VectorXd X4 = theNodes[3]->GetCoordinates();

    Eigen::VectorXd X(8);
    X << X1, X2, X3, X4;

    //Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
        //Jacobian matrix.
        Eigen::MatrixXd Jij = ComputeJacobianMatrix(xi(i,0), xi(i,1));

        //Compute Strain-Displacement Matrix at Gauss Point.
        Eigen::MatrixXd Bij = ComputeStrainDisplacementMatrix(xi(i,0), xi(i,1), Jij);

        //Gets material tangent matrix at Gauss point.
        Eigen::MatrixXd Cij = theMaterial[i]->GetTangentStiffness();

        double Gmax = theMaterial[i]->GetShearModulus();
        double rho  = theMaterial[i]->GetDensity();
        double vs   = sqrt(Gmax/rho);

        Eigen::MatrixXd Hij = ComputeShapeFunctionMatrix(xi(i,0), xi(i,1));

        Eigen::VectorXd XGauss = Hij*X;

        double z = XGauss(1);

        Eigen::VectorXd EqLinParam = ComputeGGmaxDamping(vs,z,rho);

        double GGmax = EqLinParam(0);

        //Numerical integration.
        StiffnessMatrix += GGmax*wi(i)*t*fabs(Jij.determinant())*Bij.transpose()*Cij*Bij;
    }

    return StiffnessMatrix;
}

//Compute the damping matrix of the element using gauss-integration.
Eigen::MatrixXd 
TIEQlin2DQuad4::ComputeDampingMatrix(){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Damping matrix definition.
    Eigen::MatrixXd DampingMatrix(8, 8);
    DampingMatrix.fill(0.0);

    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Quad", wi, xi);

    //Gets the element coordinates in deformed configuration. 
    Eigen::VectorXd X1 = theNodes[0]->GetCoordinates();
    Eigen::VectorXd X2 = theNodes[1]->GetCoordinates();
    Eigen::VectorXd X3 = theNodes[2]->GetCoordinates();
    Eigen::VectorXd X4 = theNodes[3]->GetCoordinates();

    Eigen::VectorXd X(8);
    X << X1, X2, X3, X4;

    //Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
        //Jacobian matrix.
        Eigen::MatrixXd Jij = ComputeJacobianMatrix(xi(i,0), xi(i,1));

        //Compute Strain-Displacement Matrix at Gauss Point.
        Eigen::MatrixXd Bij = ComputeStrainDisplacementMatrix(xi(i,0), xi(i,1), Jij);

        //Compute Interpolation Function Matrix at Gauss Point:
        Eigen::MatrixXd Hij = ComputeShapeFunctionMatrix(xi(i,0), xi(i,1));

        //Gets material tangent matrix at Gauss point.
        //----------------------------------------------
        //Eigen::MatrixXd Cij = theMaterial[i]->GetTangentStiffness();
        Eigen::MatrixXd Cij = theMaterial[i]->GetInitialTangentStiffness(); 
        //----------------------------------------------

        double Gmax = theMaterial[i]->GetShearModulus();
        double rho  = theMaterial[i]->GetDensity();
        double vs   = sqrt(Gmax/rho);

        Eigen::VectorXd XGauss = Hij*X;

        double z = XGauss(1);

        Eigen::VectorXd EqLinParam = ComputeGGmaxDamping(vs,z,rho);

		double xi_R = EqLinParam(1);
		double GGmax = EqLinParam(0);

		double d  = cf2/4.0/cf1 - cf1/4.0/cf2;
		double aM = (PI*cf2*xi_R - PI*cf1*xi_R)/d;
		double aK = (-xi_R/4.0/PI/cf2 + xi_R/4.0/PI/cf1)/d;

        //Numerical integration.
        DampingMatrix += aM*wi(i)*rho*t*fabs(Jij.determinant())*Hij.transpose()*Hij + aK*GGmax*wi(i)*t*fabs(Jij.determinant())*Bij.transpose()*Cij*Bij;
    }
    
    return DampingMatrix;
}

//Compute the PML history matrix for Perfectly-Matched Layer (PML).
Eigen::MatrixXd 
TIEQlin2DQuad4::ComputePMLMatrix(){
    Eigen::MatrixXd Kpml;
    return Kpml;
}

//Compute the internal forces acting on the element.
Eigen::VectorXd 
TIEQlin2DQuad4::ComputeInternalForces(){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Stiffness matrix definition:
    Eigen::VectorXd InternalForces(8);
    InternalForces.fill(0.0);

    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Quad", wi, xi);

    //Gets the element coordinates in deformed configuration. 
    Eigen::VectorXd X1 = theNodes[0]->GetCoordinates();
    Eigen::VectorXd X2 = theNodes[1]->GetCoordinates();
    Eigen::VectorXd X3 = theNodes[2]->GetCoordinates();
    Eigen::VectorXd X4 = theNodes[3]->GetCoordinates();

    Eigen::VectorXd X(8);
    X << X1, X2, X3, X4;

    //Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
        //Jacobian matrix.
        Eigen::MatrixXd Jij = ComputeJacobianMatrix(xi(i,0), xi(i,1));

        //Compute Strain-Displacement Matrix at Gauss Point.
        Eigen::MatrixXd Bij = ComputeStrainDisplacementMatrix(xi(i,0), xi(i,1), Jij);

        //Gets material strain at Gauss point.
        Eigen::VectorXd Stress = theMaterial[i]->GetStress();

        double Gmax = theMaterial[i]->GetShearModulus();
        double rho  = theMaterial[i]->GetDensity();
        double vs   = sqrt(Gmax/rho);

        Eigen::MatrixXd Hij = ComputeShapeFunctionMatrix(xi(i,0), xi(i,1));

        Eigen::VectorXd XGauss = Hij*X;

        double z = XGauss(1);

        Eigen::VectorXd EqLinParam = ComputeGGmaxDamping(vs,z,rho);

        double GGmax = EqLinParam(0);

        //Numerical integration.
        InternalForces += GGmax*wi(i)*t*fabs(Jij.determinant())*Bij.transpose()*Stress;
    }

    return InternalForces;
}

//Compute the elastic, inertial, and viscous forces acting on the element.
Eigen::VectorXd 
TIEQlin2DQuad4::ComputeInternalDynamicForces(){
    //The Internal dynamic force vector
    Eigen::VectorXd InternalForces;

    if( HasFixedNode(theNodes) ){
        //Allocate memory for velocity/acceleraton. 
        Eigen::VectorXd V(8); 
        Eigen::VectorXd A(8);

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
TIEQlin2DQuad4::ComputeSurfaceForces(const std::shared_ptr<Load> &surfaceLoad, unsigned int face){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Local surface load vector:
    Eigen::VectorXd surfaceForces(8);
    surfaceForces.fill(0.0);

    //Gets the surface load:
    Eigen::VectorXd qs = surfaceLoad->GetLoadVector();

    //Coordinate of Gauss points.
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Line", wi, xi);

    if(face == 1){
        //Gets the face coordinates in undeformed configuration. 
        Eigen::VectorXd x1 = theNodes[0]->GetCoordinates();
        Eigen::VectorXd x2 = theNodes[1]->GetCoordinates();

        //Numerical integration in local axis r:
        for(unsigned int i = 0; i < wi.size(); i++){
            //Jacobian matrix:
            double detJij = 0.5*(x2 - x1).norm();

            //Compute Strain-Displacement Matrix at Gauss Point:
            Eigen::MatrixXd Hij = ComputeShapeFunctionMatrix(xi(i,0), -1.0);

            //Numerical integration:
            surfaceForces += wi(i)*t*detJij*Hij.transpose()*qs;
        }
    }
    if(face == 2){
        //Gets the face coordinates in undeformed configuration. 
        Eigen::VectorXd x2 = theNodes[1]->GetCoordinates();
        Eigen::VectorXd x3 = theNodes[2]->GetCoordinates();

        //Numerical integration in local axis s:
        for(unsigned int i = 0; i < wi.size(); i++){
            //Jacobian matrix:
            double detJij = 0.5*(x3 - x2).norm();

            //Compute Strain-Displacement Matrix at Gauss Point:
            Eigen::MatrixXd Hij = ComputeShapeFunctionMatrix(1.0, xi(i,0));

            //Numerical integration:
            surfaceForces += wi(i)*t*detJij*Hij.transpose()*qs;
        }
    }
    if(face == 3){
        //Gets the face coordinates in undeformed configuration. 
        Eigen::VectorXd x3 = theNodes[2]->GetCoordinates();
        Eigen::VectorXd x4 = theNodes[3]->GetCoordinates();

        //Numerical integration in local axis r:
        for(unsigned int i = 0; i < wi.size(); i++){
            //Jacobian matrix:
            double detJij = 0.5*(x3 - x4).norm();

            //Compute Strain-Displacement Matrix at Gauss Point:
            Eigen::MatrixXd Hij = ComputeShapeFunctionMatrix(xi(i,0), 1.0);

            //Numerical integration:
            surfaceForces += wi(i)*t*detJij*Hij.transpose()*qs;
        }
    }
    if(face == 4){
        //Gets the face coordinates in undeformed configuration. 
        Eigen::VectorXd x1 = theNodes[0]->GetCoordinates();
        Eigen::VectorXd x4 = theNodes[3]->GetCoordinates();

        //Numerical integration in local axis s:
        for(unsigned int i = 0; i < wi.size(); i++){
            //Jacobian matrix:
            double detJij = 0.5*(x4 - x1).norm();

            //Compute Strain-Displacement Matrix at Gauss Point:
            Eigen::MatrixXd Hij = ComputeShapeFunctionMatrix(-1.0, xi(i,0));

            //Numerical integration:
            surfaceForces += wi(i)*t*detJij*Hij.transpose()*qs;
        }
    }

    return surfaceForces;
}

//Compute the body forces acting on the element.
Eigen::VectorXd 
TIEQlin2DQuad4::ComputeBodyForces(const std::shared_ptr<Load> &bodyLoad, unsigned int k){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Local body load vector:
    Eigen::VectorXd bodyForces(8);
    bodyForces.fill(0.0);

    //Gets the body force:
    Eigen::VectorXd qb = bodyLoad->GetLoadVector(k);
    
    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Quad", wi, xi);

    //Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
        //Gets material properties:
        double rho = theMaterial[i]->GetDensity();

        //Jacobian matrix:
        Eigen::MatrixXd Jij = ComputeJacobianMatrix(xi(i,0), xi(i,1));

        //Compute Interpolation Function Matrix at Gauss Point:
        Eigen::MatrixXd Hij = ComputeShapeFunctionMatrix(xi(i,0), xi(i,1));

        //Numerical integration:
        bodyForces += wi(i)*rho*t*fabs(Jij.determinant())*Hij.transpose()*qb;
    }

    return bodyForces;
}

//Compute the domain reduction forces acting on the element.
Eigen::VectorXd 
TIEQlin2DQuad4::ComputeDomainReductionForces(const std::shared_ptr<Load> &drm, unsigned int k){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Get the Domain-Reduction field motion.
    Eigen::VectorXd x1 = theNodes[0]->GetDomainReductionMotion(k);
    Eigen::VectorXd x2 = theNodes[1]->GetDomainReductionMotion(k);
    Eigen::VectorXd x3 = theNodes[2]->GetDomainReductionMotion(k);
    Eigen::VectorXd x4 = theNodes[3]->GetDomainReductionMotion(k);

    //Constructs the domain reduction boundary/exterior connectivity.
    std::vector<bool> DRMcond(8);
    std::vector<unsigned int> conn = GetNodes();

    for(unsigned int i = 0; i < conn.size(); i++){
        bool condition = drm->GetDRMCondition(conn[i]);
        DRMcond[2*i  ] = condition;
        DRMcond[2*i+1] = condition;
    }

    //Constructs the displacement, velocity and acceleration vectors. 
    Eigen::VectorXd Uo(8); 
    Eigen::VectorXd Vo(8);
    Eigen::VectorXd Ao(8);
 
    Uo << x1(0), x1(1), x2(0), x2(1), x3(0), x3(1), x4(0), x4(1);
    Vo << x1(2), x1(3), x2(2), x2(3), x3(2), x3(3), x4(2), x4(3);
    Ao << x1(4), x1(5), x2(4), x2(5), x3(4), x3(5), x4(4), x4(5);

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
    Eigen::VectorXd DRMForces(8);
    DRMForces = MassMatrix*Ao + DampingMatrix*Vo + StiffnessMatrix*Uo;

    return DRMForces;
}

//Update strain in the element.
Eigen::VectorXd 
TIEQlin2DQuad4::ComputeStrain(const Eigen::MatrixXd &Bij) const{
    //Gets the element coordinates in deformed configuration. 
    Eigen::VectorXd U1 = theNodes[0]->GetDisplacements() + theNodes[0]->GetIncrementalDisplacements();
    Eigen::VectorXd U2 = theNodes[1]->GetDisplacements() + theNodes[1]->GetIncrementalDisplacements();
    Eigen::VectorXd U3 = theNodes[2]->GetDisplacements() + theNodes[2]->GetIncrementalDisplacements();
    Eigen::VectorXd U4 = theNodes[3]->GetDisplacements() + theNodes[3]->GetIncrementalDisplacements();

    Eigen::VectorXd nodalDisplacement(8);
    nodalDisplacement << U1, U2, U3, U4;

    //Strain vector:
    Eigen::VectorXd Strain(3); 
    Strain = Bij*nodalDisplacement;

    return Strain;
}

//Update strain rate in the element.
Eigen::VectorXd 
TIEQlin2DQuad4::ComputeStrainRate(const Eigen::MatrixXd& UNUSED(Bij)) const{
    //TODO: Compute strain rate.
    //Strain vector definition:
    Eigen::VectorXd strainrate(3);
    strainrate.fill(0.0);

    return strainrate;
}

//Computes the jacobian of the transformation. 
Eigen::MatrixXd 
TIEQlin2DQuad4::ComputeJacobianMatrix(const double ri, const double si) const{
    //Gets the element coordinates in deformed configuration. 
    Eigen::VectorXd X1 = theNodes[0]->GetCoordinates();
    Eigen::VectorXd X2 = theNodes[1]->GetCoordinates();
    Eigen::VectorXd X3 = theNodes[2]->GetCoordinates();
    Eigen::VectorXd X4 = theNodes[3]->GetCoordinates();

    //Jacobian coefficients:
    double J11 = -1.0/4.0*(1.0 - si)*X1(0) + 1.0/4.0*(1.0 - si)*X2(0) + 1.0/4.0*(1.0 + si)*X3(0) - 1.0/4.0*(1.0 + si)*X4(0);
    double J12 = -1.0/4.0*(1.0 - si)*X1(1) + 1.0/4.0*(1.0 - si)*X2(1) + 1.0/4.0*(1.0 + si)*X3(1) - 1.0/4.0*(1.0 + si)*X4(1); 
    double J21 = -1.0/4.0*(1.0 - ri)*X1(0) - 1.0/4.0*(1.0 + ri)*X2(0) + 1.0/4.0*(1.0 + ri)*X3(0) + 1.0/4.0*(1.0 - ri)*X4(0); 
    double J22 = -1.0/4.0*(1.0 - ri)*X1(1) - 1.0/4.0*(1.0 + ri)*X2(1) + 1.0/4.0*(1.0 + ri)*X3(1) + 1.0/4.0*(1.0 - ri)*X4(1); 

    //Jacobia Matrix definition:
    Eigen::MatrixXd Jij(2,2);
    Jij << J11, J12,
           J21, J22;

    return Jij;
}

//Compute Shape Function at Gauss Point:
Eigen::MatrixXd 
TIEQlin2DQuad4::ComputeShapeFunctionMatrix(const double ri, const double si) const{
    //Shape function coefficients:
    double H11 = 1.0/4.0*(1.0 - ri)*(1.0 - si);
    double H22 = 1.0/4.0*(1.0 + ri)*(1.0 - si);
    double H33 = 1.0/4.0*(1.0 + ri)*(1.0 + si);
    double H44 = 1.0/4.0*(1.0 - ri)*(1.0 + si);

    //Shape function matrix definition:
    Eigen::MatrixXd Hij(2,8);
    Hij << H11, 0.0, H22, 0.0, H33, 0.0, H44, 0.0,
           0.0, H11, 0.0, H22, 0.0, H33, 0.0, H44;

    return Hij;
}

//Evaluates the lumped-mass matrix matrix at a given Gauss point.
Eigen::MatrixXd 
TIEQlin2DQuad4::ComputeStrainDisplacementMatrix(const double ri, const double si, const Eigen::MatrixXd &Jij) const{
    //Inverse jacobian matrix:
    Eigen::MatrixXd J = Jij.inverse();

    //Strain-displacement matrix coefficients:
    double B11 = 1.0/4.0*J(0,1)*(-1.0 + ri) + 1.0/4.0*J(0,0)*(-1.0 + si); 
    double B21 = 1.0/4.0*J(0,1)*(-1.0 - ri) + 1.0/4.0*J(0,0)*( 1.0 - si); 
    double B31 = 1.0/4.0*J(0,1)*( 1.0 + ri) + 1.0/4.0*J(0,0)*( 1.0 + si); 
    double B41 = 1.0/4.0*J(0,1)*( 1.0 - ri) + 1.0/4.0*J(0,0)*(-1.0 - si);

    double B12 = 1.0/4.0*J(1,1)*(-1.0 + ri) + 1.0/4.0*J(1,0)*(-1.0 + si); 
    double B22 = 1.0/4.0*J(1,1)*(-1.0 - ri) + 1.0/4.0*J(1,0)*( 1.0 - si);
    double B32 = 1.0/4.0*J(1,1)*( 1.0 + ri) + 1.0/4.0*J(1,0)*( 1.0 + si);
    double B42 = 1.0/4.0*J(1,1)*( 1.0 - ri) + 1.0/4.0*J(1,0)*(-1.0 - si);

    //Deformation matrix definition:
    Eigen::MatrixXd Bij(3,8);
    Bij << B11, 0.0, B21, 0.0, B31, 0.0, B41, 0.0,
           0.0, B12, 0.0, B22, 0.0, B32, 0.0, B42,
           B12, B11, B22, B21, B32, B31, B42, B41;

    return Bij;
}

//Compute GGmax value from empirical equation
Eigen::VectorXd
TIEQlin2DQuad4::ComputeGGmaxDamping(const double vs, const double z, const double rho) const{

	double GGmax, xi, Dmin, F;
	double g  = 9.81;  //gravity acceleration [m/s^2]
	double K0 = 0.50;  //coefficient of lateral earth pressure
	double PIndex;     //plasticity index
	double H = zref-z; //soil column height

    if (H < 0.0)
        H = 0.0;

    double sigmav = rho*g*H; //vertical pressure
    double sigma0 = (1.0 + 2.0*K0)*sigmav/3.0;//lateral pressure
    double ppre   = 0.106*pow(vs,1.47)*1000; // [Pa]
    double OCR    = ppre/sigmav;//overconsolidation ratio

    //emprical coefficients -- based on Table 8.7 of [1]
    Eigen::VectorXd phi(7);
    phi << 3.52E-2 , 7.07E-4, 3.69E-1 , 2.97E-1 , 9.5E-1 , 6.00E-1 , 0.0;

    //based on [2]
    if (vs <= 200.0) {
        PIndex = 10.0;
    }
    else if (vs > 200.0 && vs <= 360.0){
        PIndex = 5.0;
    }
    else{
        PIndex = 0.0;
    }

    //based on [1]
    double gammar = (phi(0) + phi(1)*PIndex*pow(OCR,phi(2)))*pow(sigma0*9.86923E-6, phi(3))/100;//reference strain

    double c1 = -1.1143*phi(4)*phi(4) + 1.8618*phi(4) + 0.2523;
    double c2 =  0.0805*phi(4)*phi(4) - 0.0710*phi(4) - 0.0095;
    double c3 = -0.0005*phi(4)*phi(4) + 0.0002*phi(4) + 0.0003;

    double Dmasinga1 = 100.0/PI*(4.0*(eref - gammar*log((eref+gammar)/gammar))/(eref*eref/(eref + gammar)) - 2.0);
    double Dmasing   = c1*Dmasinga1 + c2*pow(Dmasinga1,2) + c3*pow(Dmasinga1,3);
    double b         = phi(4) + phi(5)*log(1.0);

    Dmin = 0.01;//(0.8005+0.0129*PIndex*pow(OCR,-0.1069))*pow(sigma0,-0.2889)*(1+0.2919*log(1.0));//small strain damping

    if (strcasecmp(Type.c_str(),"DARENDELI") == 0) {
        GGmax  = 1.0/(1.0 + pow(eref/gammar,phi(4)));
        F  = b*pow(GGmax,0.1);
        xi = F*Dmasing/100 + Dmin;
        if (eref < 1.0E-8) {
            GGmax = 1.0;
            xi    = Dmin;
        }
    }
    else if (strcasecmp(Type.c_str(),"SMALLSTRAIN") == 0) {
        GGmax = 1.0;
        xi    = Dmin;
    }

    Eigen::VectorXd result(2);
    result << GGmax, xi;

    return result;
}
