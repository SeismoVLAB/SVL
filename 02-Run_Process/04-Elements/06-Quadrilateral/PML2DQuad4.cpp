#include <cmath>
#include <iostream>
#include <Eigen/LU> 
#include "Material.hpp"
#include "PML2DQuad4.hpp"
#include "GaussQuadrature.hpp"
#include "LobattoQuadrature.hpp"
#include "Definitions.hpp"

//Define VTK cell value for Paraview:
const unsigned int VTKCELL = 9;

//Overload constructor.
PML2DQuad4::PML2DQuad4(const std::vector<unsigned int> nodes, std::unique_ptr<Material> &material, const std::vector<double> parameters, const std::string quadrature, const unsigned int nGauss, bool massform) :
Element("PML2DQuad4", nodes, 20, VTKCELL), t(parameters[0]), MassForm(massform), m_pml(parameters[1]), L_pml(parameters[2]), R_pml(parameters[3]), x0_pml(parameters[4]), y0_pml(parameters[5]), nx_pml(parameters[6]), ny_pml(parameters[7]) {
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
PML2DQuad4::~PML2DQuad4(){
    //Does nothing.
}

//Save the material states in the element.
void 
PML2DQuad4::CommitState(){
    //Updates the viscous material components.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    if(theMaterial[0]->IsViscous()){
        //Gets the quadrature information.    
        Eigen::VectorXd wi;
        Eigen::MatrixXd xi;
        QuadraturePoints->GetQuadraturePoints("Quad", wi, xi);

        //Update material states.
        for(unsigned int k = 0; k < wi.size(); k++){
            //Jacobian matrix.
            Eigen::MatrixXd Jij = ComputeJacobianMatrix(xi(k,0), xi(k,1));

            //Compute Strain-Displacement Matrix at Gauss Point.
            Eigen::MatrixXd Bij = ComputeStrainDisplacementMatrix(xi(k,0), xi(k,1), Jij);

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
PML2DQuad4::UpdateState(){
    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Quad", wi, xi);

    //Update material states.
    for(unsigned int k = 0; k < wi.size(); k++){
        //Jacobian matrix.
        Eigen::MatrixXd Jij = ComputeJacobianMatrix(xi(k,0), xi(k,1));

        //Compute Strain-Displacement Matrix at Gauss Point.
        Eigen::MatrixXd Bij = ComputeStrainDisplacementMatrix(xi(k,0), xi(k,1), Jij);

        //Computes strain vector.
        Eigen::VectorXd strain = ComputeStrain(Bij);

        //Update the material state.
        theMaterial[k]->UpdateState(strain, 1);
    }
}

//Sets the finite element dependance among objects.
void 
PML2DQuad4::SetDomain(std::map<unsigned int, std::shared_ptr<Node> > &nodes){
    //Gets the global element connectivity.
    std::vector<unsigned int> conn = GetNodes();

    //Assign the element to mesh node pointer.  
    for(unsigned int i = 0; i < GetNumberOfNodes(); i++){
        theNodes[i] = nodes[conn[i]];
    }
}

//Sets the damping model.
void 
PML2DQuad4::SetDamping(const std::shared_ptr<Damping> &damping){
    //The damping model
    theDamping = damping;
}

//Gets the list of total-degree of freedom of this element.
std::vector<unsigned int> 
PML2DQuad4::GetTotalDegreeOfFreedom() const{
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
PML2DQuad4::GetStrain() const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    Eigen::MatrixXd theStrain(nPoints,3);
    for(unsigned int k = 0; k < nPoints; k++)
        theStrain.row(k) = theMaterial[k]->GetStrain();

    return theStrain;
}

//Returns the material stress at integration points.
Eigen::MatrixXd 
PML2DQuad4::GetStress() const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    Eigen::MatrixXd theStress(nPoints,3);
    for(unsigned int k = 0; k < nPoints; k++)
        theStress.row(k) = theMaterial[k]->GetTotalStress();

    return theStress;
}

//Returns the material strain-rate at integration points.
Eigen::MatrixXd 
PML2DQuad4::GetStrainRate() const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    Eigen::MatrixXd theStrainRate(nPoints,3);
    for(unsigned int k = 0; k < nPoints; k++)
        theStrainRate.row(k) = theMaterial[k]->GetStrainRate();

    return theStrainRate;
}

//Gets the material strain in section at  coordinate (x3,x2).
Eigen::MatrixXd 
PML2DQuad4::GetStrainAt(double x3, double x2) const{
    UNUNSED_PARAMETER(x3);
    UNUNSED_PARAMETER(x2);

    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    //Stress at coordinate is define within section.
    Eigen::MatrixXd theStrain(nPoints, 3);
    theStrain.fill(0.0);

    return theStrain;
}

//Gets the material stress in section at  coordinate (x3,x2).
Eigen::MatrixXd 
PML2DQuad4::GetStressAt(double x3, double x2) const{
    UNUNSED_PARAMETER(x3);
    UNUNSED_PARAMETER(x2);

    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    //Stress at coordinate is define within section.
    Eigen::MatrixXd theStress(nPoints, 3);
    theStress.fill(0.0);

    return theStress;
}

//Gets the element internal response in VTK format.
Eigen::VectorXd 
PML2DQuad4::GetVTKResponse(std::string response) const{
    UNUNSED_PARAMETER(response);

    //IMPORTANT: Since PML is a buffer for absorbing waves, we decided not to show results. 
    Eigen::VectorXd theResponse(6);
    theResponse.fill(0.0);

    return theResponse;
}

//Compute the mass matrix of the element using gauss-integration.
Eigen::MatrixXd 
PML2DQuad4::ComputeMassMatrix(){
    //Mass matrix definition:
    Eigen::MatrixXd MassMatrix(20,20);
    MassMatrix.fill(0.0);

    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Quad", wi, xi);

    //Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
        //Gets material properties:
        double rho    = theMaterial[i]->GetDensity();
        double nu     = theMaterial[i]->GetPoissonRatio();
        double mu     = theMaterial[i]->GetShearModulus();
        double E      = theMaterial[i]->GetElasticityModulus();
        double lambda = E*nu/(1.0 + nu)/(1.0 - 2.0*nu);

        Eigen::VectorXd abc_pml = ComputePMLStretchingFactors(xi(i,0), xi(i,1), rho, mu, lambda);

        double a = abc_pml(0)*abc_pml(1);

        //Compute the Shape Function Coefficients.
        double H11 = 1.0/4.0*(1.0 - xi(i,0))*(1.0 - xi(i,1));
        double H22 = 1.0/4.0*(1.0 + xi(i,0))*(1.0 - xi(i,1));
        double H33 = 1.0/4.0*(1.0 + xi(i,0))*(1.0 + xi(i,1));
        double H44 = 1.0/4.0*(1.0 - xi(i,0))*(1.0 + xi(i,1));

        Eigen::VectorXd SF(4);
        SF << H11, H22, H33, H44;

        //Jacobian matrix:
        Eigen::MatrixXd Jij = ComputeJacobianMatrix(xi(i,0), xi(i,1));
        double D = fabs(Jij.determinant());

        //Computes the PML Mass Matrix at Gauss Integration point.
        for (unsigned int j = 0; j < 4 ; j++) {
            for (unsigned int k = 0; k < 4 ; k++) {
                MassMatrix(5*j  ,5*k  ) +=  rho*a*SF(j)*SF(k)*t*wi(i)*D;
                MassMatrix(5*j+1,5*k+1) +=  rho*a*SF(j)*SF(k)*t*wi(i)*D;
                MassMatrix(5*j+2,5*k+2) += -a*(lambda + 2.0*mu)/4.0/mu/(lambda + mu)*SF(j)*SF(k)*t*wi(i)*D;
                MassMatrix(5*j+3,5*k+3) += -a*(lambda + 2.0*mu)/4.0/mu/(lambda + mu)*SF(j)*SF(k)*t*wi(i)*D;
                MassMatrix(5*j+4,5*k+4) += -a/mu*SF(j)*SF(k)*t*wi(i)*D;
                MassMatrix(5*j+2,5*k+3) +=  a*lambda/4.0/mu/(lambda + mu)*SF(j)*SF(k)*t*wi(i)*D;
                MassMatrix(5*j+3,5*k+2) +=  a*lambda/4.0/mu/(lambda + mu)*SF(j)*SF(k)*t*wi(i)*D;
            }
        }
    }

    //TODO: Lumped Mass Formulation for PML.
    if(MassForm){
        //Lumped Mass in diagonal terms.
    }
    
    return MassMatrix;
}

//Compute the stiffness matrix of the element using gauss-integration.
Eigen::MatrixXd 
PML2DQuad4::ComputeStiffnessMatrix(){
    //Stiffness matrix definition:
    Eigen::MatrixXd StiffnessMatrix(20, 20);
    StiffnessMatrix.fill(0.0);

    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Quad", wi, xi);

    //Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
        //Gets material properties:
        double rho    = theMaterial[i]->GetDensity();
        double nu     = theMaterial[i]->GetPoissonRatio();
        double mu     = theMaterial[i]->GetShearModulus();
        double E      = theMaterial[i]->GetElasticityModulus();
        double lambda = E*nu/(1.0 + nu)/(1.0 - 2.0*nu);

        Eigen::VectorXd abc_pml = ComputePMLStretchingFactors(xi(i,0), xi(i,1), rho, mu, lambda);

        double bx = abc_pml(2);
        double by = abc_pml(3);
        double c  = bx*by;

        //Shape function matrix.
        double H11 = 1.0/4.0*(1.0 - xi(i,0))*(1.0 - xi(i,1));
        double H22 = 1.0/4.0*(1.0 + xi(i,0))*(1.0 - xi(i,1));
        double H33 = 1.0/4.0*(1.0 + xi(i,0))*(1.0 + xi(i,1));
        double H44 = 1.0/4.0*(1.0 - xi(i,0))*(1.0 + xi(i,1));

        Eigen::VectorXd SF(4);
        SF << H11, H22, H33, H44;

        //Strain Displacement Matrix.
        Eigen::MatrixXd Jij = ComputeJacobianMatrix(xi(i,0), xi(i,1));
        Eigen::MatrixXd J = Jij.inverse();
        double D = fabs(Jij.determinant());

        double B11 = 1.0/4.0*J(0,1)*(-1.0 + xi(i,0)) + 1.0/4.0*J(0,0)*(-1.0 + xi(i,1)); 
        double B21 = 1.0/4.0*J(0,1)*(-1.0 - xi(i,0)) + 1.0/4.0*J(0,0)*( 1.0 - xi(i,1)); 
        double B31 = 1.0/4.0*J(0,1)*( 1.0 + xi(i,0)) + 1.0/4.0*J(0,0)*( 1.0 + xi(i,1)); 
        double B41 = 1.0/4.0*J(0,1)*( 1.0 - xi(i,0)) + 1.0/4.0*J(0,0)*(-1.0 - xi(i,1));
        double B12 = 1.0/4.0*J(1,1)*(-1.0 + xi(i,0)) + 1.0/4.0*J(1,0)*(-1.0 + xi(i,1)); 
        double B22 = 1.0/4.0*J(1,1)*(-1.0 - xi(i,0)) + 1.0/4.0*J(1,0)*( 1.0 - xi(i,1));
        double B32 = 1.0/4.0*J(1,1)*( 1.0 + xi(i,0)) + 1.0/4.0*J(1,0)*( 1.0 + xi(i,1));
        double B42 = 1.0/4.0*J(1,1)*( 1.0 - xi(i,0)) + 1.0/4.0*J(1,0)*(-1.0 - xi(i,1));

        Eigen::VectorXd dSFdx(4);
        dSFdx << B11, B21, B31, B41;

        Eigen::VectorXd dSFdy(4);
        dSFdy << B12, B22, B32, B42;

        //Computes the PML Stiffness Matrix at Gauss Integration point.
        for (unsigned int j = 0; j < 4 ; j++) {
            for (unsigned int k = 0; k < 4 ; k++) {
                StiffnessMatrix(5*j  ,5*k  ) +=  rho*c*SF(j)*SF(k)*t*wi(i)*D;
                StiffnessMatrix(5*j+1,5*k+1) +=  rho*c*SF(j)*SF(k)*t*wi(i)*D;
                StiffnessMatrix(5*j+2,5*k+2) += -c*(lambda + 2.0*mu)/4.0/mu/(lambda + mu)*SF(j)*SF(k)*t*wi(i)*D;
                StiffnessMatrix(5*j+3,5*k+3) += -c*(lambda + 2.0*mu)/4.0/mu/(lambda + mu)*SF(j)*SF(k)*t*wi(i)*D;
                StiffnessMatrix(5*j+4,5*k+4) += -c/mu*SF(j)*SF(k)*t*wi(i)*D;
                StiffnessMatrix(5*j  ,5*k+2) +=  by*dSFdx(j)*SF(k)*t*wi(i)*D;
                StiffnessMatrix(5*j+2,5*k  ) +=  by*dSFdx(k)*SF(j)*t*wi(i)*D;
                StiffnessMatrix(5*j  ,5*k+4) +=  bx*dSFdy(j)*SF(k)*t*wi(i)*D;
                StiffnessMatrix(5*j+4,5*k  ) +=  bx*dSFdy(k)*SF(j)*t*wi(i)*D;
                StiffnessMatrix(5*j+1,5*k+3) +=  bx*dSFdy(j)*SF(k)*t*wi(i)*D;
                StiffnessMatrix(5*j+3,5*k+1) +=  bx*dSFdy(k)*SF(j)*t*wi(i)*D;
                StiffnessMatrix(5*j+1,5*k+4) +=  by*dSFdx(j)*SF(k)*t*wi(i)*D;
                StiffnessMatrix(5*j+4,5*k+1) +=  by*dSFdx(k)*SF(j)*t*wi(i)*D;
                StiffnessMatrix(5*j+2,5*k+3) +=  c*lambda/4.0/mu/(lambda + mu)*SF(j)*SF(k)*t*wi(i)*D;
                StiffnessMatrix(5*j+3,5*k+2) +=  c*lambda/4.0/mu/(lambda + mu)*SF(j)*SF(k)*t*wi(i)*D;
            }
        }
    }

    return StiffnessMatrix;
}

//Compute the damping matrix of the element using gauss-integration.
Eigen::MatrixXd 
PML2DQuad4::ComputeDampingMatrix(){
    //Damping matrix definition:
    Eigen::MatrixXd DampingMatrix(20, 20);
    DampingMatrix.fill(0.0);

    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Quad", wi, xi);

    //Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
        //Gets material properties:
        double rho    = theMaterial[i]->GetDensity();
        double nu     = theMaterial[i]->GetPoissonRatio();
        double mu     = theMaterial[i]->GetShearModulus();
        double E      = theMaterial[i]->GetElasticityModulus();
        double lambda = E*nu/(1.0 + nu)/(1.0 - 2.0*nu);

        Eigen::VectorXd abc_pml = ComputePMLStretchingFactors(xi(i,0), xi(i,1), rho, mu, lambda);

        double ax = abc_pml(0);
        double ay = abc_pml(1);
        double bx = abc_pml(2);
        double by = abc_pml(3);
        
        double b  = ax*by + ay*bx;

        //Shape function matrix.
        double H11 = 1.0/4.0*(1.0 - xi(i,0))*(1.0 - xi(i,1));
        double H22 = 1.0/4.0*(1.0 + xi(i,0))*(1.0 - xi(i,1));
        double H33 = 1.0/4.0*(1.0 + xi(i,0))*(1.0 + xi(i,1));
        double H44 = 1.0/4.0*(1.0 - xi(i,0))*(1.0 + xi(i,1));

        Eigen::VectorXd SF(4);
        SF << H11, H22, H33, H44;

        //Strain Displacement Matrix.
        Eigen::MatrixXd Jij = ComputeJacobianMatrix(xi(i,0), xi(i,1));
        Eigen::MatrixXd J = Jij.inverse();
        double D = fabs(Jij.determinant());

        double B11 = 1.0/4.0*J(0,1)*(-1.0 + xi(i,0)) + 1.0/4.0*J(0,0)*(-1.0 + xi(i,1)); 
        double B21 = 1.0/4.0*J(0,1)*(-1.0 - xi(i,0)) + 1.0/4.0*J(0,0)*( 1.0 - xi(i,1)); 
        double B31 = 1.0/4.0*J(0,1)*( 1.0 + xi(i,0)) + 1.0/4.0*J(0,0)*( 1.0 + xi(i,1)); 
        double B41 = 1.0/4.0*J(0,1)*( 1.0 - xi(i,0)) + 1.0/4.0*J(0,0)*(-1.0 - xi(i,1));
        double B12 = 1.0/4.0*J(1,1)*(-1.0 + xi(i,0)) + 1.0/4.0*J(1,0)*(-1.0 + xi(i,1)); 
        double B22 = 1.0/4.0*J(1,1)*(-1.0 - xi(i,0)) + 1.0/4.0*J(1,0)*( 1.0 - xi(i,1));
        double B32 = 1.0/4.0*J(1,1)*( 1.0 + xi(i,0)) + 1.0/4.0*J(1,0)*( 1.0 + xi(i,1));
        double B42 = 1.0/4.0*J(1,1)*( 1.0 - xi(i,0)) + 1.0/4.0*J(1,0)*(-1.0 - xi(i,1));

        Eigen::VectorXd dSFdx(4);
        dSFdx << B11, B21, B31, B41;

        Eigen::VectorXd dSFdy(4);
        dSFdy << B12, B22, B32, B42;

        //Computes the PML Damping Matrix at Gauss Integration point.
        for (unsigned int j = 0; j < 4 ; j++) {
            for (unsigned int k = 0; k < 4 ; k++) {
                DampingMatrix(5*j  ,5*k  ) +=  rho*b*SF(j)*SF(k)*t*wi(i)*D;
                DampingMatrix(5*j+1,5*k+1) +=  rho*b*SF(j)*SF(k)*t*wi(i)*D;
                DampingMatrix(5*j+2,5*k+2) += -b*(lambda + 2.0*mu)/4.0/mu/(lambda + mu)*SF(j)*SF(k)*t*wi(i)*D;
                DampingMatrix(5*j+3,5*k+3) += -b*(lambda + 2.0*mu)/4.0/mu/(lambda + mu)*SF(j)*SF(k)*t*wi(i)*D;
                DampingMatrix(5*j+4,5*k+4) += -b/mu*SF(j)*SF(k)*t*wi(i)*D;
                DampingMatrix(5*j  ,5*k+2) +=  ay*dSFdx(j)*SF(k)*t*wi(i)*D;
                DampingMatrix(5*j+2,5*k  ) +=  ay*dSFdx(k)*SF(j)*t*wi(i)*D;
                DampingMatrix(5*j  ,5*k+4) +=  ax*dSFdy(j)*SF(k)*t*wi(i)*D;
                DampingMatrix(5*j+4,5*k  ) +=  ax*dSFdy(k)*SF(j)*t*wi(i)*D;
                DampingMatrix(5*j+1,5*k+3) +=  ax*dSFdy(j)*SF(k)*t*wi(i)*D;
                DampingMatrix(5*j+3,5*k+1) +=  ax*dSFdy(k)*SF(j)*t*wi(i)*D;
                DampingMatrix(5*j+1,5*k+4) +=  ay*dSFdx(j)*SF(k)*t*wi(i)*D;
                DampingMatrix(5*j+4,5*k+1) +=  ay*dSFdx(k)*SF(j)*t*wi(i)*D;
                DampingMatrix(5*j+2,5*k+3) +=  b*lambda/4.0/mu/(lambda + mu)*SF(j)*SF(k)*t*wi(i)*D;
                DampingMatrix(5*j+3,5*k+2) +=  b*lambda/4.0/mu/(lambda + mu)*SF(j)*SF(k)*t*wi(i)*D;
            }
        }
    }
    
    return DampingMatrix;
}

//Compute the PML history matrix for Perfectly-Matched Layer (PML).
Eigen::MatrixXd 
PML2DQuad4::ComputePMLMatrix(){
    Eigen::MatrixXd Kpml;
    return Kpml;
}

//Compute the internal forces acting on the element.
Eigen::VectorXd 
PML2DQuad4::ComputeInternalForces(){
    //Gets the element coordinates in deformed configuration. 
    Eigen::VectorXd U1 = theNodes[0]->GetDisplacements() + theNodes[0]->GetIncrementalDisplacements();
    Eigen::VectorXd U2 = theNodes[1]->GetDisplacements() + theNodes[1]->GetIncrementalDisplacements();
    Eigen::VectorXd U3 = theNodes[2]->GetDisplacements() + theNodes[2]->GetIncrementalDisplacements();
    Eigen::VectorXd U4 = theNodes[3]->GetDisplacements() + theNodes[3]->GetIncrementalDisplacements();

    Eigen::VectorXd nodalDisplacement(20);
    nodalDisplacement << U1, U2, U3, U4;

    //Computes the Stiffness Matrix.
    Eigen::MatrixXd K = ComputeStiffnessMatrix();

    //Computes the Internal Force Vetor.
    Eigen::VectorXd InternalForces = K*nodalDisplacement;

    return InternalForces;
}

//Compute the elastic, inertial, and vicous forces acting on the element.
Eigen::VectorXd 
PML2DQuad4::ComputeInternalDynamicForces(){
    //The Internal dynamic force vector
    Eigen::VectorXd InternalForces;

    if( HasFixedNode(theNodes) ){
        //Allocate memory for velocity/acceleraton. 
        Eigen::VectorXd V(20); 
        Eigen::VectorXd A(20);

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
PML2DQuad4::ComputeSurfaceForces(const std::shared_ptr<Load> &surface, unsigned int face){
    UNUNSED_PARAMETER(face);
    UNUNSED_PARAMETER(surface);

    //PML surface forces are not supported, i.e., makes no sense.
    Eigen::VectorXd surfaceForces(20);
    surfaceForces.fill(0.0);

    return surfaceForces;
}

//Compute the body forces acting on the element.
Eigen::VectorXd 
PML2DQuad4::ComputeBodyForces(const std::shared_ptr<Load> &bodyLoad, unsigned int k){
    UNUNSED_PARAMETER(k);
    UNUNSED_PARAMETER(bodyLoad);

    //PML body forces are not supported, i.e., makes no sense.
    Eigen::VectorXd bodyForces(20);
    bodyForces.fill(0.0);

    return bodyForces;
}

//Compute the domain reduction forces acting on the element.
Eigen::VectorXd 
PML2DQuad4::ComputeDomainReductionForces(const std::shared_ptr<Load> &drm, unsigned int k){
    UNUNSED_PARAMETER(k);
    UNUNSED_PARAMETER(drm);

    //DRM method in PML domain is not supported, i.e., makes no sense.
    Eigen::VectorXd DRMForces(20);
    DRMForces.fill(0.0);

    return DRMForces;
}

//Update strain in the element.
Eigen::VectorXd 
PML2DQuad4::ComputeStrain(const Eigen::MatrixXd &Bij) const{
    //Gets the element coordinates in deformed configuration. 
    Eigen::VectorXd U1 = theNodes[0]->GetDisplacements() + theNodes[0]->GetIncrementalDisplacements();
    Eigen::VectorXd U2 = theNodes[1]->GetDisplacements() + theNodes[1]->GetIncrementalDisplacements();
    Eigen::VectorXd U3 = theNodes[2]->GetDisplacements() + theNodes[2]->GetIncrementalDisplacements();
    Eigen::VectorXd U4 = theNodes[3]->GetDisplacements() + theNodes[3]->GetIncrementalDisplacements();

    Eigen::VectorXd nodalDisplacement(8);
    nodalDisplacement << U1(0), U1(1), U2(0), U2(1), U3(0), U3(1), U4(0), U4(1);

    //Strain vector:
    Eigen::VectorXd Strain(3); 
    Strain = Bij*nodalDisplacement;

    return Strain;
}

//Update strain rate in the element.
Eigen::VectorXd 
PML2DQuad4::ComputeStrainRate(const Eigen::MatrixXd &Bij) const{
    UNUNSED_PARAMETER(Bij);

    //TODO: Compute strain rate.
    //Strain vector definition:
    Eigen::VectorXd strainrate(3);
    strainrate.fill(0.0);

    return strainrate;
}

//Computes the jacobian of the transformation. 
Eigen::MatrixXd 
PML2DQuad4::ComputeJacobianMatrix(const double ri, const double si) const{
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
PML2DQuad4::ComputeShapeFunctionMatrix(const double ri, const double si) const{
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
PML2DQuad4::ComputeStrainDisplacementMatrix(const double ri, const double si, const Eigen::MatrixXd &Jij) const{
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

//Evaluates the stretching parameters of PML
Eigen::VectorXd
PML2DQuad4::ComputePMLStretchingFactors(const double ri, const double si, const double rho, const double mu, const double lambda) const{
    //Gets the element coordinates in deformed configuration. 
    Eigen::VectorXd X1 = theNodes[0]->GetCoordinates();
    Eigen::VectorXd X2 = theNodes[1]->GetCoordinates();
    Eigen::VectorXd X3 = theNodes[2]->GetCoordinates();
    Eigen::VectorXd X4 = theNodes[3]->GetCoordinates();

    Eigen::VectorXd X(8);
    X << X1(0), X1(1), X2(0), X2(1), X3(0), X3(1), X4(0), X4(1);

    Eigen::MatrixXd Hij = ComputeShapeFunctionMatrix(ri, si);

    Eigen::VectorXd XGauss = Hij*X;

    double x = XGauss(0);
    double y = XGauss(1);

    //P wave velocity
    double V_pml = sqrt((lambda + 2.0*mu)/rho);
    
    //characteristic length
    double b_pml = L_pml/10.0;

    //stretching parameters
    double a0 = (m_pml + 1.0)*b_pml/2.0/L_pml*log(1.0/R_pml);
    double b0 = (m_pml + 1.0)*V_pml/2.0/L_pml*log(1.0/R_pml);

    double ax = 1.0 + a0*pow((x - x0_pml)*nx_pml/L_pml,m_pml);
    double ay = 1.0 + a0*pow((y - y0_pml)*ny_pml/L_pml,m_pml);

    double bx = b0*pow((x - x0_pml)*nx_pml/L_pml,m_pml);
    double by = b0*pow((y - y0_pml)*ny_pml/L_pml,m_pml);

    Eigen::VectorXd abc_pml(4);
    abc_pml << ax, ay, bx, by;

    return abc_pml;

}
