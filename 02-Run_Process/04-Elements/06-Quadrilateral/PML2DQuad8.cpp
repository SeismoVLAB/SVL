#include <cmath>
#include <Eigen/LU> 
#include "Material.hpp"
#include "PML2DQuad8.hpp"
#include "GaussQuadrature.hpp"
#include "LobattoQuadrature.hpp"
#include "Definitions.hpp"
#include "Profiler.hpp"

//Define VTK cell value for Paraview:
const unsigned int VTKCELL = 23;

//Overload constructor.
PML2DQuad8::PML2DQuad8(const std::vector<unsigned int> nodes, std::unique_ptr<Material> &material, const std::vector<double> parameters, const std::string quadrature, const unsigned int nGauss) :
Element("PML2DQuad8", nodes, 40, VTKCELL), t(parameters[0]), m_pml(parameters[1]), L_pml(parameters[2]), R_pml(parameters[3]), x0_pml(parameters[4]), y0_pml(parameters[5]), nx_pml(parameters[6]), ny_pml(parameters[7]) {
    //The element nodes.
    theNodes.resize(8);

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
PML2DQuad8::~PML2DQuad8(){
    //Does nothing.
}

//Save the material states in the element.
void 
PML2DQuad8::CommitState(){
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

//Reverse the material states to previous converged state in this element.
void 
PML2DQuad8::ReverseState(){
    //Reverse the material components.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    for(unsigned int k = 0; k < nPoints; k++)
        theMaterial[k]->ReverseState();
}

//Brings the material state to its initial state in this element.
void 
PML2DQuad8::InitialState(){
    //Brings the material components to initial state.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    for(unsigned int k = 0; k < nPoints; k++)
        theMaterial[k]->InitialState();
}

//Update the material states in the element.
void 
PML2DQuad8::UpdateState(){
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
PML2DQuad8::SetDomain(std::map<unsigned int, std::shared_ptr<Node> > &nodes){
    //Gets the global element connectivity.
    std::vector<unsigned int> conn = GetNodes();

    //Assign the element to mesh node pointer.  
    for(unsigned int i = 0; i < GetNumberOfNodes(); i++){
        theNodes[i] = nodes[conn[i]];
    }
}

//Sets the damping model.
void 
PML2DQuad8::SetDamping(const std::shared_ptr<Damping> &damping){
    //The damping model
    theDamping = damping;
}

//Gets the list of total-degree of freedom of this element.
std::vector<unsigned int> 
PML2DQuad8::GetTotalDegreeOfFreedom() const{
    //Total number of degree-of-freedom.
    unsigned int nDofs = GetNumberOfDegreeOfFreedom();

    //Reserve memory for the element list of degree-of-freedom.
    std::vector<unsigned int> dofs(nDofs);

    //Construct the element list of degree-of-freedom for assembly.
    for(unsigned int j = 0; j < 8; j++){    
        unsigned int LengthDofs = theNodes[j]->GetNumberOfDegreeOfFreedom();
        std::vector<int> totalDofs = theNodes[j]->GetTotalDegreeOfFreedom();

        for(unsigned int i = 0; i < LengthDofs; i++)
            dofs[i + LengthDofs*j] = totalDofs[i];    
    }

    return dofs;
}

//Returns the material strain at integration points.
Eigen::MatrixXd 
PML2DQuad8::GetStrain() const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    Eigen::MatrixXd theStrain(nPoints,3);
    for(unsigned int k = 0; k < nPoints; k++)
        theStrain.row(k) = theMaterial[k]->GetStrain();

    return theStrain;
}

//Returns the material stress at integration points.
Eigen::MatrixXd 
PML2DQuad8::GetStress() const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    Eigen::MatrixXd theStress(nPoints,3);
    for(unsigned int k = 0; k < nPoints; k++)
        theStress.row(k) = theMaterial[k]->GetTotalStress();

    return theStress;
}

//Returns the material strain-rate at integration points.
Eigen::MatrixXd 
PML2DQuad8::GetStrainRate() const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    Eigen::MatrixXd theStrainRate(nPoints,3);
    for(unsigned int k = 0; k < nPoints; k++)
        theStrainRate.row(k) = theMaterial[k]->GetStrainRate();

    return theStrainRate;
}

//Gets the material strain in section at  coordinate (x3,x2).
Eigen::MatrixXd 
PML2DQuad8::GetStrainAt(double UNUSED(x3), double UNUSED(x2)) const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    //Stress at coordinate is define within section.
    Eigen::MatrixXd theStrain(nPoints, 3);
    theStrain.fill(0.0);

    return theStrain;
}

//Gets the material stress in section at  coordinate (x3,x2).
Eigen::MatrixXd 
PML2DQuad8::GetStressAt(double UNUSED(x3), double UNUSED(x2)) const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    //Stress at coordinate is define within section.
    Eigen::MatrixXd theStress(nPoints, 3);
    theStress.fill(0.0);

    return theStress;
}

//Gets the element internal response in VTK format.
Eigen::VectorXd 
PML2DQuad8::GetVTKResponse(std::string UNUSED(response)) const{
    //IMPORTANT: Since PML is a buffer for absorbing waves, we decided not to show results. 
    Eigen::VectorXd theResponse(6);
    theResponse.fill(0.0);

    return theResponse;
}

//Computes the element energy for a given deformation.
double 
PML2DQuad8::ComputeEnergy(){
    //TODO: Integrate over element volume to compute the energy
    return 0.0;
}

//Compute the mass matrix of the element using gauss-integration.
Eigen::MatrixXd 
PML2DQuad8::ComputeMassMatrix(){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Use consistent mass definition:
    Eigen::MatrixXd MassMatrix(40, 40);
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

        double ri = xi(i,0);
        double si = xi(i,1);

        Eigen::VectorXd abc_pml = ComputePMLStretchingFactors(ri, si, rho, mu, lambda);

        double a = abc_pml(0)*abc_pml(1);

        //Compute the Shape Function Coefficients.
        double H11 = -1.0/4.0*(1.0 - ri)*(1.0 - si)*(1.0 + ri + si);
        double H22 = -1.0/4.0*(1.0 + ri)*(1.0 - si)*(1.0 - ri + si);
        double H33 = -1.0/4.0*(1.0 + ri)*(1.0 + si)*(1.0 - ri - si);
        double H44 = -1.0/4.0*(1.0 - ri)*(1.0 + si)*(1.0 + ri - si);
        double H55 =  1.0/2.0*(1.0 - ri)*(1.0 + ri)*(1.0 - si);
        double H66 =  1.0/2.0*(1.0 + ri)*(1.0 - si)*(1.0 + si);
        double H77 =  1.0/2.0*(1.0 - ri)*(1.0 + ri)*(1.0 + si);
        double H88 =  1.0/2.0*(1.0 - ri)*(1.0 - si)*(1.0 + si);

        Eigen::VectorXd SF(8);
        SF << H11, H22, H33, H44, H55, H66, H77, H88;

        //Jacobian matrix:
        Eigen::MatrixXd Jij = ComputeJacobianMatrix(ri, si);
        double D = fabs(Jij.determinant());

        //Computes the PML Mass Matrix at Gauss Integration point.
        for (unsigned int j = 0; j < 8 ; j++) {
            for (unsigned int k = 0; k < 8 ; k++) {
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
    if(MassFormulation){
        //Lumped Mass in diagonal terms.
    }
    
    return MassMatrix;
}

//Compute the stiffness matrix of the element using gauss-integration.
Eigen::MatrixXd 
PML2DQuad8::ComputeStiffnessMatrix(){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Stiffness matrix definition:
    Eigen::MatrixXd StiffnessMatrix(40, 40);
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

        //Gauss-Quadrature Coordinate.
        double ri = xi(i,0);
        double si = xi(i,1);

        Eigen::VectorXd abc_pml = ComputePMLStretchingFactors(ri, si, rho, mu, lambda);

        double bx = abc_pml(2);
        double by = abc_pml(3);
        double c  = bx*by;

        //Shape function matrix.
        double H11 = -1.0/4.0*(1.0 - ri)*(1.0 - si)*(1.0 + ri + si);
        double H22 = -1.0/4.0*(1.0 + ri)*(1.0 - si)*(1.0 - ri + si);
        double H33 = -1.0/4.0*(1.0 + ri)*(1.0 + si)*(1.0 - ri - si);
        double H44 = -1.0/4.0*(1.0 - ri)*(1.0 + si)*(1.0 + ri - si);
        double H55 =  1.0/2.0*(1.0 - ri)*(1.0 + ri)*(1.0 - si);
        double H66 =  1.0/2.0*(1.0 + ri)*(1.0 - si)*(1.0 + si);
        double H77 =  1.0/2.0*(1.0 - ri)*(1.0 + ri)*(1.0 + si);
        double H88 =  1.0/2.0*(1.0 - ri)*(1.0 - si)*(1.0 + si);

        Eigen::VectorXd SF(8);
        SF << H11, H22, H33, H44, H55, H66, H77, H88;

        //Computes the Strain Displacement Matrix.
        Eigen::MatrixXd Jij = ComputeJacobianMatrix(ri, si);
        Eigen::MatrixXd J = Jij.inverse();
        double D = fabs(Jij.determinant());

        double dN11 = -((2.0*ri + si)*(si - 1.0))/4.0;
        double dN12 = -((ri + 2.0*si)*(ri - 1.0))/4.0;
        double dN21 = -((2.0*ri - si)*(si - 1.0))/4.0;
        double dN22 = -((ri - 2.0*si)*(ri + 1.0))/4.0;
        double dN31 =  ((2.0*ri + si)*(si + 1.0))/4.0;
        double dN32 =  ((ri + 2.0*si)*(ri + 1.0))/4.0;
        double dN41 =  ((2.0*ri - si)*(si + 1.0))/4.0;
        double dN42 =  ((ri - 2.0*si)*(ri - 1.0))/4.0;
        double dN51 =  ri*(si - 1.0);
        double dN52 =  ri*ri/2.0 - 1.0/2.0;
        double dN61 =  1.0/2.0 - si*si/2.0;
        double dN62 = -si*(ri + 1.0);
        double dN71 = -ri*(si + 1.0);
        double dN72 =  1.0/2.0 - ri*ri/2.0;
        double dN81 =  si*si/2.0 - 1.0/2.0;
        double dN82 =  si*(ri - 1.0);

        //Strain-displacement matrix coefficients:
        double B11 = J(0,0)*dN11 + J(0,1)*dN12;
        double B12 = J(1,0)*dN11 + J(1,1)*dN12;
        double B21 = J(0,0)*dN21 + J(0,1)*dN22;
        double B22 = J(1,0)*dN21 + J(1,1)*dN22; 
        double B31 = J(0,0)*dN31 + J(0,1)*dN32;
        double B32 = J(1,0)*dN31 + J(1,1)*dN32;
        double B41 = J(0,0)*dN41 + J(0,1)*dN42;
        double B42 = J(1,0)*dN41 + J(1,1)*dN42;
        double B51 = J(0,0)*dN51 + J(0,1)*dN52;
        double B52 = J(1,0)*dN51 + J(1,1)*dN52;
        double B61 = J(0,0)*dN61 + J(0,1)*dN62;
        double B62 = J(1,0)*dN61 + J(1,1)*dN62;
        double B71 = J(0,0)*dN71 + J(0,1)*dN72;
        double B72 = J(1,0)*dN71 + J(1,1)*dN72;
        double B81 = J(0,0)*dN81 + J(0,1)*dN82;
        double B82 = J(1,0)*dN81 + J(1,1)*dN82;

        Eigen::VectorXd dSFdx(8);
        dSFdx << B11, B21, B31, B41, B51, B61, B71, B81;

        Eigen::VectorXd dSFdy(8);
        dSFdy << B12, B22, B32, B42, B52, B62, B72, B82;

        //Computes the PML Stiffness Matrix at Gauss Integration point.
        for (unsigned int j = 0; j < 8 ; j++) {
            for (unsigned int k = 0; k < 8 ; k++) {
                StiffnessMatrix(5*j  ,5*k  ) += rho*c*SF(j)*SF(k)*t*wi(i)*D;
                StiffnessMatrix(5*j+1,5*k+1) += rho*c*SF(j)*SF(k)*t*wi(i)*D;
                StiffnessMatrix(5*j+2,5*k+2) += -c*(lambda + 2.0*mu)/4.0/mu/(lambda+mu)*SF(j)*SF(k)*t*wi(i)*D;
                StiffnessMatrix(5*j+3,5*k+3) += -c*(lambda + 2.0*mu)/4.0/mu/(lambda+mu)*SF(j)*SF(k)*t*wi(i)*D;
                StiffnessMatrix(5*j+4,5*k+4) += -c/mu*SF(j)*SF(k)*t*wi(i)*D;
                StiffnessMatrix(5*j  ,5*k+2) += by*dSFdx(j)*SF(k)*t*wi(i)*D;
                StiffnessMatrix(5*j+2,5*k  ) += by*dSFdx(k)*SF(j)*t*wi(i)*D;
                StiffnessMatrix(5*j  ,5*k+4) += bx*dSFdy(j)*SF(k)*t*wi(i)*D;
                StiffnessMatrix(5*j+4,5*k  ) += bx*dSFdy(k)*SF(j)*t*wi(i)*D;
                StiffnessMatrix(5*j+1,5*k+3) += bx*dSFdy(j)*SF(k)*t*wi(i)*D;
                StiffnessMatrix(5*j+3,5*k+1) += bx*dSFdy(k)*SF(j)*t*wi(i)*D;
                StiffnessMatrix(5*j+1,5*k+4) += by*dSFdx(j)*SF(k)*t*wi(i)*D;
                StiffnessMatrix(5*j+4,5*k+1) += by*dSFdx(k)*SF(j)*t*wi(i)*D;
                StiffnessMatrix(5*j+2,5*k+3) += c*lambda/4.0/mu/(lambda + mu)*SF(j)*SF(k)*t*wi(i)*D;
                StiffnessMatrix(5*j+3,5*k+2) += c*lambda/4.0/mu/(lambda + mu)*SF(j)*SF(k)*t*wi(i)*D;
            }
        }
    }

    return StiffnessMatrix;
}

//Compute the damping matrix of the element using gauss-integration.
Eigen::MatrixXd 
PML2DQuad8::ComputeDampingMatrix(){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Damping matrix definition:
    Eigen::MatrixXd DampingMatrix(40, 40);
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

        //Gauss-Quadrature Coordinate.
        double ri = xi(i,0);
        double si = xi(i,1);

        Eigen::VectorXd abc_pml = ComputePMLStretchingFactors(ri, si, rho, mu, lambda);

        double ax = abc_pml(0);
        double ay = abc_pml(1);
        double bx = abc_pml(2);
        double by = abc_pml(3);
        
        double b  = ax*by + ay*bx;

        //Shape function matrix.
        double H11 = -1.0/4.0*(1.0 - ri)*(1.0 - si)*(1.0 + ri + si);
        double H22 = -1.0/4.0*(1.0 + ri)*(1.0 - si)*(1.0 - ri + si);
        double H33 = -1.0/4.0*(1.0 + ri)*(1.0 + si)*(1.0 - ri - si);
        double H44 = -1.0/4.0*(1.0 - ri)*(1.0 + si)*(1.0 + ri - si);
        double H55 =  1.0/2.0*(1.0 - ri)*(1.0 + ri)*(1.0 - si);
        double H66 =  1.0/2.0*(1.0 + ri)*(1.0 - si)*(1.0 + si);
        double H77 =  1.0/2.0*(1.0 - ri)*(1.0 + ri)*(1.0 + si);
        double H88 =  1.0/2.0*(1.0 - ri)*(1.0 - si)*(1.0 + si);

        Eigen::VectorXd SF(8);
        SF << H11, H22, H33, H44, H55, H66, H77, H88;

        //Strain Displacement Matrix.
        Eigen::MatrixXd Jij = ComputeJacobianMatrix(ri, si);
        Eigen::MatrixXd J = Jij.inverse();
        double D = fabs(Jij.determinant());

        double dN11 = -((2.0*ri + si)*(si - 1.0))/4.0;
        double dN12 = -((ri + 2.0*si)*(ri - 1.0))/4.0;
        double dN21 = -((2.0*ri - si)*(si - 1.0))/4.0;
        double dN22 = -((ri - 2.0*si)*(ri + 1.0))/4.0;
        double dN31 =  ((2.0*ri + si)*(si + 1.0))/4.0;
        double dN32 =  ((ri + 2.0*si)*(ri + 1.0))/4.0;
        double dN41 =  ((2.0*ri - si)*(si + 1.0))/4.0;
        double dN42 =  ((ri - 2.0*si)*(ri - 1.0))/4.0;
        double dN51 =  ri*(si - 1.0);
        double dN52 =  ri*ri/2.0 - 1.0/2.0;
        double dN61 =  1.0/2.0 - si*si/2.0;
        double dN62 = -si*(ri + 1.0);
        double dN71 = -ri*(si + 1.0);
        double dN72 =  1.0/2.0 - ri*ri/2.0;
        double dN81 =  si*si/2.0 - 1.0/2.0;
        double dN82 =  si*(ri - 1.0);

        //Strain-displacement matrix coefficients:
        double B11 = J(0,0)*dN11 + J(0,1)*dN12;
        double B12 = J(1,0)*dN11 + J(1,1)*dN12;
        double B21 = J(0,0)*dN21 + J(0,1)*dN22;
        double B22 = J(1,0)*dN21 + J(1,1)*dN22; 
        double B31 = J(0,0)*dN31 + J(0,1)*dN32;
        double B32 = J(1,0)*dN31 + J(1,1)*dN32;
        double B41 = J(0,0)*dN41 + J(0,1)*dN42;
        double B42 = J(1,0)*dN41 + J(1,1)*dN42;
        double B51 = J(0,0)*dN51 + J(0,1)*dN52;
        double B52 = J(1,0)*dN51 + J(1,1)*dN52;
        double B61 = J(0,0)*dN61 + J(0,1)*dN62;
        double B62 = J(1,0)*dN61 + J(1,1)*dN62;
        double B71 = J(0,0)*dN71 + J(0,1)*dN72;
        double B72 = J(1,0)*dN71 + J(1,1)*dN72;
        double B81 = J(0,0)*dN81 + J(0,1)*dN82;
        double B82 = J(1,0)*dN81 + J(1,1)*dN82;

        Eigen::VectorXd dSFdx(8);
        dSFdx << B11, B21, B31, B41, B51, B61, B71, B81;

        Eigen::VectorXd dSFdy(8);
        dSFdy << B12, B22, B32, B42, B52, B62, B72, B82;

        //Computes the PML Damping Matrix at Gauss Integration point.
        for (unsigned int j = 0; j < 8 ; j++) {
            for (unsigned int k = 0; k < 8 ; k++) {           
                DampingMatrix(5*j  ,5*k  ) +=  rho*b*SF(j)*SF(k)*t*wi(i)*D;
                DampingMatrix(5*j+1,5*k+1) +=  rho*b*SF(j)*SF(k)*t*wi(i)*D;
                DampingMatrix(5*j+2,5*k+2) += -b*(lambda + 2.0*mu)/4.0/mu/(lambda+mu)*SF(j)*SF(k)*t*wi(i)*D;
                DampingMatrix(5*j+3,5*k+3) += -b*(lambda + 2.0*mu)/4.0/mu/(lambda+mu)*SF(j)*SF(k)*t*wi(i)*D;
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
PML2DQuad8::ComputePMLMatrix(){
    Eigen::MatrixXd Kpml;
    return Kpml;
}

//Compute the internal forces acting on the element.
Eigen::VectorXd 
PML2DQuad8::ComputeInternalForces(){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Gets the element coordinates in deformed configuration. 
    Eigen::VectorXd U1 = theNodes[0]->GetDisplacements() + theNodes[0]->GetIncrementalDisplacements();
    Eigen::VectorXd U2 = theNodes[1]->GetDisplacements() + theNodes[1]->GetIncrementalDisplacements();
    Eigen::VectorXd U3 = theNodes[2]->GetDisplacements() + theNodes[2]->GetIncrementalDisplacements();
    Eigen::VectorXd U4 = theNodes[3]->GetDisplacements() + theNodes[3]->GetIncrementalDisplacements();
    Eigen::VectorXd U5 = theNodes[4]->GetDisplacements() + theNodes[4]->GetIncrementalDisplacements();
    Eigen::VectorXd U6 = theNodes[5]->GetDisplacements() + theNodes[5]->GetIncrementalDisplacements();
    Eigen::VectorXd U7 = theNodes[6]->GetDisplacements() + theNodes[6]->GetIncrementalDisplacements();
    Eigen::VectorXd U8 = theNodes[7]->GetDisplacements() + theNodes[7]->GetIncrementalDisplacements();

    Eigen::VectorXd nodalDisplacement(40);
    nodalDisplacement << U1, U2, U3, U4, U5, U6, U7, U8;

    //Computes the Stiffness Matrix.
    Eigen::MatrixXd K = ComputeStiffnessMatrix();

    //Computes the Internal Force Vector.
    Eigen::VectorXd InternalForces = K*nodalDisplacement;

    return InternalForces;
}

//Compute the elastic, inertial, and viscous forces acting on the element.
Eigen::VectorXd 
PML2DQuad8::ComputeInternalDynamicForces(){
    //The Internal dynamic force vector
    Eigen::VectorXd InternalForces;

    if( HasFixedNode(theNodes) ){
        //Allocate memory for velocity/acceleraton. 
        Eigen::VectorXd V(40); 
        Eigen::VectorXd A(40);

        //Fills the response vectors with velocity/acceleraton values.
        V << theNodes[0]->GetVelocities(), theNodes[1]->GetVelocities(), theNodes[2]->GetVelocities(), theNodes[3]->GetVelocities(), 
             theNodes[4]->GetVelocities(), theNodes[5]->GetVelocities(), theNodes[6]->GetVelocities(), theNodes[7]->GetVelocities();

        A << theNodes[0]->GetAccelerations(), theNodes[1]->GetAccelerations(), theNodes[2]->GetAccelerations(), theNodes[3]->GetAccelerations(), 
             theNodes[4]->GetAccelerations(), theNodes[5]->GetAccelerations(), theNodes[6]->GetAccelerations(), theNodes[7]->GetAccelerations();

        //Compute the inertial/viscous/elastic dynamic force contribution.
        InternalForces = ComputeInternalForces() + ComputeDampingMatrix()*V + ComputeMassMatrix()*A;
    }

    return InternalForces;
}

//Compute the surface forces acting on the element.
Eigen::VectorXd 
PML2DQuad8::ComputeSurfaceForces(const std::shared_ptr<Load>& UNUSED(surface), unsigned int UNUSED(face)){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //PML surface forces are not supported, i.e., makes no sense.
    Eigen::VectorXd surfaceForces(40);
    surfaceForces.fill(0.0);

    return surfaceForces;
}

//Compute the body forces acting on the element.
Eigen::VectorXd 
PML2DQuad8::ComputeBodyForces(const std::shared_ptr<Load>& UNUSED(bodyLoad), unsigned int UNUSED(k)){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //PML body forces are not supported, i.e., makes no sense.
    Eigen::VectorXd bodyForces(40);
    bodyForces.fill(0.0);

    return bodyForces;
}

//Compute the domain reduction forces acting on the element.
Eigen::VectorXd 
PML2DQuad8::ComputeDomainReductionForces(const std::shared_ptr<Load>& UNUSED(drm), unsigned int UNUSED(k)){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //DRM method in PML domain is not supported, i.e., makes no sense.
    Eigen::VectorXd DRMForces(40);
    DRMForces.fill(0.0);

    return DRMForces;
}

//Update strain in the element.
Eigen::VectorXd 
PML2DQuad8::ComputeStrain(const Eigen::MatrixXd &Bij) const{
    //Gets the element coordinates in deformed configuration. 
    Eigen::VectorXd U1 = theNodes[0]->GetDisplacements() + theNodes[0]->GetIncrementalDisplacements();
    Eigen::VectorXd U2 = theNodes[1]->GetDisplacements() + theNodes[1]->GetIncrementalDisplacements();
    Eigen::VectorXd U3 = theNodes[2]->GetDisplacements() + theNodes[2]->GetIncrementalDisplacements();
    Eigen::VectorXd U4 = theNodes[3]->GetDisplacements() + theNodes[3]->GetIncrementalDisplacements();
    Eigen::VectorXd U5 = theNodes[4]->GetDisplacements() + theNodes[4]->GetIncrementalDisplacements();
    Eigen::VectorXd U6 = theNodes[5]->GetDisplacements() + theNodes[5]->GetIncrementalDisplacements();
    Eigen::VectorXd U7 = theNodes[6]->GetDisplacements() + theNodes[6]->GetIncrementalDisplacements();
    Eigen::VectorXd U8 = theNodes[7]->GetDisplacements() + theNodes[7]->GetIncrementalDisplacements();

    Eigen::VectorXd nodalDisplacement(16);
    nodalDisplacement << U1(0), U1(1), U2(0), U2(1), U3(0), U3(1), U4(0), U4(1), U5(0), U5(1), U6(0), U6(1), U7(0), U7(1), U8(0), U8(1);

    //Strain vector:
    Eigen::VectorXd Strain(3); 
    Strain = Bij*nodalDisplacement;

    return Strain;
}

//Update strain rate in the element.
Eigen::VectorXd 
PML2DQuad8::ComputeStrainRate(const Eigen::MatrixXd& UNUSED(Bij)) const{
    //TODO: Compute strain rate.
    //Strain vector definition:
    Eigen::VectorXd strainrate(3);
    strainrate.fill(0.0);

    return strainrate;
}

//Computes the jacobian of the transformation. 
Eigen::MatrixXd 
PML2DQuad8::ComputeJacobianMatrix(const double ri, const double si) const{
    //Gets the element coordinates in deformed configuration. 
    Eigen::VectorXd X1 = theNodes[0]->GetCoordinates();
    Eigen::VectorXd X2 = theNodes[1]->GetCoordinates();
    Eigen::VectorXd X3 = theNodes[2]->GetCoordinates();
    Eigen::VectorXd X4 = theNodes[3]->GetCoordinates();
    Eigen::VectorXd X5 = theNodes[4]->GetCoordinates();
    Eigen::VectorXd X6 = theNodes[5]->GetCoordinates();
    Eigen::VectorXd X7 = theNodes[6]->GetCoordinates();
    Eigen::VectorXd X8 = theNodes[7]->GetCoordinates();

    double dN11 = -((2.0*ri + si)*(si - 1.0))/4.0;
    double dN12 = -((ri + 2.0*si)*(ri - 1.0))/4.0;
    double dN21 = -((2.0*ri - si)*(si - 1.0))/4.0;
    double dN22 = -((ri - 2.0*si)*(ri + 1.0))/4.0;
    double dN31 = ((2.0*ri + si)*(si + 1.0))/4.0;
    double dN32 = ((ri + 2.0*si)*(ri + 1.0))/4.0;
    double dN41 = ((2.0*ri - si)*(si + 1.0))/4.0;
    double dN42 = ((ri - 2.0*si)*(ri - 1.0))/4.0;
    double dN51 = ri*(si - 1.0);
    double dN52 = ri*ri/2.0 - 1.0/2.0;
    double dN61 = 1.0/2.0 - si*si/2.0;
    double dN62 = -si*(ri + 1.0);
    double dN71 = -ri*(si + 1.0);
    double dN72 = 1.0/2.0 - ri*ri/2.0;
    double dN81 = si*si/2.0 - 1.0/2.0;
    double dN82 = si*(ri - 1.0);

    //Jacobian coefficients:
    double J11 = dN11*X1(0) + dN21*X2(0) + dN31*X3(0) + dN41*X4(0) + dN51*X5(0) + dN61*X6(0) + dN71*X7(0) + dN81*X8(0);
    double J12 = dN11*X1(1) + dN21*X2(1) + dN31*X3(1) + dN41*X4(1) + dN51*X5(1) + dN61*X6(1) + dN71*X7(1) + dN81*X8(1);
    double J21 = dN12*X1(0) + dN22*X2(0) + dN32*X3(0) + dN42*X4(0) + dN52*X5(0) + dN62*X6(0) + dN72*X7(0) + dN82*X8(0);
    double J22 = dN12*X1(1) + dN22*X2(1) + dN32*X3(1) + dN42*X4(1) + dN52*X5(1) + dN62*X6(1) + dN72*X7(1) + dN82*X8(1);

    //Jacobian Matrix definition:
    Eigen::MatrixXd Jij(2,2);
    Jij << J11, J12,
           J21, J22;

    return Jij;
}

//Compute Shape Function at Gauss Point:
Eigen::MatrixXd 
PML2DQuad8::ComputeShapeFunctionMatrix(const double ri, const double si) const{
    //Shape function coefficients:
    double H11 = -1.0/4.0*(1.0 - ri)*(1.0 - si)*(1.0 + ri + si);
    double H22 = -1.0/4.0*(1.0 + ri)*(1.0 - si)*(1.0 - ri + si);
    double H33 = -1.0/4.0*(1.0 + ri)*(1.0 + si)*(1.0 - ri - si);
    double H44 = -1.0/4.0*(1.0 - ri)*(1.0 + si)*(1.0 + ri - si);
    double H55 =  1.0/2.0*(1.0 - ri)*(1.0 + ri)*(1.0 - si);
    double H66 =  1.0/2.0*(1.0 + ri)*(1.0 - si)*(1.0 + si);
    double H77 =  1.0/2.0*(1.0 - ri)*(1.0 + ri)*(1.0 + si);
    double H88 =  1.0/2.0*(1.0 - ri)*(1.0 - si)*(1.0 + si);

    //Shape function matrix definition:
    Eigen::MatrixXd Hij(2,16);
    Hij << H11, 0.0, H22, 0.0, H33, 0.0, H44, 0.0, H55, 0.0, H66, 0.0, H77, 0.0, H88, 0.0,
           0.0, H11, 0.0, H22, 0.0, H33, 0.0, H44, 0.0, H55, 0.0, H66, 0.0, H77, 0.0, H88;

    return Hij;
}

//Evaluates the lumped-mass matrix matrix at a given Gauss point.
Eigen::MatrixXd 
PML2DQuad8::ComputeStrainDisplacementMatrix(const double ri, const double si, const Eigen::MatrixXd &Jij) const{
    //Inverse jacobian matrix:
    Eigen::MatrixXd J = Jij.inverse();

    //Strain-displacement matrix coefficients:
    double dN11 = -((2.0*ri + si)*(si - 1.0))/4.0;
    double dN12 = -((ri + 2.0*si)*(ri - 1.0))/4.0;
    double dN21 = -((2.0*ri - si)*(si - 1.0))/4.0;
    double dN22 = -((ri - 2.0*si)*(ri + 1.0))/4.0;
    double dN31 =  ((2.0*ri + si)*(si + 1.0))/4.0;
    double dN32 =  ((ri + 2.0*si)*(ri + 1.0))/4.0;
    double dN41 =  ((2.0*ri - si)*(si + 1.0))/4.0;
    double dN42 =  ((ri - 2.0*si)*(ri - 1.0))/4.0;
    double dN51 =  ri*(si - 1.0);
    double dN52 =  ri*ri/2.0 - 1.0/2.0;
    double dN61 =  1.0/2.0 - si*si/2.0;
    double dN62 = -si*(ri + 1.0);
    double dN71 = -ri*(si + 1.0);
    double dN72 =  1.0/2.0 - ri*ri/2.0;
    double dN81 =  si*si/2.0 - 1.0/2.0;
    double dN82 =  si*(ri - 1.0);

    //Strain-displacement matrix coefficients:
    double B11 = J(0,0)*dN11 + J(0,1)*dN12;
    double B12 = J(1,0)*dN11 + J(1,1)*dN12;
    double B21 = J(0,0)*dN21 + J(0,1)*dN22;
    double B22 = J(1,0)*dN21 + J(1,1)*dN22; 
    double B31 = J(0,0)*dN31 + J(0,1)*dN32;
    double B32 = J(1,0)*dN31 + J(1,1)*dN32;
    double B41 = J(0,0)*dN41 + J(0,1)*dN42;
    double B42 = J(1,0)*dN41 + J(1,1)*dN42;
    double B51 = J(0,0)*dN51 + J(0,1)*dN52;
    double B52 = J(1,0)*dN51 + J(1,1)*dN52;
    double B61 = J(0,0)*dN61 + J(0,1)*dN62;
    double B62 = J(1,0)*dN61 + J(1,1)*dN62;
    double B71 = J(0,0)*dN71 + J(0,1)*dN72;
    double B72 = J(1,0)*dN71 + J(1,1)*dN72;
    double B81 = J(0,0)*dN81 + J(0,1)*dN82;
    double B82 = J(1,0)*dN81 + J(1,1)*dN82;

    //Deformation matrix definition:
    Eigen::MatrixXd Bij(3,16);
    Bij << B11, 0.0, B21, 0.0, B31, 0.0, B41, 0.0, B51, 0.0, B61, 0.0, B71, 0.0, B81, 0.0,
           0.0, B12, 0.0, B22, 0.0, B32, 0.0, B42, 0.0, B52, 0.0, B62, 0.0, B72, 0.0, B82,
           B12, B11, B22, B21, B32, B31, B42, B41, B52, B51, B62, B61, B72, B71, B82, B81;

    return Bij;
}

//Evaluates the stretching parameters of PML
Eigen::VectorXd
PML2DQuad8::ComputePMLStretchingFactors(const double ri, const double si, const double rho, const double mu, const double lambda) const {
    //Gets the element coordinates in deformed configuration. 
    Eigen::VectorXd X1 = theNodes[0]->GetCoordinates();
    Eigen::VectorXd X2 = theNodes[1]->GetCoordinates();
    Eigen::VectorXd X3 = theNodes[2]->GetCoordinates();
    Eigen::VectorXd X4 = theNodes[3]->GetCoordinates();
    Eigen::VectorXd X5 = theNodes[4]->GetCoordinates();
    Eigen::VectorXd X6 = theNodes[5]->GetCoordinates();
    Eigen::VectorXd X7 = theNodes[6]->GetCoordinates();
    Eigen::VectorXd X8 = theNodes[7]->GetCoordinates();

    Eigen::VectorXd X(16);
    X << X1(0),X1(1), X2(0), X2(1), X3(0), X3(1), X4(0), X4(1), X5(0), X5(1), X6(0), X6(1), X7(0), X7(1), X8(0), X8(1);

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
