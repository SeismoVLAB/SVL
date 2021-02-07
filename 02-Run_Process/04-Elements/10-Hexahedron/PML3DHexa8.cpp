#include <cmath>
#include <iostream>
#include <Eigen/LU> 

#include "Material.hpp"
#include "PML3DHexa8.hpp"
#include "GaussQuadrature.hpp"
#include "LobattoQuadrature.hpp"
#include "Definitions.hpp"
#include "Profiler.hpp"

//THIS IS THE VTK NUMBER FOR HEXA IN PARAVIEW
const unsigned int VTKCELL = 12;

//Overload constructor.
PML3DHexa8::PML3DHexa8(const std::vector<unsigned int> nodes, std::unique_ptr<Material> &material, const std::vector<double> parameters, const std::string quadrature, const unsigned int nGauss) :
Element("PML3DHexa8", nodes, 72, VTKCELL), m_pml(parameters[0]), L_pml(parameters[1]), R_pml(parameters[2]), x0_pml(parameters[3]), y0_pml(parameters[4]), z0_pml(parameters[5]), nx_pml(parameters[6]), ny_pml(parameters[7]), nz_pml(parameters[8]) {
    //The element nodes.
    theNodes.resize(8);

    //Numerical integration rule.
    if(strcasecmp(quadrature.c_str(),"GAUSS") == 0)
        QuadraturePoints = std::make_unique<GaussQuadrature>("Hexa", nGauss);
    else if(strcasecmp(quadrature.c_str(),"LOBATTO") == 0)
        QuadraturePoints = std::make_unique<LobattoQuadrature>("Hexa", nGauss);

	//The element material.
    theMaterial.resize(nGauss);
    for(unsigned int i = 0; i < nGauss; i++)
        theMaterial[i] = material->CopyMaterial(); 
}

//Destructor.
PML3DHexa8::~PML3DHexa8(){
    //Does nothing.
}

//Save the material states in the element.
void 
PML3DHexa8::CommitState(){
	//Updates the viscous material components.
	unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    //We are assuming the material is not viscous, right?

	for(unsigned int k = 0; k < nPoints; k++)
		theMaterial[k]->CommitState();
}

//Update the material states in the element.
void 
PML3DHexa8::UpdateState(){
	//Gets the quadrature information.
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Hexa", wi, xi);

	//Update material states.
	for(unsigned int k = 0; k < wi.size(); k++){
		//Jacobian matrix.
		Eigen::MatrixXd Jij = ComputeJacobianMatrix(xi(k,0), xi(k,1), xi(k,2));

		//Compute Strain-Displacement Matrix at Gauss Point.
		Eigen::MatrixXd Bij = ComputeStrainDisplacementMatrix(xi(k,0), xi(k,1), xi(k,2), Jij);

		//Computes strain vector.
		Eigen::VectorXd strain = ComputeStrain(Bij);

		//Update the material state.
		theMaterial[k]->UpdateState(strain, 1);
	}
}

//Sets the finite element dependance among objects.
void 
PML3DHexa8::SetDomain(std::map<unsigned int, std::shared_ptr<Node> > &nodes){
	//Gets the global element connectivity.
	std::vector<unsigned int> conn = GetNodes();

	//Assign the element to mesh node pointer.  
	for(unsigned int i = 0; i < GetNumberOfNodes(); i++){
		theNodes[i] = nodes[conn[i]];
	}
}

//Sets the damping model.
void 
PML3DHexa8::SetDamping(const std::shared_ptr<Damping> &damping){
	//The damping model
	theDamping = damping;
}

//Gets the list of total-degree of freedom of this element.
std::vector<unsigned int> 
PML3DHexa8::GetTotalDegreeOfFreedom() const{
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
PML3DHexa8::GetStrain() const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    Eigen::MatrixXd theStrain(nPoints,6);
    for(unsigned int k = 0; k < nPoints; k++)
        theStrain.row(k) = theMaterial[k]->GetStrain();

    return theStrain;
}

//Returns the material stress at integration points.
Eigen::MatrixXd 
PML3DHexa8::GetStress() const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    Eigen::MatrixXd theStress(nPoints,6);
    for(unsigned int k = 0; k < nPoints; k++)
        theStress.row(k) = theMaterial[k]->GetTotalStress();

    return theStress;
}

//Returns the material strain-rate at integration points.
Eigen::MatrixXd 
PML3DHexa8::GetStrainRate() const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    Eigen::MatrixXd theStrainRate(nPoints,6);
    for(unsigned int k = 0; k < nPoints; k++)
        theStrainRate.row(k) = theMaterial[k]->GetStrainRate();

    return theStrainRate;
}

//Gets the material strain in section at  coordinate (x3,x2).
Eigen::MatrixXd 
PML3DHexa8::GetStrainAt(double UNUSED(x3), double UNUSED(x2)) const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    //Stress at coordinate is define within section.
    Eigen::MatrixXd theStrain(nPoints, 6);
    theStrain.fill(0.0);

    return theStrain;
}

//Gets the material stress in section at  coordinate (x3,x2).
Eigen::MatrixXd 
PML3DHexa8::GetStressAt(double UNUSED(x3), double UNUSED(x2)) const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    //Stress at coordinate is define within section.
    Eigen::MatrixXd theStress(nPoints, 6);
    theStress.fill(0.0);

    return theStress;
}

//Gets the element internal response in VTK format.
Eigen::VectorXd 
PML3DHexa8::GetVTKResponse(std::string UNUSED(response)) const{
    //IMPORTANT: Since PML is a buffer for absorbing waves, we decided not to show results. 
    Eigen::VectorXd theResponse(6);
    theResponse.fill(0.0);

    return theResponse;
}

//Computes the element energy for a given deformation.
double 
PML3DHexa8::ComputeEnergy(){
    //TODO: Integrate over element volume to compute the energy
    return 0.0;
}

//Compute the Mass matrix of the element using gauss-integration.
Eigen::MatrixXd 
PML3DHexa8::ComputeMassMatrix(){
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

	//ux, uy, uz, sxx, syy, szz, sxy, syz, sxz
	Eigen::MatrixXd MassMatrix(72,72);
	MassMatrix.fill(0.0);

	//Gets the quadrature information.
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Hexa", wi, xi);

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
		double ti = xi(i,2);

		Eigen::VectorXd abc_pml = ComputePMLStretchingFactors(ri, si, ti, rho, mu, lambda);

		double a = abc_pml(0)*abc_pml(1)*abc_pml(2);

		double H11 = 1.0/8.0*(1.0 - ri)*(1.0 - si)*(1.0 - ti);
		double H22 = 1.0/8.0*(1.0 + ri)*(1.0 - si)*(1.0 - ti);
		double H33 = 1.0/8.0*(1.0 + ri)*(1.0 + si)*(1.0 - ti);
		double H44 = 1.0/8.0*(1.0 - ri)*(1.0 + si)*(1.0 - ti);
		double H55 = 1.0/8.0*(1.0 - ri)*(1.0 - si)*(1.0 + ti);
		double H66 = 1.0/8.0*(1.0 + ri)*(1.0 - si)*(1.0 + ti);
		double H77 = 1.0/8.0*(1.0 + ri)*(1.0 + si)*(1.0 + ti);
		double H88 = 1.0/8.0*(1.0 - ri)*(1.0 + si)*(1.0 + ti);

		Eigen::VectorXd SF(8);
		SF << H11, H22, H33, H44, H55, H66, H77, H88;

		//Jacobian matrix:
		Eigen::MatrixXd Jij = ComputeJacobianMatrix(ri, si, ti);
		double D = fabs(Jij.determinant());

		for (unsigned int j = 0; j < 8 ; j++) {
			for (unsigned int k = 0; k < 8 ; k++) {
				MassMatrix(9*j  ,9*k  ) += rho*a*SF(j)*SF(k)*wi(i)*D;
				MassMatrix(9*j+1,9*k+1) += rho*a*SF(j)*SF(k)*wi(i)*D;
				MassMatrix(9*j+2,9*k+2) += rho*a*SF(j)*SF(k)*wi(i)*D;
				MassMatrix(9*j+3,9*k+3) += -a*(lambda+mu)/mu/(3*lambda+2*mu)*SF(j)*SF(k)*wi(i)*D;
				MassMatrix(9*j+4,9*k+4) += -a*(lambda+mu)/mu/(3*lambda+2*mu)*SF(j)*SF(k)*wi(i)*D;
				MassMatrix(9*j+5,9*k+5) += -a*(lambda+mu)/mu/(3*lambda+2*mu)*SF(j)*SF(k)*wi(i)*D;
				MassMatrix(9*j+6,9*k+6) += -a/mu*SF(j)*SF(k)*wi(i)*D;
				MassMatrix(9*j+7,9*k+7) += -a/mu*SF(j)*SF(k)*wi(i)*D;
				MassMatrix(9*j+8,9*k+8) += -a/mu*SF(j)*SF(k)*wi(i)*D;
				
				MassMatrix(9*j+3,9*k+4) += a*lambda/2/mu/(3*lambda+2*mu)*SF(j)*SF(k)*wi(i)*D;
				MassMatrix(9*j+3,9*k+5) += a*lambda/2/mu/(3*lambda+2*mu)*SF(j)*SF(k)*wi(i)*D;
				MassMatrix(9*j+4,9*k+5) += a*lambda/2/mu/(3*lambda+2*mu)*SF(j)*SF(k)*wi(i)*D;

				MassMatrix(9*j+4,9*k+3) += a*lambda/2/mu/(3*lambda+2*mu)*SF(j)*SF(k)*wi(i)*D;
				MassMatrix(9*j+5,9*k+3) += a*lambda/2/mu/(3*lambda+2*mu)*SF(j)*SF(k)*wi(i)*D;
				MassMatrix(9*j+5,9*k+4) += a*lambda/2/mu/(3*lambda+2*mu)*SF(j)*SF(k)*wi(i)*D;
			}
		}
	}
	
	return MassMatrix;
}

//Compute the stiffness matrix of the element using gauss-integration.
Eigen::MatrixXd 
PML3DHexa8::ComputeStiffnessMatrix(){
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

	//Stiffness matrix definition:
	Eigen::MatrixXd StiffnessMatrix(72,72);
	StiffnessMatrix.fill(0.0);

	//Gets the quadrature information.
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Hexa", wi, xi);

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
		double ti = xi(i,2);

		Eigen::VectorXd abc_pml = ComputePMLStretchingFactors(ri, si, ti, rho, mu, lambda);

		double ax = abc_pml(0);
		double ay = abc_pml(1);
		double az = abc_pml(2);
		double bx = abc_pml(3);
		double by = abc_pml(4);
		double bz = abc_pml(5);

		double c  = ax*by*bz + bx*by*az + bx*ay*bz;

		double H11 = 1.0/8.0*(1.0 - ri)*(1.0 - si)*(1.0 - ti);
		double H22 = 1.0/8.0*(1.0 + ri)*(1.0 - si)*(1.0 - ti);
		double H33 = 1.0/8.0*(1.0 + ri)*(1.0 + si)*(1.0 - ti);
		double H44 = 1.0/8.0*(1.0 - ri)*(1.0 + si)*(1.0 - ti);
		double H55 = 1.0/8.0*(1.0 - ri)*(1.0 - si)*(1.0 + ti);
		double H66 = 1.0/8.0*(1.0 + ri)*(1.0 - si)*(1.0 + ti);
		double H77 = 1.0/8.0*(1.0 + ri)*(1.0 + si)*(1.0 + ti);
		double H88 = 1.0/8.0*(1.0 - ri)*(1.0 + si)*(1.0 + ti);

		//Jacobian matrix:
		Eigen::MatrixXd Jij = ComputeJacobianMatrix(ri, si, ti);
		
		double D = fabs(Jij.determinant());

		Eigen::MatrixXd J = Jij.inverse();

		//Strain-displacement matrix coefficients:
		double B11 = -1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 - ti);
		double B21 = -1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 - ti);
		double B31 = -1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 - ti);
		double B41 = -1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 - ti);
		double B51 =  1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 + ti);
		double B61 =  1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 + ti);
		double B71 =  1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 + ti);
		double B81 =  1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 + ti);

		double B12 = -1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 - ti);
		double B22 = -1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 - ti);
		double B32 = -1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 - ti);
		double B42 = -1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 - ti);
		double B52 =  1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 + ti);
		double B62 =  1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 + ti);
		double B72 =  1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 + ti);
		double B82 =  1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 + ti);

		double B13 = -1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 - ti);
		double B23 = -1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 - ti);
		double B33 = -1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 - ti);
		double B43 = -1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 - ti);
		double B53 =  1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 + ti);
		double B63 =  1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 + ti);
		double B73 =  1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 + ti);
		double B83 =  1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 + ti);

		Eigen::VectorXd SF(8);
		SF << H11, H22, H33, H44, H55, H66, H77, H88;

		Eigen::VectorXd dSFdx(8);
		dSFdx << B11, B21, B31, B41, B51, B61, B71, B81;

		Eigen::VectorXd dSFdy(8);
		dSFdy << B12, B22, B32, B42, B52, B62, B72, B82;

		Eigen::VectorXd dSFdz(8);
		dSFdz << B13, B23, B33, B43, B53, B63, B73, B83;

		for (unsigned int j = 0; j < 8 ; j++) {
			for (unsigned int k = 0; k < 8 ; k++) {
				//The state vector: U = [u1, u2, u3, s11, s22, s33, s12, s23, s13]
				StiffnessMatrix(9*j  ,9*k  ) += rho*c*SF(j)*SF(k)*wi(i)*D;
				StiffnessMatrix(9*j+1,9*k+1) += rho*c*SF(j)*SF(k)*wi(i)*D;
				StiffnessMatrix(9*j+2,9*k+2) += rho*c*SF(j)*SF(k)*wi(i)*D;
				StiffnessMatrix(9*j+3,9*k+3) += -c*(lambda + mu)/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				StiffnessMatrix(9*j+4,9*k+4) += -c*(lambda + mu)/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				StiffnessMatrix(9*j+5,9*k+5) += -c*(lambda + mu)/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				StiffnessMatrix(9*j+6,9*k+6) += -c/mu*SF(j)*SF(k)*wi(i)*D;
				StiffnessMatrix(9*j+7,9*k+7) += -c/mu*SF(j)*SF(k)*wi(i)*D;
				StiffnessMatrix(9*j+8,9*k+8) += -c/mu*SF(j)*SF(k)*wi(i)*D;
				
				StiffnessMatrix(9*j+3,9*k+4) += c*lambda/2.0/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				StiffnessMatrix(9*j+3,9*k+5) += c*lambda/2.0/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				StiffnessMatrix(9*j+4,9*k+5) += c*lambda/2.0/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;

				StiffnessMatrix(9*j+4,9*k+3) += c*lambda/2.0/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				StiffnessMatrix(9*j+5,9*k+3) += c*lambda/2.0/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				StiffnessMatrix(9*j+5,9*k+4) += c*lambda/2.0/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;

				StiffnessMatrix(9*j  ,9*k+3) += dSFdx(j)*SF(k)*(ay*bz+az*by)*wi(i)*D; 
				StiffnessMatrix(9*j  ,9*k+6) += dSFdy(j)*SF(k)*(ax*bz+az*bx)*wi(i)*D;
				StiffnessMatrix(9*j  ,9*k+8) += dSFdz(j)*SF(k)*(ax*by+ay*bx)*wi(i)*D;
				StiffnessMatrix(9*j+1,9*k+4) += dSFdy(j)*SF(k)*(ax*bz+az*bx)*wi(i)*D;
				StiffnessMatrix(9*j+1,9*k+6) += dSFdx(j)*SF(k)*(ay*bz+az*by)*wi(i)*D;
				StiffnessMatrix(9*j+1,9*k+7) += dSFdz(j)*SF(k)*(ax*by+ay*bx)*wi(i)*D;
				StiffnessMatrix(9*j+2,9*k+5) += dSFdz(j)*SF(k)*(ax*by+ay*bx)*wi(i)*D;
				StiffnessMatrix(9*j+2,9*k+8) += dSFdx(j)*SF(k)*(ay*bz+az*by)*wi(i)*D;
				StiffnessMatrix(9*j+2,9*k+7) += dSFdy(j)*SF(k)*(ax*bz+az*bx)*wi(i)*D;

				StiffnessMatrix(9*j+3,9*k  ) += dSFdx(k)*SF(j)*(ay*bz+az*by)*wi(i)*D; 
				StiffnessMatrix(9*j+6,9*k  ) += dSFdy(k)*SF(j)*(ax*bz+az*bx)*wi(i)*D;
				StiffnessMatrix(9*j+8,9*k  ) += dSFdz(k)*SF(j)*(ax*by+ay*bx)*wi(i)*D;
				StiffnessMatrix(9*j+4,9*k+1) += dSFdy(k)*SF(j)*(ax*bz+az*bx)*wi(i)*D;
				StiffnessMatrix(9*j+6,9*k+1) += dSFdx(k)*SF(j)*(ay*bz+az*by)*wi(i)*D;
				StiffnessMatrix(9*j+7,9*k+1) += dSFdz(k)*SF(j)*(ax*by+ay*bx)*wi(i)*D;
				StiffnessMatrix(9*j+5,9*k+2) += dSFdz(k)*SF(j)*(ax*by+ay*bx)*wi(i)*D;
				StiffnessMatrix(9*j+8,9*k+2) += dSFdx(k)*SF(j)*(ay*bz+az*by)*wi(i)*D;
				StiffnessMatrix(9*j+7,9*k+2) += dSFdy(k)*SF(j)*(ax*bz+az*bx)*wi(i)*D;
			}
		}
	}

	return StiffnessMatrix;
}

//Compute the Damping matrix of the element using gauss-integration.
Eigen::MatrixXd 
PML3DHexa8::ComputeDampingMatrix(){
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

	//Stiffness matrix definition:
	Eigen::MatrixXd DampingMatrix(72,72);
	DampingMatrix.fill(0.0);

	//Gets the quadrature information.
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Hexa", wi, xi);

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
		double ti = xi(i,2);

		Eigen::VectorXd abc_pml = ComputePMLStretchingFactors(ri, si, ti, rho, mu, lambda);

		double ax = abc_pml(0);
		double ay = abc_pml(1);
		double az = abc_pml(2);
		double bx = abc_pml(3);
		double by = abc_pml(4);
		double bz = abc_pml(5);

		double b  = ax*ay*bz + ax*by*az + bx*ay*az;

		double H11 = 1.0/8.0*(1.0 - ri)*(1.0 - si)*(1.0 - ti);
		double H22 = 1.0/8.0*(1.0 + ri)*(1.0 - si)*(1.0 - ti);
		double H33 = 1.0/8.0*(1.0 + ri)*(1.0 + si)*(1.0 - ti);
		double H44 = 1.0/8.0*(1.0 - ri)*(1.0 + si)*(1.0 - ti);
		double H55 = 1.0/8.0*(1.0 - ri)*(1.0 - si)*(1.0 + ti);
		double H66 = 1.0/8.0*(1.0 + ri)*(1.0 - si)*(1.0 + ti);
		double H77 = 1.0/8.0*(1.0 + ri)*(1.0 + si)*(1.0 + ti);
		double H88 = 1.0/8.0*(1.0 - ri)*(1.0 + si)*(1.0 + ti);

		//Jacobian matrix:
		Eigen::MatrixXd Jij = ComputeJacobianMatrix(ri, si, ti);
		
		double D = fabs(Jij.determinant());

		Eigen::MatrixXd J = Jij.inverse();

		//Strain-displacement matrix coefficients:
		double B11 = -1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 - ti);
		double B21 = -1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 - ti);
		double B31 = -1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 - ti);
		double B41 = -1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 - ti);
		double B51 =  1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 + ti);
		double B61 =  1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 + ti);
		double B71 =  1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 + ti);
		double B81 =  1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 + ti);

		double B12 = -1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 - ti);
		double B22 = -1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 - ti);
		double B32 = -1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 - ti);
		double B42 = -1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 - ti);
		double B52 =  1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 + ti);
		double B62 =  1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 + ti);
		double B72 =  1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 + ti);
		double B82 =  1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 + ti);

		double B13 = -1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 - ti);
		double B23 = -1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 - ti);
		double B33 = -1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 - ti);
		double B43 = -1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 - ti);
		double B53 =  1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 + ti);
		double B63 =  1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 + ti);
		double B73 =  1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 + ti);
		double B83 =  1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 + ti);

		Eigen::VectorXd SF(8);
		SF << H11, H22, H33, H44, H55, H66, H77, H88;

		Eigen::VectorXd dSFdx(8);
		dSFdx << B11, B21, B31, B41, B51, B61, B71, B81;

		Eigen::VectorXd dSFdy(8);
		dSFdy << B12, B22, B32, B42, B52, B62, B72, B82;

		Eigen::VectorXd dSFdz(8);
		dSFdz << B13, B23, B33, B43, B53, B63, B73, B83;

		for (unsigned int j = 0; j < 8 ; j++) {
			for (unsigned int k = 0; k < 8 ; k++) {
				//The state vector: U = [u1, u2, u3, s11, s22, s33, s12, s23, s13]
				DampingMatrix(9*j  ,9*k  ) += rho*b*SF(j)*SF(k)*wi(i)*D;
				DampingMatrix(9*j+1,9*k+1) += rho*b*SF(j)*SF(k)*wi(i)*D;
				DampingMatrix(9*j+2,9*k+2) += rho*b*SF(j)*SF(k)*wi(i)*D;
				DampingMatrix(9*j+3,9*k+3) += -b*(lambda + mu)/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				DampingMatrix(9*j+4,9*k+4) += -b*(lambda + mu)/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				DampingMatrix(9*j+5,9*k+5) += -b*(lambda + mu)/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				DampingMatrix(9*j+6,9*k+6) += -b/mu*SF(j)*SF(k)*wi(i)*D;
				DampingMatrix(9*j+7,9*k+7) += -b/mu*SF(j)*SF(k)*wi(i)*D;
				DampingMatrix(9*j+8,9*k+8) += -b/mu*SF(j)*SF(k)*wi(i)*D;
				
				DampingMatrix(9*j+3,9*k+4) += b*lambda/2.0/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				DampingMatrix(9*j+3,9*k+5) += b*lambda/2.0/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				DampingMatrix(9*j+4,9*k+5) += b*lambda/2.0/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;

				DampingMatrix(9*j+4,9*k+3) += b*lambda/2.0/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				DampingMatrix(9*j+5,9*k+3) += b*lambda/2.0/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				DampingMatrix(9*j+5,9*k+4) += b*lambda/2.0/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;

				DampingMatrix(9*j  ,9*k+3) += dSFdx(j)*SF(k)*(ay*az)*wi(i)*D; 
				DampingMatrix(9*j  ,9*k+6) += dSFdy(j)*SF(k)*(ax*az)*wi(i)*D;
				DampingMatrix(9*j  ,9*k+8) += dSFdz(j)*SF(k)*(ax*ay)*wi(i)*D;
				DampingMatrix(9*j+1,9*k+4) += dSFdy(j)*SF(k)*(ax*az)*wi(i)*D;
				DampingMatrix(9*j+1,9*k+6) += dSFdx(j)*SF(k)*(ay*az)*wi(i)*D;
				DampingMatrix(9*j+1,9*k+7) += dSFdz(j)*SF(k)*(ax*ay)*wi(i)*D;
				DampingMatrix(9*j+2,9*k+5) += dSFdz(j)*SF(k)*(ax*ay)*wi(i)*D;
				DampingMatrix(9*j+2,9*k+8) += dSFdx(j)*SF(k)*(ay*az)*wi(i)*D;
				DampingMatrix(9*j+2,9*k+7) += dSFdy(j)*SF(k)*(ax*az)*wi(i)*D;

				DampingMatrix(9*j+3,9*k  ) += dSFdx(k)*SF(j)*(ay*az)*wi(i)*D; 
				DampingMatrix(9*j+6,9*k  ) += dSFdy(k)*SF(j)*(ax*az)*wi(i)*D;
				DampingMatrix(9*j+8,9*k  ) += dSFdz(k)*SF(j)*(ax*ay)*wi(i)*D;
				DampingMatrix(9*j+4,9*k+1) += dSFdy(k)*SF(j)*(ax*az)*wi(i)*D;
				DampingMatrix(9*j+6,9*k+1) += dSFdx(k)*SF(j)*(ay*az)*wi(i)*D;
				DampingMatrix(9*j+7,9*k+1) += dSFdz(k)*SF(j)*(ax*ay)*wi(i)*D;
				DampingMatrix(9*j+5,9*k+2) += dSFdz(k)*SF(j)*(ax*ay)*wi(i)*D;
				DampingMatrix(9*j+8,9*k+2) += dSFdx(k)*SF(j)*(ay*az)*wi(i)*D;
				DampingMatrix(9*j+7,9*k+2) += dSFdy(k)*SF(j)*(ax*az)*wi(i)*D;
			}
		}
	}
	
	return DampingMatrix;
}

//Compute the PML history matrix for Perfectly-Matched Layer (PML).
Eigen::MatrixXd 
PML3DHexa8::ComputePMLMatrix(){
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

	//Impedance buffer matrix definition:
    Eigen::MatrixXd ImpedanceMatrix(72,72);
	ImpedanceMatrix.fill(0.0);

	//Gets the quadrature information.
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Hexa", wi, xi);

	//Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
		//Gets material properties:
        double rho    = theMaterial[i]->GetDensity();
        double nu     = theMaterial[i]->GetPoissonRatio();
        double mu     = theMaterial[i]->GetShearModulus();
        double E      = theMaterial[i]->GetElasticityModulus();
        double lambda = E*nu/(1.0 + nu)/(1.0 - 2.0*nu);

        //Gets the integration point coordinates.
		double ri = xi(i,0);
		double si = xi(i,1);
		double ti = xi(i,2);

		Eigen::VectorXd abc_pml = ComputePMLStretchingFactors(ri, si, ti, rho, mu, lambda);

		//double ax = abc_pml(0);
		//double ay = abc_pml(1);
		//double az = abc_pml(2);
		double bx = abc_pml(3);
		double by = abc_pml(4);
		double bz = abc_pml(5);

		double d  = bx*by*bz;

		double H11 = 1.0/8.0*(1.0 - ri)*(1.0 - si)*(1.0 - ti);
		double H22 = 1.0/8.0*(1.0 + ri)*(1.0 - si)*(1.0 - ti);
		double H33 = 1.0/8.0*(1.0 + ri)*(1.0 + si)*(1.0 - ti);
		double H44 = 1.0/8.0*(1.0 - ri)*(1.0 + si)*(1.0 - ti);
		double H55 = 1.0/8.0*(1.0 - ri)*(1.0 - si)*(1.0 + ti);
		double H66 = 1.0/8.0*(1.0 + ri)*(1.0 - si)*(1.0 + ti);
		double H77 = 1.0/8.0*(1.0 + ri)*(1.0 + si)*(1.0 + ti);
		double H88 = 1.0/8.0*(1.0 - ri)*(1.0 + si)*(1.0 + ti);

		//Jacobian matrix:
		Eigen::MatrixXd Jij = ComputeJacobianMatrix(ri, si, ti);
		
		double D = fabs(Jij.determinant());

		Eigen::MatrixXd J = Jij.inverse();

		//Strain-displacement matrix coefficients:
		double B11 = -1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 - ti);
		double B21 = -1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 - ti);
		double B31 = -1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 - ti);
		double B41 = -1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 - ti);
		double B51 =  1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 + ti);
		double B61 =  1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 + ti);
		double B71 =  1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 + ti);
		double B81 =  1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 + ti);

		double B12 = -1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 - ti);
		double B22 = -1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 - ti);
		double B32 = -1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 - ti);
		double B42 = -1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 - ti);
		double B52 =  1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 + ti);
		double B62 =  1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 + ti);
		double B72 =  1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 + ti);
		double B82 =  1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 + ti);

		double B13 = -1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 - ti);
		double B23 = -1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 - ti);
		double B33 = -1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 - ti);
		double B43 = -1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 - ti);
		double B53 =  1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 + ti);
		double B63 =  1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 + ti);
		double B73 =  1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 + ti);
		double B83 =  1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 + ti);

		Eigen::VectorXd SF(8);
		SF << H11, H22, H33, H44, H55, H66, H77, H88;

		Eigen::VectorXd dSFdx(8);
		dSFdx << B11, B21, B31, B41, B51, B61, B71, B81;

		Eigen::VectorXd dSFdy(8);
		dSFdy << B12, B22, B32, B42, B52, B62, B72, B82;

		Eigen::VectorXd dSFdz(8);
		dSFdz << B13, B23, B33, B43, B53, B63, B73, B83;

		for (unsigned int j = 0; j < 8 ; j++) {
			for (unsigned int k = 0; k < 8 ; k++) {
				//The state vector: U = [u1, u2, u3, s11, s22, s33, s12, s23, s13]
				ImpedanceMatrix(9*j  ,9*k  ) += rho*d*SF(j)*SF(k)*wi(i)*D;
				ImpedanceMatrix(9*j+1,9*k+1) += rho*d*SF(j)*SF(k)*wi(i)*D;
				ImpedanceMatrix(9*j+2,9*k+2) += rho*d*SF(j)*SF(k)*wi(i)*D;

				ImpedanceMatrix(9*j+3,9*k+3) += -d*(lambda + mu)/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				ImpedanceMatrix(9*j+4,9*k+4) += -d*(lambda + mu)/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				ImpedanceMatrix(9*j+5,9*k+5) += -d*(lambda + mu)/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				ImpedanceMatrix(9*j+6,9*k+6) += -d/mu*SF(j)*SF(k)*wi(i)*D;
				ImpedanceMatrix(9*j+7,9*k+7) += -d/mu*SF(j)*SF(k)*wi(i)*D;
				ImpedanceMatrix(9*j+8,9*k+8) += -d/mu*SF(j)*SF(k)*wi(i)*D;
				
				ImpedanceMatrix(9*j+3,9*k+4) += d*lambda/2.0/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				ImpedanceMatrix(9*j+3,9*k+5) += d*lambda/2.0/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				ImpedanceMatrix(9*j+4,9*k+5) += d*lambda/2.0/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;

				ImpedanceMatrix(9*j+4,9*k+3) += d*lambda/2.0/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				ImpedanceMatrix(9*j+5,9*k+3) += d*lambda/2.0/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				ImpedanceMatrix(9*j+5,9*k+4) += d*lambda/2.0/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;

				ImpedanceMatrix(9*j  ,9*k+3) += dSFdx(j)*SF(k)*(by*bz)*wi(i)*D; 
				ImpedanceMatrix(9*j  ,9*k+6) += dSFdy(j)*SF(k)*(bx*bz)*wi(i)*D;
				ImpedanceMatrix(9*j  ,9*k+8) += dSFdz(j)*SF(k)*(bx*by)*wi(i)*D;
				ImpedanceMatrix(9*j+1,9*k+4) += dSFdy(j)*SF(k)*(bx*bz)*wi(i)*D;
				ImpedanceMatrix(9*j+1,9*k+6) += dSFdx(j)*SF(k)*(by*bz)*wi(i)*D;
				ImpedanceMatrix(9*j+1,9*k+7) += dSFdz(j)*SF(k)*(bx*by)*wi(i)*D;
				ImpedanceMatrix(9*j+2,9*k+5) += dSFdz(j)*SF(k)*(bx*by)*wi(i)*D;
				ImpedanceMatrix(9*j+2,9*k+8) += dSFdx(j)*SF(k)*(by*bz)*wi(i)*D;
				ImpedanceMatrix(9*j+2,9*k+7) += dSFdy(j)*SF(k)*(bx*bz)*wi(i)*D;

				ImpedanceMatrix(9*j+3,9*k  ) += dSFdx(k)*SF(j)*(by*bz)*wi(i)*D; 
				ImpedanceMatrix(9*j+6,9*k  ) += dSFdy(k)*SF(j)*(bx*bz)*wi(i)*D;
				ImpedanceMatrix(9*j+8,9*k  ) += dSFdz(k)*SF(j)*(bx*by)*wi(i)*D;
				ImpedanceMatrix(9*j+4,9*k+1) += dSFdy(k)*SF(j)*(bx*bz)*wi(i)*D;
				ImpedanceMatrix(9*j+6,9*k+1) += dSFdx(k)*SF(j)*(by*bz)*wi(i)*D;
				ImpedanceMatrix(9*j+7,9*k+1) += dSFdz(k)*SF(j)*(bx*by)*wi(i)*D;
				ImpedanceMatrix(9*j+5,9*k+2) += dSFdz(k)*SF(j)*(bx*by)*wi(i)*D;
				ImpedanceMatrix(9*j+8,9*k+2) += dSFdx(k)*SF(j)*(by*bz)*wi(i)*D;
				ImpedanceMatrix(9*j+7,9*k+2) += dSFdy(k)*SF(j)*(bx*bz)*wi(i)*D;
			}
		}
	}

    return ImpedanceMatrix;
}

//Compute the internal forces acting on the element.
Eigen::VectorXd 
PML3DHexa8::ComputeInternalForces(){
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

	//Stiffness matrix definition:
	Eigen::VectorXd InternalForces(72);
	InternalForces.fill(0.0) ;

	//Gets the element coordinates in deformed configuration. 
	Eigen::VectorXd U1 = theNodes[0]->GetDisplacements() + theNodes[0]->GetIncrementalDisplacements();
	Eigen::VectorXd U2 = theNodes[1]->GetDisplacements() + theNodes[1]->GetIncrementalDisplacements();
	Eigen::VectorXd U3 = theNodes[2]->GetDisplacements() + theNodes[2]->GetIncrementalDisplacements();
	Eigen::VectorXd U4 = theNodes[3]->GetDisplacements() + theNodes[3]->GetIncrementalDisplacements();
	Eigen::VectorXd U5 = theNodes[4]->GetDisplacements() + theNodes[4]->GetIncrementalDisplacements();
	Eigen::VectorXd U6 = theNodes[5]->GetDisplacements() + theNodes[5]->GetIncrementalDisplacements();
	Eigen::VectorXd U7 = theNodes[6]->GetDisplacements() + theNodes[6]->GetIncrementalDisplacements();
	Eigen::VectorXd U8 = theNodes[7]->GetDisplacements() + theNodes[7]->GetIncrementalDisplacements();

	Eigen::VectorXd nodalDisplacement(72);
	nodalDisplacement << U1, U2, U3, U4, U5, U6, U7, U8;

	Eigen::MatrixXd K = ComputeStiffnessMatrix();

	InternalForces = K*nodalDisplacement;

	return InternalForces;
}

//Compute the elastic, inertial, and vicous forces acting on the element.
Eigen::VectorXd 
PML3DHexa8::ComputeInternalDynamicForces(){
    //The Internal dynamic force vector
    Eigen::VectorXd InternalForces;

    if( HasFixedNode(theNodes) ){
        //Allocate memory for velocity/acceleraton. 
        Eigen::VectorXd V(72); 
        Eigen::VectorXd A(72);

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

//Compute the body forces acting on the element.
Eigen::VectorXd 
PML3DHexa8::ComputeBodyForces(const std::shared_ptr<Load>& UNUSED(bodyLoad), unsigned int UNUSED(k)){
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

	//Local body load vector:
	Eigen::VectorXd bodyForces(72);
	bodyForces.fill(0.0);

	return bodyForces;
}

//Compute the surface forces acting on the element.
Eigen::VectorXd 
PML3DHexa8::ComputeSurfaceForces(const std::shared_ptr<Load>& UNUSED(surface), unsigned int UNUSED(face)){
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

    //PML surface forces are not supported, i.e., makes no sense.
    Eigen::VectorXd surfaceForces(72);
    surfaceForces.fill(0.0);

    return surfaceForces;
}

//Compute the domain reduction forces acting on the element.
Eigen::VectorXd 
PML3DHexa8::ComputeDomainReductionForces(const std::shared_ptr<Load>& UNUSED(drm), unsigned int UNUSED(k)){
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

    //PML is not supposed to generate DRM forces.
	Eigen::VectorXd DRMForces(72);
	DRMForces.fill(0.0);

	return DRMForces;
}

//Update strain in the element.
Eigen::VectorXd 
PML3DHexa8::ComputeStrain(const Eigen::MatrixXd &Bij) const{
	//Gets the element coordinates in deformed configuration. 
	Eigen::VectorXd U1 = theNodes[0]->GetDisplacements() + theNodes[0]->GetIncrementalDisplacements();
	Eigen::VectorXd U2 = theNodes[1]->GetDisplacements() + theNodes[1]->GetIncrementalDisplacements();
	Eigen::VectorXd U3 = theNodes[2]->GetDisplacements() + theNodes[2]->GetIncrementalDisplacements();
	Eigen::VectorXd U4 = theNodes[3]->GetDisplacements() + theNodes[3]->GetIncrementalDisplacements();
	Eigen::VectorXd U5 = theNodes[4]->GetDisplacements() + theNodes[4]->GetIncrementalDisplacements();
	Eigen::VectorXd U6 = theNodes[5]->GetDisplacements() + theNodes[5]->GetIncrementalDisplacements();
	Eigen::VectorXd U7 = theNodes[6]->GetDisplacements() + theNodes[6]->GetIncrementalDisplacements();
	Eigen::VectorXd U8 = theNodes[7]->GetDisplacements() + theNodes[7]->GetIncrementalDisplacements();

	Eigen::VectorXd nodalDisplacement(24);
	nodalDisplacement << U1(0), U1(1), U1(2),
						 U2(0), U2(1), U2(2),
						 U3(0), U3(1), U3(2),
						 U4(0), U4(1), U4(2),
						 U5(0), U5(1), U5(2),
						 U6(0), U6(1), U6(2),
						 U7(0), U7(1), U7(2),
						 U8(0), U8(1), U8(2);

    //Strain vector:
	Eigen::VectorXd Strain(6); 
    Strain = Bij*nodalDisplacement;

	return Strain;
}

//Update strain rate in the element.
Eigen::VectorXd 
PML3DHexa8::ComputeStrainRate(const Eigen::MatrixXd& UNUSED(Bij)) const{
	//TODO: Compute strain rate.
	//Strain vector definition:
    Eigen::VectorXd strainrate(6);
	strainrate.fill(0.0);

	return strainrate;
}

//Computes the jacobian of the transformation. 
Eigen::MatrixXd 
PML3DHexa8::ComputeJacobianMatrix(const double ri, const double si, const double ti) const{
	//Gets the element coordinates in deformed configuration. 
	Eigen::VectorXd X1 = theNodes[0]->GetCoordinates();
	Eigen::VectorXd X2 = theNodes[1]->GetCoordinates();
	Eigen::VectorXd X3 = theNodes[2]->GetCoordinates();
	Eigen::VectorXd X4 = theNodes[3]->GetCoordinates();
	Eigen::VectorXd X5 = theNodes[4]->GetCoordinates();
	Eigen::VectorXd X6 = theNodes[5]->GetCoordinates();
	Eigen::VectorXd X7 = theNodes[6]->GetCoordinates();
	Eigen::VectorXd X8 = theNodes[7]->GetCoordinates();

	//Jacobian coefficients:
	double J11 = -1.0/8.0*(1.0 - si)*(1.0 - ti)*X1(0) + 1.0/8.0*(1.0 - si)*(1.0 - ti)*X2(0) + 1.0/8.0*(1.0 + si)*(1.0 - ti)*X3(0) - 1.0/8.0*(1.0 + si)*(1.0 - ti)*X4(0) - 1.0/8.0*(1.0 - si)*(1.0 + ti)*X5(0) + 1.0/8.0*(1.0 - si)*(1.0 + ti)*X6(0) + 1.0/8.0*(1.0 + si)*(1.0 + ti)*X7(0) - 1.0/8.0*(1.0 + si)*(1.0 + ti)*X8(0);
	double J12 = -1.0/8.0*(1.0 - si)*(1.0 - ti)*X1(1) + 1.0/8.0*(1.0 - si)*(1.0 - ti)*X2(1) + 1.0/8.0*(1.0 + si)*(1.0 - ti)*X3(1) - 1.0/8.0*(1.0 + si)*(1.0 - ti)*X4(1) - 1.0/8.0*(1.0 - si)*(1.0 + ti)*X5(1) + 1.0/8.0*(1.0 - si)*(1.0 + ti)*X6(1) + 1.0/8.0*(1.0 + si)*(1.0 + ti)*X7(1) - 1.0/8.0*(1.0 + si)*(1.0 + ti)*X8(1);
	double J13 = -1.0/8.0*(1.0 - si)*(1.0 - ti)*X1(2) + 1.0/8.0*(1.0 - si)*(1.0 - ti)*X2(2) + 1.0/8.0*(1.0 + si)*(1.0 - ti)*X3(2) - 1.0/8.0*(1.0 + si)*(1.0 - ti)*X4(2) - 1.0/8.0*(1.0 - si)*(1.0 + ti)*X5(2) + 1.0/8.0*(1.0 - si)*(1.0 + ti)*X6(2) + 1.0/8.0*(1.0 + si)*(1.0 + ti)*X7(2) - 1.0/8.0*(1.0 + si)*(1.0 + ti)*X8(2); 
	double J21 = -1.0/8.0*(1.0 - ri)*(1.0 - ti)*X1(0) - 1.0/8.0*(1.0 + ri)*(1.0 - ti)*X2(0) + 1.0/8.0*(1.0 + ri)*(1.0 - ti)*X3(0) + 1.0/8.0*(1.0 - ri)*(1.0 - ti)*X4(0) - 1.0/8.0*(1.0 - ri)*(1.0 + ti)*X5(0) - 1.0/8.0*(1.0 + ri)*(1.0 + ti)*X6(0) + 1.0/8.0*(1.0 + ri)*(1.0 + ti)*X7(0) + 1.0/8.0*(1.0 - ri)*(1.0 + ti)*X8(0);
	double J22 = -1.0/8.0*(1.0 - ri)*(1.0 - ti)*X1(1) - 1.0/8.0*(1.0 + ri)*(1.0 - ti)*X2(1) + 1.0/8.0*(1.0 + ri)*(1.0 - ti)*X3(1) + 1.0/8.0*(1.0 - ri)*(1.0 - ti)*X4(1) - 1.0/8.0*(1.0 - ri)*(1.0 + ti)*X5(1) - 1.0/8.0*(1.0 + ri)*(1.0 + ti)*X6(1) + 1.0/8.0*(1.0 + ri)*(1.0 + ti)*X7(1) + 1.0/8.0*(1.0 - ri)*(1.0 + ti)*X8(1);
	double J23 = -1.0/8.0*(1.0 - ri)*(1.0 - ti)*X1(2) - 1.0/8.0*(1.0 + ri)*(1.0 - ti)*X2(2) + 1.0/8.0*(1.0 + ri)*(1.0 - ti)*X3(2) + 1.0/8.0*(1.0 - ri)*(1.0 - ti)*X4(2) - 1.0/8.0*(1.0 - ri)*(1.0 + ti)*X5(2) - 1.0/8.0*(1.0 + ri)*(1.0 + ti)*X6(2) + 1.0/8.0*(1.0 + ri)*(1.0 + ti)*X7(2) + 1.0/8.0*(1.0 - ri)*(1.0 + ti)*X8(2);
	double J31 = -1.0/8.0*(1.0 - ri)*(1.0 - si)*X1(0) - 1.0/8.0*(1.0 + ri)*(1.0 - si)*X2(0) - 1.0/8.0*(1.0 + ri)*(1.0 + si)*X3(0) - 1.0/8.0*(1.0 - ri)*(1.0 + si)*X4(0) + 1.0/8.0*(1.0 - ri)*(1.0 - si)*X5(0) + 1.0/8.0*(1.0 + ri)*(1.0 - si)*X6(0) + 1.0/8.0*(1.0 + ri)*(1.0 + si)*X7(0) + 1.0/8.0*(1.0 - ri)*(1.0 + si)*X8(0);
	double J32 = -1.0/8.0*(1.0 - ri)*(1.0 - si)*X1(1) - 1.0/8.0*(1.0 + ri)*(1.0 - si)*X2(1) - 1.0/8.0*(1.0 + ri)*(1.0 + si)*X3(1) - 1.0/8.0*(1.0 - ri)*(1.0 + si)*X4(1) + 1.0/8.0*(1.0 - ri)*(1.0 - si)*X5(1) + 1.0/8.0*(1.0 + ri)*(1.0 - si)*X6(1) + 1.0/8.0*(1.0 + ri)*(1.0 + si)*X7(1) + 1.0/8.0*(1.0 - ri)*(1.0 + si)*X8(1);
	double J33 = -1.0/8.0*(1.0 - ri)*(1.0 - si)*X1(2) - 1.0/8.0*(1.0 + ri)*(1.0 - si)*X2(2) - 1.0/8.0*(1.0 + ri)*(1.0 + si)*X3(2) - 1.0/8.0*(1.0 - ri)*(1.0 + si)*X4(2) + 1.0/8.0*(1.0 - ri)*(1.0 - si)*X5(2) + 1.0/8.0*(1.0 + ri)*(1.0 - si)*X6(2) + 1.0/8.0*(1.0 + ri)*(1.0 + si)*X7(2) + 1.0/8.0*(1.0 - ri)*(1.0 + si)*X8(2);

	//Jacobia Matrix definition:
	Eigen::MatrixXd Jij(3,3);
	Jij << J11, J12, J13,
		   J21, J22, J23,
		   J31, J32, J33;

	return Jij;
}

//Compute Shape Function at Gauss Point:
Eigen::MatrixXd 
PML3DHexa8::ComputeShapeFunctionMatrix(const double ri, const double si, const double ti) const{
	//Shape function coefficients:
	double H11 = 1.0/8.0*(1.0 - ri)*(1.0 - si)*(1.0 - ti);
	double H22 = 1.0/8.0*(1.0 + ri)*(1.0 - si)*(1.0 - ti);
	double H33 = 1.0/8.0*(1.0 + ri)*(1.0 + si)*(1.0 - ti);
	double H44 = 1.0/8.0*(1.0 - ri)*(1.0 + si)*(1.0 - ti);
	double H55 = 1.0/8.0*(1.0 - ri)*(1.0 - si)*(1.0 + ti);
	double H66 = 1.0/8.0*(1.0 + ri)*(1.0 - si)*(1.0 + ti);
	double H77 = 1.0/8.0*(1.0 + ri)*(1.0 + si)*(1.0 + ti);
	double H88 = 1.0/8.0*(1.0 - ri)*(1.0 + si)*(1.0 + ti);

	//Shape function matrix:
	Eigen::MatrixXd Hij(3,24);
	Hij << H11, 0.0, 0.0, H22, 0.0, 0.0, H33, 0.0, 0.0, H44, 0.0, 0.0, H55, 0.0, 0.0, H66, 0.0, 0.0, H77, 0.0, 0.0, H88, 0.0, 0.0,
		   0.0, H11, 0.0, 0.0, H22, 0.0, 0.0, H33, 0.0, 0.0, H44, 0.0, 0.0, H55, 0.0, 0.0, H66, 0.0, 0.0, H77, 0.0, 0.0, H88, 0.0,
		   0.0, 0.0, H11, 0.0, 0.0, H22, 0.0, 0.0, H33, 0.0, 0.0, H44, 0.0, 0.0, H55, 0.0, 0.0, H66, 0.0, 0.0, H77, 0.0, 0.0, H88;

	return Hij;
}

//Evaluates the lumped-mass matrix matrix at a given Gauss point.
Eigen::MatrixXd 
PML3DHexa8::ComputeStrainDisplacementMatrix(const double ri, const double si, const double ti, const Eigen::MatrixXd &Jij) const{
	//Inverse jacobian matrix:
	Eigen::MatrixXd J = Jij.inverse();

	//Strain-displacement matrix coefficients:
	double B11 = -1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 - ti);
	double B21 = -1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 - ti);
	double B31 = -1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 - ti);
	double B41 = -1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 - ti);
	double B51 =  1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 + ti);
	double B61 =  1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 + ti);
	double B71 =  1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 + ti);
	double B81 =  1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 + ti);

	double B12 = -1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 - ti);
	double B22 = -1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 - ti);
	double B32 = -1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 - ti);
	double B42 = -1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 - ti);
	double B52 =  1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 + ti);
	double B62 =  1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 + ti);
	double B72 =  1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 + ti);
	double B82 =  1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 + ti);

	double B13 = -1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 - ti);
	double B23 = -1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 - ti);
	double B33 = -1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 - ti);
	double B43 = -1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 - ti);
	double B53 =  1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 + ti);
	double B63 =  1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 + ti);
	double B73 =  1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 + ti);
	double B83 =  1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 + ti);

	//Deformation matrix definition:
	Eigen::MatrixXd Bij(6,24);
	Bij <<  B11, 0.0, 0.0, B21, 0.0, 0.0, B31, 0.0, 0.0, B41, 0.0, 0.0, B51, 0.0, 0.0, B61, 0.0, 0.0, B71, 0.0, 0.0, B81, 0.0, 0.0,
		    0.0, B12, 0.0, 0.0, B22, 0.0, 0.0, B32, 0.0, 0.0, B42, 0.0, 0.0, B52, 0.0, 0.0, B62, 0.0, 0.0, B72, 0.0, 0.0, B82, 0.0,
		    0.0, 0.0, B13, 0.0, 0.0, B23, 0.0, 0.0, B33, 0.0, 0.0, B43, 0.0, 0.0, B53, 0.0, 0.0, B63, 0.0, 0.0, B73, 0.0, 0.0, B83,
		    B12, B11, 0.0, B22, B21, 0.0, B32, B31, 0.0, B42, B41, 0.0, B52, B51, 0.0, B62, B61, 0.0, B72, B71, 0.0, B82, B81, 0.0,
		    0.0, B13, B12, 0.0, B23, B22, 0.0, B33, B32, 0.0, B43, B42, 0.0, B53, B52, 0.0, B63, B62, 0.0, B73, B72, 0.0, B83, B82,
		    B13, 0.0, B11, B23, 0.0, B21, B33, 0.0, B31, B43, 0.0, B41, B53, 0.0, B51, B63, 0.0, B61, B73, 0.0, B71, B83, 0.0, B81;

	return Bij;
}

//Evaluates the stretching parameters of PML
Eigen::VectorXd
PML3DHexa8::ComputePMLStretchingFactors(const double ri, const double si, const double ti, const double rho, const double mu, const double lambda) const {
	//Gets the element coordinates in deformed configuration. 
	Eigen::VectorXd X1 = theNodes[0]->GetCoordinates();
	Eigen::VectorXd X2 = theNodes[1]->GetCoordinates();
	Eigen::VectorXd X3 = theNodes[2]->GetCoordinates();
	Eigen::VectorXd X4 = theNodes[3]->GetCoordinates();
	Eigen::VectorXd X5 = theNodes[4]->GetCoordinates();
	Eigen::VectorXd X6 = theNodes[5]->GetCoordinates();
	Eigen::VectorXd X7 = theNodes[6]->GetCoordinates();
	Eigen::VectorXd X8 = theNodes[7]->GetCoordinates();

	Eigen::VectorXd X(24);
	X << X1, X2, X3, X4, X5, X6, X7, X8;

	Eigen::MatrixXd Hij = ComputeShapeFunctionMatrix(ri, si, ti);

	Eigen::VectorXd XGauss = Hij*X;

	double x = XGauss(0);
	double y = XGauss(1);
	double z = XGauss(2);

	//P wave velocity
	double V_pml = sqrt((lambda + 2.0*mu)/rho);
	
	//characteristic length
	double b_pml = L_pml/10.0;

	//stretching parameters
	double a0 = (m_pml+1.0)*b_pml/2.0/L_pml*log(1.0/R_pml);
	double b0 = (m_pml+1.0)*V_pml/2.0/L_pml*log(1.0/R_pml);

	double ax = 1.0+a0*pow((x-x0_pml)*nx_pml/L_pml,m_pml);
	double ay = 1.0+a0*pow((y-y0_pml)*ny_pml/L_pml,m_pml);
	double az = 1.0+a0*pow((z-z0_pml)*nz_pml/L_pml,m_pml);

	double bx = b0*pow((x-x0_pml)*nx_pml/L_pml,m_pml);
	double by = b0*pow((y-y0_pml)*ny_pml/L_pml,m_pml);
	double bz = b0*pow((z-z0_pml)*nz_pml/L_pml,m_pml);

	Eigen::VectorXd abc_pml(6);
	abc_pml << ax, ay, az, bx, by, bz;

	return abc_pml;
}
