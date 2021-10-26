#include <cmath>
#include <Eigen/LU> 
#include "Material.hpp"
#include "kin2DQuad4.hpp"
#include "GaussQuadrature.hpp"
#include "LobattoQuadrature.hpp"
#include "Definitions.hpp"
#include "Profiler.hpp"

//Overload constructor.
kin2DQuad4::kin2DQuad4(const std::vector<unsigned int> nodes, std::unique_ptr<Material> &material, const double th, const std::string quadrature, const unsigned int nGauss) :
Element("kin2DQuad4", nodes, 8, VTK_LINEAR_QUAD, GROUP_ELEMENT_QUAD), t(th){
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
kin2DQuad4::~kin2DQuad4(){
    //Does nothing.
}

//Save the material states in the element.
void 
kin2DQuad4::CommitState(){
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

            //Compute Linear Strain-Displacement Matrix at Gauss Point.
            Eigen::MatrixXd Bij = ComputeLinearStrainDisplacementMatrix(xi(k,0), xi(k,1), Jij);

            //Computes Strain Rate vector.
            Eigen::VectorXd strainrate = ComputeStrainRate(xi(k,0), xi(k,1), Bij);

            //Update the material state.
            theMaterial[k]->UpdateState(strainrate, 2);
        }
    }

    for(unsigned int k = 0; k < nPoints; k++)
        theMaterial[k]->CommitState();
}

//Reverse the material states to previous converged state in this element.
void 
kin2DQuad4::ReverseState(){
    //Reverse the material components.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    for(unsigned int k = 0; k < nPoints; k++)
        theMaterial[k]->ReverseState();
}

//Brings the material state to its initial state in this element.
void 
kin2DQuad4::InitialState(){
    //Brings the material components to initial state.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    for(unsigned int k = 0; k < nPoints; k++)
        theMaterial[k]->InitialState();
}

//Update the material states in the element.
void 
kin2DQuad4::UpdateState(){
    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Quad", wi, xi);

    //Update material states.
    for(unsigned int k = 0; k < wi.size(); k++){
        //Jacobian matrix.
        Eigen::MatrixXd Jij = ComputeJacobianMatrix(xi(k,0), xi(k,1));

        //Computes strain vector.
        Eigen::VectorXd strain = ComputeStrain(xi(k,0), xi(k,1), Jij);

        //Update the material state.
        theMaterial[k]->UpdateState(strain, 1);
    }
}

//Sets the finite element dependance among objects.
void 
kin2DQuad4::SetDomain(std::map<unsigned int, std::shared_ptr<Node> > &nodes){
    //Gets the global element connectivity.
    std::vector<unsigned int> conn = GetNodes();

    //Assign the element to mesh node pointer.  
    for(unsigned int i = 0; i < GetNumberOfNodes(); i++){
        theNodes[i] = nodes[conn[i]];
    }
}

//Sets the damping model.
void 
kin2DQuad4::SetDamping(const std::shared_ptr<Damping> &damping){
    //The damping model
    theDamping = damping;
}

//Gets the list of total-degree of freedom of this element.
std::vector<unsigned int> 
kin2DQuad4::GetTotalDegreeOfFreedom() const{
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
kin2DQuad4::GetStrain() const{
    //Number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    //Strain values at each integration point.
    Eigen::MatrixXd theStrain(nPoints,3);
    for(unsigned int k = 0; k < nPoints; k++)
        theStrain.row(k) = theMaterial[k]->GetStrain();

    return theStrain;
}

//Returns the material stress at integration points.
Eigen::MatrixXd 
kin2DQuad4::GetStress() const{
    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Quad", wi, xi);

    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    //Stress values at each integration point.
    Eigen::MatrixXd theStress(nPoints,3);
    for(unsigned int k = 0; k < nPoints; k++){
        //Jacobian matrix.
        Eigen::MatrixXd Jij = ComputeJacobianMatrix(xi(k,0), xi(k,1));

        //Deformation gradient matrix.
        Eigen::MatrixXd Fij = ComputeDeformationGradientMatrix(xi(k,0), xi(k,1), Jij);

        //Second Piola-Kirchhoff Stress vector.
        Eigen::VectorXd Stress = theMaterial[k]->GetTotalStress();
        Eigen::MatrixXd Sij = TransformVectorToTensor(Stress);

        Sij = 1.0/Fij.determinant()*Fij*Sij*Fij.transpose();

        //Cauchy Stress Tensor.
        theStress.row(k) = TransformTensorToVector(Sij); 
    }

    return theStress;
}

//Returns the material strain-rate at integration points.
Eigen::MatrixXd 
kin2DQuad4::GetStrainRate() const{
    //TODO: Compute strain rate using large-deformation.
    //Total number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    Eigen::MatrixXd theStrainRate(nPoints,3);
    for(unsigned int k = 0; k < nPoints; k++)
        theStrainRate.row(k) = theMaterial[k]->GetStrainRate();

    return theStrainRate;
}

//Gets the material strain in section at  coordinate (x3,x2).
Eigen::MatrixXd 
kin2DQuad4::GetStrainAt(double UNUSED(x3), double UNUSED(x2)) const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    //Stress at coordinate is define within section.
    Eigen::MatrixXd theStrain(nPoints, 3);
    theStrain.fill(0.0);

    return theStrain;
}

//Gets the material stress in section at  coordinate (x3,x2).
Eigen::MatrixXd 
kin2DQuad4::GetStressAt(double UNUSED(x3), double UNUSED(x2)) const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    //Stress at coordinate is define within section.
    Eigen::MatrixXd theStress(nPoints, 3);
    theStress.fill(0.0);

    return theStress;
}

//Gets the element internal response in VTK format.
Eigen::VectorXd 
kin2DQuad4::GetVTKResponse(std::string response) const{
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
kin2DQuad4::ComputeEnergy(){
    //TODO: Integrate over element volume to compute the energy
    return 0.0;
}

//Compute the mass matrix of the element using gauss-integration.
Eigen::MatrixXd 
kin2DQuad4::ComputeMassMatrix(){
    PROFILE_FUNCTION();

    //Use consistent mass definition:
    Eigen::MatrixXd MassMatrix(8,8);
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
kin2DQuad4::ComputeStiffnessMatrix(){
    PROFILE_FUNCTION();

    //Stiffness matrix definition:
    Eigen::MatrixXd StiffnessMatrix(8,8);
    StiffnessMatrix.fill(0.0);

    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Quad", wi, xi);

    //Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
        //Jacobian matrix.
        Eigen::MatrixXd Jij = ComputeJacobianMatrix(xi(i,0), xi(i,1));

        //Gets material tangent matrix at Gauss point.
        Eigen::MatrixXd Cij = theMaterial[i]->GetTangentStiffness();

        //Compute Linear Strain-Displacement Matrix at Gauss Point.
        Eigen::MatrixXd Bij = ComputeLinearStrainDisplacementMatrix(xi(i,0), xi(i,1), Jij);

        //Gets material strain at Gauss point.
        Eigen::VectorXd Stress = theMaterial[i]->GetStress();

        //Computes Second Piola-Kirchhoff stress tensor at Gauss point.
        Eigen::MatrixXd Sij = ComputeSecondPiolaKirchhoffMatrix(Stress);

        //Compute Linear Strain-Displacement Matrix at Gauss Point.
        Eigen::MatrixXd Gij = ComputeNonLinearStrainDisplacementMatrix(xi(i,0), xi(i,1), Jij);

        //Numerical integration.
        StiffnessMatrix += wi(i)*t*fabs(Jij.determinant())*Bij.transpose()*Cij*Bij;
        StiffnessMatrix += wi(i)*t*fabs(Jij.determinant())*Gij.transpose()*Sij*Gij;
    }

    return StiffnessMatrix;
}

//Compute the damping matrix of the element using gauss-integration.
Eigen::MatrixXd 
kin2DQuad4::ComputeDampingMatrix(){
    PROFILE_FUNCTION();

    //Damping matrix definition.
    Eigen::MatrixXd DampingMatrix(8,8);
    DampingMatrix.fill(0.0);

    //Material damping contribution.
    if(theMaterial[0]->IsViscous()){
        //Gets the quadrature information.
        Eigen::VectorXd wi;
        Eigen::MatrixXd xi;
        QuadraturePoints->GetQuadraturePoints("Quad", wi, xi);
        
        //Numerical integration.
        for(unsigned int i = 0; i < wi.size(); i++){
            //Jacobian matrix.
            Eigen::MatrixXd Jij = ComputeJacobianMatrix(xi(i,0), xi(i,1));

            //Compute Linear Strain-Displacement Matrix at Gauss Point.
            Eigen::MatrixXd Bij = ComputeLinearStrainDisplacementMatrix(xi(i,0), xi(i,1), Jij);
            
            //Gets material damping matrix at Gauss point.
            Eigen::MatrixXd Dij = theMaterial[i]->GetDamping();
            
            //Numerical integration.
            DampingMatrix += wi(i)*t*fabs(Jij.determinant())*Bij.transpose()*Dij*Bij;
        }
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
kin2DQuad4::ComputePMLMatrix(){
    Eigen::MatrixXd Kpml;
    return Kpml;
}

//Compute the internal forces acting on the element.
Eigen::VectorXd 
kin2DQuad4::ComputeInternalForces(){
    PROFILE_FUNCTION();

    //Stiffness matrix definition:
    Eigen::VectorXd InternalForces(8);
    InternalForces.fill(0.0);

    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Quad", wi, xi);

    //Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
        //Jacobian matrix.
        Eigen::MatrixXd Jij = ComputeJacobianMatrix(xi(i,0), xi(i,1));

        //Compute Linear Strain-Displacement Matrix at Gauss Point.
        Eigen::MatrixXd Bij = ComputeLinearStrainDisplacementMatrix(xi(i,0), xi(i,1), Jij);

        //Gets material strain at Gauss point.
        Eigen::VectorXd Stress = theMaterial[i]->GetStress();

        //Numerical integration.
        InternalForces += wi(i)*t*fabs(Jij.determinant())*Bij.transpose()*Stress;
    }

    return InternalForces;
}

//Compute the elastic, inertial, and viscous forces acting on the element.
Eigen::VectorXd 
kin2DQuad4::ComputeInternalDynamicForces(){
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
kin2DQuad4::ComputeSurfaceForces(const std::shared_ptr<Load> &surfaceLoad, unsigned int face){
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
kin2DQuad4::ComputeBodyForces(const std::shared_ptr<Load> &bodyLoad, unsigned int k){
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
kin2DQuad4::ComputeDomainReductionForces(const std::shared_ptr<Load> &drm, unsigned int k){
    PROFILE_FUNCTION();

    //Get the Domain-Reduction field motion.
    Eigen::VectorXd x1 = theNodes[0]->GetDomainReductionMotion(k);
    Eigen::VectorXd x2 = theNodes[1]->GetDomainReductionMotion(k);
    Eigen::VectorXd x3 = theNodes[2]->GetDomainReductionMotion(k);
    Eigen::VectorXd x4 = theNodes[3]->GetDomainReductionMotion(k);

    //Constructs the domain reduction boundary/exterior connectivity.
    std::vector<bool> DRMcond(8);
    std::vector<unsigned int> conn = GetNodes();

    for(unsigned int i = 0; i < conn.size(); i ++){
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

//Transform tensor components into vector.
Eigen::VectorXd 
kin2DQuad4::TransformTensorToVector(const Eigen::MatrixXd &Tensor) const{
    //Vector definition.
    Eigen::VectorXd Vector(3);
    Vector << Tensor(0,0), Tensor(1,1), Tensor(0,1);

    return Tensor;
}

//Transform vector components into tensor.
Eigen::MatrixXd 
kin2DQuad4::TransformVectorToTensor(const Eigen::VectorXd &Vector) const{
    //Tensor definition.
    Eigen::MatrixXd Tensor(2,2);
    Tensor << Vector(0), Vector(2),
              Vector(2), Vector(1);

    return Tensor;
}

//Update strain in the element.
Eigen::VectorXd 
kin2DQuad4::ComputeStrain(const double ri, const double si, const Eigen::MatrixXd &Jij) const{
    //Inverse jacobian matrix.
    Eigen::MatrixXd J = Jij.inverse();

    //Gets the element coordinates in deformed configuration. 
    Eigen::VectorXd U1 = theNodes[0]->GetDisplacements() + theNodes[0]->GetIncrementalDisplacements();
    Eigen::VectorXd U2 = theNodes[1]->GetDisplacements() + theNodes[1]->GetIncrementalDisplacements();
    Eigen::VectorXd U3 = theNodes[2]->GetDisplacements() + theNodes[2]->GetIncrementalDisplacements();
    Eigen::VectorXd U4 = theNodes[3]->GetDisplacements() + theNodes[3]->GetIncrementalDisplacements();

    //Shape functions derivatives with respect to undeformed configuration. 
    double N11 = 1.0/4.0*J(0,1)*(-1.0 + ri) + 1.0/4.0*J(0,0)*(-1.0 + si); 
    double N21 = 1.0/4.0*J(0,1)*(-1.0 - ri) + 1.0/4.0*J(0,0)*( 1.0 - si); 
    double N31 = 1.0/4.0*J(0,1)*( 1.0 + ri) + 1.0/4.0*J(0,0)*( 1.0 + si); 
    double N41 = 1.0/4.0*J(0,1)*( 1.0 - ri) + 1.0/4.0*J(0,0)*(-1.0 - si);

    double N12 = 1.0/4.0*J(1,1)*(-1.0 + ri) + 1.0/4.0*J(1,0)*(-1.0 + si); 
    double N22 = 1.0/4.0*J(1,1)*(-1.0 - ri) + 1.0/4.0*J(1,0)*( 1.0 - si);
    double N32 = 1.0/4.0*J(1,1)*( 1.0 + ri) + 1.0/4.0*J(1,0)*( 1.0 + si);
    double N42 = 1.0/4.0*J(1,1)*( 1.0 - ri) + 1.0/4.0*J(1,0)*(-1.0 - si);

    //Displacement derivatives (strains) with respect to undeformed configuration.
    double U11 = N11*U1(0) + N21*U2(0) + N31*U3(0) + N41*U4(0);
    double U21 = N11*U1(1) + N21*U2(1) + N31*U3(1) + N41*U4(1);
    double U12 = N12*U1(0) + N22*U2(0) + N32*U3(0) + N42*U4(0);
    double U22 = N12*U1(1) + N22*U2(1) + N32*U3(1) + N42*U4(1);

    //Green-Lagrange strain components.
    double E11 = U11 + 1.0/2.0*(U11*U11 + U21*U21);
    double E22 = U22 + 1.0/2.0*(U12*U12 + U22*U22);
    double E12 = U12 + U21 + U11*U12 + U21*U22;

    //Green-Lagrange strain vector.
    Eigen::VectorXd Strain(3); 
    Strain << E11, E22, E12;

    return Strain;
}

//Update strain rate in the element.
Eigen::VectorXd 
kin2DQuad4::ComputeStrainRate(const double UNUSED(ri), const double UNUSED(si), const Eigen::MatrixXd& UNUSED(Bij)) const{
    //TODO: Compute strain rate.
    //Strain vector definition:
    Eigen::VectorXd strainrate(3);
    strainrate.fill(0.0);

    return strainrate;
}

//Computes the jacobian of the transformation. 
Eigen::MatrixXd 
kin2DQuad4::ComputeJacobianMatrix(const double ri, const double si) const{
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

    //Jacobian matrix:
    Eigen::MatrixXd Jij(2,2);
    Jij << J11, J12,
           J21, J22;

    return Jij;
}

//Computes the deformation gradient. 
Eigen::MatrixXd 
kin2DQuad4::ComputeDeformationGradientMatrix(const double ri, const double si, const Eigen::MatrixXd &Jij) const{
    //Inverse jacobian matrix:
    Eigen::MatrixXd J = Jij.inverse();

    //Gets the element coordinates in deformed configuration. 
    Eigen::VectorXd x1 = theNodes[0]->GetCoordinates() + theNodes[0]->GetDisplacements() + theNodes[0]->GetIncrementalDisplacements();
    Eigen::VectorXd x2 = theNodes[1]->GetCoordinates() + theNodes[1]->GetDisplacements() + theNodes[1]->GetIncrementalDisplacements();
    Eigen::VectorXd x3 = theNodes[2]->GetCoordinates() + theNodes[2]->GetDisplacements() + theNodes[2]->GetIncrementalDisplacements();
    Eigen::VectorXd x4 = theNodes[3]->GetCoordinates() + theNodes[3]->GetDisplacements() + theNodes[3]->GetIncrementalDisplacements();

    //Coordinates derivatives with respect to isoparametric coordinates.
    double N11 = -1.0/4.0*(1.0 - si)*x1(0) + 1.0/4.0*(1.0 - si)*x2(0) + 1.0/4.0*(1.0 + si)*x3(0) - 1.0/4.0*(1.0 + si)*x4(0);
    double N12 = -1.0/4.0*(1.0 - si)*x1(1) + 1.0/4.0*(1.0 - si)*x2(1) + 1.0/4.0*(1.0 + si)*x3(1) - 1.0/4.0*(1.0 + si)*x4(1); 
    double N21 = -1.0/4.0*(1.0 - ri)*x1(0) - 1.0/4.0*(1.0 + ri)*x2(0) + 1.0/4.0*(1.0 + ri)*x3(0) + 1.0/4.0*(1.0 - ri)*x4(0); 
    double N22 = -1.0/4.0*(1.0 - ri)*x1(1) - 1.0/4.0*(1.0 + ri)*x2(1) + 1.0/4.0*(1.0 + ri)*x3(1) + 1.0/4.0*(1.0 - ri)*x4(1); 

    //Deformation-Gradient coefficients:
    double F11 = N11*J(0,0) + N21*J(0,1);
    double F12 = N11*J(1,0) + N21*J(1,1);
    double F21 = N12*J(0,0) + N22*J(0,1); 
    double F22 = N12*J(1,0) + N22*J(1,1);  

    //Deformation-Gradient Matrix.
    Eigen::MatrixXd Fij(2,2);

    Fij << F11, F12,
           F21, F22;

    return Fij;
}

//Compute the Second Piola-Kirchhoff Stress Tensor.
Eigen::MatrixXd 
kin2DQuad4::ComputeSecondPiolaKirchhoffMatrix(const Eigen::VectorXd &stress) const{
    //Second Piola-Kirchhoff tensor.
    Eigen::MatrixXd Stress(4,4);

    Stress << stress(0), stress(2), 0.0      , 0.0      ,
              stress(2), stress(1), 0.0      , 0.0      ,
              0.0      , 0.0      , stress(0), stress(2),
              0.0      , 0.0      , stress(2), stress(1);

    return Stress;
}

//Compute Shape Function at Gauss Point.
Eigen::MatrixXd 
kin2DQuad4::ComputeShapeFunctionMatrix(const double ri, const double si) const{
    //Shape function coefficients:
    double H11 = 1.0/4.0*(1.0 - ri)*(1.0 - si);
    double H22 = 1.0/4.0*(1.0 + ri)*(1.0 - si);
    double H33 = 1.0/4.0*(1.0 + ri)*(1.0 + si);
    double H44 = 1.0/4.0*(1.0 - ri)*(1.0 + si);

    //Shape function matrix:
    Eigen::MatrixXd Hij(2,8);
    Hij << H11, 0.0, H22, 0.0, H33, 0.0, H44, 0.0,
           0.0, H11, 0.0, H22, 0.0, H33, 0.0, H44;

    return Hij;
}

//Evaluates the Linear (material) strain-displacemnt matrix at a given Gauss point.
Eigen::MatrixXd 
kin2DQuad4::ComputeLinearStrainDisplacementMatrix(const double ri, const double si, const Eigen::MatrixXd &Jij) const{
    //Inverse jacobian matrix:
    Eigen::MatrixXd J = Jij.inverse();

    //Gets the element coordinates in deformed configuration. 
    Eigen::VectorXd U1 = theNodes[0]->GetDisplacements() + theNodes[0]->GetIncrementalDisplacements();
    Eigen::VectorXd U2 = theNodes[1]->GetDisplacements() + theNodes[1]->GetIncrementalDisplacements();
    Eigen::VectorXd U3 = theNodes[2]->GetDisplacements() + theNodes[2]->GetIncrementalDisplacements();
    Eigen::VectorXd U4 = theNodes[3]->GetDisplacements() + theNodes[3]->GetIncrementalDisplacements();

    //Shape functions derivatives with respect to undeformed configuration. 
    double N11 = 1.0/4.0*J(0,1)*(-1.0 + ri) + 1.0/4.0*J(0,0)*(-1.0 + si); 
    double N21 = 1.0/4.0*J(0,1)*(-1.0 - ri) + 1.0/4.0*J(0,0)*( 1.0 - si); 
    double N31 = 1.0/4.0*J(0,1)*( 1.0 + ri) + 1.0/4.0*J(0,0)*( 1.0 + si); 
    double N41 = 1.0/4.0*J(0,1)*( 1.0 - ri) + 1.0/4.0*J(0,0)*(-1.0 - si);

    double N12 = 1.0/4.0*J(1,1)*(-1.0 + ri) + 1.0/4.0*J(1,0)*(-1.0 + si); 
    double N22 = 1.0/4.0*J(1,1)*(-1.0 - ri) + 1.0/4.0*J(1,0)*( 1.0 - si);
    double N32 = 1.0/4.0*J(1,1)*( 1.0 + ri) + 1.0/4.0*J(1,0)*( 1.0 + si);
    double N42 = 1.0/4.0*J(1,1)*( 1.0 - ri) + 1.0/4.0*J(1,0)*(-1.0 - si);

    //Displacement derivatives (strains) with respect to undeformed configuration.
    double U11 = N11*U1(0) + N21*U2(0) + N31*U3(0) + N41*U4(0);
    double U21 = N11*U1(1) + N21*U2(1) + N31*U3(1) + N41*U4(1);
    double U12 = N12*U1(0) + N22*U2(0) + N32*U3(0) + N42*U4(0);
    double U22 = N12*U1(1) + N22*U2(1) + N32*U3(1) + N42*U4(1);

    //Linear Strain-displacement matrix coefficients.
    double B11 = N11 + U11*N11; double B21 =       U12*N12; double B31 = N12 + U11*N12 + U12*N11;
    double B12 =       U21*N11; double B22 = N12 + U22*N12; double B32 = N11 + U21*N12 + U22*N11;

    double B13 = N21 + U11*N21; double B23 =       U12*N22; double B33 = N22 + U11*N22 + U12*N21;
    double B14 =       U21*N21; double B24 = N22 + U22*N22; double B34 = N21 + U21*N22 + U22*N21;

    double B15 = N31 + U11*N31; double B25 =       U12*N32; double B35 = N32 + U11*N32 + U12*N31;
    double B16 =       U21*N31; double B26 = N32 + U22*N32; double B36 = N31 + U21*N32 + U22*N31;

    double B17 = N41 + U11*N41; double B27 =       U12*N42; double B37 = N42 + U11*N42 + U12*N41;
    double B18 =       U21*N41; double B28 = N42 + U22*N42; double B38 = N41 + U21*N42 + U22*N41;

    //Linear Strain-displacement matrix.
    Eigen::MatrixXd Bij(3,8);
    Bij << B11, B12, B13, B14, B15, B16, B17, B18,
           B21, B22, B23, B24, B25, B26, B27, B28,        
           B31, B32, B33, B34, B35, B36, B37, B38;

    return Bij;
}

//Evaluates the Non-Linear (geometric) strain-displacemnt matrix at a given Gauss point.
Eigen::MatrixXd 
kin2DQuad4::ComputeNonLinearStrainDisplacementMatrix(const double ri, const double si, const Eigen::MatrixXd &Jij) const{
    //Inverse jacobian matrix.
    Eigen::MatrixXd J = Jij.inverse();

    //Shape functions derivatives with respect to undeformed configuration. 
    double N11 = 1.0/4.0*J(0,1)*(-1.0 + ri) + 1.0/4.0*J(0,0)*(-1.0 + si); 
    double N21 = 1.0/4.0*J(0,1)*(-1.0 - ri) + 1.0/4.0*J(0,0)*( 1.0 - si); 
    double N31 = 1.0/4.0*J(0,1)*( 1.0 + ri) + 1.0/4.0*J(0,0)*( 1.0 + si); 
    double N41 = 1.0/4.0*J(0,1)*( 1.0 - ri) + 1.0/4.0*J(0,0)*(-1.0 - si);

    double N12 = 1.0/4.0*J(1,1)*(-1.0 + ri) + 1.0/4.0*J(1,0)*(-1.0 + si); 
    double N22 = 1.0/4.0*J(1,1)*(-1.0 - ri) + 1.0/4.0*J(1,0)*( 1.0 - si);
    double N32 = 1.0/4.0*J(1,1)*( 1.0 + ri) + 1.0/4.0*J(1,0)*( 1.0 + si);
    double N42 = 1.0/4.0*J(1,1)*( 1.0 - ri) + 1.0/4.0*J(1,0)*(-1.0 - si);

    //Geometric (non-linear) strain-displacement matrix.
    Eigen::MatrixXd Gij(4,8);
    Gij << N11, 0.0, N21, 0.0, N31, 0.0, N41, 0.0,
           N12, 0.0, N22, 0.0, N32, 0.0, N42, 0.0,
           0.0, N11, 0.0, N21, 0.0, N31, 0.0, N41,
           0.0, N12, 0.0, N22, 0.0, N32, 0.0, N42;

    return Gij;
}

//Evaluates the linearized deformation matrix matrix at a given Gauss point.
Eigen::MatrixXd 
kin2DQuad4::ComputeLinearizedStrainDisplacementMatrix(const double ri, const double si, const Eigen::MatrixXd &Jij) const{
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

    //Shape function matrix:
    Eigen::MatrixXd Bij(3,8);
    Bij << B11, 0.0, B21, 0.0, B31, 0.0, B41, 0.0,
           0.0, B12, 0.0, B22, 0.0, B32, 0.0, B42,
           B12, B11, B22, B21, B32, B31, B42, B41;

    return Bij;
}

//Compute the initial stiffness matrix of the element using gauss-integration.
Eigen::MatrixXd 
kin2DQuad4::ComputeInitialStiffnessMatrix() const{
    //Stiffness matrix definition:
    Eigen::MatrixXd StiffnessMatrix(8,8);
    StiffnessMatrix.fill(0.0);

    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Quad", wi, xi);

    //Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
        //Jacobian matrix.
        Eigen::MatrixXd Jij = ComputeJacobianMatrix(xi(i,0), xi(i,1));

        //Gets material tangent matrix at Gauss point.
        Eigen::MatrixXd Cij = theMaterial[i]->GetInitialTangentStiffness();

        //Compute Linear Strain-Displacement Matrix at Gauss Point.
        Eigen::MatrixXd Bij = ComputeLinearizedStrainDisplacementMatrix(xi(i,0), xi(i,1), Jij);

        //Numerical integration.
        StiffnessMatrix += wi(i)*t*fabs(Jij.determinant())*Bij.transpose()*Cij*Bij;
    }

    return StiffnessMatrix;
}
