#include <cmath>
#include <Eigen/LU> 
#include "Material.hpp"
#include "kin3DHexa8.hpp"
#include "GaussQuadrature.hpp"
#include "LobattoQuadrature.hpp"
#include "Definitions.hpp"
#include "Profiler.hpp"

//Define VTK cell value for Paraview:
const unsigned int VTKCELL = 12;

//Overload constructor.
kin3DHexa8::kin3DHexa8(const std::vector<unsigned int> nodes, std::unique_ptr<Material> &material, const std::string quadrature, const unsigned int nGauss) :
Element("kin3DHexa8", nodes, 24, VTKCELL){
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
kin3DHexa8::~kin3DHexa8(){
    //Does nothing
}

//Save the material states in the element.
void 
kin3DHexa8::CommitState(){
    //Updates the viscous material components.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    if(theMaterial[0]->IsViscous()){
        //Gets the quadrature information.    
        Eigen::VectorXd wi;
        Eigen::MatrixXd xi;
        QuadraturePoints->GetQuadraturePoints("Hexa", wi, xi);

        //Update material states.
        for(unsigned int k = 0; k < wi.size(); k++){
            //Jacobian matrix.
            Eigen::MatrixXd Jij = ComputeJacobianMatrix(xi(k,0), xi(k,1), xi(k,2));

            //Compute Strain-Displacement Matrix at Gauss Point.
            Eigen::MatrixXd Bij = ComputeLinearStrainDisplacementMatrix(xi(k,0), xi(k,1), xi(k,2), Jij);

            //Computes Strain Rate vector.
            Eigen::VectorXd strainrate = ComputeStrainRate(xi(k,0), xi(k,1), xi(k,2), Bij);

            //Update the material state.
            theMaterial[k]->UpdateState(strainrate, 2);
        }
    }

    for(unsigned int k = 0; k < nPoints; k++)
        theMaterial[k]->CommitState();
}

//Reverse the material states to previous converged state in this element.
void 
kin3DHexa8::ReverseState(){
    //Reverse the material components.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    for(unsigned int k = 0; k < nPoints; k++)
        theMaterial[k]->ReverseState();
}

//Brings the material state to its initial state in this element.
void 
kin3DHexa8::InitialState(){
    //Brings the material components to initial state.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    for(unsigned int k = 0; k < nPoints; k++)
        theMaterial[k]->InitialState();
}

//Update the material states in the element.
void 
kin3DHexa8::UpdateState(){
    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Hexa", wi, xi);

    //Updates the material components.    
    for(unsigned int k = 0; k < wi.size(); k++){
        //Jacobian matrix.
        Eigen::MatrixXd Jij = ComputeJacobianMatrix(xi(k,0), xi(k,1), xi(k,2));

        //Computes strain vector.
        Eigen::VectorXd strain = ComputeStrain(xi(k,0), xi(k,1), xi(k,2), Jij);

        //Update the material state.
        theMaterial[k]->UpdateState(strain, 1);
    }
}

//Sets the finite element dependance among objects.
void 
kin3DHexa8::SetDomain(std::map<unsigned int, std::shared_ptr<Node> > &nodes){
    //Gets the global element connectivity.
    std::vector<unsigned int> conn = GetNodes();

    //Assign the element to mesh node pointer.  
    for(unsigned int i = 0; i < GetNumberOfNodes(); i++){
        theNodes[i] = nodes[conn[i]];
    }
}

//Sets the damping model.
void 
kin3DHexa8::SetDamping(const std::shared_ptr<Damping> &damping){
    //The damping model
    theDamping = damping;
}

//Gets the list of total-degree of freedom of this element.
std::vector<unsigned int> 
kin3DHexa8::GetTotalDegreeOfFreedom() const{
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
kin3DHexa8::GetStrain() const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    Eigen::MatrixXd theStrain(nPoints,6);
    for(unsigned int k = 0; k < nPoints; k++)
        theStrain.row(k) = theMaterial[k]->GetStrain();

    return theStrain;
}

//Returns the material stress at integration points.
Eigen::MatrixXd 
kin3DHexa8::GetStress() const{
    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Hexa", wi, xi);
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    //Stress values at each integration point.
    Eigen::MatrixXd theStress(nPoints,6);
    for(unsigned int k = 0; k < wi.size(); k++){
        //Jacobian matrix.
        Eigen::MatrixXd Jij = ComputeJacobianMatrix(xi(k,0), xi(k,1), xi(k,2));

        //Deformation gradient matrix.
        Eigen::MatrixXd Fij = ComputeDeformationGradientMatrix(xi(k,0), xi(k,1), xi(k,2), Jij);

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
kin3DHexa8::GetStrainRate() const{
    //TODO: Compute strain rate using large-deformation.
    //Total number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    Eigen::MatrixXd theStrainRate(nPoints,6);
    for(unsigned int k = 0; k < nPoints; k++)
        theStrainRate.row(k) = theMaterial[k]->GetStrainRate();

    return theStrainRate;
}

//Gets the material strain in section at  coordinate (x3,x2).
Eigen::MatrixXd 
kin3DHexa8::GetStrainAt(double UNUSED(x3), double UNUSED(x2)) const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    //Stress at coordinate is define within section.
    Eigen::MatrixXd theStrain(nPoints, 6);
    theStrain.fill(0.0);

    return theStrain;
}

//Gets the material stress in section at  coordinate (x3,x2).
Eigen::MatrixXd 
kin3DHexa8::GetStressAt(double UNUSED(x3), double UNUSED(x2)) const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    //Stress at coordinate is define within section.
    Eigen::MatrixXd theStress(nPoints, 6);
    theStress.fill(0.0);

    return theStress;
}

//Gets the element internal response in VTK format.
Eigen::VectorXd 
kin3DHexa8::GetVTKResponse(std::string response) const{
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
kin3DHexa8::ComputeEnergy(){
    //TODO: Integrate over element volume to compute the energy
    return 0.0;
}

//Compute the mass matrix of the element using gauss-integration.
Eigen::MatrixXd 
kin3DHexa8::ComputeMassMatrix(){
    PROFILE_FUNCTION();

    //Use consistent mass definition:
    Eigen::MatrixXd MassMatrix(24,24);
    MassMatrix.fill(0.0);

    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Hexa", wi, xi);

    //Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
        //Gets material properties:
        double rho = theMaterial[i]->GetDensity();

        //Jacobian matrix:
        Eigen::MatrixXd Jij = ComputeJacobianMatrix(xi(i,0), xi(i,1), xi(i,2));

        //Compute Interpolation Function Matrix at Gauss Point:
        Eigen::MatrixXd Hij = ComputeShapeFunctionMatrix(xi(i,0), xi(i,1), xi(i,2));

        //Numerical integration:
        MassMatrix += wi(i)*rho*fabs(Jij.determinant())*Hij.transpose()*Hij;
    }

    //Lumped Mass Formulation.
    if(MassFormulation){
        //Lumped Mass in diagonal terms.
        for (unsigned int i = 0; i < 24; i++){
            for (unsigned int j = 0; j < 24; j++){
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
kin3DHexa8::ComputeStiffnessMatrix(){
    PROFILE_FUNCTION();

    //Stiffness matrix definition:
    Eigen::MatrixXd StiffnessMatrix(24,24);
    StiffnessMatrix.fill(0.0);

    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Hexa", wi, xi);

    //Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
        //Jacobian matrix.
        Eigen::MatrixXd Jij = ComputeJacobianMatrix(xi(i,0), xi(i,1), xi(i,2));

        //Gets material tangent matrix at Gauss point.
        Eigen::MatrixXd Cij = theMaterial[i]->GetTangentStiffness();

        //Compute Strain-Displacement Matrix at Gauss Point.
        Eigen::MatrixXd Bij = ComputeLinearStrainDisplacementMatrix(xi(i,0), xi(i,1), xi(i,2), Jij);

        //Gets material strain at Gauss point.
        Eigen::VectorXd Stress = theMaterial[i]->GetStress();

        //Computes Second Piola-Kirchhoff stress tensor at Gauss point.
        Eigen::MatrixXd Sij = ComputeSecondPiolaKirchhoffMatrix(Stress);

        //Compute Linear Strain-Displacement Matrix at Gauss Point.
        Eigen::MatrixXd Gij = ComputeNonLinearStrainDisplacementMatrix(xi(i,0), xi(i,1), xi(i,2), Jij);

        //Numerical integration.
        StiffnessMatrix += wi(i)*fabs(Jij.determinant())*Bij.transpose()*Cij*Bij;
        StiffnessMatrix += wi(i)*fabs(Jij.determinant())*Gij.transpose()*Sij*Gij;
    }

    return StiffnessMatrix;
}

//Compute the damping matrix of the element using gauss-integration.
Eigen::MatrixXd 
kin3DHexa8::ComputeDampingMatrix(){
    PROFILE_FUNCTION();

    //Damping matrix definition
    Eigen::MatrixXd DampingMatrix(24,24);
    DampingMatrix.fill(0.0);
    
    //Material damping contribution.
    if(theMaterial[0]->IsViscous()){
        //Gets the quadrature information.    
        Eigen::VectorXd wi;
        Eigen::MatrixXd xi;
        QuadraturePoints->GetQuadraturePoints("Hexa", wi, xi);
        
        //Numerical integration.
        for(unsigned int i = 0; i < wi.size(); i++){
            //Jacobian matrix.
            Eigen::MatrixXd Jij = ComputeJacobianMatrix(xi(i,0), xi(i,1), xi(i,2));
        
            //Compute Strain-Displacement Matrix at Gauss Point.
            Eigen::MatrixXd Bij = ComputeLinearStrainDisplacementMatrix(xi(i,0), xi(i,1), xi(i,2), Jij);
            
            //Gets material damping matrix at Gauss point.
            Eigen::MatrixXd Dij = theMaterial[i]->GetDamping();
            
            //Numerical integration.
            DampingMatrix += wi(i)*fabs(Jij.determinant())*Bij.transpose()*Dij*Bij;
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
kin3DHexa8::ComputePMLMatrix(){
    Eigen::MatrixXd Kpml;
    return Kpml;
}

//Compute the internal forces acting on the element.
Eigen::VectorXd 
kin3DHexa8::ComputeInternalForces(){
    PROFILE_FUNCTION();

    //Stiffness matrix definition.
    Eigen::VectorXd InternalForces(24);
    InternalForces.fill(0.0);

    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Hexa", wi, xi);

    //Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
        //Jacobian matrix.
        Eigen::MatrixXd Jij = ComputeJacobianMatrix(xi(i,0), xi(i,1), xi(i,2));

        //Compute Strain-Displacement Matrix at Gauss Point.
        Eigen::MatrixXd Bij = ComputeLinearStrainDisplacementMatrix(xi(i,0), xi(i,1), xi(i,2), Jij);

        //Gets material strain at Gauss point.
        Eigen::VectorXd Stress = theMaterial[i]->GetStress();

        //Numerical integration.
        InternalForces += wi(i)*fabs(Jij.determinant())*Bij.transpose()*Stress;
    }
    
    return InternalForces;
}

//Compute the elastic, inertial, and viscous forces acting on the element.
Eigen::VectorXd 
kin3DHexa8::ComputeInternalDynamicForces(){
    //The Internal dynamic force vector
    Eigen::VectorXd InternalForces;

    if( HasFixedNode(theNodes) ){
        //Allocate memory for velocity/acceleraton. 
        Eigen::VectorXd V(24); 
        Eigen::VectorXd A(24);

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
kin3DHexa8::ComputeSurfaceForces(const std::shared_ptr<Load> &surfaceLoad, unsigned int face){
    PROFILE_FUNCTION();

    //Local surface load vector:
    Eigen::VectorXd surfaceForces(24);
    surfaceForces.fill(0.0);

    //Gets the surface load:
    Eigen::VectorXd qs = surfaceLoad->GetLoadVector();

    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Quad", wi, xi);

    if(face == 1){
        //Gets the face coordinates in undeformed configuration. 
        Eigen::VectorXd x1 = theNodes[0]->GetCoordinates();
        Eigen::VectorXd x2 = theNodes[1]->GetCoordinates();
        Eigen::VectorXd x3 = theNodes[2]->GetCoordinates();
        Eigen::VectorXd x4 = theNodes[3]->GetCoordinates();

        //Numerical integration in local axis r-s:
        for(unsigned int i = 0; i < wi.size(); i++){
            double r = xi(i,0);
            double s = xi(i,1);

            //vectors along s and t axes.
            Eigen::Vector3d v1, v2;
            v1 << 1.0/4.0*(1.0 - s)*(x2 - x1) + 1.0/4.0*(1.0 + s)*(x3 - x4);
            v2 << 1.0/4.0*(1.0 - r)*(x4 - x1) + 1.0/4.0*(1.0 + r)*(x3 - x2);

            //Jacobian matrix:
            double detJij = v1.cross(v2).norm();

            //Compute Strain-Displacement Matrix at Gauss Point:
            Eigen::MatrixXd Hij = ComputeShapeFunctionMatrix(r, s, -1.0);

            //Numerical integration:
            surfaceForces += wi(i)*detJij*Hij.transpose()*qs;
        }
    }
    if(face == 2){
        //Gets the face coordinates in undeformed configuration. 
        Eigen::VectorXd x1 = theNodes[0]->GetCoordinates();
        Eigen::VectorXd x2 = theNodes[1]->GetCoordinates();
        Eigen::VectorXd x6 = theNodes[5]->GetCoordinates();
        Eigen::VectorXd x5 = theNodes[4]->GetCoordinates();

        //Numerical integration in local axis r-t:
        for(unsigned int i = 0; i < wi.size(); i++){
            double r = xi(i,0);
            double t = xi(i,1);

            //vectors along s and t axes.
            Eigen::Vector3d v1, v2;
            v1 << 1.0/4.0*(1.0 - t)*(x2 - x1) + 1.0/4.0*(1.0 + t)*(x6 - x5);
            v2 << 1.0/4.0*(1.0 - r)*(x5 - x1) + 1.0/4.0*(1.0 + r)*(x6 - x2);

            //Jacobian matrix:
            double detJij = v1.cross(v2).norm();

            //Compute Strain-Displacement Matrix at Gauss Point:
            Eigen::MatrixXd Hij = ComputeShapeFunctionMatrix(r, -1.0, t);

            //Numerical integration:
            surfaceForces += wi(i)*detJij*Hij.transpose()*qs;
        }
    }
    if(face == 3){
        //Gets the face coordinates in undeformed configuration. 
        Eigen::VectorXd x2 = theNodes[1]->GetCoordinates();
        Eigen::VectorXd x3 = theNodes[2]->GetCoordinates();
        Eigen::VectorXd x7 = theNodes[6]->GetCoordinates();
        Eigen::VectorXd x6 = theNodes[5]->GetCoordinates();

        //Numerical integration in local axis s-t:
        for(unsigned int i = 0; i < wi.size(); i++){
            double s = xi(i,0);
            double t = xi(i,1);

            //vectors along s and t axes.
            Eigen::Vector3d v1, v2;
            v1 << 1.0/4.0*(1.0 - t)*(x3 - x2) + 1.0/4.0*(1.0 + t)*(x7 - x6);
            v2 << 1.0/4.0*(1.0 - s)*(x6 - x2) + 1.0/4.0*(1.0 + s)*(x7 - x3);

            //Jacobian matrix:
            double detJij = v1.cross(v2).norm();

            //Compute Strain-Displacement Matrix at Gauss Point:
            Eigen::MatrixXd Hij = ComputeShapeFunctionMatrix(1.0, s, t);

            //Numerical integration:
            surfaceForces += wi(i)*detJij*Hij.transpose()*qs;
        }
    }
    if(face == 4){
        //Gets the face coordinates in undeformed configuration. 
        Eigen::VectorXd x4 = theNodes[3]->GetCoordinates();
        Eigen::VectorXd x3 = theNodes[2]->GetCoordinates();
        Eigen::VectorXd x7 = theNodes[6]->GetCoordinates();
        Eigen::VectorXd x8 = theNodes[7]->GetCoordinates();

        //Numerical integration in local axis r-t:
        for(unsigned int i = 0; i < wi.size(); i++){
            double r = xi(i,0);
            double t = xi(i,1);

            //vectors along r and t axes.
            Eigen::Vector3d v1, v2;
            v1 << 1.0/4.0*(1.0 - t)*(x3 - x4) + 1.0/4.0*(1.0 + t)*(x7 - x8);
            v2 << 1.0/4.0*(1.0 - r)*(x8 - x4) + 1.0/4.0*(1.0 + r)*(x7 - x3);

            //Jacobian matrix:
            double detJij = v1.cross(v2).norm();

            //Compute Strain-Displacement Matrix at Gauss Point:
            Eigen::MatrixXd Hij = ComputeShapeFunctionMatrix(r, 1.0, t);

            //Numerical integration:
            surfaceForces += wi(i)*detJij*Hij.transpose()*qs;
        }
    }
    if(face == 5){
        //Gets the face coordinates in undeformed configuration. 
        Eigen::VectorXd x1 = theNodes[0]->GetCoordinates();
        Eigen::VectorXd x4 = theNodes[3]->GetCoordinates();
        Eigen::VectorXd x8 = theNodes[7]->GetCoordinates();
        Eigen::VectorXd x5 = theNodes[4]->GetCoordinates();

        //Numerical integration in local axis s-t:
        for(unsigned int i = 0; i < wi.size(); i++){
            double s = xi(i,0);
            double t = xi(i,1);

            //vectors along s and t axes.
            Eigen::Vector3d v1, v2;
            v1 << 1.0/4.0*(1.0 - t)*(x4 - x1) + 1.0/4.0*(1.0 + t)*(x8 - x5);
            v2 << 1.0/4.0*(1.0 - s)*(x5 - x1) + 1.0/4.0*(1.0 + s)*(x8 - x4);

            //Jacobian matrix:
            double detJij = v1.cross(v2).norm();

            //Compute Strain-Displacement Matrix at Gauss Point:
            Eigen::MatrixXd Hij = ComputeShapeFunctionMatrix(-1.0, s, t);

            //Numerical integration:
            surfaceForces += wi(i)*detJij*Hij.transpose()*qs;
        }
    }
    if(face == 6){
        //Gets the face coordinates in undeformed configuration. 
        Eigen::VectorXd x5 = theNodes[4]->GetCoordinates();
        Eigen::VectorXd x6 = theNodes[5]->GetCoordinates();
        Eigen::VectorXd x7 = theNodes[6]->GetCoordinates();
        Eigen::VectorXd x8 = theNodes[7]->GetCoordinates();

        //Numerical integration in local axis r-s:
        for(unsigned int i = 0; i < wi.size(); i++){
            double r = xi(i,0);
            double s = xi(i,1);

            //vectors along r and s axes.
            Eigen::Vector3d v1, v2;
            v1 << 1.0/4.0*(1.0 - s)*(x6 - x5) + 1.0/4.0*(1.0 + s)*(x7 - x8);
            v2 << 1.0/4.0*(1.0 - r)*(x8 - x5) + 1.0/4.0*(1.0 + r)*(x7 - x6);

            //Jacobian matrix:
            double detJij = v1.cross(v2).norm();

            //Compute Strain-Displacement Matrix at Gauss Point:
            Eigen::MatrixXd Hij = ComputeShapeFunctionMatrix(r, s, 1.0);

            //Numerical integration:
            surfaceForces += wi(i)*detJij*Hij.transpose()*qs;
        }
    }

    return surfaceForces;
}

//Compute the body forces acting on the element.
Eigen::VectorXd 
kin3DHexa8::ComputeBodyForces(const std::shared_ptr<Load> &bodyLoad, unsigned int k){
    PROFILE_FUNCTION();

    //Local body load vector:
    Eigen::VectorXd bodyForces(24);
    bodyForces.fill(0.0);

    //Gets the body force:
    Eigen::VectorXd qb = bodyLoad->GetLoadVector(k);
    
    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Hexa", wi, xi);

    //Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
        //Gets material properties:
        double rho = theMaterial[i]->GetDensity();

        //Jacobian matrix:
        Eigen::MatrixXd Jij = ComputeJacobianMatrix(xi(i,0), xi(i,1), xi(i,2));

        //Compute Interpolation Function Matrix at Gauss Point:
        Eigen::MatrixXd Hij = ComputeShapeFunctionMatrix(xi(i,0), xi(i,1), xi(i,2));

        //Numerical integration:
        bodyForces += wi(i)*rho*fabs(Jij.determinant())*Hij.transpose()*qb;
    }
    
    return bodyForces;
}

//Compute the domain reduction forces acting on the element.
Eigen::VectorXd 
kin3DHexa8::ComputeDomainReductionForces(const std::shared_ptr<Load> &drm, unsigned int k){
    PROFILE_FUNCTION();

    //Local domain-reduction load vector.
    Eigen::VectorXd DRMForces(24);

    //Get the Domain-Reduction field motion.
    Eigen::VectorXd x1 = theNodes[0]->GetDomainReductionMotion(k);
    Eigen::VectorXd x2 = theNodes[1]->GetDomainReductionMotion(k);
    Eigen::VectorXd x3 = theNodes[2]->GetDomainReductionMotion(k);
    Eigen::VectorXd x4 = theNodes[3]->GetDomainReductionMotion(k);
    Eigen::VectorXd x5 = theNodes[4]->GetDomainReductionMotion(k);
    Eigen::VectorXd x6 = theNodes[5]->GetDomainReductionMotion(k);
    Eigen::VectorXd x7 = theNodes[6]->GetDomainReductionMotion(k);
    Eigen::VectorXd x8 = theNodes[7]->GetDomainReductionMotion(k);

    //Constructs the domain reduction boundary/exterior connectivity.
    std::vector<bool> DRMcond(24);
    std::vector<unsigned int> conn = GetNodes();

    for(unsigned int i = 0; i < conn.size(); i ++){
        bool condition = drm->GetDRMCondition(conn[i]);
        DRMcond[3*i  ] = condition;
        DRMcond[3*i+1] = condition;
        DRMcond[3*i+2] = condition;
    }

    //Constructs the displacement, velocity and acceleration vectors. 
    Eigen::VectorXd Uo(24); 
    Eigen::VectorXd Vo(24);
    Eigen::VectorXd Ao(24);
 
    Uo << x1(0), x1(1), x1(2), x2(0), x2(1), x2(2), x3(0), x3(1), x3(2), x4(0), x4(1), x4(2), x5(0), x5(1), x5(2), x6(0), x6(1), x6(2), x7(0), x7(1), x7(2), x8(0), x8(1), x8(2);
    Vo << x1(3), x1(4), x1(5), x2(3), x2(4), x2(5), x3(3), x3(4), x3(5), x4(3), x4(4), x4(5), x5(3), x5(4), x5(5), x6(3), x6(4), x6(5), x7(3), x7(4), x7(5), x8(3), x8(4), x8(5);
    Ao << x1(6), x1(7), x1(8), x2(6), x2(7), x2(8), x3(6), x3(7), x3(8), x4(6), x4(7), x4(8), x5(6), x5(7), x5(8), x6(6), x6(7), x6(8), x7(6), x7(7), x7(8), x8(6), x8(7), x8(8);

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

//Transform tensor components into vector.
Eigen::VectorXd 
kin3DHexa8::TransformTensorToVector(const Eigen::MatrixXd &Tensor) const{
    //Vector definiition.
    Eigen::VectorXd Vector(6);

    Vector << Tensor(0,0), Tensor(1,1), Tensor(2,2), Tensor(0,1), Tensor(1,2), Tensor(0,2);

    return Tensor;
}

//Transform vector components into tensor.
Eigen::MatrixXd 
kin3DHexa8::TransformVectorToTensor(const Eigen::VectorXd &Vector) const{
    //Tensor definiition.
    Eigen::MatrixXd Tensor(3,3);

    Tensor << Vector(0), Vector(3), Vector(5),
              Vector(3), Vector(1), Vector(4),
              Vector(5), Vector(4), Vector(2);

    return Tensor;
}

//Update strain in the element.
Eigen::VectorXd 
kin3DHexa8::ComputeStrain(const double ri, const double si, const double ti, const Eigen::MatrixXd &Jij) const{
    //Inverse jacobian matrix.
    Eigen::MatrixXd J = Jij.inverse();

    //Gets the element coordinates in deformed configuration. 
    Eigen::VectorXd U1 = theNodes[0]->GetDisplacements() + theNodes[0]->GetIncrementalDisplacements();
    Eigen::VectorXd U2 = theNodes[1]->GetDisplacements() + theNodes[1]->GetIncrementalDisplacements();
    Eigen::VectorXd U3 = theNodes[2]->GetDisplacements() + theNodes[2]->GetIncrementalDisplacements();
    Eigen::VectorXd U4 = theNodes[3]->GetDisplacements() + theNodes[3]->GetIncrementalDisplacements();
    Eigen::VectorXd U5 = theNodes[4]->GetDisplacements() + theNodes[4]->GetIncrementalDisplacements();
    Eigen::VectorXd U6 = theNodes[5]->GetDisplacements() + theNodes[5]->GetIncrementalDisplacements();
    Eigen::VectorXd U7 = theNodes[6]->GetDisplacements() + theNodes[6]->GetIncrementalDisplacements();
    Eigen::VectorXd U8 = theNodes[7]->GetDisplacements() + theNodes[7]->GetIncrementalDisplacements();

    //Shape functions derivatives with respect to undeformed configuration (dNI/dXj).
    double N11 = -1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 - ti);
    double N21 = -1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 - ti);
    double N31 = -1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 - ti);
    double N41 = -1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 - ti);
    double N51 =  1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 + ti);
    double N61 =  1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 + ti);
    double N71 =  1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 + ti);
    double N81 =  1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 + ti);

    double N12 = -1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 - ti);
    double N22 = -1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 - ti);
    double N32 = -1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 - ti);
    double N42 = -1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 - ti);
    double N52 =  1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 + ti);
    double N62 =  1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 + ti);
    double N72 =  1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 + ti);
    double N82 =  1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 + ti);

    double N13 = -1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 - ti);
    double N23 = -1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 - ti);
    double N33 = -1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 - ti);
    double N43 = -1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 - ti);
    double N53 =  1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 + ti);
    double N63 =  1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 + ti);
    double N73 =  1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 + ti);
    double N83 =  1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 + ti);

    //Displacement derivatives (strains) with respect to undeformed configuration dUi/dXj.
    double U11 = N11*U1(0) + N21*U2(0) + N31*U3(0) + N41*U4(0) + N51*U5(0) + N61*U6(0) + N71*U7(0) + N81*U8(0);
    double U21 = N11*U1(1) + N21*U2(1) + N31*U3(1) + N41*U4(1) + N51*U5(1) + N61*U6(1) + N71*U7(1) + N81*U8(1);
    double U31 = N11*U1(2) + N21*U2(2) + N31*U3(2) + N41*U4(2) + N51*U5(2) + N61*U6(2) + N71*U7(2) + N81*U8(2);

    double U12 = N12*U1(0) + N22*U2(0) + N32*U3(0) + N42*U4(0) + N52*U5(0) + N62*U6(0) + N72*U7(0) + N82*U8(0);
    double U22 = N12*U1(1) + N22*U2(1) + N32*U3(1) + N42*U4(1) + N52*U5(1) + N62*U6(1) + N72*U7(1) + N82*U8(1);
    double U32 = N12*U1(2) + N22*U2(2) + N32*U3(2) + N42*U4(2) + N52*U5(2) + N62*U6(2) + N72*U7(2) + N82*U8(2);

    double U13 = N13*U1(0) + N23*U2(0) + N33*U3(0) + N43*U4(0) + N53*U5(0) + N63*U6(0) + N73*U7(0) + N83*U8(0);
    double U23 = N13*U1(1) + N23*U2(1) + N33*U3(1) + N43*U4(1) + N53*U5(1) + N63*U6(1) + N73*U7(1) + N83*U8(1);
    double U33 = N13*U1(2) + N23*U2(2) + N33*U3(2) + N43*U4(2) + N53*U5(2) + N63*U6(2) + N73*U7(2) + N83*U8(2);

    //Green-Lagrange strain components.
    double E11 = U11 + 1.0/2.0*(U11*U11 + U21*U21 + U31*U31);
    double E22 = U22 + 1.0/2.0*(U12*U12 + U22*U22 + U32*U32);
    double E33 = U33 + 1.0/2.0*(U13*U13 + U23*U23 + U33*U33);
    double E12 = U12 + U21 + U11*U12 + U21*U22 + U31*U32;
    double E23 = U23 + U32 + U12*U13 + U22*U23 + U32*U33;
    double E31 = U13 + U31 + U11*U13 + U21*U23 + U31*U33;

    //Strain vector:
    Eigen::VectorXd Strain(6);
    Strain << E11, E22, E33, E12, E23, E31; 

    return Strain;
}

//Update strain rate in the element.
Eigen::VectorXd 
kin3DHexa8::ComputeStrainRate(const double UNUSED(ri), const double UNUSED(si), const double UNUSED(ti), const Eigen::MatrixXd& UNUSED(Bij)) const{
    //TODO: Compute strain rate.
    //Strain vector definition:
    Eigen::VectorXd strainrate(6);
    strainrate.fill(0.0);

    return strainrate;
}

//Computes the jacobian of the transformation. 
Eigen::MatrixXd 
kin3DHexa8::ComputeJacobianMatrix(const double ri, const double si, const double ti) const{
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

    //Jacobian matrix:
    Eigen::MatrixXd Jij(3,3);
    Jij << J11, J12, J13,
           J21, J22, J23,
           J31, J32, J33;

    return Jij;
}

//Computes the jacobian of the transformation. 
Eigen::MatrixXd 
kin3DHexa8::ComputeDeformationGradientMatrix(const double ri, const double si, const double ti, const Eigen::MatrixXd &Jij) const{
    //Inverse jacobian matrix:
    Eigen::MatrixXd J = Jij.inverse();

    //Gets the element coordinates in deformed configuration. 
    Eigen::VectorXd x1 = theNodes[0]->GetCoordinates() + theNodes[0]->GetDisplacements() + theNodes[0]->GetIncrementalDisplacements();
    Eigen::VectorXd x2 = theNodes[1]->GetCoordinates() + theNodes[1]->GetDisplacements() + theNodes[1]->GetIncrementalDisplacements();
    Eigen::VectorXd x3 = theNodes[2]->GetCoordinates() + theNodes[2]->GetDisplacements() + theNodes[2]->GetIncrementalDisplacements();
    Eigen::VectorXd x4 = theNodes[3]->GetCoordinates() + theNodes[3]->GetDisplacements() + theNodes[3]->GetIncrementalDisplacements();
    Eigen::VectorXd x5 = theNodes[4]->GetCoordinates() + theNodes[4]->GetDisplacements() + theNodes[4]->GetIncrementalDisplacements();
    Eigen::VectorXd x6 = theNodes[5]->GetCoordinates() + theNodes[5]->GetDisplacements() + theNodes[5]->GetIncrementalDisplacements();
    Eigen::VectorXd x7 = theNodes[6]->GetCoordinates() + theNodes[6]->GetDisplacements() + theNodes[6]->GetIncrementalDisplacements();
    Eigen::VectorXd x8 = theNodes[7]->GetCoordinates() + theNodes[7]->GetDisplacements() + theNodes[7]->GetIncrementalDisplacements();

    //Coordinates derivatives with respect to isoparametric coordinates dxi/drk.
    double N11 = -1.0/8.0*(1.0 - si)*(1.0 - ti)*x1(0) + 1.0/8.0*(1.0 - si)*(1.0 - ti)*x2(0) + 1.0/8.0*(1.0 + si)*(1.0 - ti)*x3(0) - 1.0/8.0*(1.0 + si)*(1.0 - ti)*x4(0) - 1.0/8.0*(1.0 - si)*(1.0 + ti)*x5(0) + 1.0/8.0*(1.0 - si)*(1.0 + ti)*x6(0) + 1.0/8.0*(1.0 + si)*(1.0 + ti)*x7(0) - 1.0/8.0*(1.0 + si)*(1.0 + ti)*x8(0);
    double N12 = -1.0/8.0*(1.0 - si)*(1.0 - ti)*x1(1) + 1.0/8.0*(1.0 - si)*(1.0 - ti)*x2(1) + 1.0/8.0*(1.0 + si)*(1.0 - ti)*x3(1) - 1.0/8.0*(1.0 + si)*(1.0 - ti)*x4(1) - 1.0/8.0*(1.0 - si)*(1.0 + ti)*x5(1) + 1.0/8.0*(1.0 - si)*(1.0 + ti)*x6(1) + 1.0/8.0*(1.0 + si)*(1.0 + ti)*x7(1) - 1.0/8.0*(1.0 + si)*(1.0 + ti)*x8(1);
    double N13 = -1.0/8.0*(1.0 - si)*(1.0 - ti)*x1(2) + 1.0/8.0*(1.0 - si)*(1.0 - ti)*x2(2) + 1.0/8.0*(1.0 + si)*(1.0 - ti)*x3(2) - 1.0/8.0*(1.0 + si)*(1.0 - ti)*x4(2) - 1.0/8.0*(1.0 - si)*(1.0 + ti)*x5(2) + 1.0/8.0*(1.0 - si)*(1.0 + ti)*x6(2) + 1.0/8.0*(1.0 + si)*(1.0 + ti)*x7(2) - 1.0/8.0*(1.0 + si)*(1.0 + ti)*x8(2); 
    double N21 = -1.0/8.0*(1.0 - ri)*(1.0 - ti)*x1(0) - 1.0/8.0*(1.0 + ri)*(1.0 - ti)*x2(0) + 1.0/8.0*(1.0 + ri)*(1.0 - ti)*x3(0) + 1.0/8.0*(1.0 - ri)*(1.0 - ti)*x4(0) - 1.0/8.0*(1.0 - ri)*(1.0 + ti)*x5(0) - 1.0/8.0*(1.0 + ri)*(1.0 + ti)*x6(0) + 1.0/8.0*(1.0 + ri)*(1.0 + ti)*x7(0) + 1.0/8.0*(1.0 - ri)*(1.0 + ti)*x8(0);
    double N22 = -1.0/8.0*(1.0 - ri)*(1.0 - ti)*x1(1) - 1.0/8.0*(1.0 + ri)*(1.0 - ti)*x2(1) + 1.0/8.0*(1.0 + ri)*(1.0 - ti)*x3(1) + 1.0/8.0*(1.0 - ri)*(1.0 - ti)*x4(1) - 1.0/8.0*(1.0 - ri)*(1.0 + ti)*x5(1) - 1.0/8.0*(1.0 + ri)*(1.0 + ti)*x6(1) + 1.0/8.0*(1.0 + ri)*(1.0 + ti)*x7(1) + 1.0/8.0*(1.0 - ri)*(1.0 + ti)*x8(1);
    double N23 = -1.0/8.0*(1.0 - ri)*(1.0 - ti)*x1(2) - 1.0/8.0*(1.0 + ri)*(1.0 - ti)*x2(2) + 1.0/8.0*(1.0 + ri)*(1.0 - ti)*x3(2) + 1.0/8.0*(1.0 - ri)*(1.0 - ti)*x4(2) - 1.0/8.0*(1.0 - ri)*(1.0 + ti)*x5(2) - 1.0/8.0*(1.0 + ri)*(1.0 + ti)*x6(2) + 1.0/8.0*(1.0 + ri)*(1.0 + ti)*x7(2) + 1.0/8.0*(1.0 - ri)*(1.0 + ti)*x8(2);
    double N31 = -1.0/8.0*(1.0 - ri)*(1.0 - si)*x1(0) - 1.0/8.0*(1.0 + ri)*(1.0 - si)*x2(0) - 1.0/8.0*(1.0 + ri)*(1.0 + si)*x3(0) - 1.0/8.0*(1.0 - ri)*(1.0 + si)*x4(0) + 1.0/8.0*(1.0 - ri)*(1.0 - si)*x5(0) + 1.0/8.0*(1.0 + ri)*(1.0 - si)*x6(0) + 1.0/8.0*(1.0 + ri)*(1.0 + si)*x7(0) + 1.0/8.0*(1.0 - ri)*(1.0 + si)*x8(0);
    double N32 = -1.0/8.0*(1.0 - ri)*(1.0 - si)*x1(1) - 1.0/8.0*(1.0 + ri)*(1.0 - si)*x2(1) - 1.0/8.0*(1.0 + ri)*(1.0 + si)*x3(1) - 1.0/8.0*(1.0 - ri)*(1.0 + si)*x4(1) + 1.0/8.0*(1.0 - ri)*(1.0 - si)*x5(1) + 1.0/8.0*(1.0 + ri)*(1.0 - si)*x6(1) + 1.0/8.0*(1.0 + ri)*(1.0 + si)*x7(1) + 1.0/8.0*(1.0 - ri)*(1.0 + si)*x8(1);
    double N33 = -1.0/8.0*(1.0 - ri)*(1.0 - si)*x1(2) - 1.0/8.0*(1.0 + ri)*(1.0 - si)*x2(2) - 1.0/8.0*(1.0 + ri)*(1.0 + si)*x3(2) - 1.0/8.0*(1.0 - ri)*(1.0 + si)*x4(2) + 1.0/8.0*(1.0 - ri)*(1.0 - si)*x5(2) + 1.0/8.0*(1.0 + ri)*(1.0 - si)*x6(2) + 1.0/8.0*(1.0 + ri)*(1.0 + si)*x7(2) + 1.0/8.0*(1.0 - ri)*(1.0 + si)*x8(2);

    //Deformation-Gradient coefficients:
    double F11 = N11*J(0,0) + N21*J(0,1) + N31*J(0,2);
    double F12 = N11*J(1,0) + N21*J(1,1) + N31*J(1,2);
    double F13 = N11*J(2,0) + N21*J(2,1) + N31*J(2,2);
    double F21 = N12*J(0,0) + N22*J(0,1) + N32*J(0,2);
    double F22 = N12*J(1,0) + N22*J(1,1) + N32*J(1,2);
    double F23 = N12*J(2,0) + N22*J(2,1) + N32*J(2,2);
    double F31 = N13*J(0,0) + N23*J(0,1) + N33*J(0,2);
    double F32 = N13*J(1,0) + N23*J(1,1) + N33*J(1,2);
    double F33 = N13*J(2,0) + N23*J(2,1) + N33*J(2,2); 

    //Deformation-Gradient Matrix.
    Eigen::MatrixXd Fij(3,3);

    Fij << F11, F12, F13,
           F21, F22, F23,
           F31, F32, F33;

    return Fij;
}

//Compute the Second Piola-Kirchhoff Stress Tensor.
Eigen::MatrixXd 
kin3DHexa8::ComputeSecondPiolaKirchhoffMatrix(const Eigen::VectorXd &stress) const{
    //Second Piola-Kirchhoff tensor.
    Eigen::MatrixXd Stress(9,9);

    Stress << stress(0), stress(3), stress(5), 0.0      , 0.0      , 0.0      , 0.0      , 0.0      , 0.0      ,
              stress(3), stress(1), stress(4), 0.0      , 0.0      , 0.0      , 0.0      , 0.0      , 0.0      ,
              stress(5), stress(4), stress(2), 0.0      , 0.0      , 0.0      , 0.0      , 0.0      , 0.0      ,
              0.0      , 0.0      , 0.0      , stress(0), stress(3), stress(5), 0.0      , 0.0      , 0.0      , 
              0.0      , 0.0      , 0.0      , stress(3), stress(1), stress(4), 0.0      , 0.0      , 0.0      ,
              0.0      , 0.0      , 0.0      , stress(5), stress(4), stress(2), 0.0      , 0.0      , 0.0      ,
              0.0      , 0.0      , 0.0      , 0.0      , 0.0      , 0.0      , stress(0), stress(3), stress(5), 
              0.0      , 0.0      , 0.0      , 0.0      , 0.0      , 0.0      , stress(3), stress(1), stress(4),
              0.0      , 0.0      , 0.0      , 0.0      , 0.0      , 0.0      , stress(5), stress(4), stress(2);

    return Stress;
}

//Compute Shape Function at Gauss Point.
Eigen::MatrixXd 
kin3DHexa8::ComputeShapeFunctionMatrix(const double ri, const double si, const double ti) const{
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

//Evaluates the deformation matrix at a given Gauss point.
Eigen::MatrixXd 
kin3DHexa8::ComputeLinearizedStrainDisplacementMatrix(const double ri, const double si, const double ti, const Eigen::MatrixXd &Jij) const{
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

//Evaluates the Linear (material) deformation matrix at a given Gauss point.
Eigen::MatrixXd 
kin3DHexa8::ComputeLinearStrainDisplacementMatrix(const double ri, const double si, const double ti, const Eigen::MatrixXd &Jij) const{
    //Inverse jacobian matrix:
    Eigen::MatrixXd J = Jij.inverse();

    //Gets the element coordinates in deformed configuration.
    Eigen::VectorXd U1 = theNodes[0]->GetDisplacements() + theNodes[0]->GetIncrementalDisplacements();
    Eigen::VectorXd U2 = theNodes[1]->GetDisplacements() + theNodes[1]->GetIncrementalDisplacements();
    Eigen::VectorXd U3 = theNodes[2]->GetDisplacements() + theNodes[2]->GetIncrementalDisplacements();
    Eigen::VectorXd U4 = theNodes[3]->GetDisplacements() + theNodes[3]->GetIncrementalDisplacements();
    Eigen::VectorXd U5 = theNodes[4]->GetDisplacements() + theNodes[4]->GetIncrementalDisplacements();
    Eigen::VectorXd U6 = theNodes[5]->GetDisplacements() + theNodes[5]->GetIncrementalDisplacements();
    Eigen::VectorXd U7 = theNodes[6]->GetDisplacements() + theNodes[6]->GetIncrementalDisplacements();
    Eigen::VectorXd U8 = theNodes[7]->GetDisplacements() + theNodes[7]->GetIncrementalDisplacements();

    //Strain-displacement matrix coefficients:
    double N11 = -1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 - ti);
    double N21 = -1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 - ti);
    double N31 = -1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 - ti);
    double N41 = -1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 - ti);
    double N51 =  1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 + ti);
    double N61 =  1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 + ti);
    double N71 =  1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 + ti);
    double N81 =  1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 + ti);

    double N12 = -1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 - ti);
    double N22 = -1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 - ti);
    double N32 = -1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 - ti);
    double N42 = -1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 - ti);
    double N52 =  1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 + ti);
    double N62 =  1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 + ti);
    double N72 =  1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 + ti);
    double N82 =  1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 + ti);

    double N13 = -1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 - ti);
    double N23 = -1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 - ti);
    double N33 = -1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 - ti);
    double N43 = -1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 - ti);
    double N53 =  1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 + ti);
    double N63 =  1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 + ti);
    double N73 =  1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 + ti);
    double N83 =  1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 + ti);

    //Displacement derivatives (strains) with respect to undeformed configuration.
    double U11 = N11*U1(0) + N21*U2(0) + N31*U3(0) + N41*U4(0) + N51*U5(0) + N61*U6(0) + N71*U7(0) + N81*U8(0);
    double U21 = N11*U1(1) + N21*U2(1) + N31*U3(1) + N41*U4(1) + N51*U5(1) + N61*U6(1) + N71*U7(1) + N81*U8(1);
    double U31 = N11*U1(2) + N21*U2(2) + N31*U3(2) + N41*U4(2) + N51*U5(2) + N61*U6(2) + N71*U7(2) + N81*U8(2);

    double U12 = N12*U1(0) + N22*U2(0) + N32*U3(0) + N42*U4(0) + N52*U5(0) + N62*U6(0) + N72*U7(0) + N82*U8(0);
    double U22 = N12*U1(1) + N22*U2(1) + N32*U3(1) + N42*U4(1) + N52*U5(1) + N62*U6(1) + N72*U7(1) + N82*U8(1);
    double U32 = N12*U1(2) + N22*U2(2) + N32*U3(2) + N42*U4(2) + N52*U5(2) + N62*U6(2) + N72*U7(2) + N82*U8(2);

    double U13 = N13*U1(0) + N23*U2(0) + N33*U3(0) + N43*U4(0) + N53*U5(0) + N63*U6(0) + N73*U7(0) + N83*U8(0);
    double U23 = N13*U1(1) + N23*U2(1) + N33*U3(1) + N43*U4(1) + N53*U5(1) + N63*U6(1) + N73*U7(1) + N83*U8(1);
    double U33 = N13*U1(2) + N23*U2(2) + N33*U3(2) + N43*U4(2) + N53*U5(2) + N63*U6(2) + N73*U7(2) + N83*U8(2);

    //Linear Strain-displacement matrix coefficients.
    double B11 = N11 + U11*N11;  double B21 =       U12*N12;  double B31 =       U13*N13;  double B41 = N12 + U11*N12 + U12*N11;  double B51 =       U12*N13 + U13*N12;  double B61 = N13 + U11*N13 + U13*N11;
    double B12 =       U21*N11;  double B22 = N12 + U22*N12;  double B32 =       U23*N13;  double B42 = N11 + U21*N12 + U22*N11;  double B52 = N13 + U22*N13 + U23*N12;  double B62 =       U21*N13 + U23*N11;
    double B13 =       U31*N11;  double B23 =       U32*N12;  double B33 = N13 + U33*N13;  double B43 =       U31*N12 + U32*N11;  double B53 = N12 + U32*N13 + U33*N12;  double B63 = N11 + U31*N13 + U33*N11;

    double B14 = N21 + U11*N21;  double B24 =       U12*N22;  double B34 =       U13*N23;  double B44 = N22 + U11*N22 + U12*N21;  double B54 =       U12*N23 + U13*N22;  double B64 = N23 + U11*N23 + U13*N21;
    double B15 =       U21*N21;  double B25 = N22 + U22*N22;  double B35 =       U23*N23;  double B45 = N21 + U21*N22 + U22*N21;  double B55 = N23 + U22*N23 + U23*N22;  double B65 =       U21*N23 + U23*N21;
    double B16 =       U31*N21;  double B26 =       U32*N22;  double B36 = N23 + U33*N23;  double B46 =       U31*N22 + U32*N21;  double B56 = N22 + U32*N23 + U33*N22;  double B66 = N21 + U31*N23 + U33*N21;

    double B17 = N31 + U11*N31;  double B27 =       U12*N32;  double B37 =       U13*N33;  double B47 = N32 + U11*N32 + U12*N31;  double B57 =       U12*N33 + U13*N32;  double B67 = N33 + U11*N33 + U13*N31;
    double B18 =       U21*N31;  double B28 = N32 + U22*N32;  double B38 =       U23*N33;  double B48 = N31 + U21*N32 + U22*N31;  double B58 = N33 + U22*N33 + U23*N32;  double B68 =       U21*N33 + U23*N31;
    double B19 =       U31*N31;  double B29 =       U32*N32;  double B39 = N33 + U33*N33;  double B49 =       U31*N32 + U32*N31;  double B59 = N32 + U32*N33 + U33*N32;  double B69 = N31 + U31*N33 + U33*N31;

    double B110 = N41 + U11*N41; double B210 =       U12*N42; double B310 =       U13*N43; double B410 = N42 + U11*N42 + U12*N41; double B510 =       U12*N43 + U13*N42; double B610 = N43 + U11*N43 + U13*N41;
    double B111 =       U21*N41; double B211 = N42 + U22*N42; double B311 =       U23*N43; double B411 = N41 + U21*N42 + U22*N41; double B511 = N43 + U22*N43 + U23*N42; double B611 =       U21*N43 + U23*N41;
    double B112 =       U31*N41; double B212 =       U32*N42; double B312 = N43 + U33*N43; double B412 =       U31*N42 + U32*N41; double B512 = N42 + U32*N43 + U33*N42; double B612 = N41 + U31*N43 + U33*N41;

    double B113 = N51 + U11*N51; double B213 =       U12*N52; double B313 =       U13*N53; double B413 = N52 + U11*N52 + U12*N51; double B513 =       U12*N53 + U13*N52; double B613 = N53 + U11*N53 + U13*N51;
    double B114 =       U21*N51; double B214 = N52 + U22*N52; double B314 =       U23*N53; double B414 = N51 + U21*N52 + U22*N51; double B514 = N53 + U22*N53 + U23*N52; double B614 =       U21*N53 + U23*N51;
    double B115 =       U31*N51; double B215 =       U32*N52; double B315 = N53 + U33*N53; double B415 =       U31*N52 + U32*N51; double B515 = N52 + U32*N53 + U33*N52; double B615 = N51 + U31*N53 + U33*N51;

    double B116 = N61 + U11*N61; double B216 =       U12*N62; double B316 =       U13*N63; double B416 = N62 + U11*N62 + U12*N61; double B516 =       U12*N63 + U13*N62; double B616 = N63 + U11*N63 + U13*N61;
    double B117 =       U21*N61; double B217 = N62 + U22*N62; double B317 =       U23*N63; double B417 = N61 + U21*N62 + U22*N61; double B517 = N63 + U22*N63 + U23*N62; double B617 =       U21*N63 + U23*N61;
    double B118 =       U31*N61; double B218 =       U32*N62; double B318 = N63 + U33*N63; double B418 =       U31*N62 + U32*N61; double B518 = N62 + U32*N63 + U33*N62; double B618 = N61 + U31*N63 + U33*N61;

    double B119 = N71 + U11*N71; double B219 =       U12*N72; double B319 =       U13*N73; double B419 = N72 + U11*N72 + U12*N71; double B519 =       U12*N73 + U13*N72; double B619 = N73 + U11*N73 + U13*N71;
    double B120 =       U21*N71; double B220 = N72 + U22*N72; double B320 =       U23*N73; double B420 = N71 + U21*N72 + U22*N71; double B520 = N73 + U22*N73 + U23*N72; double B620 =       U21*N73 + U23*N71;
    double B121 =       U31*N71; double B221 =       U32*N72; double B321 = N73 + U33*N73; double B421 =       U31*N72 + U32*N71; double B521 = N72 + U32*N73 + U33*N72; double B621 = N71 + U31*N73 + U33*N71;

    double B122 = N81 + U11*N81; double B222 =       U12*N82; double B322 =       U13*N83; double B422 = N82 + U11*N82 + U12*N81; double B522 =       U12*N83 + U13*N82; double B622 = N83 + U11*N83 + U13*N81;
    double B123 =       U21*N81; double B223 = N82 + U22*N82; double B323 =       U23*N83; double B423 = N81 + U21*N82 + U22*N81; double B523 = N83 + U22*N83 + U23*N82; double B623 =       U21*N83 + U23*N81;
    double B124 =       U31*N81; double B224 =       U32*N82; double B324 = N83 + U33*N83; double B424 =       U31*N82 + U32*N81; double B524 = N82 + U32*N83 + U33*N82; double B624 = N81 + U31*N83 + U33*N81;


    //Shape function matrix:
    Eigen::MatrixXd Bij(6,24);
    Bij <<  B11, B12, B13, B14, B15, B16, B17, B18, B19, B110, B111, B112, B113, B114, B115, B116, B117, B118, B119, B120, B121, B122, B123, B124, 
            B21, B22, B23, B24, B25, B26, B27, B28, B29, B210, B211, B212, B213, B214, B215, B216, B217, B218, B219, B220, B221, B222, B223, B224, 
            B31, B32, B33, B34, B35, B36, B37, B38, B39, B310, B311, B312, B313, B314, B315, B316, B317, B318, B319, B320, B321, B322, B323, B324,
            B41, B42, B43, B44, B45, B46, B47, B48, B49, B410, B411, B412, B413, B414, B415, B416, B417, B418, B419, B420, B421, B422, B423, B424,
            B51, B52, B53, B54, B55, B56, B57, B58, B59, B510, B511, B512, B513, B514, B515, B516, B517, B518, B519, B520, B521, B522, B523, B524,
            B61, B62, B63, B64, B65, B66, B67, B68, B69, B610, B611, B612, B613, B614, B615, B616, B617, B618, B619, B620, B621, B622, B623, B624; 

    return Bij;
}

//Evaluates the Non-Linear (geometric) strain-displacemnt matrix at a given Gauss point.
Eigen::MatrixXd 
kin3DHexa8::ComputeNonLinearStrainDisplacementMatrix(const double ri, const double si, const double ti, const Eigen::MatrixXd &Jij) const{
    //Inverse jacobian matrix.
    Eigen::MatrixXd J = Jij.inverse();

    //Shape functions derivatives with respect to undeformed configuration. 
    double N11 = -1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 - ti);
    double N21 = -1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 - ti);
    double N31 = -1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 - ti);
    double N41 = -1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 - ti);
    double N51 =  1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 + ti);
    double N61 =  1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 + ti);
    double N71 =  1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 + ti);
    double N81 =  1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 + ti);

    double N12 = -1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 - ti);
    double N22 = -1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 - ti);
    double N32 = -1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 - ti);
    double N42 = -1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 - ti);
    double N52 =  1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 + ti);
    double N62 =  1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 + ti);
    double N72 =  1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 + ti);
    double N82 =  1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 + ti);

    double N13 = -1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 - ti);
    double N23 = -1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 - ti);
    double N33 = -1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 - ti);
    double N43 = -1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 - ti);
    double N53 =  1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 + ti);
    double N63 =  1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 + ti);
    double N73 =  1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 + ti);
    double N83 =  1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 + ti);

    //Geometric (non-linear) strain-displacement matrix.
    Eigen::MatrixXd Gij(9,24);
    Gij << N11, 0.0, 0.0, N21, 0.0, 0.0, N31, 0.0, 0.0, N41, 0.0, 0.0, N51, 0.0, 0.0, N61, 0.0, 0.0, N71, 0.0, 0.0, N81, 0.0, 0.0,
           N12, 0.0, 0.0, N22, 0.0, 0.0, N32, 0.0, 0.0, N42, 0.0, 0.0, N52, 0.0, 0.0, N62, 0.0, 0.0, N72, 0.0, 0.0, N82, 0.0, 0.0,
           N13, 0.0, 0.0, N23, 0.0, 0.0, N33, 0.0, 0.0, N43, 0.0, 0.0, N53, 0.0, 0.0, N63, 0.0, 0.0, N73, 0.0, 0.0, N83, 0.0, 0.0,
           0.0, N11, 0.0, 0.0, N21, 0.0, 0.0, N31, 0.0, 0.0, N41, 0.0, 0.0, N51, 0.0, 0.0, N61, 0.0, 0.0, N71, 0.0, 0.0, N81, 0.0,
           0.0, N12, 0.0, 0.0, N22, 0.0, 0.0, N32, 0.0, 0.0, N42, 0.0, 0.0, N52, 0.0, 0.0, N62, 0.0, 0.0, N72, 0.0, 0.0, N82, 0.0,
           0.0, N13, 0.0, 0.0, N23, 0.0, 0.0, N33, 0.0, 0.0, N43, 0.0, 0.0, N53, 0.0, 0.0, N63, 0.0, 0.0, N73, 0.0, 0.0, N83, 0.0,
           0.0, 0.0, N11, 0.0, 0.0, N21, 0.0, 0.0, N31, 0.0, 0.0, N41, 0.0, 0.0, N51, 0.0, 0.0, N61, 0.0, 0.0, N71, 0.0, 0.0, N81,
           0.0, 0.0, N12, 0.0, 0.0, N22, 0.0, 0.0, N32, 0.0, 0.0, N42, 0.0, 0.0, N52, 0.0, 0.0, N62, 0.0, 0.0, N72, 0.0, 0.0, N82,
           0.0, 0.0, N13, 0.0, 0.0, N23, 0.0, 0.0, N33, 0.0, 0.0, N43, 0.0, 0.0, N53, 0.0, 0.0, N63, 0.0, 0.0, N73, 0.0, 0.0, N83;

    return Gij;
}

//Compute the initial stiffness matrix of the element using gauss-integration.
Eigen::MatrixXd 
kin3DHexa8::ComputeInitialStiffnessMatrix() const{
    //Stiffness matrix definition:
    Eigen::MatrixXd StiffnessMatrix(24,24);
    StiffnessMatrix.fill(0.0);

    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Hexa", wi, xi);

    //Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
        //Jacobian matrix.
        Eigen::MatrixXd Jij = ComputeJacobianMatrix(xi(i,0), xi(i,1), xi(i,2));

        //Compute Strain-Displacement Matrix at Gauss Point.
        Eigen::MatrixXd Bij = ComputeLinearizedStrainDisplacementMatrix(xi(i,0), xi(i,1), xi(i,2), Jij);

        //Gets material tangent matrix at Gauss point.
        Eigen::MatrixXd Cij = theMaterial[i]->GetInitialTangentStiffness();

        //Numerical integration.
        StiffnessMatrix += wi(i)*fabs(Jij.determinant())*Bij.transpose()*Cij*Bij;
    }

    return StiffnessMatrix;
}
