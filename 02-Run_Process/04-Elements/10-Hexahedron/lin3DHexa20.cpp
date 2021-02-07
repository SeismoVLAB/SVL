#include <cmath>
#include <iostream>
#include <Eigen/LU> 
#include "Material.hpp"
#include "lin3DHexa20.hpp"
#include "GaussQuadrature.hpp"
#include "LobattoQuadrature.hpp"
#include "Definitions.hpp"
#include "Profiler.hpp"

//Define VTK cell value for Paraview:
const unsigned int VTKCELL = 25;

//Overload constructor.
lin3DHexa20::lin3DHexa20(const std::vector<unsigned int> nodes, std::unique_ptr<Material> &material, const std::string quadrature, const unsigned int nGauss) :
Element("lin3DHexa20", nodes, 60, VTKCELL){
    //The element nodes.
    theNodes.resize(20);

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
lin3DHexa20::~lin3DHexa20(){
    //Does nothing.
}

//Save the material states in the element.
void 
lin3DHexa20::CommitState(){
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
            Eigen::MatrixXd Bij = ComputeStrainDisplacementMatrix(xi(k,0), xi(k,1), xi(k,2), Jij);

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
lin3DHexa20::UpdateState(){
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
lin3DHexa20::SetDomain(std::map<unsigned int, std::shared_ptr<Node> > &nodes){
    //Gets the global element connectivity.
    std::vector<unsigned int> conn = GetNodes();

    //Assign the element to mesh node pointer.  
    for(unsigned int i = 0; i < GetNumberOfNodes(); i++){
        theNodes[i] = nodes[conn[i]];
    }
}

//Sets the damping model.
void 
lin3DHexa20::SetDamping(const std::shared_ptr<Damping> &damping){
    //The damping model
    theDamping = damping;
}

//Gets the list of total-degree of freedom of this element.
std::vector<unsigned int> 
lin3DHexa20::GetTotalDegreeOfFreedom() const{
    //Total number of degree-of-freedom.
    unsigned int nDofs = GetNumberOfDegreeOfFreedom();

    //Reserve memory for the element list of degree-of-freedom.
    std::vector<unsigned int> dofs(nDofs);

    //Construct the element list of degree-of-freedom for assembly.
    for(unsigned int j = 0; j < 20; j++){    
        unsigned int LengthDofs = theNodes[j]->GetNumberOfDegreeOfFreedom();
        std::vector<int> totalDofs = theNodes[j]->GetTotalDegreeOfFreedom();

        for(unsigned int i = 0; i < LengthDofs; i++)
            dofs[i + LengthDofs*j] = totalDofs[i];    
    }

    return dofs;
}

//Returns the material strain at integration points.
Eigen::MatrixXd 
lin3DHexa20::GetStrain() const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    Eigen::MatrixXd theStrain(nPoints,6);
    for(unsigned int k = 0; k < nPoints; k++)
        theStrain.row(k) = theMaterial[k]->GetStrain();

    return theStrain;
}

//Returns the material stress at integration points.
Eigen::MatrixXd 
lin3DHexa20::GetStress() const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    Eigen::MatrixXd theStress(nPoints,6);
    for(unsigned int k = 0; k < nPoints; k++)
        theStress.row(k) = theMaterial[k]->GetTotalStress();

    return theStress;
}

//Returns the material strain-rate at integration points.
Eigen::MatrixXd 
lin3DHexa20::GetStrainRate() const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    Eigen::MatrixXd theStrainRate(nPoints,6);
    for(unsigned int k = 0; k < nPoints; k++)
        theStrainRate.row(k) = theMaterial[k]->GetStrainRate();

    return theStrainRate;
}

//Gets the material strain in section at  coordinate (x3,x2).
Eigen::MatrixXd 
lin3DHexa20::GetStrainAt(double UNUSED(x3), double UNUSED(x2)) const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    //Stress at coordinate is define within section.
    Eigen::MatrixXd theStrain(nPoints, 6);
    theStrain.fill(0.0);

    return theStrain;
}

//Gets the material stress in section at  coordinate (x3,x2).
Eigen::MatrixXd 
lin3DHexa20::GetStressAt(double UNUSED(x3), double UNUSED(x2)) const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    //Stress at coordinate is define within section.
    Eigen::MatrixXd theStress(nPoints, 6);
    theStress.fill(0.0);

    return theStress;
}

//Gets the element internal response in VTK format.
Eigen::VectorXd 
lin3DHexa20::GetVTKResponse(std::string response) const{
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
lin3DHexa20::ComputeEnergy(){
    //TODO: Integrate over element volume to compute the energy
    return 0.0;
}

//Compute the mass matrix of the element using gauss-integration.
Eigen::MatrixXd 
lin3DHexa20::ComputeMassMatrix(){
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

    //Use consistent mass definition:
    Eigen::MatrixXd MassMatrix(60, 60);
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

    //Lumped Mass Formulation
    if(MassFormulation){
        //Lumped Mass in diagonal terms.
        double TotalMass = MassMatrix.sum()/3.00;
        double TraceMass = MassMatrix.trace()/3.00;
        for (unsigned int i = 0; i < 60; i++){
            for (unsigned int j = 0; j < 60; j++){
                if(i == j)
                    MassMatrix(i,i) = MassMatrix(i,j)*TotalMass/TraceMass;
                else
                    MassMatrix(i,j) = 0.0;
            }
        }
    }

    return MassMatrix;
}

//Compute the stiffness matrix of the element using gauss-integration.
Eigen::MatrixXd 
lin3DHexa20::ComputeStiffnessMatrix(){
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

    //Stiffness matrix definition:
    Eigen::MatrixXd StiffnessMatrix(60, 60);
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
        Eigen::MatrixXd Bij = ComputeStrainDisplacementMatrix(xi(i,0), xi(i,1), xi(i,2), Jij);

        //Gets material tangent matrix at Gauss point.
        Eigen::MatrixXd Cij = theMaterial[i]->GetTangentStiffness();

        //Numerical integration.
        StiffnessMatrix += wi(i)*fabs(Jij.determinant())*Bij.transpose()*Cij*Bij;
    }

    return StiffnessMatrix;
}

//Compute the damping matrix of the element using gauss-integration.
Eigen::MatrixXd 
lin3DHexa20::ComputeDampingMatrix(){
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

    //Damping matrix definition
    Eigen::MatrixXd DampingMatrix(60,60);
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
            Eigen::MatrixXd Bij = ComputeStrainDisplacementMatrix(xi(i,0), xi(i,1), xi(i,2), Jij);
            
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
lin3DHexa20::ComputePMLMatrix(){
    Eigen::MatrixXd Kpml;
    return Kpml;
}

//Compute the internal forces acting on the element.
Eigen::VectorXd 
lin3DHexa20::ComputeInternalForces(){
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

    //Stiffness matrix definition.
    Eigen::VectorXd InternalForces(60);
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
        Eigen::MatrixXd Bij = ComputeStrainDisplacementMatrix(xi(i,0), xi(i,1), xi(i,2), Jij);

        //Gets material strain at Gauss point.
        Eigen::VectorXd Stress = theMaterial[i]->GetStress();

        //Numerical integration.
        InternalForces += wi(i)*fabs(Jij.determinant())*Bij.transpose()*Stress;
    }
    
    return InternalForces;
}

//Compute the elastic, inertial, and vicous forces acting on the element.
Eigen::VectorXd 
lin3DHexa20::ComputeInternalDynamicForces(){
    //The Internal dynamic force vector
    Eigen::VectorXd InternalForces;

    if( HasFixedNode(theNodes) ){
        //Allocate memory for velocity/acceleraton. 
        Eigen::VectorXd V(60); 
        Eigen::VectorXd A(60);

        //Fills the response vectors with velocity/acceleraton values.
        V << theNodes[ 0]->GetVelocities(), theNodes[ 1]->GetVelocities(), theNodes[ 2]->GetVelocities(), theNodes[ 3]->GetVelocities(), theNodes[ 4]->GetVelocities(),
             theNodes[ 5]->GetVelocities(), theNodes[ 6]->GetVelocities(), theNodes[ 7]->GetVelocities(), theNodes[ 8]->GetVelocities(), theNodes[ 9]->GetVelocities(),
             theNodes[10]->GetVelocities(), theNodes[11]->GetVelocities(), theNodes[12]->GetVelocities(), theNodes[13]->GetVelocities(), theNodes[14]->GetVelocities(),
             theNodes[15]->GetVelocities(), theNodes[16]->GetVelocities(), theNodes[17]->GetVelocities(), theNodes[18]->GetVelocities(), theNodes[19]->GetVelocities();

        A << theNodes[ 0]->GetAccelerations(), theNodes[ 1]->GetAccelerations(), theNodes[ 2]->GetAccelerations(), theNodes[ 3]->GetAccelerations(), theNodes[ 4]->GetAccelerations(),
             theNodes[ 5]->GetAccelerations(), theNodes[ 6]->GetAccelerations(), theNodes[ 7]->GetAccelerations(), theNodes[ 8]->GetAccelerations(), theNodes[ 9]->GetAccelerations(),
             theNodes[10]->GetAccelerations(), theNodes[11]->GetAccelerations(), theNodes[12]->GetAccelerations(), theNodes[13]->GetAccelerations(), theNodes[14]->GetAccelerations(),
             theNodes[15]->GetAccelerations(), theNodes[16]->GetAccelerations(), theNodes[17]->GetAccelerations(), theNodes[18]->GetAccelerations(), theNodes[19]->GetAccelerations();

        //Compute the inertial/viscous/elastic dynamic force contribution.
        InternalForces = ComputeInternalForces() + ComputeDampingMatrix()*V + ComputeMassMatrix()*A;
    }

    return InternalForces;
}

//Compute the surface forces acting on the element.
Eigen::VectorXd 
lin3DHexa20::ComputeSurfaceForces(const std::shared_ptr<Load> &surfaceLoad, unsigned int face){
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

    //Local surface load vector:
    Eigen::VectorXd surfaceForces(60);
    surfaceForces.fill(0.0);

    //Gets the surface load:
    Eigen::VectorXd qs = surfaceLoad->GetLoadVector();

    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Quad", wi, xi);

    if(face == 1){
        //Gets the face coordinates in undeformed configuration. 
        Eigen::VectorXd x1  = theNodes[0]->GetCoordinates();
        Eigen::VectorXd x2  = theNodes[1]->GetCoordinates();
        Eigen::VectorXd x3  = theNodes[2]->GetCoordinates();
        Eigen::VectorXd x4  = theNodes[3]->GetCoordinates();
        Eigen::VectorXd x9  = theNodes[8]->GetCoordinates();
        Eigen::VectorXd x10 = theNodes[9]->GetCoordinates();
        Eigen::VectorXd x11 = theNodes[10]->GetCoordinates();
        Eigen::VectorXd x12 = theNodes[11]->GetCoordinates();

        //Numerical integration in local axis r-s:
        for(unsigned int i = 0; i < wi.size(); i++){
            double r = xi(i,0);
            double s = xi(i,1);

            //vectors along s and t axes.
            Eigen::Vector3d v1, v2;
            v1 << 1.0/4.0*(1.0 - s)*(2.0*r + s)*x1 + 1.0/4.0*(1.0 - s)*(2.0*r - s)*x2 + 1.0/4.0*(1.0 + s)*(2.0*r + s)*x3 + 1.0/4.0*(1.0 + s)*(2.0*r - s)*x4 - r*(1.0 - s)*x9 + 1.0/2.0*(1.0 - s*s)*x10 - r*(1.0 + s)*x11 - 1.0/2.0*(1.0 - s*s)*x12;
            v2 << 1.0/4.0*(1.0 - r)*(r + 2.0*s)*x1 - 1.0/4.0*(1.0 + r)*(r - 2.0*s)*x2 + 1.0/4.0*(1.0 + r)*(r + 2.0*s)*x3 - 1.0/4.0*(1.0 - r)*(r - 2.0*s)*x4 - 1.0/2.0*(1.0 - r*r)*x9 - (1.0 + r)*s*x10 + 1.0/2.0*(1.0 - r*r)*x11 - (1.0 - r)*s*x12;

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
        Eigen::VectorXd x1  = theNodes[0]->GetCoordinates();
        Eigen::VectorXd x2  = theNodes[1]->GetCoordinates();
        Eigen::VectorXd x6  = theNodes[5]->GetCoordinates();
        Eigen::VectorXd x5  = theNodes[4]->GetCoordinates();
        Eigen::VectorXd x9  = theNodes[8]->GetCoordinates();
        Eigen::VectorXd x18 = theNodes[17]->GetCoordinates();
        Eigen::VectorXd x13 = theNodes[12]->GetCoordinates();
        Eigen::VectorXd x17 = theNodes[16]->GetCoordinates();

        //Numerical integration in local axis r-t:
        for(unsigned int i = 0; i < wi.size(); i++){
            double r = xi(i,0);
            double t = xi(i,1);

            //vectors along s and t axes.
            Eigen::Vector3d v1, v2;
            v1 << 1.0/4.0*(1.0 - t)*(2.0*r + t)*x1 + 1.0/4.0*(1.0 - t)*(2.0*r - t)*x2 + 1.0/4.0*(1.0 + t)*(2.0*r - t)*x5 + 1.0/4.0*(1.0 + t)*(2.0*r + t)*x6 - r*(1.0 - t)*x9 - r*(1.0 + t)*x13 - 1.0/2.0*(1.0 - t*t)*x17 + 1.0/2.0*(1.0 - t*t)*x18;
            v2 << 1.0/4.0*(1.0 - r)*(r + 2.0*t)*x1 - 1.0/4.0*(1.0 + r)*(r - 2.0*t)*x2 - 1.0/4.0*(1.0 - r)*(r - 2.0*t)*x5 + 1.0/4.0*(1.0 + r)*(r + 2.0*t)*x6 - 1.0/2.0*(1.0 - r*r)*x9 + 1.0/2.0*(1.0 - r*r)*x13 - (1.0 - r)*t*x17 - (1.0 + r)*t*x18;

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
        Eigen::VectorXd x2  = theNodes[1]->GetCoordinates();
        Eigen::VectorXd x3  = theNodes[2]->GetCoordinates();
        Eigen::VectorXd x7  = theNodes[6]->GetCoordinates();
        Eigen::VectorXd x6  = theNodes[5]->GetCoordinates();
        Eigen::VectorXd x10 = theNodes[9]->GetCoordinates();
        Eigen::VectorXd x19 = theNodes[18]->GetCoordinates();
        Eigen::VectorXd x14 = theNodes[13]->GetCoordinates();
        Eigen::VectorXd x18 = theNodes[17]->GetCoordinates();

        //Numerical integration in local axis s-t:
        for(unsigned int i = 0; i < wi.size(); i++){
            double s = xi(i,0);
            double t = xi(i,1);

            //vectors along s and t axes.
            Eigen::Vector3d v1, v2;
            v1 << 1.0/4.0*(1.0 - t)*(2.0*s + t)*x2 + 1.0/4.0*(1.0 - t)*(2.0*s - t)*x3 + 1.0/4.0*(1.0 + t)*(2.0*s - t)*x6 + 1.0/4.0*(1.0 + t)*(2.0*s + t)*x7 - s*(1.0 - t)*x10 - s*(1.0 + t)*x14 - 1.0/2.0*(1.0 - t*t)*x18 + 1.0/2.0*(1.0 - t*t)*x19;
            v2 << 1.0/4.0*(1.0 - s)*(s + 2.0*t)*x2 - 1.0/4.0*(1.0 + s)*(s - 2.0*t)*x3 - 1.0/4.0*(1.0 - s)*(s - 2.0*t)*x6 + 1.0/4.0*(1.0 + s)*(s + 2.0*t)*x7 - 1.0/2.0*(1.0 - s*s)*x10 + 1.0/2.0*(1.0 - s*s)*x14 - (1.0 - s)*t*x18 - (1.0 + s)*t*x19;

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
        Eigen::VectorXd x4  = theNodes[3]->GetCoordinates();
        Eigen::VectorXd x3  = theNodes[2]->GetCoordinates();
        Eigen::VectorXd x7  = theNodes[6]->GetCoordinates();
        Eigen::VectorXd x8  = theNodes[7]->GetCoordinates();
        Eigen::VectorXd x11 = theNodes[10]->GetCoordinates();
        Eigen::VectorXd x19 = theNodes[18]->GetCoordinates();
        Eigen::VectorXd x15 = theNodes[14]->GetCoordinates();
        Eigen::VectorXd x20 = theNodes[19]->GetCoordinates();

        //Numerical integration in local axis r-t:
        for(unsigned int i = 0; i < wi.size(); i++){
            double r = xi(i,0);
            double t = xi(i,1);

            //vectors along r and t axes.
            Eigen::Vector3d v1, v2;
            v1 <<  1.0/4.0*(1.0 - t)*(2.0*r - t)*x3 + 1.0/4.0*(1.0 - t)*(2.0*r + t)*x4 + 1.0/4.0*(1.0 + t)*(2.0*r + t)*x7 + 1.0/4.0*(1.0 + t)*(2.0*r - t)*x8 - r*(1.0 - t)*x11 - r*(1.0 + t)*x15 + 1.0/2.0*(1.0 - t*t)*x19 - 1.0/2.0*(1.0 - t*t)*x20;
            v2 << -1.0/4.0*(1.0 + r)*(r - 2.0*t)*x3 + 1.0/4.0*(1.0 - r)*(r + 2.0*t)*x4 + 1.0/4.0*(1.0 + r)*(r + 2.0*t)*x7 - 1.0/4.0*(1.0 - r)*(r - 2.0*t)*x8 - 1.0/2.0*(1.0 - r*r)*x11 + 1.0/2.0*(1.0 - r*r)*x15 - (1.0 + r)*t*x19 - (1.0 - r)*t*x20;

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
        Eigen::VectorXd x1  = theNodes[0]->GetCoordinates();
        Eigen::VectorXd x4  = theNodes[3]->GetCoordinates();
        Eigen::VectorXd x8  = theNodes[7]->GetCoordinates();
        Eigen::VectorXd x5  = theNodes[4]->GetCoordinates();
        Eigen::VectorXd x12 = theNodes[11]->GetCoordinates();
        Eigen::VectorXd x20 = theNodes[19]->GetCoordinates();
        Eigen::VectorXd x16 = theNodes[15]->GetCoordinates();
        Eigen::VectorXd x17 = theNodes[16]->GetCoordinates();

        //Numerical integration in local axis s-t:
        for(unsigned int i = 0; i < wi.size(); i++){
            double s = xi(i,0);
            double t = xi(i,1);

            //vectors along s and t axes.
            Eigen::Vector3d v1, v2;
            v1 << 1.0/4.0*(1.0 - t)*(2.0*s + t)*x1 + 1.0/4.0*(1.0 - t)*(2.0*s - t)*x4 + 1.0/4.0*(1.0 + t)*(2.0*s - t)*x5 + 1.0/4.0*(1.0 + t)*(2.0*s + t)*x8 - s*(1.0 - t)*x12 - s*(1.0 + t)*x16 - 1.0/2.0*(1.0 - t*t)*x17 + 1.0/2.0*(1.0 - t*t)*x20;
            v2 << 1.0/4.0*(1.0 - s)*(s + 2.0*t)*x1 - 1.0/4.0*(1.0 + s)*(s - 2.0*t)*x4 - 1.0/4.0*(1.0 - s)*(s - 2.0*t)*x5 + 1.0/4.0*(1.0 + s)*(s + 2.0*t)*x8 - 1.0/2.0*(1.0 - s*s)*x12 + 1.0/2.0*(1.0 - s*s)*x16 - (1.0 - s)*t*x17 - (1.0 + s)*t*x20;

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
        Eigen::VectorXd x5  = theNodes[4]->GetCoordinates();
        Eigen::VectorXd x6  = theNodes[5]->GetCoordinates();
        Eigen::VectorXd x7  = theNodes[6]->GetCoordinates();
        Eigen::VectorXd x8  = theNodes[7]->GetCoordinates();
        Eigen::VectorXd x13 = theNodes[12]->GetCoordinates();
        Eigen::VectorXd x14 = theNodes[13]->GetCoordinates();
        Eigen::VectorXd x15 = theNodes[14]->GetCoordinates();
        Eigen::VectorXd x16 = theNodes[15]->GetCoordinates();

        //Numerical integration in local axis r-s:
        for(unsigned int i = 0; i < wi.size(); i++){
            double r = xi(i,0);
            double s = xi(i,1);

            //vectors along r and s axes.
            Eigen::Vector3d v1, v2;
            v1 << 1.0/4.0*(1.0 - s)*(2.0*r + s)*x5 + 1.0/4.0*(1.0 - s)*(2.0*r - s)*x6 + 1.0/4.0*(1.0 + s)*(2.0*r + s)*x7 + 1.0/4.0*(1.0 + s)*(2.0*r - s)*x8 - r*(1.0 - s)*x13 + 1.0/2.0*(1.0 - s*s)*x14 - r*(1.0 + s)*x15 - 1.0/2.0*(1.0 - s*s)*x16;
            v2 << 1.0/4.0*(1.0 - r)*(r + 2.0*s)*x5 - 1.0/4.0*(1.0 + r)*(r - 2.0*s)*x6 + 1.0/4.0*(1.0 + r)*(r + 2.0*s)*x7 - 1.0/4.0*(1.0 - r)*(r - 2.0*s)*x8 - 1.0/2.0*(1.0 - r*r)*x13 - (1.0 + r)*s*x14 + 1.0/2.0*(1.0 - r*r)*x15 - (1.0 - r)*s*x16;

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
lin3DHexa20::ComputeBodyForces(const std::shared_ptr<Load> &bodyLoad, unsigned int k){
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

    //Local body load vector.
    Eigen::VectorXd bodyForces(60);
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
lin3DHexa20::ComputeDomainReductionForces(const std::shared_ptr<Load> &drm, unsigned int k){
    //Starts profiling this funtion.
    PROFILE_FUNCTION();

    //Local domain-reduction load vector.
    Eigen::VectorXd DRMForces(60);

    //Get the Domain-Reduction field motion.
    Eigen::VectorXd x1  = theNodes[0]->GetDomainReductionMotion(k);
    Eigen::VectorXd x2  = theNodes[1]->GetDomainReductionMotion(k);
    Eigen::VectorXd x3  = theNodes[2]->GetDomainReductionMotion(k);
    Eigen::VectorXd x4  = theNodes[3]->GetDomainReductionMotion(k);
    Eigen::VectorXd x5  = theNodes[4]->GetDomainReductionMotion(k);
    Eigen::VectorXd x6  = theNodes[5]->GetDomainReductionMotion(k);
    Eigen::VectorXd x7  = theNodes[6]->GetDomainReductionMotion(k);
    Eigen::VectorXd x8  = theNodes[7]->GetDomainReductionMotion(k);
    Eigen::VectorXd x9  = theNodes[8]->GetDomainReductionMotion(k);
    Eigen::VectorXd x10 = theNodes[9]->GetDomainReductionMotion(k);
    Eigen::VectorXd x11 = theNodes[10]->GetDomainReductionMotion(k);
    Eigen::VectorXd x12 = theNodes[11]->GetDomainReductionMotion(k);
    Eigen::VectorXd x13 = theNodes[12]->GetDomainReductionMotion(k);
    Eigen::VectorXd x14 = theNodes[13]->GetDomainReductionMotion(k);
    Eigen::VectorXd x15 = theNodes[14]->GetDomainReductionMotion(k);
    Eigen::VectorXd x16 = theNodes[15]->GetDomainReductionMotion(k);
    Eigen::VectorXd x17 = theNodes[16]->GetDomainReductionMotion(k);
    Eigen::VectorXd x18 = theNodes[17]->GetDomainReductionMotion(k);
    Eigen::VectorXd x19 = theNodes[18]->GetDomainReductionMotion(k);
    Eigen::VectorXd x20 = theNodes[19]->GetDomainReductionMotion(k);

    //Constructs the domain reduction boundary/exterior connectivity.
    std::vector<bool> DRMcond(60);
    std::vector<unsigned int> conn = GetNodes();

    for(unsigned int i = 0; i < conn.size(); i ++){
        bool condition = drm->GetDRMCondition(conn[i]);
        DRMcond[3*i  ] = condition;
        DRMcond[3*i+1] = condition;
        DRMcond[3*i+2] = condition;
    }

    //Constructs the displacement, velocity and acceleration vectors. 
    Eigen::VectorXd Uo(60); 
    Eigen::VectorXd Vo(60);
    Eigen::VectorXd Ao(60);
 
    Uo << x1(0), x1(1), x1(2), x2(0), x2(1), x2(2), x3(0), x3(1), x3(2), x4(0), x4(1), x4(2), x5(0), x5(1), x5(2), x6(0), x6(1), x6(2), x7(0), x7(1), x7(2), x8(0), x8(1), x8(2), x9(0), x9(1), x9(2), x10(0), x10(1), x10(2), x11(0), x11(1), x11(2), x12(0), x12(1), x12(2), x13(0), x13(1), x13(2), x14(0), x14(1), x14(2), x15(0), x15(1), x15(2), x16(0), x16(1), x16(2), x17(0), x17(1), x17(2), x18(0), x18(1), x18(2), x19(0), x19(1), x19(2), x20(0), x20(1), x20(2);
    Vo << x1(3), x1(4), x1(5), x2(3), x2(4), x2(5), x3(3), x3(4), x3(5), x4(3), x4(4), x4(5), x5(3), x5(4), x5(5), x6(3), x6(4), x6(5), x7(3), x7(4), x7(5), x8(3), x8(4), x8(5), x9(3), x9(4), x9(5), x10(3), x10(4), x10(5), x11(3), x11(4), x11(5), x12(3), x12(4), x12(5), x13(3), x13(4), x13(5), x14(3), x14(4), x14(5), x15(3), x15(4), x15(5), x16(3), x16(4), x16(5), x17(3), x17(4), x17(5), x18(3), x18(4), x18(5), x19(3), x19(4), x19(5), x20(3), x20(4), x20(5);
    Ao << x1(6), x1(7), x1(8), x2(6), x2(7), x2(8), x3(6), x3(7), x3(8), x4(6), x4(7), x4(8), x5(6), x5(7), x5(8), x6(6), x6(7), x6(8), x7(6), x7(7), x7(8), x8(6), x8(7), x8(8), x9(6), x9(7), x9(8), x10(6), x10(7), x10(8), x11(6), x11(7), x11(8), x12(6), x12(7), x12(8), x13(6), x13(7), x13(8), x14(6), x14(7), x14(8), x15(6), x15(7), x15(8), x16(6), x16(7), x16(8), x17(6), x17(7), x17(8), x18(6), x18(7), x18(8), x19(6), x19(7), x19(8), x20(6), x20(7), x20(8);

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
lin3DHexa20::ComputeStrain(const Eigen::MatrixXd &Bij) const{
    //Gets the element coordinates in deformed configuration. 
    Eigen::VectorXd U1  = theNodes[0]->GetDisplacements()  + theNodes[0]->GetIncrementalDisplacements();
    Eigen::VectorXd U2  = theNodes[1]->GetDisplacements()  + theNodes[1]->GetIncrementalDisplacements();
    Eigen::VectorXd U3  = theNodes[2]->GetDisplacements()  + theNodes[2]->GetIncrementalDisplacements();
    Eigen::VectorXd U4  = theNodes[3]->GetDisplacements()  + theNodes[3]->GetIncrementalDisplacements();
    Eigen::VectorXd U5  = theNodes[4]->GetDisplacements()  + theNodes[4]->GetIncrementalDisplacements();
    Eigen::VectorXd U6  = theNodes[5]->GetDisplacements()  + theNodes[5]->GetIncrementalDisplacements();
    Eigen::VectorXd U7  = theNodes[6]->GetDisplacements()  + theNodes[6]->GetIncrementalDisplacements();
    Eigen::VectorXd U8  = theNodes[7]->GetDisplacements()  + theNodes[7]->GetIncrementalDisplacements();
    Eigen::VectorXd U9  = theNodes[8]->GetDisplacements()  + theNodes[8]->GetIncrementalDisplacements();
    Eigen::VectorXd U10 = theNodes[9]->GetDisplacements()  + theNodes[9]->GetIncrementalDisplacements();
    Eigen::VectorXd U11 = theNodes[10]->GetDisplacements() + theNodes[10]->GetIncrementalDisplacements();
    Eigen::VectorXd U12 = theNodes[11]->GetDisplacements() + theNodes[11]->GetIncrementalDisplacements();
    Eigen::VectorXd U13 = theNodes[12]->GetDisplacements() + theNodes[12]->GetIncrementalDisplacements();
    Eigen::VectorXd U14 = theNodes[13]->GetDisplacements() + theNodes[13]->GetIncrementalDisplacements();
    Eigen::VectorXd U15 = theNodes[14]->GetDisplacements() + theNodes[14]->GetIncrementalDisplacements();
    Eigen::VectorXd U16 = theNodes[15]->GetDisplacements() + theNodes[15]->GetIncrementalDisplacements();
    Eigen::VectorXd U17 = theNodes[16]->GetDisplacements() + theNodes[16]->GetIncrementalDisplacements();
    Eigen::VectorXd U18 = theNodes[17]->GetDisplacements() + theNodes[17]->GetIncrementalDisplacements();
    Eigen::VectorXd U19 = theNodes[18]->GetDisplacements() + theNodes[18]->GetIncrementalDisplacements();
    Eigen::VectorXd U20 = theNodes[19]->GetDisplacements() + theNodes[19]->GetIncrementalDisplacements();

    Eigen::VectorXd nodalDisplacement(60);
    nodalDisplacement << U1, U2, U3, U4, U5, U6, U7, U8, U9, U10, U11, U12, U13, U14, U15, U16, U17, U18, U19, U20;

    //Strain vector:
    Eigen::VectorXd Strain(6); 
    Strain = Bij*nodalDisplacement;

    return Strain;
}

//Update strain rate in the element.
Eigen::VectorXd 
lin3DHexa20::ComputeStrainRate(const Eigen::MatrixXd& UNUSED(Bij)) const{
     //TODO: Compute strain rate.
    //Strain vector definition:
    Eigen::VectorXd strainrate(6);
    strainrate.fill(0.0);

    return strainrate;
}

//Computes the jacobian of the transformation. 
Eigen::MatrixXd 
lin3DHexa20::ComputeJacobianMatrix(const double ri, const double si, const double ti) const{
    //Gets the element coordinates in deformed configuration. 
    Eigen::VectorXd X1  = theNodes[0]->GetCoordinates();
    Eigen::VectorXd X2  = theNodes[1]->GetCoordinates();
    Eigen::VectorXd X3  = theNodes[2]->GetCoordinates();
    Eigen::VectorXd X4  = theNodes[3]->GetCoordinates();
    Eigen::VectorXd X5  = theNodes[4]->GetCoordinates();
    Eigen::VectorXd X6  = theNodes[5]->GetCoordinates();
    Eigen::VectorXd X7  = theNodes[6]->GetCoordinates();
    Eigen::VectorXd X8  = theNodes[7]->GetCoordinates();
    Eigen::VectorXd X9  = theNodes[8]->GetCoordinates();
    Eigen::VectorXd X10 = theNodes[9]->GetCoordinates();
    Eigen::VectorXd X11 = theNodes[10]->GetCoordinates();
    Eigen::VectorXd X12 = theNodes[11]->GetCoordinates();
    Eigen::VectorXd X13 = theNodes[12]->GetCoordinates();
    Eigen::VectorXd X14 = theNodes[13]->GetCoordinates();
    Eigen::VectorXd X15 = theNodes[14]->GetCoordinates();
    Eigen::VectorXd X16 = theNodes[15]->GetCoordinates();
    Eigen::VectorXd X17 = theNodes[16]->GetCoordinates();
    Eigen::VectorXd X18 = theNodes[17]->GetCoordinates();
    Eigen::VectorXd X19 = theNodes[18]->GetCoordinates();
    Eigen::VectorXd X20 = theNodes[19]->GetCoordinates();
 
    //Shape function derivatives:
    //(a) Corner Nodes.
    double dN1r =  1.0/8.0*(-1.0 + si)*(-1.0 + ti)*( 1.0 + 2.0*ri + si + ti);
    double dN1s =  1.0/8.0*(-1.0 + ri)*(-1.0 + ti)*( 1.0 + ri + 2.0*si + ti);
    double dN1t =  1.0/8.0*(-1.0 + ri)*(-1.0 + si)*( 1.0 + ri + si + 2.0*ti);
    double dN2r = -1.0/8.0*(-1.0 + si)*(-1.0 + ti)*( 1.0 - 2.0*ri + si + ti);
    double dN2s =  1.0/8.0*( 1.0 + ri)*(-1.0 + ti)*(-1.0 + ri - 2.0*si - ti);
    double dN2t =  1.0/8.0*( 1.0 + ri)*(-1.0 + si)*(-1.0 + ri - si - 2.0*ti);
    double dN3r = -1.0/8.0*(1.0 + si)*(-1.0 + ti)*(-1.0 + 2.0*ri + si - ti);
    double dN3s = -1.0/8.0*(1.0 + ri)*(-1.0 + ti)*(-1.0 + ri + 2.0*si - ti);
    double dN3t = -1.0/8.0*(1.0 + ri)*( 1.0 + si)*(-1.0 + ri + si - 2.0*ti);
    double dN4r =  1.0/8.0*( 1.0 + si)*(-1.0 + ti)*(-1.0 - 2.0*ri + si - ti);
    double dN4s = -1.0/8.0*(-1.0 + ri)*(-1.0 + ti)*( 1.0 + ri - 2.0*si + ti);
    double dN4t = -1.0/8.0*(-1.0 + ri)*( 1.0 + si)*( 1.0 + ri - si + 2.0*ti);
    double dN5r = -1.0/8.0*(-1.0 + si)*( 1.0 + ti)*(1.0 + 2.0*ri + si - ti);
    double dN5s = -1.0/8.0*(-1.0 + ri)*( 1.0 + ti)*(1.0 + ri + 2.0*si - ti);
    double dN5t = -1.0/8.0*(-1.0 + ri)*(-1.0 + si)*(1.0 + ri + si - 2.0*ti);
    double dN6r =  1.0/8.0*(-1.0 + si)*( 1.0 + ti)*( 1.0 - 2.0*ri + si - ti);
    double dN6s = -1.0/8.0*( 1.0 + ri)*( 1.0 + ti)*(-1.0 + ri - 2.0*si + ti);
    double dN6t = -1.0/8.0*( 1.0 + ri)*(-1.0 + si)*(-1.0 + ri - si + 2.0*ti);
    double dN7r =  1.0/8.0*(1.0 + si)*(1.0 + ti)*(-1.0 + 2.0*ri + si + ti);
    double dN7s =  1.0/8.0*(1.0 + ri)*(1.0 + ti)*(-1.0 + ri + 2.0*si + ti);
    double dN7t =  1.0/8.0*(1.0 + ri)*(1.0 + si)*(-1.0 + ri + si + 2.0*ti);
    double dN8r = -1.0/8.0*(1.0 + si)*(1.0 + ti)*(-1.0 - 2.0*ri + si + ti);
    double dN8s =  1.0/8.0*(-1.0 + ri)*(1.0 + ti)*( 1.0 + ri - 2.0*si - ti);
    double dN8t =  1.0/8.0*(-1.0 + ri)*(1.0 + si)*( 1.0 + ri - si - 2.0*ti);

    //(b) Side Nodes.
    double dN9r = -1.0/2.0*ri*(-1.0 + si)*(-1.0 + ti);
    double dN9s = -1.0/4.0*(-1.0 + ri*ri)*(-1.0 + ti);
    double dN9t = -1.0/4.0*(-1.0 + ri*ri)*(-1.0 + si);
    double dN10r =  1.0/4.0*(-1.0 + si*si)*(-1.0 + ti);
    double dN10s =  1.0/2.0*( 1.0 + ri)*si*(-1.0 + ti);
    double dN10t =  1.0/4.0*( 1.0 + ri)*(-1.0 + si*si);
    double dN11r =  1.0/2.0*ri*(1.0 + si)*(-1.0 + ti);
    double dN11s =  1.0/4.0*(-1.0 + ri*ri)*(-1.0 + ti);
    double dN11t =  1.0/4.0*(-1.0 + ri*ri)*(1.0 + si);
    double dN12r = -1.0/4.0*(-1.0 + si*si)*(-1.0 + ti);
    double dN12s = -1.0/2.0*(-1.0 + ri)*si*(-1.0 + ti);
    double dN12t = -1.0/4.0*(-1.0 + ri)*(-1.0 + si*si);
    double dN13r =  1.0/2.0*ri*(-1.0 + si)*(1.0 + ti);
    double dN13s =  1.0/4.0*(-1.0 + ri*ri)*(1.0 + ti);
    double dN13t =  1.0/4.0*(-1.0 + ri*ri)*(-1.0 + si);
    double dN14r = -1.0/4.0*(-1.0 + si*si)*(1.0 + ti);
    double dN14s = -1.0/2.0*(1.0 + ri)*si*(1.0 + ti);
    double dN14t = -1.0/4.0*(1.0 + ri)*(-1.0 + si*si);
    double dN15r = -1.0/2.0*ri*(1.0 + si)*(1.0 + ti);
    double dN15s = -1.0/4.0*(-1.0 + ri*ri)*(1.0 + ti);
    double dN15t = -1.0/4.0*(-1.0 + ri*ri)*(1.0 + si);
    double dN16r =  1.0/4.0*(-1.0 + si*si)*(1.0 + ti);
    double dN16s =  1.0/2.0*(-1.0 + ri)*si*(1.0 + ti);
    double dN16t =  1.0/4.0*(-1.0 + ri)*(-1.0 + si*si);
    double dN17r = -1.0/4.0*(-1.0 + si)*(-1.0 + ti*ti);
    double dN17s = -1.0/4.0*(-1.0 + ri)*(-1.0 + ti*ti);
    double dN17t = -1.0/2.0*(-1.0 + ri)*(-1.0 + si)*ti;
    double dN18r =  1.0/4.0*(-1.0 + si)*(-1.0 + ti*ti);
    double dN18s =  1.0/4.0*(1.0 + ri)*(-1.0 + ti*ti);
    double dN18t =  1.0/2.0*(1.0 + ri)*(-1.0 + si)*ti;
    double dN19r = -1.0/4.0*(1.0 + si)*(-1.0 + ti*ti);
    double dN19s = -1.0/4.0*(1.0 + ri)*(-1.0 + ti*ti);
    double dN19t = -1.0/2.0*(1.0 + ri)*(1.0 + si)*ti;
    double dN20r =  1.0/4.0*(1.0 + si)*(-1.0 + ti*ti);
    double dN20s =  1.0/4.0*(-1.0 + ri)*(-1.0 + ti*ti);
    double dN20t =  1.0/2.0*(-1.0 + ri)*(1.0 + si)*ti;

    //Jacobian coefficients:
    double J11 = dN1r*X1(0) + dN2r*X2(0) + dN3r*X3(0) + dN4r*X4(0) + dN5r*X5(0) + dN6r*X6(0) + dN7r*X7(0) + dN8r*X8(0) + dN9r*X9(0) + dN10r*X10(0) + dN11r*X11(0) + dN12r*X12(0) + dN13r*X13(0) + dN14r*X14(0) + dN15r*X15(0) + dN16r*X16(0) + dN17r*X17(0) + dN18r*X18(0) + dN19r*X19(0) + dN20r*X20(0); 
    double J12 = dN1r*X1(1) + dN2r*X2(1) + dN3r*X3(1) + dN4r*X4(1) + dN5r*X5(1) + dN6r*X6(1) + dN7r*X7(1) + dN8r*X8(1) + dN9r*X9(1) + dN10r*X10(1) + dN11r*X11(1) + dN12r*X12(1) + dN13r*X13(1) + dN14r*X14(1) + dN15r*X15(1) + dN16r*X16(1) + dN17r*X17(1) + dN18r*X18(1) + dN19r*X19(1) + dN20r*X20(1); 
    double J13 = dN1r*X1(2) + dN2r*X2(2) + dN3r*X3(2) + dN4r*X4(2) + dN5r*X5(2) + dN6r*X6(2) + dN7r*X7(2) + dN8r*X8(2) + dN9r*X9(2) + dN10r*X10(2) + dN11r*X11(2) + dN12r*X12(2) + dN13r*X13(2) + dN14r*X14(2) + dN15r*X15(2) + dN16r*X16(2) + dN17r*X17(2) + dN18r*X18(2) + dN19r*X19(2) + dN20r*X20(2); 
    double J21 = dN1s*X1(0) + dN2s*X2(0) + dN3s*X3(0) + dN4s*X4(0) + dN5s*X5(0) + dN6s*X6(0) + dN7s*X7(0) + dN8s*X8(0) + dN9s*X9(0) + dN10s*X10(0) + dN11s*X11(0) + dN12s*X12(0) + dN13s*X13(0) + dN14s*X14(0) + dN15s*X15(0) + dN16s*X16(0) + dN17s*X17(0) + dN18s*X18(0) + dN19s*X19(0) + dN20s*X20(0); 
    double J22 = dN1s*X1(1) + dN2s*X2(1) + dN3s*X3(1) + dN4s*X4(1) + dN5s*X5(1) + dN6s*X6(1) + dN7s*X7(1) + dN8s*X8(1) + dN9s*X9(1) + dN10s*X10(1) + dN11s*X11(1) + dN12s*X12(1) + dN13s*X13(1) + dN14s*X14(1) + dN15s*X15(1) + dN16s*X16(1) + dN17s*X17(1) + dN18s*X18(1) + dN19s*X19(1) + dN20s*X20(1);
    double J23 = dN1s*X1(2) + dN2s*X2(2) + dN3s*X3(2) + dN4s*X4(2) + dN5s*X5(2) + dN6s*X6(2) + dN7s*X7(2) + dN8s*X8(2) + dN9s*X9(2) + dN10s*X10(2) + dN11s*X11(2) + dN12s*X12(2) + dN13s*X13(2) + dN14s*X14(2) + dN15s*X15(2) + dN16s*X16(2) + dN17s*X17(2) + dN18s*X18(2) + dN19s*X19(2) + dN20s*X20(2);
    double J31 = dN1t*X1(0) + dN2t*X2(0) + dN3t*X3(0) + dN4t*X4(0) + dN5t*X5(0) + dN6t*X6(0) + dN7t*X7(0) + dN8t*X8(0) + dN9t*X9(0) + dN10t*X10(0) + dN11t*X11(0) + dN12t*X12(0) + dN13t*X13(0) + dN14t*X14(0) + dN15t*X15(0) + dN16t*X16(0) + dN17t*X17(0) + dN18t*X18(0) + dN19t*X19(0) + dN20t*X20(0); 
    double J32 = dN1t*X1(1) + dN2t*X2(1) + dN3t*X3(1) + dN4t*X4(1) + dN5t*X5(1) + dN6t*X6(1) + dN7t*X7(1) + dN8t*X8(1) + dN9t*X9(1) + dN10t*X10(1) + dN11t*X11(1) + dN12t*X12(1) + dN13t*X13(1) + dN14t*X14(1) + dN15t*X15(1) + dN16t*X16(1) + dN17t*X17(1) + dN18t*X18(1) + dN19t*X19(1) + dN20t*X20(1); 
    double J33 = dN1t*X1(2) + dN2t*X2(2) + dN3t*X3(2) + dN4t*X4(2) + dN5t*X5(2) + dN6t*X6(2) + dN7t*X7(2) + dN8t*X8(2) + dN9t*X9(2) + dN10t*X10(2) + dN11t*X11(2) + dN12t*X12(2) + dN13t*X13(2) + dN14t*X14(2) + dN15t*X15(2) + dN16t*X16(2) + dN17t*X17(2) + dN18t*X18(2) + dN19t*X19(2) + dN20t*X20(2); 

    //Jacobian matrix:
    Eigen::MatrixXd Jij(3,3);
    Jij << J11, J12, J13,
           J21, J22, J23,
           J31, J32, J33;

    return Jij;
}

//Compute Shape Function at Gauss Point:
Eigen::MatrixXd 
lin3DHexa20::ComputeShapeFunctionMatrix(const double ri, const double si, const double ti) const{
    //Shape function coefficients:
    double H01 = -1.0/8.0*(1.0 - ri)*(1.0 - si)*(1.0 - ti)*(2 + ri + si + ti);
    double H02 = -1.0/8.0*(1.0 + ri)*(1.0 - si)*(1.0 - ti)*(2 - ri + si + ti);
    double H03 = -1.0/8.0*(1.0 + ri)*(1.0 + si)*(1.0 - ti)*(2 - ri - si + ti);
    double H04 = -1.0/8.0*(1.0 - ri)*(1.0 + si)*(1.0 - ti)*(2 + ri - si + ti);
    double H05 = -1.0/8.0*(1.0 - ri)*(1.0 - si)*(2 + ri + si - ti)*(1.0 + ti);
    double H06 = -1.0/8.0*(1.0 + ri)*(1.0 - si)*(2 - ri + si - ti)*(1.0 + ti);
    double H07 = -1.0/8.0*(1.0 + ri)*(1.0 + si)*(2 - ri - si - ti)*(1.0 + ti);
    double H08 = -1.0/8.0*(1.0 - ri)*(1.0 + si)*(2 + ri - si - ti)*(1.0 + ti);
    double H09 = 1.0/4.0*(1.0 - ri*ri)*(1.0 - si)*(1.0 - ti);
    double H10 = 1.0/4.0*(1.0 + ri)*(1.0 - si*si)*(1.0 - ti);
    double H11 = 1.0/4.0*(1.0 - ri*ri)*(1.0 + si)*(1.0 - ti);
    double H12 = 1.0/4.0*(1.0 - ri)*(1.0 - si*si)*(1.0 - ti);
    double H13 = 1.0/4.0*(1.0 - ri*ri)*(1.0 - si)*(1.0 + ti);
    double H14 = 1.0/4.0*(1.0 + ri)*(1.0 - si*si)*(1.0 + ti);
    double H15 = 1.0/4.0*(1.0 - ri*ri)*(1.0 + si)*(1.0 + ti);
    double H16 = 1.0/4.0*(1.0 - ri)*(1.0 - si*si)*(1.0 + ti);
    double H17 = 1.0/4.0*(1.0 - ri)*(1.0 - si)*(1.0 - ti*ti);
    double H18 = 1.0/4.0*(1.0 + ri)*(1.0 - si)*(1.0 - ti*ti);
    double H19 = 1.0/4.0*(1.0 + ri)*(1.0 + si)*(1.0 - ti*ti);
    double H20 = 1.0/4.0*(1.0 - ri)*(1.0 + si)*(1.0 - ti*ti);

    //Shape function matrix:
    Eigen::MatrixXd Hij(3,60);
    Hij << H01, 0.0, 0.0, H02, 0.0, 0.0, H03, 0.0, 0.0, H04, 0.0, 0.0, H05, 0.0, 0.0, H06, 0.0, 0.0, H07, 0.0, 0.0, H08, 0.0, 0.0, H09, 0.0, 0.0, H10, 0.0, 0.0, H11, 0.0, 0.0, H12, 0.0, 0.0, H13, 0.0, 0.0, H14, 0.0, 0.0, H15, 0.0, 0.0, H16, 0.0, 0.0, H17, 0.0, 0.0, H18, 0.0, 0.0, H19, 0.0, 0.0, H20, 0.0, 0.0,
           0.0, H01, 0.0, 0.0, H02, 0.0, 0.0, H03, 0.0, 0.0, H04, 0.0, 0.0, H05, 0.0, 0.0, H06, 0.0, 0.0, H07, 0.0, 0.0, H08, 0.0, 0.0, H09, 0.0, 0.0, H10, 0.0, 0.0, H11, 0.0, 0.0, H12, 0.0, 0.0, H13, 0.0, 0.0, H14, 0.0, 0.0, H15, 0.0, 0.0, H16, 0.0, 0.0, H17, 0.0, 0.0, H18, 0.0, 0.0, H19, 0.0, 0.0, H20, 0.0, 
           0.0, 0.0, H01, 0.0, 0.0, H02, 0.0, 0.0, H03, 0.0, 0.0, H04, 0.0, 0.0, H05, 0.0, 0.0, H06, 0.0, 0.0, H07, 0.0, 0.0, H08, 0.0, 0.0, H09, 0.0, 0.0, H10, 0.0, 0.0, H11, 0.0, 0.0, H12, 0.0, 0.0, H13, 0.0, 0.0, H14, 0.0, 0.0, H15, 0.0, 0.0, H16, 0.0, 0.0, H17, 0.0, 0.0, H18, 0.0, 0.0, H19, 0.0, 0.0, H20;

    return Hij;
}

//Evaluates the deformation matrix at a given Gauss point.
Eigen::MatrixXd 
lin3DHexa20::ComputeStrainDisplacementMatrix(const double ri, const double si, const double ti, const Eigen::MatrixXd &Jij) const{
    //Inverse jacobian matrix:
    Eigen::MatrixXd J = Jij.inverse();

    //Strain-displacement matrix coefficients:
    double B11  =  1.0/8.0*(-1.0 + si)*(-1.0 + ti)*( 1.0 + 2.0*ri + si + ti)*J(0,0) + 1.0/8.0*(-1.0 + ri)*(-1.0 + ti)*( 1.0 + ri + 2.0*si + ti)*J(0,1) + 1.0/8.0*(-1.0 + ri)*(-1.0 + si)*( 1.0 + ri + si + 2.0*ti)*J(0,2);
    double B21  = -1.0/8.0*(-1.0 + si)*(-1.0 + ti)*( 1.0 - 2.0*ri + si + ti)*J(0,0) + 1.0/8.0*( 1.0 + ri)*(-1.0 + ti)*(-1.0 + ri - 2.0*si - ti)*J(0,1) + 1.0/8.0*( 1.0 + ri)*(-1.0 + si)*(-1.0 + ri - si - 2.0*ti)*J(0,2);
    double B31  = -1.0/8.0*(1.0 + si)*(-1.0 + ti)*(-1.0 + 2.0*ri + si - ti)*J(0,0) - 1.0/8.0*(1.0 + ri)*(-1.0 + ti)*(-1.0 + ri + 2.0*si - ti)*J(0,1) - 1.0/8.0*(1.0 + ri)*( 1.0 + si)*(-1.0 + ri + si - 2.0*ti)*J(0,2);
    double B41  =  1.0/8.0*( 1.0 + si)*(-1.0 + ti)*(-1.0 - 2.0*ri + si - ti)*J(0,0) - 1.0/8.0*(-1.0 + ri)*(-1.0 + ti)*( 1.0 + ri - 2.0*si + ti)*J(0,1) - 1.0/8.0*(-1.0 + ri)*( 1.0 + si)*( 1.0 + ri - si + 2.0*ti)*J(0,2);
    double B51  = -1.0/8.0*(-1.0 + si)*( 1.0 + ti)*(1.0 + 2.0*ri + si - ti)*J(0,0) - 1.0/8.0*(-1.0 + ri)*( 1.0 + ti)*(1.0 + ri + 2.0*si - ti)*J(0,1) - 1.0/8.0*(-1.0 + ri)*(-1.0 + si)*(1.0 + ri + si - 2.0*ti)*J(0,2);
    double B61  =  1.0/8.0*(-1.0 + si)*( 1.0 + ti)*( 1.0 - 2.0*ri + si - ti)*J(0,0) - 1.0/8.0*( 1.0 + ri)*( 1.0 + ti)*(-1.0 + ri - 2.0*si + ti)*J(0,1) - 1.0/8.0*( 1.0 + ri)*(-1.0 + si)*(-1.0 + ri - si + 2.0*ti)*J(0,2);
    double B71  = 1.0/8.0*(1.0 + si)*(1.0 + ti)*(-1.0 + 2.0*ri + si + ti)*J(0,0) + 1.0/8.0*(1.0 + ri)*(1.0 + ti)*(-1.0 + ri + 2.0*si + ti)*J(0,1) + 1.0/8.0*(1.0 + ri)*(1.0 + si)*(-1.0 + ri + si + 2.0*ti)*J(0,2);
    double B81  = -1.0/8.0*(1.0 + si)*(1.0 + ti)*(-1.0 - 2.0*ri + si + ti)*J(0,0) + 1.0/8.0*(-1.0 + ri)*(1.0 + ti)*( 1.0 + ri - 2.0*si - ti)*J(0,1) + 1.0/8.0*(-1.0 + ri)*(1.0 + si)*( 1.0 + ri - si - 2.0*ti)*J(0,2);
    double B91  = -1.0/2.0*ri*(-1.0 + si)*(-1.0 + ti)*J(0,0) - 1.0/4.0*(-1.0 + ri*ri)*(-1.0 + ti)*J(0,1) - 1.0/4.0*(-1.0 + ri*ri)*(-1.0 + si)*J(0,2);
    double B101  = 1.0/4.0*(-1.0 + si*si)*(-1.0 + ti)*J(0,0) + 1.0/2.0*( 1.0 + ri)*si*(-1.0 + ti)*J(0,1) + 1.0/4.0*( 1.0 + ri)*(-1.0 + si*si)*J(0,2);
    double B111  = 1.0/2.0*ri*(1.0 + si)*(-1.0 + ti)*J(0,0) + 1.0/4.0*(-1.0 + ri*ri)*(-1.0 + ti)*J(0,1) + 1.0/4.0*(-1.0 + ri*ri)*(1.0 + si)*J(0,2);
    double B121  = -1.0/4.0*(-1.0 + si*si)*(-1.0 + ti)*J(0,0) - 1.0/2.0*(-1.0 + ri)*si*(-1.0 + ti)*J(0,1) - 1.0/4.0*(-1.0 + ri)*(-1.0 + si*si)*J(0,2);
    double B131  = 1.0/2.0*ri*(-1.0 + si)*(1.0 + ti)*J(0,0) + 1.0/4.0*(-1.0 + ri*ri)*(1.0 + ti)*J(0,1) + 1.0/4.0*(-1.0 + ri*ri)*(-1.0 + si)*J(0,2);
    double B141  = -1.0/4.0*(-1.0 + si*si)*(1.0 + ti)*J(0,0) - 1.0/2.0*(1.0 + ri)*si*(1.0 + ti)*J(0,1) - 1.0/4.0*(1.0 + ri)*(-1.0 + si*si)*J(0,2);
    double B151  = -1.0/2.0*ri*(1.0 + si)*(1.0 + ti)*J(0,0) - 1.0/4.0*(-1.0 + ri*ri)*(1.0 + ti)*J(0,1) - 1.0/4.0*(-1.0 + ri*ri)*(1.0 + si)*J(0,2);
    double B161  = 1.0/4.0*(-1.0 + si*si)*(1.0 + ti)*J(0,0) + 1.0/2.0*(-1.0 + ri)*si*(1.0 + ti)*J(0,1) + 1.0/4.0*(-1.0 + ri)*(-1.0 + si*si)*J(0,2);
    double B171  = -1.0/4.0*(-1.0 + si)*(-1.0 + ti*ti)*J(0,0) - 1.0/4.0*(-1.0 + ri)*(-1.0 + ti*ti)*J(0,1) - 1.0/2.0*(-1.0 + ri)*(-1.0 + si)*ti*J(0,2);
    double B181  = 1.0/4.0*(-1.0 + si)*(-1.0 + ti*ti)*J(0,0) + 1.0/4.0*(1.0 + ri)*(-1.0 + ti*ti)*J(0,1) + 1.0/2.0*(1.0 + ri)*(-1.0 + si)*ti*J(0,2);
    double B191  = -1.0/4.0*(1.0 + si)*(-1.0 + ti*ti)*J(0,0) - 1.0/4.0*(1.0 + ri)*(-1.0 + ti*ti)*J(0,1) - 1.0/2.0*(1.0 + ri)*(1.0 + si)*ti*J(0,2);
    double B201  = 1.0/4.0*(1.0 + si)*(-1.0 + ti*ti)*J(0,0) + 1.0/4.0*(-1.0 + ri)*(-1.0 + ti*ti)*J(0,1) + 1.0/2.0*(-1.0 + ri)*(1.0 + si)*ti*J(0,2);

    double B12  =  1.0/8.0*(-1.0 + si)*(-1.0 + ti)*( 1.0 + 2.0*ri + si + ti)*J(1,0) + 1.0/8.0*(-1.0 + ri)*(-1.0 + ti)*( 1.0 + ri + 2.0*si + ti)*J(1,1) + 1.0/8.0*(-1.0 + ri)*(-1.0 + si)*( 1.0 + ri + si + 2.0*ti)*J(1,2);
    double B22  = -1.0/8.0*(-1.0 + si)*(-1.0 + ti)*( 1.0 - 2.0*ri + si + ti)*J(1,0) + 1.0/8.0*( 1.0 + ri)*(-1.0 + ti)*(-1.0 + ri - 2.0*si - ti)*J(1,1) + 1.0/8.0*( 1.0 + ri)*(-1.0 + si)*(-1.0 + ri - si - 2.0*ti)*J(1,2);
    double B32  = -1.0/8.0*(1.0 + si)*(-1.0 + ti)*(-1.0 + 2.0*ri + si - ti)*J(1,0) - 1.0/8.0*(1.0 + ri)*(-1.0 + ti)*(-1.0 + ri + 2.0*si - ti)*J(1,1) - 1.0/8.0*(1.0 + ri)*( 1.0 + si)*(-1.0 + ri + si - 2.0*ti)*J(1,2);
    double B42  =  1.0/8.0*( 1.0 + si)*(-1.0 + ti)*(-1.0 - 2.0*ri + si - ti)*J(1,0) - 1.0/8.0*(-1.0 + ri)*(-1.0 + ti)*( 1.0 + ri - 2.0*si + ti)*J(1,1) - 1.0/8.0*(-1.0 + ri)*( 1.0 + si)*( 1.0 + ri - si + 2.0*ti)*J(1,2);
    double B52  = -1.0/8.0*(-1.0 + si)*( 1.0 + ti)*(1.0 + 2.0*ri + si - ti)*J(1,0) - 1.0/8.0*(-1.0 + ri)*( 1.0 + ti)*(1.0 + ri + 2.0*si - ti)*J(1,1) - 1.0/8.0*(-1.0 + ri)*(-1.0 + si)*(1.0 + ri + si - 2.0*ti)*J(1,2);
    double B62  =  1.0/8.0*(-1.0 + si)*( 1.0 + ti)*( 1.0 - 2.0*ri + si - ti)*J(1,0) - 1.0/8.0*( 1.0 + ri)*( 1.0 + ti)*(-1.0 + ri - 2.0*si + ti)*J(1,1) - 1.0/8.0*( 1.0 + ri)*(-1.0 + si)*(-1.0 + ri - si + 2.0*ti)*J(1,2);
    double B72  = 1.0/8.0*(1.0 + si)*(1.0 + ti)*(-1.0 + 2.0*ri + si + ti)*J(1,0) + 1.0/8.0*(1.0 + ri)*(1.0 + ti)*(-1.0 + ri + 2.0*si + ti)*J(1,1) + 1.0/8.0*(1.0 + ri)*(1.0 + si)*(-1.0 + ri + si + 2.0*ti)*J(1,2);
    double B82  = -1.0/8.0*(1.0 + si)*(1.0 + ti)*(-1.0 - 2.0*ri + si + ti)*J(1,0) + 1.0/8.0*(-1.0 + ri)*(1.0 + ti)*( 1.0 + ri - 2.0*si - ti)*J(1,1) + 1.0/8.0*(-1.0 + ri)*(1.0 + si)*( 1.0 + ri - si - 2.0*ti)*J(1,2);
    double B92  = -1.0/2.0*ri*(-1.0 + si)*(-1.0 + ti)*J(1,0) - 1.0/4.0*(-1.0 + ri*ri)*(-1.0 + ti)*J(1,1) - 1.0/4.0*(-1.0 + ri*ri)*(-1.0 + si)*J(1,2);
    double B102  = 1.0/4.0*(-1.0 + si*si)*(-1.0 + ti)*J(1,0) + 1.0/2.0*( 1.0 + ri)*si*(-1.0 + ti)*J(1,1) + 1.0/4.0*( 1.0 + ri)*(-1.0 + si*si)*J(1,2);
    double B112  = 1.0/2.0*ri*(1.0 + si)*(-1.0 + ti)*J(1,0) + 1.0/4.0*(-1.0 + ri*ri)*(-1.0 + ti)*J(1,1) + 1.0/4.0*(-1.0 + ri*ri)*(1.0 + si)*J(1,2);
    double B122  = -1.0/4.0*(-1.0 + si*si)*(-1.0 + ti)*J(1,0) - 1.0/2.0*(-1.0 + ri)*si*(-1.0 + ti)*J(1,1) - 1.0/4.0*(-1.0 + ri)*(-1.0 + si*si)*J(1,2);
    double B132  = 1.0/2.0*ri*(-1.0 + si)*(1.0 + ti)*J(1,0) + 1.0/4.0*(-1.0 + ri*ri)*(1.0 + ti)*J(1,1) + 1.0/4.0*(-1.0 + ri*ri)*(-1.0 + si)*J(1,2);
    double B142  = -1.0/4.0*(-1.0 + si*si)*(1.0 + ti)*J(1,0) - 1.0/2.0*(1.0 + ri)*si*(1.0 + ti)*J(1,1) - 1.0/4.0*(1.0 + ri)*(-1.0 + si*si)*J(1,2);
    double B152  = -1.0/2.0*ri*(1.0 + si)*(1.0 + ti)*J(1,0) - 1.0/4.0*(-1.0 + ri*ri)*(1.0 + ti)*J(1,1) - 1.0/4.0*(-1.0 + ri*ri)*(1.0 + si)*J(1,2);
    double B162  = 1.0/4.0*(-1.0 + si*si)*(1.0 + ti)*J(1,0) + 1.0/2.0*(-1.0 + ri)*si*(1.0 + ti)*J(1,1) + 1.0/4.0*(-1.0 + ri)*(-1.0 + si*si)*J(1,2);
    double B172  = -1.0/4.0*(-1.0 + si)*(-1.0 + ti*ti)*J(1,0) - 1.0/4.0*(-1.0 + ri)*(-1.0 + ti*ti)*J(1,1) - 1.0/2.0*(-1.0 + ri)*(-1.0 + si)*ti*J(1,2);
    double B182  = 1.0/4.0*(-1.0 + si)*(-1.0 + ti*ti)*J(1,0) + 1.0/4.0*(1.0 + ri)*(-1.0 + ti*ti)*J(1,1) + 1.0/2.0*(1.0 + ri)*(-1.0 + si)*ti*J(1,2);
    double B192  = -1.0/4.0*(1.0 + si)*(-1.0 + ti*ti)*J(1,0) - 1.0/4.0*(1.0 + ri)*(-1.0 + ti*ti)*J(1,1) - 1.0/2.0*(1.0 + ri)*(1.0 + si)*ti*J(1,2);
    double B202  = 1.0/4.0*(1.0 + si)*(-1.0 + ti*ti)*J(1,0) + 1.0/4.0*(-1.0 + ri)*(-1.0 + ti*ti)*J(1,1) + 1.0/2.0*(-1.0 + ri)*(1.0 + si)*ti*J(1,2);

    double B13  =  1.0/8.0*(-1.0 + si)*(-1.0 + ti)*( 1.0 + 2.0*ri + si + ti)*J(2,0) + 1.0/8.0*(-1.0 + ri)*(-1.0 + ti)*( 1.0 + ri + 2.0*si + ti)*J(2,1) + 1.0/8.0*(-1.0 + ri)*(-1.0 + si)*( 1.0 + ri + si + 2.0*ti)*J(2,2);
    double B23  = -1.0/8.0*(-1.0 + si)*(-1.0 + ti)*( 1.0 - 2.0*ri + si + ti)*J(2,0) + 1.0/8.0*( 1.0 + ri)*(-1.0 + ti)*(-1.0 + ri - 2.0*si - ti)*J(2,1) + 1.0/8.0*( 1.0 + ri)*(-1.0 + si)*(-1.0 + ri - si - 2.0*ti)*J(2,2);
    double B33  = -1.0/8.0*(1.0 + si)*(-1.0 + ti)*(-1.0 + 2.0*ri + si - ti)*J(2,0) - 1.0/8.0*(1.0 + ri)*(-1.0 + ti)*(-1.0 + ri + 2.0*si - ti)*J(2,1) - 1.0/8.0*(1.0 + ri)*( 1.0 + si)*(-1.0 + ri + si - 2.0*ti)*J(2,2);
    double B43  =  1.0/8.0*( 1.0 + si)*(-1.0 + ti)*(-1.0 - 2.0*ri + si - ti)*J(2,0) - 1.0/8.0*(-1.0 + ri)*(-1.0 + ti)*( 1.0 + ri - 2.0*si + ti)*J(2,1) - 1.0/8.0*(-1.0 + ri)*( 1.0 + si)*( 1.0 + ri - si + 2.0*ti)*J(2,2);
    double B53  = -1.0/8.0*(-1.0 + si)*( 1.0 + ti)*(1.0 + 2.0*ri + si - ti)*J(2,0) - 1.0/8.0*(-1.0 + ri)*( 1.0 + ti)*(1.0 + ri + 2.0*si - ti)*J(2,1) - 1.0/8.0*(-1.0 + ri)*(-1.0 + si)*(1.0 + ri + si - 2.0*ti)*J(2,2);
    double B63  =  1.0/8.0*(-1.0 + si)*( 1.0 + ti)*( 1.0 - 2.0*ri + si - ti)*J(2,0) - 1.0/8.0*( 1.0 + ri)*( 1.0 + ti)*(-1.0 + ri - 2.0*si + ti)*J(2,1) - 1.0/8.0*( 1.0 + ri)*(-1.0 + si)*(-1.0 + ri - si + 2.0*ti)*J(2,2);
    double B73  = 1.0/8.0*(1.0 + si)*(1.0 + ti)*(-1.0 + 2.0*ri + si + ti)*J(2,0) + 1.0/8.0*(1.0 + ri)*(1.0 + ti)*(-1.0 + ri + 2.0*si + ti)*J(2,1) + 1.0/8.0*(1.0 + ri)*(1.0 + si)*(-1.0 + ri + si + 2.0*ti)*J(2,2);
    double B83  = -1.0/8.0*(1.0 + si)*(1.0 + ti)*(-1.0 - 2.0*ri + si + ti)*J(2,0) + 1.0/8.0*(-1.0 + ri)*(1.0 + ti)*( 1.0 + ri - 2.0*si - ti)*J(2,1) + 1.0/8.0*(-1.0 + ri)*(1.0 + si)*( 1.0 + ri - si - 2.0*ti)*J(2,2);
    double B93  = -1.0/2.0*ri*(-1.0 + si)*(-1.0 + ti)*J(2,0) - 1.0/4.0*(-1.0 + ri*ri)*(-1.0 + ti)*J(2,1) - 1.0/4.0*(-1.0 + ri*ri)*(-1.0 + si)*J(2,2);
    double B103  = 1.0/4.0*(-1.0 + si*si)*(-1.0 + ti)*J(2,0) + 1.0/2.0*( 1.0 + ri)*si*(-1.0 + ti)*J(2,1) + 1.0/4.0*( 1.0 + ri)*(-1.0 + si*si)*J(2,2);
    double B113  = 1.0/2.0*ri*(1.0 + si)*(-1.0 + ti)*J(2,0) + 1.0/4.0*(-1.0 + ri*ri)*(-1.0 + ti)*J(2,1) + 1.0/4.0*(-1.0 + ri*ri)*(1.0 + si)*J(2,2);
    double B123  = -1.0/4.0*(-1.0 + si*si)*(-1.0 + ti)*J(2,0) - 1.0/2.0*(-1.0 + ri)*si*(-1.0 + ti)*J(2,1) - 1.0/4.0*(-1.0 + ri)*(-1.0 + si*si)*J(2,2);
    double B133  = 1.0/2.0*ri*(-1.0 + si)*(1.0 + ti)*J(2,0) + 1.0/4.0*(-1.0 + ri*ri)*(1.0 + ti)*J(2,1) + 1.0/4.0*(-1.0 + ri*ri)*(-1.0 + si)*J(2,2);
    double B143  = -1.0/4.0*(-1.0 + si*si)*(1.0 + ti)*J(2,0) - 1.0/2.0*(1.0 + ri)*si*(1.0 + ti)*J(2,1) - 1.0/4.0*(1.0 + ri)*(-1.0 + si*si)*J(2,2);
    double B153  = -1.0/2.0*ri*(1.0 + si)*(1.0 + ti)*J(2,0) - 1.0/4.0*(-1.0 + ri*ri)*(1.0 + ti)*J(2,1) - 1.0/4.0*(-1.0 + ri*ri)*(1.0 + si)*J(2,2);
    double B163  = 1.0/4.0*(-1.0 + si*si)*(1.0 + ti)*J(2,0) + 1.0/2.0*(-1.0 + ri)*si*(1.0 + ti)*J(2,1) + 1.0/4.0*(-1.0 + ri)*(-1.0 + si*si)*J(2,2);
    double B173  = -1.0/4.0*(-1.0 + si)*(-1.0 + ti*ti)*J(2,0) - 1.0/4.0*(-1.0 + ri)*(-1.0 + ti*ti)*J(2,1) - 1.0/2.0*(-1.0 + ri)*(-1.0 + si)*ti*J(2,2);
    double B183  = 1.0/4.0*(-1.0 + si)*(-1.0 + ti*ti)*J(2,0) + 1.0/4.0*(1.0 + ri)*(-1.0 + ti*ti)*J(2,1) + 1.0/2.0*(1.0 + ri)*(-1.0 + si)*ti*J(2,2);
    double B193  = -1.0/4.0*(1.0 + si)*(-1.0 + ti*ti)*J(2,0) - 1.0/4.0*(1.0 + ri)*(-1.0 + ti*ti)*J(2,1) - 1.0/2.0*(1.0 + ri)*(1.0 + si)*ti*J(2,2);
    double B203  = 1.0/4.0*(1.0 + si)*(-1.0 + ti*ti)*J(2,0) + 1.0/4.0*(-1.0 + ri)*(-1.0 + ti*ti)*J(2,1) + 1.0/2.0*(-1.0 + ri)*(1.0 + si)*ti*J(2,2);

    //Deformation matrix definition:
    Eigen::MatrixXd Bij(6,60);
    Bij << B11, 0.0, 0.0, B21, 0.0, 0.0, B31, 0.0, 0.0, B41, 0.0, 0.0, B51, 0.0, 0.0, B61, 0.0, 0.0, B71, 0.0, 0.0, B81, 0.0, 0.0, B91, 0.0, 0.0, B101,  0.0,  0.0, B111,  0.0,  0.0, B121,  0.0,  0.0, B131,  0.0,  0.0, B141,  0.0,  0.0, B151,  0.0,  0.0, B161,  0.0,  0.0, B171,  0.0,  0.0, B181,  0.0,  0.0, B191,  0.0,  0.0, B201,  0.0,  0.0,
           0.0, B12, 0.0, 0.0, B22, 0.0, 0.0, B32, 0.0, 0.0, B42, 0.0, 0.0, B52, 0.0, 0.0, B62, 0.0, 0.0, B72, 0.0, 0.0, B82, 0.0, 0.0, B92, 0.0,  0.0, B102,  0.0,  0.0, B112,  0.0,  0.0, B122,  0.0,  0.0, B132,  0.0,  0.0, B142,  0.0,  0.0, B152,  0.0,  0.0, B162,  0.0,  0.0, B172,  0.0,  0.0, B182,  0.0,  0.0, B192,  0.0,  0.0, B202,  0.0,
           0.0, 0.0, B13, 0.0, 0.0, B23, 0.0, 0.0, B33, 0.0, 0.0, B43, 0.0, 0.0, B53, 0.0, 0.0, B63, 0.0, 0.0, B73, 0.0, 0.0, B83, 0.0, 0.0, B93,  0.0,  0.0, B103,  0.0,  0.0, B113,  0.0,  0.0, B123,  0.0,  0.0, B133,  0.0,  0.0, B143,  0.0,  0.0, B153,  0.0,  0.0, B163,  0.0,  0.0, B173,  0.0,  0.0, B183,  0.0,  0.0, B193,  0.0,  0.0, B203,
           B12, B11, 0.0, B22, B21, 0.0, B32, B31, 0.0, B42, B41, 0.0, B52, B51, 0.0, B62, B61, 0.0, B72, B71, 0.0, B82, B81, 0.0, B92, B91, 0.0, B102, B101,  0.0, B112, B111,  0.0, B122, B121,  0.0, B132, B131,  0.0, B142, B141,  0.0, B152, B151,  0.0, B162, B161,  0.0, B172, B171,  0.0, B182, B181,  0.0, B192, B191,  0.0, B202, B201,  0.0,
           0.0, B13, B12, 0.0, B23, B22, 0.0, B33, B32, 0.0, B43, B42, 0.0, B53, B52, 0.0, B63, B62, 0.0, B73, B72, 0.0, B83, B82, 0.0, B93, B92,  0.0, B103, B102,  0.0, B113, B112,  0.0, B123, B122,  0.0, B133, B132,  0.0, B143, B142,  0.0, B153, B152,  0.0, B163, B162,  0.0, B173, B172,  0.0, B183, B182,  0.0, B193, B192,  0.0, B203, B202,
           B13, 0.0, B11, B23, 0.0, B21, B33, 0.0, B31, B43, 0.0, B41, B53, 0.0, B51, B63, 0.0, B61, B73, 0.0, B71, B83, 0.0, B81, B93, 0.0, B91, B103,  0.0, B101, B113,  0.0, B111, B123,  0.0, B121, B133,  0.0, B131, B143,  0.0, B141, B153,  0.0, B151, B163,  0.0, B161, B173,  0.0, B171, B183,  0.0, B181, B193,  0.0, B191, B203,  0.0, B201;

    return Bij;
}

//Compute the initial stiffness matrix of the element using gauss-integration.
Eigen::MatrixXd 
lin3DHexa20::ComputeInitialStiffnessMatrix() const{
    //Stiffness matrix definition:
    Eigen::MatrixXd StiffnessMatrix(60, 60);
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
        Eigen::MatrixXd Bij = ComputeStrainDisplacementMatrix(xi(i,0), xi(i,1), xi(i,2), Jij);

        //Gets material tangent matrix at Gauss point.
        Eigen::MatrixXd Cij = theMaterial[i]->GetInitialTangentStiffness();

        //Numerical integration.
        StiffnessMatrix += wi(i)*fabs(Jij.determinant())*Bij.transpose()*Cij*Bij;
    }

    return StiffnessMatrix;
}
