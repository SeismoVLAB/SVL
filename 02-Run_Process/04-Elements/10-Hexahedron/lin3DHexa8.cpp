#include <cmath>
#include <Eigen/LU> 
#include "Material.hpp"
#include "lin3DHexa8.hpp"
#include "GaussQuadrature.hpp"
#include "LobattoQuadrature.hpp"
#include "Definitions.hpp"
#include "Profiler.hpp"

//Overload constructor.
lin3DHexa8::lin3DHexa8(const std::vector<unsigned int> nodes, std::unique_ptr<Material> &material, const std::string quadrature, const unsigned int nGauss) :
Element("lin3DHexa8", nodes, 24, VTK_LINEAR_HEXA, GROUP_ELEMENT_HEXA){
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
lin3DHexa8::~lin3DHexa8(){
    //Does nothing.
}

//Save the material states in the element.
void 
lin3DHexa8::CommitState(){
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

//Reverse the material states to previous converged state in this element.
void 
lin3DHexa8::ReverseState(){
    //Reverse the material components.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    for(unsigned int k = 0; k < nPoints; k++)
        theMaterial[k]->ReverseState();
}

//Brings the material state to its initial state in this element.
void 
lin3DHexa8::InitialState(){
    //Brings the material components to initial state.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    for(unsigned int k = 0; k < nPoints; k++)
        theMaterial[k]->InitialState();
}

//Update the material states in the element.
void 
lin3DHexa8::UpdateState(){
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
lin3DHexa8::SetDomain(std::map<unsigned int, std::shared_ptr<Node> > &nodes){
    //Gets the global element connectivity.
    std::vector<unsigned int> conn = GetNodes();

    //Assign the element to mesh node pointer.  
    for(unsigned int i = 0; i < GetNumberOfNodes(); i++){
        theNodes[i] = nodes[conn[i]];
    }
}

//Sets the damping model.
void 
lin3DHexa8::SetDamping(const std::shared_ptr<Damping> &damping){
    //The damping model
    theDamping = damping;
}

//Gets the list of total-degree of freedom of this element.
std::vector<unsigned int> 
lin3DHexa8::GetTotalDegreeOfFreedom() const{
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
lin3DHexa8::GetStrain() const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    Eigen::MatrixXd theStrain(nPoints,6);
    for(unsigned int k = 0; k < nPoints; k++)
        theStrain.row(k) = theMaterial[k]->GetStrain();

    return theStrain;
}

//Returns the material stress at integration points.
Eigen::MatrixXd 
lin3DHexa8::GetStress() const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    Eigen::MatrixXd theStress(nPoints,6);
    for(unsigned int k = 0; k < nPoints; k++)
        theStress.row(k) = theMaterial[k]->GetTotalStress();

    return theStress;
}

//Returns the material strain-rate at integration points.
Eigen::MatrixXd 
lin3DHexa8::GetStrainRate() const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    Eigen::MatrixXd theStrainRate(nPoints,6);
    for(unsigned int k = 0; k < nPoints; k++)
        theStrainRate.row(k) = theMaterial[k]->GetStrainRate();

    return theStrainRate;
}

//Gets the material strain in section at  coordinate (x3,x2).
Eigen::MatrixXd 
lin3DHexa8::GetStrainAt(double UNUSED(x3), double UNUSED(x2)) const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    //Stress at coordinate is define within section.
    Eigen::MatrixXd theStrain(nPoints, 6);
    theStrain.fill(0.0);

    return theStrain;
}

//Gets the material stress in section at  coordinate (x3,x2).
Eigen::MatrixXd 
lin3DHexa8::GetStressAt(double UNUSED(x3), double UNUSED(x2)) const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    //Stress at coordinate is define within section.
    Eigen::MatrixXd theStress(nPoints, 6);
    theStress.fill(0.0);

    return theStress;
}

//Gets the element internal response in VTK format.
Eigen::VectorXd 
lin3DHexa8::GetVTKResponse(std::string response) const{
    //The VTK response vector.
    Eigen::VectorXd theResponse(18);

    if (strcasecmp(response.c_str(),"Strain") == 0){
        Eigen::MatrixXd strain = GetStrain();
        Eigen::VectorXd Strain = strain.colwise().mean();
        theResponse << Strain(0), Strain(1), Strain(2), Strain(3), Strain(4), Strain(5), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    }
    else if(strcasecmp(response.c_str(),"Stress") == 0){
        Eigen::MatrixXd stress = GetStress();
        Eigen::VectorXd Stress = stress.colwise().mean();
        theResponse << Stress(0), Stress(1), Stress(2), Stress(3), Stress(4), Stress(5), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    }

    return theResponse;
}

//Computes the element energy for a given deformation.
double 
lin3DHexa8::ComputeEnergy(){
    //TODO: Integrate over element volume to compute the energy
    return 0.0;
}

//Compute the mass matrix of the element using gauss-integration.
Eigen::MatrixXd 
lin3DHexa8::ComputeMassMatrix(){
    //Starts profiling this function.
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

    //Lumped Mass Formulation
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
lin3DHexa8::ComputeStiffnessMatrix(){
    //Starts profiling this function.
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
lin3DHexa8::ComputeDampingMatrix(){
    //Starts profiling this function.
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
lin3DHexa8::ComputePMLMatrix(){
    Eigen::MatrixXd Kpml;
    return Kpml;
}

//Compute the internal forces acting on the element.
Eigen::VectorXd 
lin3DHexa8::ComputeInternalForces(){
    //Starts profiling this function.
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
        Eigen::MatrixXd Bij = ComputeStrainDisplacementMatrix(xi(i,0), xi(i,1), xi(i,2), Jij);

        //Gets material strain at Gauss point.
        Eigen::VectorXd Stress = theMaterial[i]->GetStress();

        //Numerical integration.
        InternalForces += wi(i)*fabs(Jij.determinant())*Bij.transpose()*Stress;
    }
    
    return InternalForces;
}

//Compute the elastic, inertial, and viscous forces acting on the element.
Eigen::VectorXd 
lin3DHexa8::ComputeInternalDynamicForces(){
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
lin3DHexa8::ComputeSurfaceForces(const std::shared_ptr<Load> &surfaceLoad, unsigned int face){
    //Starts profiling this function.
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
lin3DHexa8::ComputeBodyForces(const std::shared_ptr<Load> &bodyLoad, unsigned int k){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Local body load vector.
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
lin3DHexa8::ComputeDomainReductionForces(const std::shared_ptr<Load> &drm, unsigned int k){
    //Starts profiling this function.
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

//Update strain in the element.
Eigen::VectorXd 
lin3DHexa8::ComputeStrain(const Eigen::MatrixXd &Bij) const{
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
    nodalDisplacement << U1, U2, U3, U4, U5, U6, U7, U8;

    //Strain vector:
    Eigen::VectorXd Strain(6); 
    Strain = Bij*nodalDisplacement;

    return Strain;
}

//Update strain rate in the element.
Eigen::VectorXd 
lin3DHexa8::ComputeStrainRate(const Eigen::MatrixXd& UNUSED(Bij)) const{
    //TODO: Compute strain rate.
    //Strain vector definition:
    Eigen::VectorXd strainrate(6);
    strainrate.fill(0.0);

    return strainrate;
}

//Computes the jacobian of the transformation. 
Eigen::MatrixXd 
lin3DHexa8::ComputeJacobianMatrix(const double ri, const double si, const double ti) const{
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

//Compute Shape Function at Gauss Point:
Eigen::MatrixXd 
lin3DHexa8::ComputeShapeFunctionMatrix(const double ri, const double si, const double ti) const{
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
lin3DHexa8::ComputeStrainDisplacementMatrix(const double ri, const double si, const double ti, const Eigen::MatrixXd &Jij) const{
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

//Compute the initial stiffness matrix of the element using gauss-integration.
Eigen::MatrixXd 
lin3DHexa8::ComputeInitialStiffnessMatrix() const{
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
        Eigen::MatrixXd Bij = ComputeStrainDisplacementMatrix(xi(i,0), xi(i,1), xi(i,2), Jij);

        //Gets material tangent matrix at Gauss point.
        Eigen::MatrixXd Cij = theMaterial[i]->GetInitialTangentStiffness();

        //Numerical integration.
        StiffnessMatrix += wi(i)*fabs(Jij.determinant())*Bij.transpose()*Cij*Bij;
    }

    return StiffnessMatrix;
}
