#include <cmath>
#include "lin3DFrame2.hpp"
#include "GaussQuadrature.hpp"
#include "LobattoQuadrature.hpp"
#include "Definitions.hpp"
#include "Profiler.hpp"

//Define constant tolerance value:
const double TOL = 0.9999995;

//Overload constructor.
lin3DFrame2::lin3DFrame2(const std::vector<unsigned int> nodes, std::unique_ptr<Section> &section, bool formulation, const std::string quadrature, unsigned int nGauss) :
Element("lin3DFrame2", nodes, 12, VTK_LINEAR_LINE, GROUP_ELEMENT_FRAME), Formulation(formulation), Phiy(0.0), Phiz(0.0){
    //The element nodes.
    theNodes.resize(2);

    //Numerical integration rule.
    if(strcasecmp(quadrature.c_str(),"GAUSS") == 0)
        QuadraturePoints = std::make_unique<GaussQuadrature>("Line", nGauss);
    else if(strcasecmp(quadrature.c_str(),"LOBATTO") == 0)
        QuadraturePoints = std::make_unique<LobattoQuadrature>("Line", nGauss);

    //The element material. 
    theSection.resize(nGauss);
    for(unsigned int i = 0; i < nGauss; i++)
        theSection[i] = section->CopySection();

    //Section Shear Factor.
    if(Formulation){
        Eigen::MatrixXd Cs = section->GetInitialTangentStiffness();
        Phiy = Cs(3,3)/Cs(4,4);
        Phiz = Cs(2,2)/Cs(5,5);
    }
}

//Destructor.
lin3DFrame2::~lin3DFrame2(){
    //Does nothing.
}

//Save the material states in the element.
void 
lin3DFrame2::CommitState(){
    //It considers only linear elastic material, thus viscous material is not allowed.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    for(unsigned int k = 0; k < nPoints; k++)
        theSection[k]->CommitState();
}

//Reverse the section states to previous converged state in this element.
void 
lin3DFrame2::ReverseState(){
    //Reverse the section components.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    for(unsigned int k = 0; k < nPoints; k++)
        theSection[k]->ReverseState();
}

//Brings the section state to its initial state in this element.
void 
lin3DFrame2::InitialState(){
    //Brings the section components to initial state.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    for(unsigned int k = 0; k < nPoints; k++)
        theSection[k]->InitialState();
}

//Update the section states in the element.
void 
lin3DFrame2::UpdateState(){
    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Line", wi, xi);

    //Update material states.
    for(unsigned int k = 0; k < wi.size(); k++){
        //Compute Strain-Displacement Matrix at Gauss Point.
        Eigen::MatrixXd Bij = ComputeStrainDisplacementMatrix(xi(k,0));

        //Computes strain vector.
        Eigen::VectorXd strain = ComputeStrain(Bij);

        //Update the material state.
        theSection[k]->UpdateState(strain, 1);
    }
}

//Sets the finite element dependance among objects.
void 
lin3DFrame2::SetDomain(std::map<unsigned int, std::shared_ptr<Node> > &nodes){
    //Gets the global element connectivity.
    std::vector<unsigned int> conn = GetNodes();

    //Assign the element to mesh node pointer.  
    for(unsigned int i = 0; i < GetNumberOfNodes(); i++){
        theNodes[i] = nodes[conn[i]];
    }

    //Computes length of element. 
    Eigen::VectorXd Xi = theNodes[0]->GetCoordinates();
    Eigen::VectorXd Xj = theNodes[1]->GetCoordinates();

    L = (Xj - Xi).norm();

    //Section Shear Factor.
    if(Formulation){
        Phiy = 12.0*Phiy/L/L;
        Phiz = 12.0*Phiz/L/L;
    }
}

//Sets the damping model.
void 
lin3DFrame2::SetDamping(const std::shared_ptr<Damping> &damping){
    //The damping model
    theDamping = damping;
}

//Gets the list of total-degree of freedom of this element.
std::vector<unsigned int> 
lin3DFrame2::GetTotalDegreeOfFreedom() const{
    //Total number of degree-of-freedom.
    unsigned int nDofs = GetNumberOfDegreeOfFreedom();

    //Reserve memory for the element list of degree-of-freedom.
    std::vector<unsigned int> dofs(nDofs);

    //Construct the element list of degree-of-freedom for assembly.
    for(unsigned int j = 0; j < 2; j++){    
        unsigned int LengthDofs = theNodes[j]->GetNumberOfDegreeOfFreedom();
        std::vector<int> totalDofs = theNodes[j]->GetTotalDegreeOfFreedom();

        for(unsigned int i = 0; i < LengthDofs; i++)
            dofs[i + LengthDofs*j] = totalDofs[i];    
    }

    return dofs;
}

//Returns the section generalised strain at integration point.
Eigen::MatrixXd 
lin3DFrame2::GetStrain() const{
    //Number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    Eigen::MatrixXd theStrain(nPoints,6);
    for(unsigned int k = 0; k < nPoints; k++)
        theStrain.row(k) = theSection[k]->GetStrain();

    return theStrain;
}

//Returns the section generalised stress at integration point.
Eigen::MatrixXd 
lin3DFrame2::GetStress() const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    Eigen::MatrixXd theStress(nPoints,6);
    for(unsigned int k = 0; k < nPoints; k++)
        theStress.row(k) = theSection[k]->GetStress();

    return theStress;
}

//Returns the section generalised strain rate at integration point.
Eigen::MatrixXd 
lin3DFrame2::GetStrainRate() const{
    //TODO: No strain-rate is employed for frame section.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    Eigen::MatrixXd theStrainRate(nPoints,6);
    theStrainRate.fill(0.0);

    return theStrainRate;
}

//Gets the material strain in section at  coordinate (x3,x2).
Eigen::MatrixXd 
lin3DFrame2::GetStrainAt(double x3, double x2) const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    //Stress at coordinate is define within section.
    Eigen::MatrixXd theStrain(nPoints, 6);
    for(unsigned int k = 0; k < nPoints; k++)
        theStrain.row(k) = theSection[k]->GetStrainAt(x3, x2);

    return theStrain;
}

//Gets the material stress in section at  coordinate (x3,x2).
Eigen::MatrixXd 
lin3DFrame2::GetStressAt(double x3, double x2) const{
    //number of integration points.
    unsigned int nPoints = QuadraturePoints->GetNumberOfQuadraturePoints();

    //Stress at coordinate is define within section.
    Eigen::MatrixXd theStress(nPoints, 6);
    for(unsigned int k = 0; k < nPoints; k++)
        theStress.row(k) = theSection[k]->GetStressAt(x3, x2);

    return theStress;
}

//Gets the element internal response in VTK format.
Eigen::VectorXd 
lin3DFrame2::GetVTKResponse(std::string response) const{   
    //The VTK response vector.
    Eigen::VectorXd theResponse(18);

    if (strcasecmp(response.c_str(),"Strain") == 0){
        Eigen::MatrixXd strain = GetStrain();
        Eigen::VectorXd Strain = strain.colwise().mean();

        //[exx, phi, kyy, kzz, txy, txz] = [u,x, phi,x, v,xx, w,xx, u,y+v,x, u,z+w,x]
        theResponse << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Strain(0), Strain(4), Strain(5), Strain(1), Strain(2), Strain(3);
    }
    else if(strcasecmp(response.c_str(),"Stress") == 0){
        Eigen::MatrixXd stress = GetStress();
        Eigen::VectorXd Stress = stress.colwise().mean();

        //[Fx, Tx, Myy, Mzz, Qy, Qz] = [EA u,x, GJ phi,x, EIyy v,xx, EIzz w,xx, GAy (u,y+v,x), GAz (u,z+w,x)]
        theResponse << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Stress(0), Stress(4), Stress(5), Stress(1), Stress(2), Stress(3);
    }

    return theResponse;
}

//Computes the element energy for a given deformation.
double 
lin3DFrame2::ComputeEnergy(){
    //TODO: Integrate over element volume to compute the energy
    return 0.0;
}

//Compute the mass matrix of the element.
Eigen::MatrixXd 
lin3DFrame2::ComputeMassMatrix(){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Line", wi, xi);

    //Consistent mass definition.
    Eigen::MatrixXd MassMatrix(12,12);
    MassMatrix.fill(0.0);

    //Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
        //Gets material properties:
        Eigen::MatrixXd Rho = theSection[i]->GetDensity();

        //Compute Interpolation Function Matrix at Gauss Point:
        Eigen::MatrixXd Hij = ComputeShapeFunctionMatrix(xi(i,0));

        //Numerical integration:
        MassMatrix += wi(i)*L/2.0*Hij.transpose()*Rho*Hij;
    }

    //Lumped Mass Formulation
    if(MassFormulation){
        //Lumped Mass in diagonal terms.
        double m11 = MassMatrix(0,0) + MassMatrix(0,6);
        double m22 = MassMatrix(1,1) + MassMatrix(1,7);
        double m33 = MassMatrix(2,2) + MassMatrix(2,8);
        double m44 = MassMatrix(3,3) + MassMatrix(3,9);
        double m55 = MassMatrix(4,4);
        double m66 = MassMatrix(5,5);

        MassMatrix <<  m11,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 
                       0.0,  m22,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
                       0.0,  0.0,  m33,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
                       0.0,  0.0,  0.0,  m44,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
                       0.0,  0.0,  0.0,  0.0,  m55,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 
                       0.0,  0.0,  0.0,  0.0,  0.0,  m66,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
                       0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  m11,  0.0,  0.0,  0.0,  0.0,  0.0,
                       0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  m22,  0.0,  0.0,  0.0,  0.0,
                       0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  m33,  0.0,  0.0,  0.0, 
                       0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  m44,  0.0,  0.0,
                       0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  m55,  0.0,
                       0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  m66;
    }

    //Computes the transformation matrix.
    Eigen::MatrixXd localAxes = ComputeLocalAxes();

    //Transform Mass matrix into Global Coordinates.
    MassMatrix = localAxes.transpose()*MassMatrix*localAxes;

    return MassMatrix;
}

//Compute the stiffness matrix of the element using gauss-integration.
Eigen::MatrixXd 
lin3DFrame2::ComputeStiffnessMatrix(){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Line", wi, xi);

    //Stiffness matrix definition:
    Eigen::MatrixXd StiffnessMatrix(12, 12);
    StiffnessMatrix.fill(0.0);

    //Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
        //Compute Strain-Displacement Matrix at Gauss Point.
        Eigen::MatrixXd Bij = ComputeStrainDisplacementMatrix(xi(i,0));

        //Gets material tangent matrix at Gauss point.
        Eigen::MatrixXd Cij = theSection[i]->GetTangentStiffness();

        //Numerical integration.
        StiffnessMatrix += wi(i)*L/2.0*Bij.transpose()*Cij*Bij;
    }

    //The global axes transformation: 
    Eigen::MatrixXd localAxes = ComputeLocalAxes();

    //Stiffness matrix definition.
    StiffnessMatrix = localAxes.transpose()*StiffnessMatrix*localAxes;

    return StiffnessMatrix;
}

//Compute the damping matrix of the element.
Eigen::MatrixXd 
lin3DFrame2::ComputeDampingMatrix(){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Damping matrix definition.
    Eigen::MatrixXd DampingMatrix(12, 12);
    DampingMatrix.fill(0.0);
    
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

    //No material damping contribution is allowed.
    return DampingMatrix;
}

//Compute the PML history matrix for Perfectly-Matched Layer (PML).
Eigen::MatrixXd 
lin3DFrame2::ComputePMLMatrix(){
    Eigen::MatrixXd Kpml;
    return Kpml;
}

//Compute the element the internal forces acting on the element.
Eigen::VectorXd 
lin3DFrame2::ComputeInternalForces(){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Line", wi, xi);

    //Stiffness matrix definition:
    Eigen::VectorXd InternalForces(12);
    InternalForces.fill(0.0);

    //Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
        //Compute Strain-Displacement Matrix at Gauss Point.
        Eigen::MatrixXd Bij = ComputeStrainDisplacementMatrix(xi(i,0));

        //Gets material strain at Gauss point.
        Eigen::VectorXd Stress = theSection[i]->GetStress();

        //Numerical integration.
        InternalForces += wi(i)*L/2.0*Bij.transpose()*Stress;
    }

    //The global axes transformation.
    Eigen::MatrixXd localAxes = ComputeLocalAxes(); 

    //Internal load vector in global coordinates.
    InternalForces = localAxes.transpose()*InternalForces;

    return InternalForces;
}

//Compute the elastic, inertial, and vicous forces acting on the element.
Eigen::VectorXd 
lin3DFrame2::ComputeInternalDynamicForces(){
    //The Internal dynamic force vector
    Eigen::VectorXd InternalForces;

    if( HasFixedNode(theNodes) ){
        //Allocate memory for velocity/acceleraton. 
        Eigen::VectorXd V(12); 
        Eigen::VectorXd A(12);

        //Fills the response vectors with velocity/acceleraton values.
        V << theNodes[0]->GetVelocities(), theNodes[1]->GetVelocities();
        A << theNodes[0]->GetAccelerations(), theNodes[1]->GetAccelerations();

        //Compute the inertial/viscous/elastic dynamic force contribution.
        InternalForces = ComputeInternalForces() + ComputeDampingMatrix()*V + ComputeMassMatrix()*A;
    }

    return InternalForces;
}

//Compute the surface forces acting on the element.
Eigen::VectorXd 
lin3DFrame2::ComputeSurfaceForces(const std::shared_ptr<Load> &surfaceLoad, unsigned int UNUSED(face)){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Local surface load vector:
    Eigen::VectorXd surfaceForces(12);

    //Gets the surface force:
    Eigen::VectorXd qs = surfaceLoad->GetLoadVector();

    //Transformation matrix to local coordinates.
    Eigen::MatrixXd localAxes = ComputeLocalAxes();
    Eigen::MatrixXd rotationMatrix = localAxes.block(0,0,3,3);

    //Transform load from global to local coordinates.
    qs = rotationMatrix*qs;

    surfaceForces << qs(0)*L/2, 
                     qs(1)*L/2, 
                     qs(2)*L/2,
                     0.0,
                    -qs(2)*L*L/12.0,
                     qs(1)*L*L/12.0, 
                     qs(0)*L/2, 
                     qs(1)*L/2, 
                     qs(2)*L/2,
                     0.0,
                     qs(2)*L*L/12.0,
                    -qs(1)*L*L/12.0;

    //Node load vector in global coordinates.
    surfaceForces = localAxes.transpose()*surfaceForces;

    return surfaceForces;
}

//Compute the body forces acting on the element.
Eigen::VectorXd 
lin3DFrame2::ComputeBodyForces(const std::shared_ptr<Load> &bodyLoad, unsigned int k){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //Local body load vector:
    Eigen::VectorXd bodyForces(12);

    //Gets section density .
    Eigen::MatrixXd rho = theSection[0]->GetDensity();

    //Gets the body force:
    Eigen::VectorXd qb = rho(0,0)*bodyLoad->GetLoadVector(k);

    //Transformation matrix to local coordinates.
    Eigen::MatrixXd localAxes = ComputeLocalAxes();
    Eigen::MatrixXd rotationMatrix = localAxes.block(0,0,3,3);

    //Transform load into local coordinates.
    qb = rotationMatrix*qb;

    bodyForces << qb(0)*L/2, 
                  qb(1)*L/2, 
                  qb(2)*L/2,
                  0.0,
                 -qb(2)*L*L/12.0,
                  qb(1)*L*L/12.0, 
                  qb(0)*L/2, 
                  qb(1)*L/2, 
                  qb(2)*L/2,
                  0.0,
                  qb(2)*L*L/12.0,
                 -qb(1)*L*L/12.0;

    //Node load vector in global coordinates.
    bodyForces = localAxes.transpose()*bodyForces;

    return bodyForces;
}

//Compute the domain reduction forces acting on the element.
Eigen::VectorXd 
lin3DFrame2::ComputeDomainReductionForces(const std::shared_ptr<Load>& UNUSED(drm), unsigned int UNUSED(k)){
    //Starts profiling this function.
    PROFILE_FUNCTION();

    //TODO: Domain reduction forces not implemented for frame.
    Eigen::VectorXd DRMForces(12);
    DRMForces.fill(0.0);

    return DRMForces;
}

//Compute the current length of the element.
double
lin3DFrame2::ComputeLength() const{
    //Gets the element coordinates in undeformed configuration.  
    Eigen::VectorXd Xi = theNodes[0]->GetCoordinates();
    Eigen::VectorXd Xj = theNodes[1]->GetCoordinates();

    //Gets the element coordinates in deformed configuration. 
    Eigen::VectorXd xi = Xi + (theNodes[0]->GetDisplacements() + theNodes[0]->GetIncrementalDisplacements()).head(3);
    Eigen::VectorXd xj = Xj + (theNodes[1]->GetDisplacements() + theNodes[1]->GetIncrementalDisplacements()).head(3);

    //Local axis-1.
    Eigen::Vector3d v1;
    v1 = Xj - Xi;
    v1 = v1/v1.norm();

    //Current length of element:
    double l = v1.dot(xj - xi); 

    return l;
}

//Compute/update the local axis of the element.
Eigen::MatrixXd
lin3DFrame2::ComputeLocalAxes() const{
    //Gets the element coordinates in undeformed configuration.  
    Eigen::VectorXd Xi = theNodes[0]->GetCoordinates();
    Eigen::VectorXd Xj = theNodes[1]->GetCoordinates();

    //Local axis definition.
    Eigen::Vector3d v1;
    Eigen::Vector3d v2;
    Eigen::Vector3d v3;

    //Local axis 1.
    v1 = Xj - Xi;
    v1 = v1/v1.norm();

    //Local Axis 3.
    if(fabs(v1(2)) > TOL){
        v3 << 0.0, v1(2), -v1(1);
        v3 = v3/v3.norm();
    }
    else{
        v3 << v1(1), -v1(0), 0.0;
        v3 = v3/v3.norm();
    }

    //Local Axis 2.
    v2 = v3.cross(v1);
    v2 = v2/v2.norm();

    //The global axes transformation. 
    Eigen::MatrixXd localAxes(12, 12);

    localAxes << v1(0), v1(1), v1(2),  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
                 v2(0), v2(1), v2(2),  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
                 v3(0), v3(1), v3(2),  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
                   0.0,   0.0,   0.0, v1(0), v1(1), v1(2), 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
                   0.0,   0.0,   0.0, v2(0), v2(1), v2(2), 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
                   0.0,   0.0,   0.0, v3(0), v3(1), v3(2), 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
                   0.0,   0.0,   0.0,  0.0,   0.0,   0.0, v1(0), v1(1), v1(2), 0.0,   0.0,   0.0,
                   0.0,   0.0,   0.0,  0.0,   0.0,   0.0, v2(0), v2(1), v2(2), 0.0,   0.0,   0.0,
                   0.0,   0.0,   0.0,  0.0,   0.0,   0.0, v3(0), v3(1), v3(2), 0.0,   0.0,   0.0,
                   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0, v1(0), v1(1), v1(2),
                   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0, v2(0), v2(1), v2(2),
                   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0, v3(0), v3(1), v3(2);

    return localAxes;
}

//Update strain in the element.
Eigen::VectorXd 
lin3DFrame2::ComputeStrain(const Eigen::MatrixXd &Bij) const{
    //The global axes transformation.
    Eigen::MatrixXd localAxes    = ComputeLocalAxes(); 
    Eigen::MatrixXd rotationAxes = localAxes.block(0,0,6,6);
    
    //Gets the element displacements in local coordinates. 
    Eigen::VectorXd U1 = rotationAxes*(theNodes[0]->GetDisplacements() + theNodes[0]->GetIncrementalDisplacements());
    Eigen::VectorXd U2 = rotationAxes*(theNodes[1]->GetDisplacements() + theNodes[1]->GetIncrementalDisplacements());

    Eigen::VectorXd nodalDisplacement(12);
    nodalDisplacement << U1, U2;

    //Strain vector:
    //[u,x, phi, v,xx, w,xx, u,y+v,x, u,z+w,x]
    Eigen::VectorXd Strain(6); 
    Strain = Bij*nodalDisplacement;

    return Strain;
}

//Update strain rate in the element.
Eigen::VectorXd 
lin3DFrame2::ComputeStrainRate(const Eigen::MatrixXd &Bij) const{
    //The global axes transformation.
    Eigen::MatrixXd localAxes    = ComputeLocalAxes(); 
    Eigen::MatrixXd rotationAxes = localAxes.block(0,0,6,6);

    //Gets the element velocities in local coordinates. 
    Eigen::VectorXd V1 = rotationAxes*theNodes[0]->GetVelocities();
    Eigen::VectorXd V2 = rotationAxes*theNodes[1]->GetVelocities();

    Eigen::VectorXd nodalVelocities(12);
    nodalVelocities << V1, V2;

    //Strain rate vector: 
    Eigen::VectorXd StrainRate(6); 
    StrainRate = Bij*nodalVelocities;

    return StrainRate;
}

//Evaluates the shape function matrix at a given Gauss point.
Eigen::MatrixXd 
lin3DFrame2::ComputeShapeFunctionMatrix(const double ri) const{
    //Shape function coefficients:
    double N11 = 1.0/2.0*(1.0 - ri);
    double N21 = 1.0/2.0*(1.0 + ri);
    double N31 = 1.0/(4.0*(1.0 + Phiy))*(1.0 - ri)*(1.0 - ri)*(2.0 + ri) + Phiy/(2.0*(1.0 + Phiy))*(1.0 - ri);
    double N41 = L/(8.0*(1.0 + Phiy))*(1.0 - ri)*(1.0 - ri)*(ri + 1.0) + L*Phiy/(8.0*(1.0 + Phiy))*(1.0 - ri*ri);
    double N51 = 1.0/(4.0*(1.0 + Phiy))*(1.0 + ri)*(1.0 + ri)*(2.0 - ri) + Phiy/(2.0*(1.0 + Phiy))*(1.0 + ri);
    double N61 = L/(8.0*(1.0 + Phiy))*(1.0 + ri)*(1.0 + ri)*(ri - 1.0) - L*Phiy/(8.0*(1.0 + Phiy))*(1.0 - ri*ri);
    double N71 = 1.0/(4.0*(1.0 + Phiz))*(1.0 - ri)*(1.0 - ri)*(2.0 + ri) + Phiz/(2.0*(1.0 + Phiz))*(1.0 - ri);
    double N81 = L/(8.0*(1.0 + Phiz))*(1.0 - ri)*(1.0 - ri)*(ri + 1.0) + L*Phiz/(8.0*(1.0 + Phiz))*(1.0 - ri*ri);
    double N91 = 1.0/(4.0*(1.0 + Phiz))*(1.0 + ri)*(1.0 + ri)*(2.0 - ri) + Phiz/(2.0*(1.0 + Phiz))*(1.0 + ri);
    double N10 = L/(8.0*(1.0 + Phiz))*(1.0 + ri)*(1.0 + ri)*(ri - 1.0) - L*Phiz/(8.0*(1.0 + Phiz))*(1.0 - ri*ri);

    //Shape function matrix:
    Eigen::MatrixXd Nij(4,12);
    Nij << N11,  0.0,  0.0,  0.0,  0.0,  0.0,  N21,  0.0,  0.0,  0.0,  0.0,  0.0,
           0.0,  0.0,  0.0,  N11,  0.0,  0.0,  0.0,  0.0,  0.0,  N21,  0.0,  0.0,
           0.0,  N31,  0.0,  0.0,  0.0,  N41,  0.0,  N51,  0.0,  0.0,  0.0,  N61,
           0.0,  0.0,  N71,  0.0, -N81,  0.0,  0.0,  0.0,  N91,  0.0, -N10,  0.0;

    return Nij;
}

//Evaluates the strain-displacement matrix at a given Gauss point.
Eigen::MatrixXd 
lin3DFrame2::ComputeStrainDisplacementMatrix(const double ri) const{
    //Strain-displacement matrix coefficients for bending.
    double B11 = -1.0/L;
    double B21 =  1.0/L;
    double B31 =  6.0*ri/(L*L*(1.0 + Phiy));
    double B41 =  (3.0*ri - 1.0)/(L*(1.0 + Phiy)) - Phiy/(L*(1.0 + Phiy));
    double B51 = -6.0*ri/(L*L*(1.0 + Phiy));
    double B61 =  (3.0*ri + 1.0)/(L*(1.0 + Phiy)) + Phiy/(L*(1.0 + Phiy));
    double B71 =  6.0*ri/(L*L*(1.0 + Phiz));
    double B81 =  (3.0*ri - 1.0)/(L*(1.0 + Phiz)) - Phiz/(L*(1.0 + Phiz));
    double B91 = -6.0*ri/(L*L*(1.0 + Phiz));
    double B10 =  (3.0*ri + 1.0)/(L*(1.0 + Phiz)) + Phiz/(L*(1.0 + Phiz));

    //Strain-displacement matrix coefficients for shearing.
    double S11 = -Phiy/(L*(1.0 + Phiy));
    double S21 = -Phiy/(2.0*(1.0 + Phiy));
    double S31 =  Phiy/(L*(1.0 + Phiy));
    double S41 = -Phiy/(2.0*(1.0 + Phiy));
    double S51 = -Phiz/(L*(1.0 + Phiz));
    double S61 = -Phiz/(2.0*(1.0 + Phiz));
    double S71 =  Phiz/(L*(1.0 + Phiz));
    double S81 = -Phiz/(2.0*(1.0 + Phiz));

    //Shape function matrix:
    Eigen::MatrixXd Bij(6,12);
    Bij << B11,  0.0,  0.0,  0.0,  0.0,  0.0,  B21,  0.0,  0.0,  0.0,  0.0,  0.0,
           0.0,  0.0,  0.0,  B11,  0.0,  0.0,  0.0,  0.0,  0.0,  B21,  0.0,  0.0,
           0.0,  0.0, -B71,  0.0,  B81,  0.0,  0.0,  0.0, -B91,  0.0,  B10,  0.0,
           0.0, -B31,  0.0,  0.0,  0.0, -B41,  0.0, -B51,  0.0,  0.0,  0.0, -B61,
           0.0,  S11,  0.0,  0.0,  0.0,  S21,  0.0,  S31,  0.0,  0.0,  0.0,  S41,
           0.0,  0.0,  S51,  0.0,  S61,  0.0,  0.0,  0.0,  S71,  0.0,  S81,  0.0;

    return Bij;
}

//Compute the initial stiffness matrix of the element
Eigen::MatrixXd 
lin3DFrame2::ComputeInitialStiffnessMatrix() const{
    //Gets the quadrature information.    
    Eigen::VectorXd wi;
    Eigen::MatrixXd xi;
    QuadraturePoints->GetQuadraturePoints("Line", wi, xi);

    //Stiffness matrix definition:
    Eigen::MatrixXd StiffnessMatrix(12,12);
    StiffnessMatrix.fill(0.0);

    //Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
        //Compute Strain-Displacement Matrix at Gauss Point.
        Eigen::MatrixXd Bij = ComputeStrainDisplacementMatrix(xi(i,0));

        //Gets material tangent matrix at Gauss point.
        Eigen::MatrixXd Cij = theSection[i]->GetInitialTangentStiffness();

        //Numerical integration.
        StiffnessMatrix += wi(i)*L/2.0*Bij.transpose()*Cij*Bij;
    }

    //The global axes transformation: 
    Eigen::MatrixXd localAxes = ComputeLocalAxes();

    //Stiffness matrix definition.
    StiffnessMatrix = localAxes.transpose()*StiffnessMatrix*localAxes;

    return StiffnessMatrix;
}
