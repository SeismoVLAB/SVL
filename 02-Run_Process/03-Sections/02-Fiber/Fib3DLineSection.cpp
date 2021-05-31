#include <cmath>
#include "Fib3DLineSection.hpp"
#include "Definitions.hpp"

//Define constant value PI:
const double PI = 3.1415926535897932;

//Overload Constructor.
Fib3DLineSection::Fib3DLineSection(double h, double b, const std::vector<std::unique_ptr<Material> > &fibers, const std::vector<double> &zi, const std::vector<double> &yi, const std::vector<double> &Ai, double k2, double k3, double theta, unsigned int ip) : 
Section("Fib3DLineSection"), h(h), b(b), kappa2(k2), kappa3(k3), Theta(theta), InsertPoint(ip), zi(zi), yi(yi), Ai(Ai){
    //Initialize material strain.
    Strain.resize(6);
    Strain << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

    //The section fibers:
    NumberOfFibers = fibers.size();
    theFiber.resize(NumberOfFibers);
    for(unsigned int i = 0; i < NumberOfFibers; i++)
        theFiber[i] = fibers[i]->CopyMaterial();

    //Compute the section centroid.
    ComputeSectionCenter();

    //Transform the rotation angle into radians.
    Theta = PI*theta/180.0;
}

//Destructor.
Fib3DLineSection::~Fib3DLineSection(){
    //Does nothing.
}

//Clone the 'Fib3DLineSection' section.
std::unique_ptr<Section>
Fib3DLineSection::CopySection(){
    return std::make_unique<Fib3DLineSection>(h, b, theFiber, zi, yi, Ai, kappa2, kappa3, 180.0*Theta/PI, InsertPoint);
}

//Returns the section generalized strain in element coordinates.
Eigen::VectorXd
Fib3DLineSection::GetStrain(){
    return Strain;
}

//Returns the section generalized stress in element coordinates.
Eigen::VectorXd
Fib3DLineSection::GetStress(){
    //Section stress vector for flexural and shear
    Eigen::VectorXd Fm(3); Fm.fill(0.0);
    Eigen::VectorXd Fv(3); Fv.fill(0.0);
    
    for(unsigned int k = 0; k < NumberOfFibers; k++){
        Eigen::VectorXd Sk = theFiber[k]->GetStress();
        Fm(0) += Sk(0)*Ai[k];
        Fm(1) += Sk(0)*Ai[k]*(zi[k] - zcm);
        Fm(2) -= Sk(0)*Ai[k]*(yi[k] - ycm);

        double Gk = theFiber[k]->GetShearModulus();
        Fv(0) += Gk*Ai[k]*((zi[k] - zcm)*(zi[k] - zcm) + (yi[k] - ycm)*(yi[k] - ycm))*Strain(1);
        Fv(1) += kappa2*Gk*Ai[k]*Strain(4);
        Fv(2) += kappa3*Gk*Ai[k]*Strain(5);
    }

    //Compute the rotation matrix.
    Eigen::MatrixXd T = GetLineRotationMatrix(Theta);

    //Compute the translation matrix.
    Eigen::MatrixXd L = GetLineTranslationMatrix(h, b, zcm, ycm, InsertPoint);

    //Computes the generalized section force in element coordinates.
    Fv = T.transpose()*Fv;
    Fm = T.transpose()*L.transpose()*Fm;

    Eigen::VectorXd Stress(6);
    Stress << Fm(0), Fv(0), Fm(1), Fm(2), Fv(1), Fv(2);

    return Stress;
}

//Returns the section density matrix.
Eigen::MatrixXd
Fib3DLineSection::GetDensity(){
    //Section area Properties
    double A = 0.0;
    double J = 0.0;
    double Rho = 0.0;

    //Numerical integration of the section's mass properties
    for(unsigned int k = 0; k < NumberOfFibers; k++){
        double ri = theFiber[k]->GetDensity();
        double Gi = theFiber[k]->GetShearModulus();
        A += Gi*Ai[k];
        J += Gi*Ai[k]*((zi[k] - zcm)*(zi[k] - zcm) + (yi[k] - ycm)*(yi[k] - ycm));
        Rho += Ai[k]*ri;
    }

    //Returns the section density matrix for 3-dimensions.
    Eigen::MatrixXd SectionDensity(4,4);
    SectionDensity <<  Rho, 0.0, 0.0, 0.0,
                       0.0, J/A, 0.0, 0.0,
                       0.0, 0.0, Rho, 0.0,
                       0.0, 0.0, 0.0, Rho;

    return SectionDensity;
}

//Returns the section axial stiffness.
Eigen::MatrixXd
Fib3DLineSection::GetTangentStiffness(){
    //Section axial properties.
    double EA = 0.0;
    double GJ = 0.0;
    double GA2 = 0.0;
    double GA3 = 0.0;

    //Section flexural properties.
    double EI22 = 0.0;
    double EI33 = 0.0;
    double EI23 = 0.0;

    //Numerical integration of the section properties
    for(unsigned int k = 0; k < NumberOfFibers; k++){
        double Gk = theFiber[k]->GetShearModulus();
        GJ += Gk*Ai[k]*((zi[k] - zcm)*(zi[k] - zcm) + (yi[k] - ycm)*(yi[k] - ycm));
        GA2 += kappa2*Gk*Ai[k];
        GA3 += kappa3*Gk*Ai[k];

        Eigen::MatrixXd Ek = theFiber[k]->GetTangentStiffness();
        EA += Ek(0,0)*Ai[k];
        EI22 += Ek(0,0)*Ai[k]*(zi[k] - zcm)*(zi[k] - zcm);
        EI33 += Ek(0,0)*Ai[k]*(yi[k] - ycm)*(yi[k] - ycm);
        EI23 -= Ek(0,0)*Ai[k]*(zi[k] - zcm)*(yi[k] - ycm);
    }

    //Section shear Stiffness at the centroid.
    Eigen::MatrixXd Cs(3,3);
    Cs << GJ , 0.0, 0.0,
          0.0, GA2, 0.0,
          0.0, 0.0, GA3;

    //Section Flexural Stiffness at the centroid.
    Eigen::MatrixXd Cm(3,3);
    Cm << EA ,  0.0,  0.0,
          0.0, EI22, EI23,
          0.0, EI23, EI33;

    //Compute the rotation matrix.
    Eigen::MatrixXd T = GetLineRotationMatrix(Theta);

    //Compute the translation matrix.
    Eigen::MatrixXd L = GetLineTranslationMatrix(h, b, zcm, ycm, InsertPoint);

    //Computes the global section stiffness matrix.
    Cs = T.transpose()*Cs*T;
    Cm = T.transpose()*L.transpose()*Cm*L*T;

    //Returns the section stiffness for 3-dimensions.
    Eigen::MatrixXd SectionStiffness(6,6);
    SectionStiffness << Cm(0,0),   0.0  , Cm(0,1), Cm(0,2),   0.0  ,   0.0  ,
                          0.0  , Cs(0,0),   0.0  ,   0.0  , Cs(0,1), Cs(0,2),
                        Cm(0,1),   0.0  , Cm(1,1), Cm(2,1),   0.0  ,   0.0  ,
                        Cm(0,2),   0.0  , Cm(2,1), Cm(2,2),   0.0  ,   0.0  ,
                          0.0  , Cs(0,1),   0.0  ,   0.0  , Cs(1,1), Cs(2,1),
                          0.0  , Cs(0,2),   0.0  ,   0.0  , Cs(2,1), Cs(2,2);

    return SectionStiffness;
}

//Returns the section initial tangent stiffness matrix.
Eigen::MatrixXd 
Fib3DLineSection::GetInitialTangentStiffness(){
//Section axial properties.
    double EA = 0.0;
    double GJ = 0.0;
    double GA2 = 0.0;
    double GA3 = 0.0;

    //Section flexural properties.
    double EI22 = 0.0;
    double EI33 = 0.0;
    double EI23 = 0.0;

    //Numerical integration of the section properties
    for(unsigned int k = 0; k < NumberOfFibers; k++){
        double Gk = theFiber[k]->GetShearModulus();
        GJ += Gk*Ai[k]*((zi[k] - zcm)*(zi[k] - zcm) + (yi[k] - ycm)*(yi[k] - ycm));
        GA2 += kappa2*Gk*Ai[k];
        GA3 += kappa3*Gk*Ai[k];

        Eigen::MatrixXd Ek = theFiber[k]->GetInitialTangentStiffness();
        EA += Ek(0,0)*Ai[k];
        EI22 += Ek(0,0)*Ai[k]*(zi[k] - zcm)*(zi[k] - zcm);
        EI33 += Ek(0,0)*Ai[k]*(yi[k] - ycm)*(yi[k] - ycm);
        EI23 -= Ek(0,0)*Ai[k]*(zi[k] - zcm)*(yi[k] - ycm);
    }

    //Section shear Stiffness at the centroid.
    Eigen::MatrixXd Cs(3,3);
    Cs << GJ , 0.0, 0.0,
          0.0, GA2, 0.0,
          0.0, 0.0, GA3;

    //Section Flexural Stiffness at the centroid.
    Eigen::MatrixXd Cm(3,3);
    Cm << EA ,  0.0,  0.0,
          0.0, EI22, EI23,
          0.0, EI23, EI33;

    //Compute the rotation matrix.
    Eigen::MatrixXd T = GetLineRotationMatrix(Theta);

    //Compute the translation matrix.
    Eigen::MatrixXd L = GetLineTranslationMatrix(h, b, zcm, ycm, InsertPoint);

    //Computes the global section stiffness matrix.
    Cs = T.transpose()*Cs*T;
    Cm = T.transpose()*L.transpose()*Cm*L*T;

    //Returns the section stiffness for 3-dimensions.
    Eigen::MatrixXd SectionStiffness(6,6);
    SectionStiffness << Cm(0,0),   0.0  , Cm(0,1), Cm(0,2),   0.0  ,   0.0  ,
                          0.0  , Cs(0,0),   0.0  ,   0.0  , Cs(0,1), Cs(0,2),
                        Cm(0,1),   0.0  , Cm(1,1), Cm(2,1),   0.0  ,   0.0  ,
                        Cm(0,2),   0.0  , Cm(2,1), Cm(2,2),   0.0  ,   0.0  ,
                          0.0  , Cs(0,1),   0.0  ,   0.0  , Cs(1,1), Cs(2,1),
                          0.0  , Cs(0,2),   0.0  ,   0.0  , Cs(2,1), Cs(2,2);

    return SectionStiffness;
}

//Returns the section strain at given position in section local coordinates.
Eigen::VectorXd 
Fib3DLineSection::GetStrainAt(double x3, double x2){
    //Transform coordinates into section local axis.
    x3 = zcm - x3;
    x2 = x2 - ycm;

    //Finds the closer fiber to (x3,x2) coordinate
     unsigned int k = FiberIndexAt(x3, x2);

    //Gets the strain at the fiber
    Eigen::VectorXd exx = theFiber[k]->GetStrain();

    //TODO: Include shearing strains exy, exz in future
    //Epsilon = [exx, 0.0, 0.0, exy, 0.0, exz]
    Eigen::VectorXd theStrain(6);
    theStrain << exx(0), 
                0.0, 
                0.0, 
                0.0, 
                0.0,
                0.0;

    return theStrain;
}

//Returns the section stress at given position in section local coordinates.
Eigen::VectorXd 
Fib3DLineSection::GetStressAt(double x3, double x2){
    //Transform coordinates into section local axis.
    x3 = zcm - x3;
    x2 = x2 - ycm;

    //Finds the closer fiber to (x3,x2) coordinate
    unsigned int m = FiberIndexAt(x3, x2);

    //Gets the strain at the fiber
    Eigen::VectorXd Sxx = theFiber[m]->GetStress();

    //TODO: Include shearing stresses txy, txz in future
    //Sigma = [Sxx, 0.0, 0.0, txy, 0.0, txz]
    Eigen::VectorXd theStress(6);
    theStress << Sxx(0), 
                0.0, 
                0.0, 
                0.0, 
                0.0,
                0.0;

    return theStress;
}

//Perform converged section state update.
void 
Fib3DLineSection::CommitState(){
    for(unsigned int k = 0; k < NumberOfFibers; k++)
        theFiber[k]->CommitState();
}

//Reverse the section states to previous converged state.
void 
Fib3DLineSection::ReverseState(){
    for(unsigned int k = 0; k < NumberOfFibers; k++)
        theFiber[k]->ReverseState();
}

//Brings the section states to its initial state.
void 
Fib3DLineSection::InitialState(){
    Strain.fill(0.0);
    for(unsigned int k = 0; k < NumberOfFibers; k++)
        theFiber[k]->InitialState();
}

//Update the section state for this iteration.
void
Fib3DLineSection::UpdateState(Eigen::VectorXd strain, unsigned int cond){
    //Transform the strain in element coordinates into section coordinates
    Eigen::VectorXd sstrain = ComputeLineLocalAxes(h, b, zcm, ycm, Theta, InsertPoint)*strain;

    //Update the fiber behavior.
    for(unsigned int k = 0; k < NumberOfFibers; k++){
        Eigen::VectorXd FibStrain(1);
        FibStrain << (sstrain(0) + (zi[k] - zcm)*sstrain(2) - (yi[k] - ycm)*sstrain(3));
        theFiber[k]->UpdateState(FibStrain, cond);
    }

    //Update the section strain.
    Strain = strain;
}

//Gets the section centroid.
void
Fib3DLineSection::ComputeSectionCenter(){
    //Cross Section-Centroid.
    double ziAi = 0.0;
    double yiAi = 0.0;
    double Area = 0.0;
    for(unsigned int k = 0; k < NumberOfFibers; k++){
        Eigen::MatrixXd Ek = theFiber[k]->GetInitialTangentStiffness();
        Area += Ek(0,0)*Ai[k];
        ziAi += Ek(0,0)*zi[k]*Ai[k];
        yiAi += Ek(0,0)*yi[k]*Ai[k];
    }

    //Centroid coordinates
    ycm = yiAi/Area;
    zcm = ziAi/Area;
}

//Return the fiber strain closer to (x3, x2)
unsigned int 
Fib3DLineSection::FiberIndexAt(double x3, double x2){
    unsigned int ind = 0;
    std::vector<double> v(NumberOfFibers);

    //Loop over all fibers to see which one is closer to (x3, x2)
    v[0] = sqrt((x3 - (zcm - zi[0]))*(x3 - (zcm - zi[0])) + (x2 - (yi[0] - ycm))*(x2 - (yi[0] - ycm)));
    for(unsigned int k = 1; k < NumberOfFibers; k++){
        v[k] = sqrt((x3 - (zcm - zi[k]))*(x3 - (zcm - zi[k])) + (x2 - (yi[k] - ycm))*(x2 - (yi[k] - ycm)));
        if (v[k] < v[ind])
            ind = k;
    }
    
    return ind;
}