#include "Integrator.hpp"
#include "Definitions.hpp"

//Default constructor.
Integrator::Integrator(const std::shared_ptr<Mesh> &mesh){
    //Operator that enforced restrain/constraint. 
    Total2FreeMatrix = mesh->GetTotalToFreeMatrix();

    //Nodal support motion vector.
    SupportMotion.resize(numberOfTotalDofs);
    SupportMotion.fill(0.0);
}

//Virtual destructor.
Integrator::~Integrator(){
    //Does nothing.
}