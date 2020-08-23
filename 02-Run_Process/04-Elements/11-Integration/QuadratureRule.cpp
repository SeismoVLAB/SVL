#include "QuadratureRule.hpp"

//Overload Constructor.
QuadratureRule::QuadratureRule(std::string name) :
Name(name) {
    //Does nothing.
}

//Destructor.
QuadratureRule::~QuadratureRule(){
    //Does nothing.
}

//Gets the Element's Quadrature Name. 
std::string 
QuadratureRule::GetQuadratureName(){
    return Name;
}
