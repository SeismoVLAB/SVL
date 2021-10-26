#include "Damping.hpp"

//Overload Constructor.
Damping::Damping(std::string name, const std::vector<double> parameters) :
  Name(name), Parameters(parameters){
    //Does nothing.
}

//Destructor.
Damping::~Damping(){
    //Does nothing.
}

//Gets damping model's name
std::string 
Damping::GetName(){
    return Name;
}

//Gets damping parameters vector
std::vector<double>
Damping::GetParameters(){
  return Parameters;
}

//Set the name of the damping model.
void
Damping::SetName(std::string name){
    Name = name;
}

//Set the damping parameters.
void 
Damping::SetParameters(std::vector<double> param){
    Parameters = param;
}
