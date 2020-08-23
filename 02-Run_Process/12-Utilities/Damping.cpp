#include <iostream>
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
