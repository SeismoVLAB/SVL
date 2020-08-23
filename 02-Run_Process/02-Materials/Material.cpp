#include "Material.hpp"
#include "Definitions.hpp"

//Default constructor.
Material::Material(std::string name, bool model) : 
Name(name), isViscous(model){
    //Does nothing.
}

//Virtual destructor.
Material::~Material(){
    //Does nothing.
}

//Gets the Material's Name. 
std::string 
Material::GetName(){
    return Name;
}

//Gets the Material's stress model. 
bool
Material::IsViscous(){
    return isViscous;
}
