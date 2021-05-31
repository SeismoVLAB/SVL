#include "LoadCombo.hpp"
#include "Definitions.hpp"

//Default constructor:
LoadCombo::LoadCombo(std::string name, std::vector<unsigned int> loads, std::vector<double> factors) : 
Name(name), Loads(loads), Factors(factors){
    //Does nothing.
}

//Virtual destructor:
LoadCombo::~LoadCombo(){
    //Does nothing.
}

//Returns the name of the load combination. 
std::string 
LoadCombo::GetCombinationName() const{
    return Name;
}

//Returns the loads factors.
std::vector<double> 
LoadCombo::GetLoadFactors() const{
    return Factors;
}

//Returns the loads to combine.
std::vector<unsigned int> 
LoadCombo::GetLoadCombination() const{
    return Loads;
}