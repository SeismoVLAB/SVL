#include "Constraint.hpp"
#include "Definitions.hpp"

//Default constructor.
Constraint::Constraint(){
    //does nothing.
}

//Overload constructor.
Constraint::Constraint(unsigned int slave, std::vector<unsigned int> master, std::vector<double> factors):
Slave(slave), Master(master), Coefficients(factors){
    //does nothing.
}

//Default destructor.
Constraint::~Constraint(){
    //Does nothing.
}

//Set the node's location.
void 
Constraint::SetSlaveInformation(unsigned int slave){
    Slave = slave;
}

//Set the free degree of freedom number list.
void 
Constraint::SetMasterInformation(std::vector<unsigned int> master){
    Master = master;
}

//Set the global degree of freedom number list.
void 
Constraint::SetCombinationFactors(std::vector<double> factors){
    Coefficients = factors;
}

//Return the slave free degree-of-freedom of this constraint.
unsigned int 
Constraint::GetSlaveInformation() const{
    return Slave;
}

//Return the number of constraint applied to this slave degree-of-freedom.
unsigned int 
Constraint::GetNumberOfConstraints() const{
    return Coefficients.size();
}

//Return the master list of total degree-of-freedom to be combined.
const std::vector<unsigned int> 
Constraint::GetMasterInformation() const{
    return Master;
}

//Return the combinational factor list for each degree-of-freedom.
const std::vector<double> 
Constraint::GetCombinationFactors() const{
    return Coefficients;
}