#include "Load.hpp"

//Default constructor.
Load::Load(unsigned int type) : Classification(type){
    //Does nothing.
}

//Overload constructor.
Load::Load(Eigen::VectorXd dir, std::vector<double> val, unsigned int type) : 
Classification(type), ForceDirection(dir), ForceAmplitude(val){
    //Does nothing.
}

//Virtual destructor.
Load::~Load(){
    //Does nothing.
}

//Adds node that share this loaded.
void 
Load::AddNodes(std::vector<unsigned int> tags){
    Tags = tags;
}

//Adds faces that share this loaded.
void 
Load::AddFaces(std::vector<unsigned int> tags){
    Faces = tags;
}

//Adds element that share this loaded.
void 
Load::AddElements(std::vector<unsigned int> tags){
    Tags = tags;
}

//Adds the exterior/interior condition for the domain reduction node.
void 
Load::AddDRMCondition(unsigned int tag, bool cond){
    DRMConditions[tag] = cond;
}

//Returns the load classification. 
unsigned int 
Load::GetClassification() const{
    return Classification;
}

//Returns the nodes index that share this load.
std::vector<unsigned int> 
Load::GetNodes() const{
    return Tags;
}

//Returns the faces index that share this load.
std::vector<unsigned int> 
Load::GetFaces() const{
    return Faces;
}

//Returns the nodes index that share this load.
std::vector<unsigned int> 
Load::GetElements() const{
    return Tags;
}

//Returns the force vector.
Eigen::VectorXd 
Load::GetLoadVector(const unsigned int step) const{
    return ForceAmplitude[step]*ForceDirection;
}

//Returns the exterior/interior condition for the domain reduction node.
bool 
Load::GetDRMCondition(unsigned int tag){
    return DRMConditions[tag];
}
