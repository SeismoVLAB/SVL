#include "Element.hpp"
#include "Definitions.hpp"

//Default Constructor.
Element::Element(std::string name, const std::vector<unsigned int> nodes, unsigned int ndofs, unsigned int VTKcell, unsigned int SVLcell) : 
Name(name), VTKCell(VTKcell), SVLCell(SVLcell), NumberOfNodes(nodes.size()), NumberOfDegreeOfFreedom(ndofs), Nodes(nodes){
    //Does nothing.
}

//Virtual Destructor.
Element::~Element(){
    //Does nothing.
}

//Gets the Element's Name. 
std::string 
Element::GetName() const{
    return Name;
}

//Gets the Element VTK cell type.
unsigned int 
Element::GetVTKCellType() const{
    return VTKCell;
}

//Gets the Element SeismoVLAB cell type.
unsigned int 
Element::GetVTKGroupType() const{
    return SVLCell;
}

//Returns the number of nodes in element.
unsigned int 
Element::GetNumberOfNodes() const{
    return NumberOfNodes;
}

//Returns total number of degree of freedom in the element.
unsigned int 
Element::GetNumberOfDegreeOfFreedom() const{
    return NumberOfDegreeOfFreedom;
}

//Returns the Node Connectivity Indexes.
std::vector<unsigned int> 
Element::GetNodes() const{
    return Nodes;
}

//Returns if the element has fixed nodes.
bool 
Element::HasFixedNode(const std::vector<std::shared_ptr<Node> > &nodes) const{
    //Check if this element is restrained.
    bool IsFixed = false;
    for(unsigned int k = 0; k < NumberOfNodes; k++){
        if( nodes[k]->IsFixed() ){
            IsFixed = true;
            break;
        }
    }

    return IsFixed;
}