#include <vector>
#include <string>

//The processor number.
int rank;

//The number of partitions.
int size;

//The problem dimension (1D, 2D, 3D). 
unsigned int nDimensions;

//The folder where the input file is loaded.
std::string filePath;

//The input file name (.$.) to be loaded.
std::vector<std::string> fileName;

///Whether the driver (JSON) file is provided  
bool driverFile;

//The element mass formulation.
bool MassFormulation;

//The update option for internal variables in Mesh.
std::string UpdateOption;

//Maximum memory for lumped storage sparse matrix.
unsigned int LumpedStorage;

//Maximum memory for consistent storage sparse matrix.
unsigned int ConsistentStorage;

//Total number of free-degree-of-freedom.
unsigned int numberOfFreeDofs;

//Total number of total-degree-of-freedom.
unsigned int numberOfTotalDofs;

//Total number of constrained-degree-of-freedom.
unsigned int numberOfConstrainedDofs;

//Maximum memory for PML 3D storage sparse matrix.
unsigned int PMLStorage;
