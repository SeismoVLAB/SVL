#include <string>

//The processor number.
int rank;

//The number of partitions.
int size;

//The problem dimension (1D, 2D, 3D). 
unsigned int nDimensions;

//The input file name (.$.) to be loaded.
std::string fileName;

//The folder where the input file is loaded.
std::string filePath;

//The Execution for each simulation.
bool FormOfExecution;

//The element mass formulation.
bool MassFormulation;

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
