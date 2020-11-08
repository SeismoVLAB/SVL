#include <ctime>
#include <string>
#include <fstream> 
#include <sstream>  
#include <iostream>

//The processor number.
int rank;

//The number of partitions.
int size;

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
