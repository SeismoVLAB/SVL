#COMPILERS DEFINITION.  
CC    := g++
NVCC  := nvcc
MPICC := mpic++

#LIBRARIES USED TO BE INCLUDED IN COMPILATION.
MPI_DIR    = /usr/include/mpi
EIGEN_DIR  = /home/danilo/Libraries/Eigen
PETSC_DIR  = /home/danilo/Libraries/Petsc
MUMPS_DIR  = /home/danilo/Libraries/Mumps
METIS_DIR  = /usr/lib/x86_64-linux-gnu
SCOTCH_DIR = /usr/lib/x86_64-linux-gnu

#INCLUDE PETSC VARIABLES.
include ${PETSC_DIR}/lib/petsc/conf/variables

#ORDERING LIBRARIES.
LPORD    = -L$(MUMPS_DIR)/PORD -lpord
LMETIS   = -L$(METIS_DIR) -lparmetis -lmetis
LSCOTCH  = -L$(SCOTCH_DIR) -lptesmumps -lptscotch -lptscotcherr
LORDERS  = $(LMETIS) $(LPORD) $(LSCOTCH)

#MUMPS LIBRARY.
LIBMUMPS  = -L$(MUMPS_DIR) -ldmumps -lmumps_common -lpord

#OTHER LIBRARY.
LIBSPAR   = -lscalapack-openmpi -lblacs-openmpi -llapack
LIBOTHERS = -lblas -lpthread

#LIBRARY PATH OPTIONS.
LPATH = 
LIBS  = $(LORDERS) $(LIBMUMPS) $(LIBSPAR) $(LIBOTHERS) 

#COMPILATION MODE. 
ifneq (, $(filter $(DEBUG), NO FALSE DISABLE No False Disable no false disable))
	CFLAGS   = -O3 -m64 -fopenmp -std=c++14
	EIGFLAGS = -DEIGEN_NO_DEBUG
	NVFLAGS  = --ptxas-options=-O3 -m$(ARCH) -arch compute_$(NVARCH) -code sm_$(NVARCH)
else
	CFLAGS   = -Wall -Wextra -g -fopenmp -std=c++14
	EIGFLAGS = -DEIGEN_DEBUG
	NVFLAGS  = -g -G -m$(ARCH) -arch compute_$(NVARCH) -code sm_$(NVARCH)
endif

#COLORS FOR DISPLAY.
RED    := \033[1;31m
GREEN  := \033[1;32m
BLUE   := \033[1;34m
YELLOW := \033[1;33m
NC     := \033[1;0m
