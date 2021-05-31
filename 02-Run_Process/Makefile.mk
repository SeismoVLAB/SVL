#COMPILERS DEFINITION.  
CC    := g++
NVCC  := nvcc
MPICC := mpic++

#LIBRARIES USED TO BE INCLUDED IN COMPILATION.
EIGEN_DIR = /usr/include/eigen
PETSC_DIR = /usr/lib/petsc
MPI_DIR   = /usr/include/mpi

#INCLUDE PETSC VARIABLES.
include ${PETSC_DIR}/lib/petsc/conf/variables

#ORDERING LIBRARIES.
LPORD   = -lpord
LMETIS  = -lparmetis -lmetis
LSCOTCH = -lptesmumps -lptscotch -lptscotcherr
LORDERS = $(LMETIS) $(LPORD) $(LSCOTCH)

#MUMPS LIBRARY.
LIBMUMPS  = -ldmumps -lmumps_common -lpord

#OTHER LIBRARY.
LIBSPAR   = 
LIBOTHERS = 

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
