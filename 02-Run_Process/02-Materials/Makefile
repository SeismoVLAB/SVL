COND = $(DEBUG)

#INCLUDE RULE FOR COMPILATION OF LOCAL FILES.
include ../Makefile.mk

#FINDS THE FILES TO BE COMPILED:
SRC := $(shell find . -maxdepth 1 -mindepth 1 -iname  "*.cpp") $(shell find . -maxdepth 1 -mindepth 1 -iname  "*.cc") $(shell find . -maxdepth 1 -mindepth 1 -iname  "*.cu")
OBJ := $(SRC:.cpp=.o)
OBJ := $(OBJ:.cc=.o)
OBJ := $(OBJ:.cu=.o)

#COMPILE THE LOCAL FILES.
all: $(OBJ)

%.o: %.cpp
	@echo "$(GREEN)[Compiling:]$(YELLOW)" $< "$(NC)"
	$(MPICC) $(CFLAGS) $(EIGFLAGS) -I$(EIGEN_DIR) -I$(MPI_DIR) $(DIRS) -c $?

%.o: %.cu
	@echo "$(GREEN)[Compiling:]$(YELLOW)" $< "$(NC)"
	$(NVCC) $(NVFLAGS) $(EIGFLAGS) -I$(EIGEN_DIR) -I$(MPI_DIR) $(DIRS) -c $?

%.o: %.cc
	@echo "$(GREEN)[Compiling:]$(YELLOW)" $< "$(NC)"
	${PETSC_COMPILE} $(EIGFLAGS) -I$(EIGEN_DIR) -I$(PETSC_DIR)/include $(DIRS) -c $?
