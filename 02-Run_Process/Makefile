#==============================================================================
#
#                    (Seismo)s (V)irtual (Lab)oratory
#             Module for Serial and Parallel Analysis of seismic 
#         wave propagation and soil-structure interaction simulation
#         Copyright (C) 2018, The California Institute of Technology
#                           All Rights Reserved.
#
# Commercial use of this program without express permission of the California
# Institute of Technology, is strictly  prohibited. See  file "COPYRIGHT"  in
# main  directory  for  information on  usage  and  redistribution, and for a
# DISCLAIMER OF ALL WARRANTIES.
#
#==============================================================================
#
# Written by:
#   Danilo S. Kusanovic (dkusanov@caltech.edu)
#
# Version     :
# Date        : 12-12-2017
# Modified	  : 02-11-2019
# Source      : Makefile
# Keywords    : <Compilation, Make, Executable>
# Description : Compiles the Module for Serial/Parallel analysis.
#   
# Example of Usage :
#  make              (compile in debug mode)
#  make -w DEBUG=YES (compile in debug mode and no print directory)
#  make -s DEBUG=NO  (compile in optimized mode and silent mode)
#------------------------------------------------------------------------------

#DEBUGGING's CONDITION.
COND := $(DEBUG)

#PROGRAM's NAME.
EXE	 := SeismoVLab.exe

#INCLUDE RULE FOR COMPILATION.
include Makefile.mk

#LIST OF DIRECTORIES TO BE INCLUDED.
DIRS := $(patsubst %, -I%, $(dir $(shell find $$(pwd)/ -type f -iname 'Makefile')))
export DIRS:=$(DIRS)

#LIST OF SUB-DIRECTORIES TO MAKE. 
SUBDIRS := $(dir $(shell find $$(pwd)/ -mindepth 2 -type f -iname 'Makefile'))

#ALL THE SOURCE AND OBJECT FILES.
SRC := $(shell find . -iname  "*.cpp") $(shell find . -iname  "*.cc") $(shell find . -iname  "*.cu")
OBJ := $(SRC:.cpp=.o)
OBJ := $(OBJ:.cc=.o)
OBJ := $(OBJ:.cu=.o)

#COMPILE THE PROGRAM.
all: $(EXE)

$(EXE): $(SUBDIRS)
	@echo "$(BLUE)Folder:$(YELLOW)" $@ "$(BLUE)was successfully compiled.$(NC)"
	${PETSC_COMPILE} $(EIGFLAGS) -I$(EIGEN_DIR) -I$(PETSC_DIR)/include $(DIRS) -c main.cpp
	@echo "$(GREEN)Linking the $(EXE) project$(NC)"
	$(MPICC) $(CFLAGS) $(EIGFLAGS) -DAdd__ -fopenmp $(OBJ) ${PETSC_KSP_LIB} -I$(EIG_DIR) -I$(PETSC_DIR)/include $(LPATH) $(LIBS) -o $@
	@echo "$(GREEN)SUCCESSFUL LINKAGE. The executable$(RED) $(EXE) $(GREEN)was generated.$(NC)"

.PHONY: $(SUBDIRS) all wipe clean

$(SUBDIRS):
	$(MAKE) -C $@
	@echo "$(BLUE)Folder:$(YELLOW)" $@ "$(BLUE)was successfully compiled.$(NC)"

doc:
	doxygen

wipe:
	@echo "$(RED)The following list of objects were removed.$(NC)"
	rm -rf $(OBJ)

clean:
	@echo "$(RED)The executable and all objects were removed.$(NC)"
	rm -rf $(OBJ) $(EXE) ../04-Documentation/html ../04-Documentation/latex
