# Makefile for project

#################################################################
## Location of the HDF5 instalation
#################################################################
HDFLIB=/home/jmw2g14/hdf5/lib/
HDFINC=/home/jmw2g14/hdf5/include/
## IRIDIS BUILD - Using intel compilers
#HDFLIB=/local/software/gdf5/1.8.11/intel-par/lib
#HDFINC=/local/software/gdf5/1.8.11/intel-par/include

#################################################################
## Directories
#################################################################

INC_PATH = includes
LIB_PATH = lib
OBJ_PATH = obj
TEST_PATH = tests

#################################################################
## Files to use
#################################################################

TEST_FILES = $(wildcard $(TEST_PATH)/*.hpp)
NOMAIN_FILES = $(filter-out $(LIB_PATH)/main.cpp, $(wildcard $(LIB_PATH)/*.cpp))
## NON-INTEL
SOURCE_FILES = $(filter-out $(LIB_PATH)/mklrand.hpp, $(NOMAIN_FILES))
## INTEL
# SOURCE_FILES = $(filter-out $(LIB_PATH)/stdrand.hpp, $(NOMAIN_FILES))
OBJS = $(addprefix $(OBJ_PATH)/, $(notdir $(SOURCE_FILES:.cpp=.o)))

#################################################################
## Compile options
#################################################################
## INTEL
# CPPFLAGS = -std=c++11 -Ofast -qopenmp -DMKL_ILP64 -I${MKLROOT}/include -I${HDFINC} -ipo
## NON-INTEL
CPPFLAGS = -std=c++11 -Ofast -fopenmp -I${HDFINC} -Wno-narrowing -flto

#################################################################
## Link options
#################################################################
## INTEL
# LDFLAGS = -L${HDFLIB}  -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl -lhdf5 -lz -ipo
## NON-INTEL
LDFLAGS = -L${HDFLIB} -lhdf5 -lz -flto

#################################################################
## Build scripts
#################################################################

default: $(OBJS) $(OBJ_PATH)/main.o
	$(CC) $(OBJS) $(OBJ_PATH)/main.o -o run $(LDFLAGS)

$(OBJ_PATH)/%.o: $(LIB_PATH)/%.cpp
	$(CC) $(CPPFLAGS) -c -o $@ $<

$(TEST_PATH)/test.o: $(TEST_FILES) tests/test.cpp
	$(CC) tests/test.cpp $(CPPFLAGS) -c -o $(TEST_PATH)/test.o

test: $(OBJS) $(TEST_PATH)/test.o
	$(CC) $(OBJS) $(TEST_PATH)/test.o $(LDFLAGS) -lgtest -o test

clean:
	-rm -f $(OBJS)
	-rm -f $(OBJ_PATH)/main.o
	-rm -f run
	-rm -f $(TEST_PATH)/test.o
	-rm -f test
