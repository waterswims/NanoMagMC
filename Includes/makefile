# Makefile for project

OS := $(shell uname)
HOST := $(shell hostname)
INC_PATH = ../Includes

# Check if mac
ifeq ($(OS),Darwin)
CXX = mpic++
LDFLAGS =  -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -ipo
AVXCHECK := $(shell sysctl -a | grep AVX2)
else
AVXCHECK := $(shell grep avx2 /proc/cpuinfo)
# Check if ARCHER
ifneq (,$(findstring eslogin, $(HOST)))
CXX = CC
LDFLAGS = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -ipo
else
CXX = mpiicpc
LDFLAGS = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -ipo
endif
endif

# Check for systems I know have AVX2
ifneq (,$(AVXCHECK))
CPPFLAGS = -O3 -I${MKLROOT}/include -xCORE-AVX2 -ipo
else
CPPFLAGS = -O3 -I${MKLROOT}/include -ipo
endif

OBJS = functions.o cpoints.o output.o param_read.o mklrand.o array_alloc.o field_type.o hamiltonian.o state.o shape.o mpifuncs.o
TEST_FILES = cpoint_test.hpp mklrand_test.hpp ising_test.hpp heis_test.hpp fept_test.hpp test_functions.hpp skyrm_test.hpp

all: $(OBJS) main.o

mpifuncs.o: mpifuncs.cpp mpifuncs.hpp array_alloc.hpp
	$(CXX) $(CPPFLAGS) -c mpifuncs.cpp

array_alloc.o: array_alloc.cpp array_alloc.hpp
	$(CXX) $(CPPFLAGS) -c array_alloc.cpp

field_type.o: field_type.cpp field_type.hpp array_alloc.hpp mklrand.h hamiltonian.hpp
	$(CXX) $(CPPFLAGS) -c field_type.cpp

param_read.o: param_read.cpp param_read.hpp array_alloc.hpp
	$(CXX) $(CPPFLAGS) -c param_read.cpp

mklrand.o: mklrand.cpp mklrand.h
	$(CXX) $(CPPFLAGS) -c mklrand.cpp

shape.o: shape.cpp shape.hpp
	$(CXX) $(CPPFLAGS) -c shape.cpp

hamiltonian.o: hamiltonian.cpp hamiltonian.hpp field_type.hpp array_alloc.hpp functions.h
	$(CXX) $(CPPFLAGS) -c hamiltonian.cpp

state.o: state.cpp state.hpp functions.h mklrand.h hamiltonian.hpp field_type.hpp shape.hpp
	$(CXX) $(CPPFLAGS) -c state.cpp

functions.o: functions.cpp functions.h
	$(CXX) $(CPPFLAGS) -c functions.cpp

cpoints.o: cpoints.cpp cpoints.hpp functions.h
	$(CXX) $(CPPFLAGS) -c cpoints.cpp

output.o: output.cpp output.hpp functions.h
	$(CXX) $(CPPFLAGS) -c output.cpp

main.o: main.cpp state.hpp mklrand.h functions.h cpoints.hpp output.hpp array_alloc.hpp mpifuncs.hpp
	$(CXX) $(CPPFLAGS) -c main.cpp

run: $(OBJS) main.o
	$(CXX) $(OBJS) main.o -o run $(LDFLAGS)

test.o: test.cpp mklrand.h cpoints.hpp hamiltonian.hpp field_type.hpp $(TEST_FILES)
	$(CXX) $(CPPFLAGS) -c test.cpp -lgtest

test: $(OBJS) test.o
	$(CXX) $(OBJS) test.o -o test -lgtest $(LDFLAGS)

1skrm.o: 1skrm.cpp state.hpp mklrand.h functions.h output.hpp array_alloc.hpp mpifuncs.hpp
	$(CXX) $(CPPFLAGS) -c 1skrm.cpp

skyrm_print: $(OBJS) 1skrm.o
	$(CXX) $(OBJS) 1skrm.o -o skyrm_print $(LDFLAGS)

clean:
	-rm *.o
	-rm run
	-rm test
