# Makefile for project

OS := $(shell uname)

ifeq ($(OS),Darwin)
CXX = mpic++
LDFLAGS =  ${MKLROOT}/lib/libmkl_intel_lp64.a ${MKLROOT}/lib/libmkl_sequential.a ${MKLROOT}/lib/libmkl_core.a -lpthread -lm -ldl
else
CXX = mpiicpc
LDFLAGS =  -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
endif

INC_PATH = ../Includes
CPPFLAGS = -O3 -g -I$(INC_PATH) -I${MKLROOT}/include

OBJS = param_read.o mklrand.o spin_type.o field_type.o shape.o hamiltonian.o state.o functions.o cpoints.o output.o
NEW_OBJS = functions.o cpoints.o output.o param_read.o mklrand.o array_alloc.o field_type_new.o hamiltonian_new.o state_new.o shape.o
TEST_FILES = cpoint_test.hpp mklrand_test.hpp ising_test.hpp heis_test.hpp fept_test.hpp test_functions.hpp skyrm_test.hpp

all: $(OBJS) $(NEW_OBJS) main.o main_new.o field_run.o

array_alloc.o: array_alloc.cpp array_alloc.hpp
	$(CXX) $(CPPFLAGS) -c array_alloc.cpp

field_type_new.o: field_type_new.cpp field_type_new.hpp array_alloc.hpp mklrand.h hamiltonian_new.hpp
	$(CXX) $(CPPFLAGS) -c field_type_new.cpp

param_read.o: param_read.cpp param_read.hpp
	$(CXX) $(CPPFLAGS) -c param_read.cpp

mklrand.o: mklrand.cpp mklrand.h
	$(CXX) $(CPPFLAGS) -c mklrand.cpp

spin_type.o: spin_type.cpp spin_type.hpp mklrand.h
	$(CXX) $(CPPFLAGS) -c spin_type.cpp

field_type.o: field_type.cpp field_type.hpp
	$(CXX) $(CPPFLAGS) -c field_type.cpp

shape.o: shape.cpp shape.hpp
	$(CXX) $(CPPFLAGS) -c shape.cpp

hamiltonian.o: hamiltonian.cpp hamiltonian.hpp field_type.hpp spin_type.hpp functions.h
	$(CXX) $(CPPFLAGS) -c hamiltonian.cpp

hamiltonian_new.o: hamiltonian_new.cpp hamiltonian_new.hpp field_type_new.hpp array_alloc.hpp functions.h
	$(CXX) $(CPPFLAGS) -c hamiltonian_new.cpp

state.o: state.cpp state.hpp hamiltonian.hpp field_type.hpp shape.hpp
	$(CXX) $(CPPFLAGS) -c state.cpp

state_new.o: state_new.cpp state_new.hpp functions.h mklrand.h hamiltonian_new.hpp field_type_new.hpp shape.hpp
	$(CXX) $(CPPFLAGS) -c state_new.cpp

functions.o: functions.cpp functions.h
	$(CXX) $(CPPFLAGS) -c functions.cpp

cpoints.o: cpoints.cpp cpoints.hpp functions.h
	$(CXX) $(CPPFLAGS) -c cpoints.cpp

output.o: output.cpp output.hpp functions.h
	$(CXX) $(CPPFLAGS) -c output.cpp

main.o: main.cpp state.hpp spin_type.hpp mklrand.h functions.h cpoints.hpp output.hpp
	$(CXX) $(CPPFLAGS) -c main.cpp

main_new.o: main_new.cpp state_new.hpp mklrand.h functions.h cpoints.hpp output.hpp
	$(CXX) $(CPPFLAGS) -c main_new.cpp

field_run.o: field_run.cpp state.hpp spin_type.hpp mklrand.h functions.h cpoints.hpp output.hpp
	$(CXX) $(CPPFLAGS) -c field_run.cpp

run: $(OBJS) main.o
	$(CXX) $(OBJS) main.o -o run $(LDFLAGS)

run_new: $(NEW_OBJS) main_new.o
	$(CXX) $(NEW_OBJS) main_new.o -o run_new $(LDFLAGS)

f_run: $(OBJS) field_run.o
	$(CXX) $(OBJS) field_run.o -o f_run $(LDFLAGS)

test.o: test.cpp mklrand.h cpoints.hpp hamiltonian_new.hpp field_type_new.hpp $(TEST_FILES)
	$(CXX) $(CPPFLAGS) -c test.cpp -lgtest

test: $(NEW_OBJS) test.o
	$(CXX) $(NEW_OBJS) test.o -o test -lgtest $(LDFLAGS)

clean:
	-rm *.o
	-rm run
	-rm f_run
	-rm test
