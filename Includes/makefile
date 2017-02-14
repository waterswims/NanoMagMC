# Makefile for project

OS := $(shell uname)

ifeq ($(OS),Darwin)
CXX = mpic++
LDFLAGS =  ${MKLROOT}/lib/libmkl_intel_lp64.a ${MKLROOT}/lib/libmkl_sequential.a ${MKLROOT}/lib/libmkl_core.a -lpthread -lm -ldl
else
CXX = mpiicpc
LDFLAGS = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
endif

INC_PATH = ../Includes
CPPFLAGS = -O3 -I$(INC_PATH) -I${MKLROOT}/include

OBJS = param_read.o mklrand.o spin_type.o field_type.o shape.o hamiltonian.o state.o functions.o cpoints.o output.o array_alloc.o

all: $(OBJS) main.o field_run.o

array_alloc.o: array_alloc.cpp array_alloc.hpp
	$(CXX) $(SPPFLAGS) - array_alloc.cpp

param_read.o: param_read.cpp param_read.hpp
	$(CXX) $(CPPFLAGS) -c param_read.cpp

mklrand.o: mklrand.cpp mklrand.h
	$(CXX) $(CPPFLAGS) -c mklrand.cpp

spin_type.o: spin_type.cpp spin_type.hpp mklrand.h
	$(CXX) $(CPPFLAGS) -c spin_type.cpp

field_type.o: field_type.cpp field_type.hpp
	$(CXX) $(CPPFLAGS) -c field_type.cpp

# FCC.o: FCC.cpp FCC.hpp field_type.hpp
# 	$(CXX) $(CPPFLAGS) -c FCC.cpp

shape.o: shape.cpp shape.hpp
	$(CXX) $(CPPFLAGS) -c shape.cpp

hamiltonian.o: hamiltonian.cpp hamiltonian.hpp field_type.hpp spin_type.hpp functions.h
	$(CXX) $(CPPFLAGS) -c hamiltonian.cpp

state.o: state.cpp state.hpp hamiltonian.hpp field_type.hpp shape.hpp
	$(CXX) $(CPPFLAGS) -c state.cpp

functions.o: functions.cpp functions.h
	$(CXX) $(CPPFLAGS) -c functions.cpp

cpoints.o: cpoints.cpp cpoints.hpp functions.h
	$(CXX) $(CPPFLAGS) -c cpoints.cpp

output.o: output.cpp output.hpp functions.h
	$(CXX) $(CPPFLAGS) -c output.cpp

main.o: main.cpp state.hpp spin_type.hpp mklrand.h functions.h cpoints.hpp output.hpp
	$(CXX) $(CPPFLAGS) -c main.cpp

field_run.o: field_run.cpp state.hpp spin_type.hpp mklrand.h functions.h cpoints.hpp output.hpp
	$(CXX) $(CPPFLAGS) -c field_run.cpp

run: $(OBJS) main.o
	$(CXX) $(LDFLAGS) $(OBJS) main.o -o run

f_run: $(OBJS) field_run.o
	$(CXX) $(LDFLAGS) $(OBJS) field_run.o -o f_run

test.o: test.cpp mklrand.h cpoints.hpp
	$(CXX) $(CPPFLAGS) -c test.cpp -lgtest

test: $(OBJS) test.o
	$(CXX) $(LDFLAGS) $(OBJS) test.o -o test -lgtest

clean:
	-rm *.o
	-rm run
	-rm f_run
	-rm test
