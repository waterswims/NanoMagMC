# Makefile for project

CXX = mpiicpc
#CXX = mpic++
INC_PATH = ../Includes
CPPFLAGS = -g -O3 -mkl -I$(INC_PATH)

OBJS = param_read.o mklrand.o spin_type.o field_type.o shape.o hamiltonian.o state.o functions.o cpoints.o output.o

all: $(OBJS) main.o field_run.o

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
	$(CXX) $(CPPFLAGS) $(OBJS) main.o -o run

f_run: $(OBJS) field_run.o
	$(CXX) $(CPPFLAGS) $(OBJS) field_run.o -o f_run

clean:
	-rm *.o
	-rm run
	-rm f_run
