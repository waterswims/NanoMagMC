#include "state.hpp"
#include "mklrand.h"
#include "functions.h"
#include <iostream>
#include <fstream>

extern mkl_irand st_rand_int;
extern mkl_drand st_rand_double;

state::state(double size, bool isPerio, char shape_code, char ham_code, double J,
    double H, double k, double Temp, double* args)
{
    h_code = ham_code;
    s_code = shape_code;
    switch(ham_code)
    {
        case 'i':
        case 'I':
            hamil = new ham_ising(H, J);
            break;
        case 'h':
        case 'H':
            hamil = new ham_heis(H, J);
            break;
        case 'f':
        case 'F':
            hamil = new ham_FePt();
            break;
        case 'c':
        case 'C':
            hamil = new ham_cluster("Includes/Js/cluster.txt");
            break;
        default:
            cerr << "Incorrect hamiltonian, exiting..." << endl;
            exit(102);
    }
}

state::state(const state& other)
{

}

state::~state()
{

}

void state::equil(int iter)
{

}

vector<double> state::magnetisation()
{

}

vector<double> state::submag(int subnumber)
{

}

double state::energy()
{

}

int state::num_spins()
{

}

int state::sub_num(int subnumber)
{

}

void state::init_lattice()
{

}

void state::change_temp(double T)
{

}

state& state::operator=(const state& other)
{

}

void state::print_latt()
{

}

void state::ptf(string fname)
{

}
