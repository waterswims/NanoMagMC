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
        default:
            cerr << "Incorrect hamiltonian, exiting..." << endl;
            exit(102);
    }
    switch (ham_code)
    {
        case 'h':
        case 'H':
        case 'f':
        case 'F':
        switch (shape_code)
        {
            case 's':
            case 'S':
                field = new field_2d_h(int(size), isPerio);
                shape = new square;
                break;
            case 'w':
            case 'W':
                field = new field_2d_h(int(2*size+10), isPerio);
                shape = new weibull((size), args[0]);
                break;
            case 'c':
            case 'C':
                field = new field_3d_h(int(size), isPerio, args[1]);
                shape = new cube;
                break;
            case 'x':
            case 'X':
                field = new field_3d_h(int(2*size+10), isPerio);
                shape = new weibull((size), args[0]);
                break;
            default:
                cerr << "Incorrect shape code, exiting" << endl;
                exit(103);
        }
        break;

        case 'i':
        case 'I':
        switch (shape_code)
        {
            case 's':
            case 'S':
                field = new field_2d_i(int(size), isPerio);
                shape = new square;
                break;
            case 'w':
            case 'W':
                field = new field_2d_i(int(2*size+10), isPerio);
                shape = new weibull((size), args[0]);
                break;
            case 'c':
            case 'C':
                field = new field_3d_i(int(size), isPerio, args[1]);
                shape = new cube;
                break;
            case 'x':
            case 'X':
                field = new field_3d_i(int(2*size+10), isPerio);
                shape = new weibull((size), args[0]);
                break;
            default:
                cerr << "Incorrect shape code, exiting" << endl;
                exit(103);
        }
        break;
    }

    k_b = k;
    if (Temp <= 0)
    {
        cerr << "Invalid temperature, exiting" << endl;
        exit(104);
    }
    beta = 1 / (k_b * Temp);
    this->init_lattice();
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
