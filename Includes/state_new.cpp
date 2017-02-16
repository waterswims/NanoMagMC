#include "state_new.hpp"
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

    this->init_points(size, isPerio, H, J, args);

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
    E = other.E;
    M = other.M;
    beta = other.beta;
    k_b = other.k_b;
    num = other.num;
    h_code = other.h_code;
    s_code = other.s_code;
    this->copy_points(other);
}

state& state::operator=(const state& other)
{
    E = other.E;
    M = other.M;
    beta = other.beta;
    k_b = other.k_b;
    num = other.num;
    h_code = other.h_code;
    s_code = other.s_code;
    this->copy_points(other);
}

void state::init_points(double size, bool isPerio, double H, double J, double* args)
{
    switch(h_code)
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
    switch (h_code)
    {
        case 'h':
        case 'H':
        case 'f':
        case 'F':
        switch (s_code)
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
                field = new field_3d_h(int(size), isPerio);
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
        switch (s_code)
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
                field = new field_3d_i(int(size), isPerio);
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
}

void state::copy_points(const state& other)
{
    switch(h_code)
    {
        case 'i':
        case 'I':
            hamil = new ham_ising(*(other.hamil));
            break;
        case 'h':
        case 'H':
            hamil = new ham_heis(*(other.hamil));
            break;
        case 'f':
        case 'F':
            hamil = new ham_FePt(*(other.hamil));
            break;
        default:
            cerr << "Incorrect hamiltonian, exiting..." << endl;
            exit(102);
    }
    switch (h_code)
    {
        case 'h':
        case 'H':
        case 'f':
        case 'F':
        switch (s_code)
        {
            case 's':
            case 'S':
                field = new field_2d_h(*(other.field));
                shape = new square;
                break;
            case 'w':
            case 'W':
                field = new field_2d_h(*(other.field));
                shape = new weibull(*(other.shape));
                break;
            case 'c':
            case 'C':
                field = new field_3d_h(*(other.field));
                shape = new cube;
                break;
            case 'x':
            case 'X':
                field = new field_3d_h(*(other.field));
                shape = new weibull(*(other.shape));
                break;
            default:
                cerr << "Incorrect shape code, exiting" << endl;
                exit(103);
        }
        break;

        case 'i':
        case 'I':
        switch (s_code)
        {
            case 's':
            case 'S':
                field = new field_2d_i(*(other.field));
                shape = new square;
                break;
            case 'w':
            case 'W':
                field = new field_2d_i(*(other.field));
                shape = new weibull(*(other.shape));
                break;
            case 'c':
            case 'C':
                field = new field_3d_i(*(other.field));
                shape = new cube;
                break;
            case 'x':
            case 'X':
                field = new field_3d_i(*(other.field));
                shape = new weibull(*(other.shape));
                break;
            default:
                cerr << "Incorrect shape code, exiting" << endl;
                exit(103);
        }
        break;
    }
}

state::~state()
{
    delete hamil;
    delete field;
    delete shape;
}

void state::init_lattice()
{
    num = 0;
    snum = 0;
    int start = 0;
    if(!field->get_perio()){start++;}
    int dim = field->get_dim();
    vector<int> pos(dim);
    for (vector<int>::iterator it = pos.begin(); it != pos.end(); it++)
    {
        *it = start;
    }
    bool finished = false;
    while (!finished)
    {
        bool fillspin = shape->check(pos, field->get_totsize());
        int possum = sum(pos);
        if (fillspin)
        {
            field->fill_rand(pos);
            num++;
            if (possum%2 == 0){snum++;}
        }
        else{field->fill_zero(pos);}
        field->next(finished, pos);
    }
    field->fill_ghost();
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

void state::change_temp(double T)
{

}

void state::print_latt()
{

}

void state::ptf(string fname)
{

}
