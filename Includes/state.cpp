#include "state.hpp"
#include "mklrand.h"
#include "functions.h"
#include "FCC.hpp"
#include <iostream>
#include <fstream>

extern mkl_irand st_rand_int;
extern mkl_drand st_rand_double;

template <class T> state<T>::state(double size, bool isPerio, char shape_code,
    char ham_code, double J, double H, double k, double Temp, double* args)
{
    h_code = ham_code;
    s_code = shape_code;
    switch (ham_code)
    {
        case 'i':
        case 'I':
            hamil = new ham_ising<T>(H, J);
            break;
        case 'h':
        case 'H':
            hamil = new ham_heis<T>(H, J);
            break;
        case 'f':
        case 'F':
            hamil = new ham_FePt<T>();
            break;
        default:
            cout << "Incorrect Hamiltonian exiting" << endl;
            exit(102);
    }
    switch (shape_code)
    {
        case 's':
        case 'S':
            field = new field_2d<T>(int(size), isPerio);
            shape = new square;
            break;
        case 'w':
        case 'W':
            field = new field_2d<T>(int(2*size+10), isPerio);
            shape = new weibull((size), args[0]);
            break;
        case 'c':
        case 'C':
            field = new field_3d<T>(int(size), isPerio, args[1]);
            shape = new cube;
            break;
        case 'x':
        case 'X':
            field = new field_3d<T>(int(2*size+10), isPerio);
            shape = new weibull((size), args[0]);
            break;
        case 'h':
        case 'H':
            field = new hex_2d<T>(int(size), isPerio);
            shape = new square;
            break;
        case 'o':
        case 'O':
            field = new hex_3d<T>(int(size), isPerio);
            shape = new cube;
            break;
        case 'f':
        case 'F':
            field = new FC_Tetr<T>(int(size), isPerio);
            shape = new cube;
            break;
        case 'l':
        case 'L':
            field = new L10<T>(int(size), isPerio);
            shape = new cube;
            break;
        default:
            cout << "Incorrect shape code, exiting" << endl;
            exit(103);
    }
    k_b = k;
    if (Temp <= 0)
    {
        cout << "Invalid temperature, exiting" << endl;
        exit(104);
    }
    beta = 1 / (k_b * Temp);
    this->init_lattice();
    E = hamil->calc_E(field);
    M = hamil->calc_M(field);
}

template <class T> state<T>::state(const state<T>& other)
{
    E = other.E;
    M = other.M;
    beta = other.beta;
    k_b = other.k_b;
    num = other.num;
    h_code = other.h_code;
    s_code = other.s_code;

    switch (h_code)
    {
        case 'i':
        case 'I':
            hamil = new ham_ising<T>(*(other.hamil));
            break;
        case 'h':
        case 'H':
            hamil = new ham_heis<T>(*(other.hamil));
            break;
        case 'f':
        case 'F':
            hamil = new ham_FePt<T>(*(other.hamil));
            break;
        default:
            cout << "Incorrect Hamiltonian exiting" << endl;
            exit(102);
    }
    switch (s_code)
    {
        case 's':
        case 'S':
            field = new field_2d<T>(*(other.field));
            shape = new square;
            break;
        case 'w':
        case 'W':
            field = new field_2d<T>(*(other.field));
            shape = new weibull(*(other.shape));
            break;
        case 'c':
        case 'C':
            field = new field_3d<T>(*(other.field));
            shape = new cube;
            break;
        case 'x':
        case 'X':
            field = new field_3d<T>(*(other.field));
            shape = new weibull(*(other.shape));
            break;
        case 'h':
        case 'H':
            field = new hex_2d<T>(*(other.field));
            shape = new square;
            break;
        case 'o':
        case 'O':
            field = new hex_3d<T>(*(other.field));
            shape = new cube;
            break;
        case 'f':
        case 'F':
            field = new FC_Tetr<T>(*(other.field));
            shape = new cube;
            break;
        case 'l':
        case 'L':
            field = new L10<T>(*(other.field));
            shape = new cube;
            break;
        default:
            cout << "Incorrect shape code, exiting" << endl;
            exit(103);
    }
    //(*hamil) = *(other.hamil);
    //(*shape) = *(other.shape);
    //(*field) = *(other.field);
}

template <class T> state<T>::~state()
{
    delete hamil;
    delete field;
    delete shape;
}

template <class T> void state<T>::init_lattice()
{
    num = 0;
    snum = 0;
    int start = 0;
    if(!field->get_perio())
    {
        start++;
    }
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
        T* curr = &(field->next(finished, pos));
        if (fillspin)
        {
            curr->rand_spin();
            num++;
            if (possum%2 == 0)
            {
                snum++;
            }
        }
        else
        {
            curr->zero_spin();
        }
    }
    field->fill_ghost(num);
}

template <class T> vector<double> state<T>::magnetisation()
{
    M = hamil->calc_M(field);
    return M;
}

template <class T> vector<double> state<T>::submag(int subnumber)
{
    subM = hamil->calc_subM(field, subnumber);
    return subM;
}

template <class T> double state<T>::energy()
{
    E = hamil->calc_E(field);
    return E;
}

template <class T> int state<T>::num_spins()
{
    return num;
}

template <class T> int state<T>::sub_num(int subnumber)
{
    if (subnumber==0)
    {
        return snum;
    }
    else
    {
        return num-snum;
    }
}

template <class T> void state<T>::equil(int iter)
{
    // get info of the problem
    int dim = field->get_dim();
    int size = field->get_insize();
    int tsize = field->get_totsize();
    int per = field->get_perio();

    int s_start = (tsize-size)/2;

    // create some variables
    vector<int> r_choice(dim);
    double dE = 0;
    double log_eta = 0;

    // counting for checks
    // int eflips=0, neflips=0, rflips=0;
    // double dEav = 0;

    // cout << "beta = " << beta << endl;

    for (int i=0; i<iter; i++)
    {
        //choose a random spin
        for (vector<int>::iterator it=r_choice.begin(), end=r_choice.end(); it!=end; it++)
        {
            *it = int(st_rand_double.gen() * size)+s_start;
        }

        //if a zero spin try again
        if((field->access(r_choice)).is_zero())
        {
            i--;
            continue;
        }

        //check dE
        dE = hamil->dE(field, r_choice);
        //check if flip
        if(dE <= 0)
        {
            //cout << "Energy flip" << endl;
            (field->access(r_choice)).change_spin(hamil->get_test());
            // eflips++;
            //E += dE;
        }
        else
        {
            // dEav += dE;
            // neflips++;
            log_eta = log(st_rand_double.gen());
			if ((-dE * beta) > log_eta)
			{
                // rflips++;
                (field->access(r_choice)).change_spin(hamil->get_test());
				//E += dE;
			}
        }
    }
    // cout << "Number of Energy flips:" << eflips << endl;
    // cout << "Number of Random flips:" << rflips << endl;
    // cout << "Average dE:" << dEav / neflips << endl;
}

template <class T> void state<T>::change_temp(double Temp)
{
    beta = 1./(k_b*Temp);
}

template <class T> state<T>& state<T>::operator=(const state<T>& other)
{
    E = other.E;
    M = other.M;
    beta = other.beta;
    k_b = other.k_b;
    num = other.num;
    h_code = other.h_code;
    s_code = other.s_code;

    switch (h_code)
    {
        case 'i':
        case 'I':
            hamil = new ham_ising<T>(*(other.hamil));
            break;
        case 'h':
        case 'H':
            hamil = new ham_heis<T>(*(other.hamil));
            break;
        case 'f':
        case 'F':
            hamil = new ham_FePt<T>(*(other.hamil));
            break;
        default:
            cout << "Incorrect Hamiltonian exiting" << endl;
            exit(102);
    }
    switch (s_code)
    {
        case 's':
        case 'S':
            field = new field_2d<T>(*(other.field));
            shape = new square;
            break;
        case 'w':
        case 'W':
            field = new field_2d<T>(*(other.field));
            shape = new weibull(*(other.shape));
            break;
        case 'c':
        case 'C':
            field = new field_3d<T>(*(other.field));
            shape = new cube;
            break;
        case 'x':
        case 'X':
            field = new field_3d<T>(*(other.field));
            shape = new weibull(*(other.shape));
            break;
        case 'h':
        case 'H':
            field = new hex_2d<T>(*(other.field));
            shape = new square;
            break;
        case 'o':
        case 'O':
            field = new hex_3d<T>(*(other.field));
            shape = new cube;
            break;
        case 'f':
        case 'F':
            field = new FC_Tetr<T>(*(other.field));
            shape = new cube;
            break;
        case 'l':
        case 'L':
            field = new L10<T>(*(other.field));
            shape = new cube;
            break;
        default:
            cout << "Incorrect shape code, exiting" << endl;
            exit(103);
    }

    //(*hamil) = *(other.hamil);
    //(*shape) = *(other.shape);
    //(*field) = *(other.field);
    return *this;
}

template <class T> void state<T>::print_latt()
{
    field->print();
}

template <class T> void state<T>::ptf(string fname)
{
    ofstream ofile;
    ofile.open(fname.c_str());
    bool finished = false;
    int sum=0;
    int start = 0;
    if(!field->get_perio())
    {
        start++;
    }
    int dim = field->get_dim();
    vector<int> pos(dim);
    for (vector<int>::iterator it = pos.begin(); it != pos.end(); it++)
    {
        *it = start;
    }
    while (!finished)
    {
        for (vector<int>::iterator it = pos.begin(); it != pos.end(); it++)
        {
            ofile << *it << " ";
        }
        ofile << (field->next(finished, pos)).i_access() << endl;
    }
    ofile.close();
}

template class state<ising_spin>;
template class state<heis_spin>;
