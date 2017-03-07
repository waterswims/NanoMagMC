#include "param_read.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>

template <class T> T read_var(string v_name, string f_name)
{
    fstream file;
    file.open(f_name.c_str(), fstream::in);
    if(!file.is_open())
    {
        cout << "Input file not opened" << endl;
        exit(105);
    }
    string line;
    string test_v_name;
    T test_out;
    while(getline(file, line))
    {
        istringstream stream(line);
        stream >> test_v_name >> test_out;
        if (test_v_name == v_name)
        {
            file.close();
            // cout << "Read Variable: " << v_name << " = " << test_out << endl;
            return test_out;
        }
    }
    file.close();
    cout << "Variable name not found" << endl;
    exit(105);
}

void read_all_vars(string f_name, int& size, double& J, double& H, double& k,
    bool& periodic, char& shape, char& hamil, int& N_av, int& N_single,
    int& pad, double& beta, bool& distrib, double& amean, double& asd,
    string& temp_name, double& K)
{
    shape = read_var<char>("LATTSHAPE", f_name);
    switch(shape)
    {
        case 'w':
        case 'W':
        case 'x':
        case 'X':
            beta = read_var<double>("WEIBULLFACT", f_name);
            break;
        default:
            beta = 50;
    }
    hamil = read_var<char>("HAMILTONIAN", f_name);
    switch(hamil)
    {
        case 'i':
        case 'I':
        case 'h':
        case 'H':
            J = read_var<double>("EXCHANGE", f_name);
            H = read_var<double>("MAGFIELD", f_name);
            break;
        case 'f':
        case 'F':
            H = read_var<double>("MAGFIELD", f_name);
            break;
        case 's':
        case 'S':
            J = read_var<double>("EXCHANGE", f_name);
            H = read_var<double>("MAGFIELD", f_name);
            K = read_var<double>("DMISTREN", f_name);
            break;
        default:
            J = 1;
            H = 0;
            K = 0;
    }
    distrib = read_var<int>("ISDISTRIB", f_name);
    if(distrib)
    {
        amean = read_var<double>("MEANSIZE", f_name);
        asd = read_var<double>("SIZEDEV", f_name);
        size = amean;
    }
    else
    {
        size = read_var<int>("SIZE", f_name);
        amean = 0;
        asd = 0;
    }
    k = read_var<double>("BOLTZMANN", f_name);
    periodic = read_var<int>("ISPERIO", f_name);
    N_av = read_var<int>("LATTSPERPROC", f_name);
    N_single = read_var<int>("MCSWEEPS", f_name);
    pad = read_var<int>("LATTPADDING", f_name);
    temp_name = read_var<string>("TEMPNAME", f_name);
}
