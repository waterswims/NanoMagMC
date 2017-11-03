#include "../includes/param_read.hpp"
#include "../includes/array_alloc.hpp"
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
            return test_out;
        }
    }
    file.close();
    cout << "Variable name '" << v_name << "' not found" << endl;
    exit(105);
}

void read_all_vars(string f_name, double& size, double& J, double& k,
    bool& periodic, char& shape, char& hamil, int& Samp_steps, int& N_samp,
    int& Eq_steps, int& N_latts, double& beta, bool& distrib, double& amean,
    double& asd, string& temp_name, string& field_name, int& protocol,
    double& K, bool& print_latt)
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
            break;
        case 's':
        case 'S':
            J = read_var<double>("EXCHANGE", f_name);
            K = read_var<double>("DMISTREN", f_name);
            break;
        default:
            J = 1;
            K = 0;
    }
    distrib = read_var<int>("ISDISTRIB", f_name);
    if(distrib)
    {
        amean = read_var<double>("MEANSIZE", f_name);
        asd = read_var<double>("SIZEDEV", f_name);
        size = amean;
        N_latts = read_var<int>("NLATTS", f_name);
    }
    else
    {
        size = read_var<double>("SIZE", f_name);
        N_latts = 1;
        amean = 0;
        asd = 0;
    }
    k = read_var<double>("BOLTZMANN", f_name);
    periodic = read_var<int>("ISPERIO", f_name);
    Samp_steps = read_var<int>("SAMPSWEEPS", f_name);
    Eq_steps = read_var<int>("EQSWEEPS", f_name);
    N_samp = read_var<int>("NSAMPS", f_name);
    temp_name = read_var<string>("TEMPNAME", f_name);
    field_name = read_var<string>("FIELDNAME", f_name);
    protocol = read_var<int>("PROTOCOL", f_name);
    print_latt = read_var<int>("PRINT_LATT", f_name);
}

void load_Hs_Ts(string Tname,
                double* &Ts,
                int& Tnum,
                string Hname,
                double* &Hs,
                int& Hnum)
{
    // load temps
    stringstream loadstream;
    string loadname;
    loadstream << "Temps/" << Tname << endl;
    loadstream >> loadname;

    ifstream f;
    f.open(loadname.c_str());
    double curr;
    bool cont = false;
    if(f >> curr) {cont = true;}
    int i = 0;
    for (; cont; i++)
    {
        if(f >> curr) {}
        else {cont = false;}
    }
    Tnum = i;
    f.close();

    Ts = alloc_1darr<double>(Tnum);

    f.open(loadname.c_str());
    cont = false;
    if(f >> curr) {cont = true;}
    for (int i = 0; cont; i++)
    {
        Ts[i] = curr;
        if(f >> curr) {}
        else {cont = false;}
    }
    f.close();

    // load fields
    loadstream << "Fields/" << Hname << endl;
    loadstream >> loadname;

    f.open(loadname.c_str());
    cont = false;
    if(f >> curr) {cont = true;}
    i = 0;
    for (; cont; i++)
    {
        if(f >> curr) {}
        else {cont = false;}
    }
    Hnum = i;
    f.close();

    Hs = alloc_1darr<double>(Hnum);

    f.open(loadname.c_str());
    cont = false;
    if(f >> curr) {cont = true;}
    for (int i = 0; cont; i++)
    {
        Hs[i] = curr;
        if(f >> curr) {}
        else {cont = false;}
    }
    f.close();
}
