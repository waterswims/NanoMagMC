#ifndef _STATE
#define _STATE

#include "hamiltonian.hpp"
#include "field_type.hpp"
#include "shape.hpp"
#include <string.h>
#include <vector>

using namespace std;

template <class T> class state
{
private:
    double E, beta, k_b;
    vector<double> M;
    vector<double> subM;
    int num, snum;
    char s_code, h_code;
    ham_type<T>* hamil;
    field_type<T>* field;
    shape_type* shape;
public:
    state(){hamil=NULL; field=NULL; shape=NULL;}
    state(double size, bool isPerio, char shape_code, char ham_code, double J,
        double H, double k, double Temp, double* args);
    state(const state<T>& other);
    ~state();
    void equil(int iter);
    vector<double> magnetisation();
    vector<double> submag(int subnumber);
    double energy();
    int num_spins();
    int sub_num(int subnumber);
    void init_lattice();
    void change_temp(double T);
    state<T>& operator=(const state<T>& other);
    void print_latt();
    void ptf(string fname);
};

#endif
