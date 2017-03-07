#ifndef _STATE
#define _STATE

#include "hamiltonian.hpp"
#include "field_type.hpp"
#include "shape.hpp"

#include <string.h>
#include <vector>

using namespace std;

class state
{
private:
    double E, beta, k_b;
    vector<double> M;
    vector<double> subM;
    int num, snum;
    char s_code, h_code;
    ham_type* hamil;
    field_type* field;
    shape_type* shape;

public:
    state(){}
    state(double size, bool isPerio, char shape_code, char ham_code, double J,
        double H, double k, double Temp, double K, double* args);
    state(const state& other);
    ~state();
    void copy_points(const state& other);
    void init_points(double size, bool isPerio, double H, double J, double K, double* args);
    void equil(int iter);
    vector<double> magnetisation();
    vector<double> submag(int subnumber);
    double energy();
    int num_spins();
    int sub_num(int subnumber);
    void init_lattice();
    void change_temp(double T);
    state& operator=(const state& other);
    void print_latt();
    void ptf(string fname);
};

#endif
