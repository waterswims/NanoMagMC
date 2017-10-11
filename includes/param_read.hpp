#include <string>

using namespace std;

template <class T> T read_var(string v_name, string f_name);

void read_all_vars(string f_name, double& size, double& J, double& k,
    bool& periodic, char& shape, char& hamil, int& Samp_steps, int& N_samp,
    int& Eq_steps, int& pad, double& beta, bool& distrib, double& amean,
    double& asd, string& temp_name, string& field_name, int& protocol,
    double& K, bool& print_latt);

void load_Hs_Ts(string Tname,
                double* &Ts,
                int& Tnum,
                string Hname,
                double* &Hs,
                int& Hnum);
