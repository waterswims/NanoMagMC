#include <string>

using namespace std;

template <class T> T read_var(string v_name, string f_name);

void read_all_vars(string f_name, int& size, double& J, double& H, double& k,
    bool& periodic, char& shape, char& hamil, int& N_av, int& N_single,
    int& pad, double& beta, bool& distrib, double& amean, double& asd,
    string& temp_name, double& K);

void load_temps(string prefix, double* &Ts, int& num);
