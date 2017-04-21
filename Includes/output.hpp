#ifndef _OUTP
#define _OUTP

#include <vector>
#include <fstream>
#include <cstring>

using namespace std;

string main_name(double H, double J, double size, double k, char shape, char hamil);

string main_name_d(double H, double J, double amean, double asd, double k, char shape, char hamil);

string full_name(double H, double J, double size, double k, char shape, char hamil, double T);

string full_name_d(double H, double J, double amean, double asd, double k, char shape, char hamil, double T);

void create_folder(double H, double J, double size, double k, char shape, char hamil);

void create_folder_d(double H, double J, double amean, double asd, double k, char shape, char hamil);

void print_full(string datname, double T, vector<double>& allener, vector<double>& allmag,
    vector<double>& allmagx, vector<double>& allmagy, vector<double>& allmagz,
    vector<double>& allsmag, vector<double>& allsmagx, vector<double>& allsmagy,
    vector<double>& allsmagz, vector<int>& allnums);

void print_avs(string avname, vector<double>& allener, vector<double>& allmag,
    vector<double>& allmagx, vector<double>& allmagy, vector<double>& allmagz,
    vector<double>& allsmag, vector<double>& allsmagx, vector<double>& allsmagy,
    vector<double>& allsmagz, vector<int>& allnums, double H, char hamil,
    double g_lattsize, double g_slattsize, double T);

void init_avs(string avname, double H, double J, double k, double size, double Tmin,
    double Tmax);

#endif
