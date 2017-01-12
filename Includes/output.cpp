#include "output.hpp"
#include "functions.h"
#include "boost/assign.hpp"
#include <iostream>
#include <cmath>
#include <sstream>

string main_name(double H, double J, int size, double k, char shape, char hamil)
{
    stringstream avstream;
    avstream << "Output/H_" << H << "-J_" << J << "-s_" << size << "-k_" << k << "-sh_" << shape << "-ha_" << hamil << ".txt";

    string avname;
    avstream >> avname;

    return avname;
}

string main_name_d(double H, double J, double amean, double asd, double k, char shape, char hamil)
{
    stringstream avstream;
    avstream << "Output/H_" << H << "-J_" << J << "-am_" << amean << "-asd_" << asd << "-k_" << k << "-sh_" << shape << "-ha_" << hamil << ".txt";

    string avname;
    avstream >> avname;

    return avname;
}

string full_name(double H, double J, int size, double k, char shape, char hamil, double T)
{
    stringstream datstream;
    string datname;

    datstream << "Output/Fulldat/H_" << H << "-J_" << J << "-s_" << size << "-k_" << k << "-sh_" << shape << "-ha_" << hamil << "/T_" << T << ".txt";
    datstream >> datname;
    return datname;
}

string full_name_d(double H, double J, double amean, double asd, double k, char shape, char hamil, double T)
{
    stringstream datstream;
    string datname;

    datstream << "Output/Fulldat/H_" << H << "-J_" << J << "-am_" << amean << "-asd_" << asd << "-k_" << k << "-sh_" << shape << "-ha_" << hamil << "/T_" << T << ".txt";
    datstream >> datname;
    return datname;
}

void create_folder(double H, double J, int size, double k, char shape, char hamil)
{
    stringstream folderstream;
    string foldername;
    char foldcreate[1024];

    folderstream << "mkdir -p Output/Fulldat/H_" << H << "-J_" << J << "-s_" << size << "-k_" << k << "-sh_" << shape << "-ha_" << hamil << endl;
    getline(folderstream, foldername);
    strcpy(foldcreate, foldername.c_str());
    system(foldcreate) == 0;
}

void create_folder_d(double H, double J, double amean, double asd, double k, char shape, char hamil)
{
    stringstream folderstream;
    string foldername;
    char foldcreate[1024];

    folderstream << "mkdir -p Output/Fulldat/H_" << H << "-J_" << J << "-am_" << amean << "-asd_" << asd << "-k_" << k << "-sh_" << shape << "-ha_" << hamil << endl;
    getline(folderstream, foldername);
    strcpy(foldcreate, foldername.c_str());
    system(foldcreate) == 0;
}

void print_full(string datname, double T, vector<double>& allener, vector<double>& allmag,
    vector<double>& allmagx, vector<double>& allmagy, vector<double>& allmagz,
    vector<double>& allsmag, vector<double>& allsmagx, vector<double>& allsmagy,
    vector<double>& allsmagz, vector<int>& allnums)
{
    ofstream f2;
    f2.open(datname.c_str());

    cout << "T = " << T << endl;
    for (unsigned int i=0; i < allmag.size(); i++)
    {
        f2 << allmag[i] << "	" << allener[i] << " "
            << allmagx[i] << " " << allmagy[i] << " "
            << allmagz[i] << " " << allnums[i] << " "
            << allsmag[i] << " " << allsmagx[i] << " "
            << allsmagy[i] << " " << allsmagz[i] << endl;
    }

    f2.close();
    f2.clear();
}

void print_avs(string avname, vector<double>& allener, vector<double>& allmag,
    vector<double>& allmagx, vector<double>& allmagy, vector<double>& allmagz,
    vector<double>& allsmag, vector<double>& allsmagx, vector<double>& allsmagy,
    vector<double>& allsmagz, vector<int>& allnums, double H, char hamil,
    double g_lattsize, double g_slattsize, double T)
{
    fstream f;
    f.open(avname.c_str(), fstream::out | fstream::app);

    double E = sum(allener)/g_lattsize;
    double Bind = 0;
    double HC = heat_cap(allener, T, allnums);

    double Mx = sum(allmagx)/g_lattsize;
    double My = sum(allmagy)/g_lattsize;
    double Mz = sum(allmagz)/g_lattsize;
    double M, MS;
    vector<double> Mvec = boost::assign::list_of(Mx)(My)(Mz);
    double sMx = sum(allsmagx)/g_slattsize;
    double sMy = sum(allsmagy)/g_slattsize;
    double sMz = sum(allsmagz)/g_slattsize;
    double sM = sum(allsmag)/g_slattsize;

    if (H==0 || hamil=='i' || hamil=='I')
    {
        MS = mag_sus(allmag, T, allnums);
        M = sum(allmag)/g_lattsize;
        f << T << " " << HC << " " << M << " " << MS << " " << E << " "
            << Bind << " " << sM << endl;
    }
    else
    {
        MS = mag_sus(allmagx, allmagy, allmagz, T, allnums);
        M = norm(Mvec);
        f << T << " " << HC << " " << M << " " << MS << " " << E << " "
            << Bind << " " << Mx << " " << My << " " << Mz << " "
            << sM << " " << sMx << " " << sMy << " " << sMz << endl;
    }
    f.close();
}

void init_avs(string avname, double H, double J, double k, int size, double Tmin,
    double Tmax)
{
    ofstream f;
    f.open(avname.c_str());

    f << H << "	" << J << "	" << k << "	" << size << "	" << Tmin << "	" << Tmax << endl;

    f.close();
}

void load_temps(string prefix, double Ts[])
{
    stringstream loadstream;
    string loadname;
    loadstream << "Includes/Temps/" << prefix << ".txt" << endl;
    loadstream >> loadname;

    ifstream f;
    f.open(loadname.c_str());
    double curr;
    bool cont = false;
    if(f >> curr) {cont = true;}
    for (int i=0; cont; i++)
    {
        Ts[i] = curr;
        if(f >> curr) {}
        else {cont = false;}
    }
    f.close();
}
