#include "../includes/output.hpp"
#include <iostream>
#include <sstream>
#include <cstdlib>

void create_folders(const double J,
                    const int size,
                    const double k,
                    const double K,
                    const char shape,
                    const char hamil)
{
    std::stringstream datstream;
    std::string datname;
    datstream << "mkdir -p Output/" << "J_" << J << "-s_" << size << "-k_" << k;
    if (hamil == 's' && hamil == 'S')
    {
        datstream << "-K_" << K;
    }
    datstream << "-sh_" << shape << "-ha_" << hamil << "/Fulldat" << std::endl;
    getline(datstream, datname);
    std::cout << datname << std::endl;
    system(datname.c_str());

    datstream << "mkdir -p Output/" << "J_" << J << "-s_" << size << "-k_" << k;
    if (hamil == 's' && hamil == 'S')
    {
        datstream << "-K_" << K;
    }
    datstream << "-sh_" << shape << "-ha_" << hamil << "/Latt_Print";
    datstream << std::endl;
    getline(datstream, datname);
    system(datname.c_str());
}

void print_sngl_HT(const double* magx,
                   const double* magy,
                   const double* magz,
                   const double* mag,
                   const double* ener,
                   const double* smagx,
                   const double* smagy,
                   const double* smagz,
                   const double* smag,
                   const int N_samp,
                   const int protocol,
                   const double var1,
                   const double var2,
                   const double J,
                   const int size,
                   const double k,
                   const double K,
                   const char shape,
                   const char hamil)
{
    std::stringstream datstream;
    std::string datname;
    datstream << "Output/" << "J_" << J << "-s_" << size << "-k_" << k;
    if (hamil == 's' && hamil == 'S')
    {
        datstream << "-K_" << K;
    }
    datstream << "-sh_" << shape << "-ha_" << hamil << "/Fulldat/H_";
    std::cout << "\33[2K\rPrinting H = ";
    switch(protocol)
    {
        case 1:
        std::cout << var1 << ", T = " << var2 << "..." << std::flush;
        datstream << var1 <<  "-T_" << var2 << ".dat";
        break;
        case 2:
        std::cout << var2 << ", T = " << var1 << "..." << std::flush;
        datstream << var2 <<  "-T_" << var1 << ".dat";
        break;
    }
    datstream >> datname;

    std::ofstream f;
    f.open(datname.c_str());

    for (int i = 0; i < N_samp; i++)
    {
        f << ener[i] << " " << mag[i] << " " << smag[i];
        if (hamil != 'i' && hamil != 'I')
        {
            f << " " << magx[i] << " " << magy[i] << " " << magz[i] << " ";
            f << smagx[i] << " " << smagy[i] << " " << smagz[i];
        }
        f << std::endl;
    }

    f.close();
}
