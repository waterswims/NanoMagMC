#ifndef _OUTP
#define _OUTP

#include <fstream>
#include <cstring>

std::string check_h5_file(const double J,
                           const double size,
                           const double am,
                           const double asd,
                           const double k,
                           const double K,
                           const char shape,
                           const char hamil,
                           const int prot,
                           const bool distrib,
                           const int N_temps,
                           const int N_fields,
                           const int N_samps,
                           const int N_latts,
                           const double* Ts,
                           const double* Hs,
                           const int v1_size,
                           const int tc_size,
                           bool** checkp,
                           bool &file_exists);

void create_h5_file(std::string prefix,
                  const char hamil,
                  const bool distrib,
                  const int N_temps,
                  const int N_fields,
                  const int N_samps,
                  const int N_latts,
                  const double* Ts,
                  const double* Hs,
                  const int v1_size,
                  const int tc_size);

void print_TD_h5(const double* magx,
             const double* magy,
             const double* magz,
             const double* mag,
             const double* ener,
             const double* smagx,
             const double* smagy,
             const double* smagz,
             const double* smag,
             double** tcs,
             const int N_samp,
             const int protocol,
             const int var1,
             const int var2,
             const int v2max,
             const int tc_size,
             const std::string prefix,
             const char hamil,
             const bool distrib,
             const int latt_num);

#endif
