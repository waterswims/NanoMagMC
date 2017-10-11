#ifndef _OUTP
#define _OUTP

#include <fstream>
#include <cstring>

void create_folders(const double J,
                    const int size,
                    const double k,
                    const double K,
                    const char shape,
                    const char hamil);

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
                   const char hamil);

#endif
