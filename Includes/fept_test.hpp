#ifndef _FEPTTEST
#define _FEPTTEST

#include "field_type_new.hpp"
#include "hamiltonian_new.hpp"

#ifndef _TESTRANDGEN
#define _TESTRANDGEN
mkl_irand st_rand_int(1e5, 1);
mkl_drand st_rand_double(1e5, 2);
mkl_lnrand rand_ln(0, 0.25, 1e5, 3);
#endif

#endif
