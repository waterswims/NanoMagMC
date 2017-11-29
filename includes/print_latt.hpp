#ifndef _PLATT
#define _PLATT

#include "field_type.hpp"
#include <cstring>

field_type* set_sum_latt(double size,
                         bool periodic,
                         char shape,
                         char hamil);

std::string av_latt_name(const int protocol,
                 const int var1,
                 const int var2);

std::string sing_latt_name(const int protocol,
              const int var1,
              const int var2);

#endif
