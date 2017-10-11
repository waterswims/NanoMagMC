#ifndef _PLATT
#define _PLATT

#include "field_type.hpp"
#include <cstring>

field_type* set_sum_latt(double size,
                         bool periodic,
                         char shape,
                         char hamil);

void print_field(field_type* field,
                int protocol,
                const double var1,
                const double var2,
                const std::string prefix);

#endif
