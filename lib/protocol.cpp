#include "../includes/protocol.hpp"

void set_protocol(
    int proto_code,
    double*& var1_list,
    double*& var2_list,
    int& var1_size,
    int& var2_size,
    double* Hs,
    double* Ts,
    int H_size,
    int T_size)
{
    switch(proto_code)
    {
        case 1:
        var1_list = Hs;
        var2_list = Ts;
        var1_size = H_size;
        var2_size = T_size;
        break;
        case 2:
        var1_list = Ts;
        var2_list = Hs;
        var1_size = T_size;
        var2_size = H_size;
        break;
        default:
        std::cout << "Invalid protocol, exiting..." << std::endl;
    }
}
