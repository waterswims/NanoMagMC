#include "../includes/protocol.hpp"

void set_protocol(
    const int proto_code,
    float*& var1_list,
    float*& var2_list,
    int& var1_size,
    int& var2_size,
    int& var1_begin,
    int& var2_begin,
    int& var1_end,
    int& var2_end,
    int& var1_final,
    float* Hs,
    float* Ts,
    const int H_size,
    const int T_size)
{
    switch(proto_code)
    {
        case 1:
        var1_list = Hs;
        var2_list = Ts;
        var1_size = H_size;
        var2_size = T_size;
        var1_begin = var2_begin = 0;
        var1_end = H_size;
        var2_end = T_size;
        var1_final = H_size - 1;
        break;
        case 2:
        var1_list = Ts;
        var2_list = Hs;
        var1_size = T_size;
        var2_size = H_size;
        var1_begin = var2_begin = 0;
        var1_end = T_size;
        var2_end = H_size;
        var1_final = T_size - 1;
        break;
        case 4:
        var1_list = Ts;
        var2_list = Hs;
        var1_size = T_size;
        var2_size = H_size;
        var1_begin = 0;
        var2_end = -1;
        var1_end = T_size;
        var2_begin = H_size - 1;
        var1_final = T_size - 1;
        break;
        default:
        std::cout << "Invalid protocol, exiting..." << std::endl;
        exit(207);
    }
}

void incr_v1(
    const int proto_code,
    int& var1_curr)
{
    switch(proto_code)
    {
        case 1:
        case 2:
        case 4:
        var1_curr++;
        break;
        default:
        std::cout << "Invalid protocol, exiting..." << std::endl;
        exit(207);
    }
}

void incr_v2(
    const int proto_code,
    int& var2_curr)
{
    switch(proto_code)
    {
        case 1:
        case 2:
        var2_curr++;
        break;
        case 4:
        var2_curr--;
        break;
        default:
        std::cout << "Invalid protocol, exiting..." << std::endl;
        exit(207);
    }
}

bool check_rank_run(
    const int proto_code,
    const int i,
    const int comm_size,
    const int rank,
    const int var1_size)
{
    switch(proto_code)
    {
        case 1:
        case 2:
        case 4:
        return rank != i%comm_size;
        default:
        std::cout << "Invalid protocol, exiting..." << std::endl;
        exit(207);
    }
}
