#include "../includes/print_latt.hpp"
#include <cstdlib>
#include <iostream>
#include <sstream>

field_type* set_sum_latt(double size,
                         bool periodic,
                         char shape,
                         char hamil)
{
    field_type* field;
    int pad = 1;
    switch (hamil)
    {
        case 'f':
        case 'F':
            double sizeab = size / 0.26233661582;
            double sizec = size / 0.385;
            switch (shape)
            {
                case 'x':
                case 'X':
                    field = new field_3d_h(int(2*sizeab+10), periodic, pad);
                    break;
                default:
                    std::cerr << "Incorrect shape code, FePt only works with Weibull, exiting" << std::endl;
                    exit(103);
            }
            break;
        case 's':
        case 'S':
            pad = 2;
        case 'h':
        case 'H':
        switch (shape)
        {
            case 's':
            case 'S':
                field = new field_2d_h(int(size), periodic, pad);
                break;
            case 'w':
            case 'W':
                field = new field_2d_h(int(2*size+10), periodic, pad);
                break;
            case 'c':
            case 'C':
                field = new field_3d_h(int(size), periodic, pad);
                break;
            case 'x':
            case 'X':
                field = new field_3d_h(int(2*size+10), periodic, pad);
                break;
            default:
                std::cerr << "Incorrect shape code, exiting" << std::endl;
                exit(103);
        }
        break;

        case 'i':
        case 'I':
        switch (shape)
        {
            case 's':
            case 'S':
                field = new field_2d_i(int(size), periodic);
                break;
            case 'w':
            case 'W':
                field = new field_2d_i(int(2*size+10), periodic);
                break;
            case 'c':
            case 'C':
                field = new field_3d_i(int(size), periodic);
                break;
            case 'x':
            case 'X':
                field = new field_3d_i(int(2*size+10), periodic);
                break;
            default:
                std::cerr << "Incorrect shape code, exiting" << std::endl;
                exit(103);
        }
        break;
    }
    return field;
}

void print_field(field_type* field,
                int protocol,
                const double var1,
                const double var2,
                const std::string prefix)
{
    std::stringstream datstream, foldstream;
    std::string foldname, datname;
    foldstream << prefix << "/Latt_Print/H_";
    switch(protocol)
    {
        case 1:
        foldstream << var1 <<  "-T_" << var2;
        break;
        case 2:
        foldstream << var2 <<  "-T_" << var1;
        break;
    }
    foldstream >> foldname;

    datstream << "mkdir -p " << foldname << std::endl;
    getline(datstream, datname);
    system(datname.c_str());

    datstream << foldname << "/row" << std::endl;
    getline(datstream, datname);

    field->print(datname);
}
