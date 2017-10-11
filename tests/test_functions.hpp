#ifndef _TESTFUNCS
#define _TESTFUNCS

#include "../includes/field_type.hpp"
#include "../includes/hamiltonian.hpp"
#include <vector>

using namespace std;

mkl_irand st_rand_int(1e5, 1);
mkl_drand st_rand_double(1e5, 2);
mkl_lnrand rand_ln(0, 0.25, 1e5, 3);

field_2d_i gen_2d_ising_fm()
{
    field_2d_i field(10, false);
    vector<int> pos(2);
    for(int i = 1; i < 11; i++)
    {
        pos[0] = i;
        for(int j = 1; j < 11; j++)
        {
            pos[1] = j;
            field.fill_val_i(pos, 1);
        }
    }
    field.fill_ghost(1);
    return field;
}

field_2d_i gen_2d_ising_afm()
{
    field_2d_i field(10, false);
    vector<int> pos(2);
    for(int i = 1; i < 11; i++)
    {
        pos[0] = i;
        for(int j = 1; j < 11; j++)
        {
            pos[1] = j;
            field.fill_val_i(pos, ((i+j)%2)*2-1);
        }
    }
    field.fill_ghost(1);
    return field;
}

field_3d_i gen_3d_ising_fm()
{
    field_3d_i field(10, false);
    vector<int> pos(3);
    for(int i = 1; i < 11; i++)
    {
        pos[0] = i;
        for(int j = 1; j < 11; j++)
        {
            pos[1] = j;
            for(int k = 1; k < 11; k++)
            {
                pos[2] = k;
                field.fill_val_i(pos, 1);
            }
        }
    }
    field.fill_ghost(1);
    return field;
}

field_3d_i gen_3d_ising_afm()
{
    field_3d_i field(10, false);
    vector<int> pos(3);
    for(int i = 1; i < 11; i++)
    {
        pos[0] = i;
        for(int j = 1; j < 11; j++)
        {
            pos[1] = j;
            for(int k = 1; k < 11; k++)
            {
                pos[2] = k;
                field.fill_val_i(pos, ((i+j+k)%2)*2-1);
            }
        }
    }
    field.fill_ghost(1);
    return field;
}

field_2d_h gen_2d_heis_fm(double x, double y, double z)
{
    field_2d_h field(10, false);
    vector<int> pos(2);
    for(int i = 1; i < 11; i++)
    {
        pos[0] = i;
        for(int j = 1; j < 11; j++)
        {
            pos[1] = j;
            field.fill_val_h(pos, x, y, z);
        }
    }
    field.fill_ghost(1);
    return field;
}

field_2d_h gen_2d_heis_afm(double x, double y, double z)
{
    field_2d_h field(10, false);
    vector<int> pos(2);
    for(int i = 1; i < 11; i++)
    {
        pos[0] = i;
        for(int j = 1; j < 11; j++)
        {
            pos[1] = j;
            field.fill_val_h(pos, (((i+j)%2)*2-1)*x, (((i+j)%2)*2-1)*y,
                             (((i+j)%2)*2-1)*z);
        }
    }
    field.fill_ghost(1);
    return field;
}

field_3d_h gen_3d_heis_fm(double x, double y, double z)
{
    field_3d_h field(10, false);
    vector<int> pos(3);
    for(int i = 1; i < 11; i++)
    {
        pos[0] = i;
        for(int j = 1; j < 11; j++)
        {
            pos[1] = j;
            for(int k = 1; k < 11; k++)
            {
                pos[2] = k;
                field.fill_val_h(pos, x, y, z);
            }
        }
    }
    field.fill_ghost(1);
    return field;
}

field_3d_h gen_3d_heis_afm(double x, double y, double z)
{
    field_3d_h field(10, false);
    vector<int> pos(3);
    for(int i = 1; i < 11; i++)
    {
        pos[0] = i;
        for(int j = 1; j < 11; j++)
        {
            pos[1] = j;
            for(int k = 1; k < 11; k++)
            {
                pos[2] = k;
                field.fill_val_h(pos, (((i+j+k)%2)*2-1)*x, (((i+j+k)%2)*2-1)*y,
                                 (((i+j+k)%2)*2-1)*z);
            }
        }
    }
    field.fill_ghost(1);
    return field;
}

#endif
