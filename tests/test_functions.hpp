#ifndef _TESTFUNCS
#define _TESTFUNCS

#ifndef _PI
#define _PI 3.141592653589793
#endif

#include "../includes/field_type.hpp"
#include "../includes/hamiltonian.hpp"

#ifdef __INTEL_COMPILER
#include "../includes/mklrand.hpp"
#define IRANDTYPE mklrand::mkl_irand
#define DRANDTYPE mklrand::mkl_drand
#define LNRANDTYPE mklrand::mkl_lnrand
#else
#include "../includes/stdrand.hpp"
#define IRANDTYPE stdrand::std_i_unirand
#define DRANDTYPE stdrand::std_d_unirand
#define LNRANDTYPE stdrand::std_lognormrand
#endif

#include <vector>
#include <iostream>
#include <cmath>

IRANDTYPE st_rand_int(1e5, 1);
DRANDTYPE st_rand_double(1e5, 2);
LNRANDTYPE rand_ln(0, 0.25, 1e5, 3);

field_2d_i gen_2d_ising_fm()
{
    field_2d_i field(10, false);
    std::vector<int> pos(2);
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
    std::vector<int> pos(2);
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
    std::vector<int> pos(3);
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
    std::vector<int> pos(3);
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
    std::vector<int> pos(2);
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
    std::vector<int> pos(2);
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

field_2d_h gen_2d_skyrm()
{
    int tsize = 50;
    float mid = 24.5, sk_R = 20;
    bool isPerio = true;

    mid += (!isPerio);

    float k = _PI / sk_R;

    field_2d_h field(tsize, isPerio);
    std::vector<int> pos(2);
    for(int i = (!isPerio); i < (tsize+(!isPerio)); i++)
    {
        pos[0] = i;
        float dx = i-mid;
        for(int j = (!isPerio); j < (tsize+(!isPerio)); j++)
        {
            pos[1] = j;
            float dy = j - mid;
            float r2 = pow(dx, 2) + pow(dy, 2);
            float r = pow(r2, 0.5);
            if(r >= sk_R)
            {
                field.fill_val_h(pos, 0, 0, -1);
            }
            else
            {
                float phi = atan2(i-mid, j-mid);
                float x = -sin(k*r)*cos(phi);
                float y = sin(k*r)*sin(phi);
                float z = cos(k*r);
                field.fill_val_h(pos, x, y, z);
            }
        }
    }
    field.fill_ghost((!isPerio));

    return field;
}

field_3d_h gen_3d_heis_fm(double x, double y, double z)
{
    field_3d_h field(10, false);
    std::vector<int> pos(3);
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
    std::vector<int> pos(3);
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

field_3d_h gen_3d_skyrm()
{
    int tsize = 50;
    float mid = 24.5, sk_R = 20;
    bool isPerio = true;

    mid += (!isPerio);

    float k = _PI / sk_R;

    field_3d_h field(tsize, isPerio);
    std::vector<int> pos(3);
    for(int k2 = (!isPerio); k2 < (tsize+(!isPerio)); k2++)
    {
        pos[2] = k2;
        for(int i = (!isPerio); i < (tsize+(!isPerio)); i++)
        {
            pos[0] = i;
            float dx = i-mid;
            for(int j = (!isPerio); j < (tsize+(!isPerio)); j++)
            {
                pos[1] = j;
                float dy = j - mid;
                float r2 = pow(dx, 2) + pow(dy, 2);
                float r = pow(r2, 0.5);
                if(r >= sk_R)
                {
                    field.fill_val_h(pos, 0, 0, -1);
                }
                else
                {
                    float phi = atan2(i-mid, j-mid);
                    float x = -sin(k*r)*cos(phi);
                    float y = sin(k*r)*sin(phi);
                    float z = cos(k*r);
                    field.fill_val_h(pos, x, y, z);
                }
            }
        }
    }
    field.fill_ghost((!isPerio));

    return field;
}

#endif
