#ifndef _HEISTEST
#define _HEISTEST
#include "field_type_new.hpp"
#include "hamiltonian_new.hpp"

#ifndef _TESTRANDGEN
#define _TESTRANDGEN
mkl_irand st_rand_int(1e5, 1);
mkl_drand st_rand_double(1e5, 2);
mkl_lnrand rand_ln(0, 0.25, 1e5, 3);
#endif

///////////////////////////////////////////////////////
// Necessary functions
///////////////////////////////////////////////////////

field_2d_h gen_2d_heis_fm_up()
{
    field_2d_h field(10, false);
    vector<int> pos(2);
    for(int i = 1; i < 11; i++)
    {
        pos[0] = i;
        for(int j = 1; j < 11; j++)
        {
            pos[1] = j;
            field.fill_val_h(pos, 1, 0, 0);
        }
    }
    field.fill_ghost();
    return field;
}

///////////////////////////////////////////////////////
// Heisenberg model tests - 3D
///////////////////////////////////////////////////////

TEST(Heis_model, 3d_ferromagnetic_energy_zero_field)
{
    field_2d_h field = gen_2d_heis_fm_up();
    ham_heis hamil(0, 1);
    hamil.init_dim(&field);
    EXPECT_EQ(-180, hamil.calc_E(&field));
}

#endif
