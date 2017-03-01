#ifndef _FEPTTEST
#define _FEPTTEST

#include "test_functions.hpp"

///////////////////////////////////////////////////////
// FePt model tests - 3D
///////////////////////////////////////////////////////

TEST(FePt, 3d_energy_zero_field)
{
    field_3d_h field = gen_3d_heis_fm(1, 0, 0);
    ham_FePt hamil;
    hamil.init_dim(&field);
    EXPECT_EQ(-2700, hamil.calc_E(&field));

    field = gen_3d_heis_fm(0, 1, 0);
    EXPECT_EQ(-2700, hamil.calc_E(&field));

    field = gen_3d_heis_fm(0, 0, 1);
    EXPECT_EQ(-2700, hamil.calc_E(&field));
}

#endif
