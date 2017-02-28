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
    field.fill_ghost();
    return field;
}

///////////////////////////////////////////////////////
// Heisenberg model tests - 2D
///////////////////////////////////////////////////////

TEST(Heis_model, 2d_ferromagnetic_energy_zero_field)
{
    field_2d_h field = gen_2d_heis_fm(1, 0, 0);
    ham_heis hamil(0, 1);
    hamil.init_dim(&field);
    EXPECT_EQ(-180, hamil.calc_E(&field));

    field = gen_2d_heis_fm(-1, 0, 0);
    EXPECT_EQ(-180, hamil.calc_E(&field));

    field = gen_2d_heis_fm(0, 1, 0);
    EXPECT_EQ(-180, hamil.calc_E(&field));

    field = gen_2d_heis_fm(0, -1, 0);
    EXPECT_EQ(-180, hamil.calc_E(&field));

    field = gen_2d_heis_fm(0, 0, 1);
    EXPECT_EQ(-180, hamil.calc_E(&field));

    field = gen_2d_heis_fm(0, 0, -1);
    EXPECT_EQ(-180, hamil.calc_E(&field));
}

TEST(Heis_model, 2d_ferromagnetic_energy_ext_field)
{
    field_2d_h field = gen_2d_heis_fm(1, 0, 0);
    ham_heis hamil(0.1, 1);
    hamil.init_dim(&field);
    EXPECT_EQ(-180, hamil.calc_E(&field));

    field = gen_2d_heis_fm(-1, 0, 0);
    EXPECT_EQ(-180, hamil.calc_E(&field));

    field = gen_2d_heis_fm(0, 1, 0);
    EXPECT_EQ(-180, hamil.calc_E(&field));

    field = gen_2d_heis_fm(0, -1, 0);
    EXPECT_EQ(-180, hamil.calc_E(&field));

    field = gen_2d_heis_fm(0, 0, 1);
    EXPECT_EQ(-190, hamil.calc_E(&field));

    field = gen_2d_heis_fm(0, 0, -1);
    EXPECT_EQ(-170, hamil.calc_E(&field));
}

TEST(Heis_model, 2d_ferromagnetic_mag)
{
    field_2d_h field = gen_2d_heis_fm(1, 0, 0);
    ham_heis hamil(0, 1);
    hamil.init_dim(&field);
    vector<double> mag = hamil.calc_M(&field);
    EXPECT_EQ(100, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_2d_heis_fm(-1, 0, 0);
    mag = hamil.calc_M(&field);
    EXPECT_EQ(-100, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_2d_heis_fm(0, 1, 0);
    mag = hamil.calc_M(&field);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(100, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_2d_heis_fm(0, -1, 0);
    mag = hamil.calc_M(&field);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(-100, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_2d_heis_fm(0, 0, 1);
    mag = hamil.calc_M(&field);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(100, mag[2]);

    field = gen_2d_heis_fm(0, 0, -1);
    mag = hamil.calc_M(&field);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(-100, mag[2]);
}

TEST(Heis_model, 2d_ferromagnetic_submag)
{
    field_2d_h field = gen_2d_heis_fm(1, 0, 0);
    ham_heis hamil(0, 1);
    hamil.init_dim(&field);
    vector<double> mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(50, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_2d_heis_fm(-1, 0, 0);
    mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(-50, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_2d_heis_fm(0, 1, 0);
    mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(50, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_2d_heis_fm(0, -1, 0);
    mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(-50, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_2d_heis_fm(0, 0, 1);
    mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(50, mag[2]);

    field = gen_2d_heis_fm(0, 0, -1);
    mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(-50, mag[2]);
}

TEST(Heis_model, 2d_dE_consist)
{
    field_2d_h field = gen_2d_heis_fm(1, 0, 0);
    ham_heis hamil(0, 1);
    hamil.init_dim(&field);
    vector<int> pos(2);
    double old_E = hamil.calc_E(&field);
    for(int i = 0; i < 1000; i++)
    {
        pos[0] = int(st_rand_double.gen()*10 + 1);
        pos[1] = int(st_rand_double.gen()*10 + 1);
        double dE = hamil.dE(&field, pos);
        field.change_to_test(pos, &hamil);
        double new_E = hamil.calc_E(&field);
        // EXPECT_EQ(old_E + dE, new_E);
        EXPECT_NEAR(old_E + dE, new_E, abs(new_E*1e-10));
        old_E = new_E;
    }
}

#endif
