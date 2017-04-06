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
    EXPECT_NEAR(-94277.38689422168, hamil.calc_E(&field), 94277.38689422168 * 1e-10);

    field = gen_3d_heis_fm(0, 1, 0);
    EXPECT_NEAR(-94277.38689422168, hamil.calc_E(&field), 94277.38689422168 * 1e-10);

    field = gen_3d_heis_fm(0, 0, 1);
    EXPECT_NEAR(-95458.5206541323, hamil.calc_E(&field), 95458.5206541323 * 1e-10);

    field = gen_3d_heis_fm(-1, 0, 0);
    EXPECT_NEAR(-94277.38689422168, hamil.calc_E(&field), 94277.38689422168 * 1e-10);

    field = gen_3d_heis_fm(0, -1, 0);
    EXPECT_NEAR(-94277.38689422168, hamil.calc_E(&field), 94277.38689422168 * 1e-10);

    field = gen_3d_heis_fm(0, 0, -1);
    EXPECT_NEAR(-95458.5206541323, hamil.calc_E(&field), 95458.5206541323 * 1e-10);

    // Anti
    field = gen_3d_heis_afm(1, 0, 0);
    EXPECT_NEAR(-93205.65884639179, hamil.calc_E(&field), 93205.65884639179 * 1e-10);

    field = gen_3d_heis_afm(0, 1, 0);
    EXPECT_NEAR(-93205.65884639179, hamil.calc_E(&field), 93205.65884639179 * 1e-10);

    field = gen_3d_heis_afm(0, 0, 1);
    EXPECT_NEAR(-93074.44198019, hamil.calc_E(&field), 93074.44198019 * 1e-10);

    field = gen_3d_heis_afm(-1, 0, 0);
    EXPECT_NEAR(-93205.65884639179, hamil.calc_E(&field), 93205.65884639179 * 1e-10);

    field = gen_3d_heis_afm(0, -1, 0);
    EXPECT_NEAR(-93205.65884639179, hamil.calc_E(&field), 93205.65884639179 * 1e-10);

    field = gen_3d_heis_afm(0, 0, -1);
    EXPECT_NEAR(-93074.44198019, hamil.calc_E(&field), 93074.44198019 * 1e-10);
}

TEST(FePt, 3d_energy_ext_field)
{
    field_3d_h field = gen_3d_heis_fm(1, 0, 0);
    ham_FePt hamil(1);
    hamil.init_dim(&field);
    EXPECT_NEAR(-94277.38689422168, hamil.calc_E(&field), 94277.38689422168 * 1e-10);

    field = gen_3d_heis_fm(0, 1, 0);
    EXPECT_NEAR(-94277.38689422168, hamil.calc_E(&field), 94277.38689422168 * 1e-10);

    field = gen_3d_heis_fm(0, 0, 1);
    EXPECT_NEAR(-96458.5206541323, hamil.calc_E(&field), 95458.5206541323 * 1e-10);

    field = gen_3d_heis_fm(-1, 0, 0);
    EXPECT_NEAR(-94277.38689422168, hamil.calc_E(&field), 94277.38689422168 * 1e-10);

    field = gen_3d_heis_fm(0, -1, 0);
    EXPECT_NEAR(-94277.38689422168, hamil.calc_E(&field), 94277.38689422168 * 1e-10);

    field = gen_3d_heis_fm(0, 0, -1);
    EXPECT_NEAR(-94458.5206541323, hamil.calc_E(&field), 95458.5206541323 * 1e-10);

    // Anti
    field = gen_3d_heis_afm(1, 0, 0);
    EXPECT_NEAR(-93205.65884639179, hamil.calc_E(&field), 93205.65884639179 * 1e-10);

    field = gen_3d_heis_afm(0, 1, 0);
    EXPECT_NEAR(-93205.65884639179, hamil.calc_E(&field), 93205.65884639179 * 1e-10);

    field = gen_3d_heis_afm(0, 0, 1);
    EXPECT_NEAR(-93074.44198019, hamil.calc_E(&field), 93074.44198019 * 1e-10);

    field = gen_3d_heis_afm(-1, 0, 0);
    EXPECT_NEAR(-93205.65884639179, hamil.calc_E(&field), 93205.65884639179 * 1e-10);

    field = gen_3d_heis_afm(0, -1, 0);
    EXPECT_NEAR(-93205.65884639179, hamil.calc_E(&field), 93205.65884639179 * 1e-10);

    field = gen_3d_heis_afm(0, 0, -1);
    EXPECT_NEAR(-93074.44198019, hamil.calc_E(&field), 93074.44198019 * 1e-10);
}

TEST(FePt, 3d_mag)
{
    field_3d_h field = gen_3d_heis_fm(1, 0, 0);
    ham_FePt hamil;
    hamil.init_dim(&field);
    vector<double> mag = hamil.calc_M(&field);
    EXPECT_EQ(1000, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_3d_heis_fm(-1, 0, 0);
    mag = hamil.calc_M(&field);
    EXPECT_EQ(-1000, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_3d_heis_fm(0, 1, 0);
    mag = hamil.calc_M(&field);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(1000, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_3d_heis_fm(0, -1, 0);
    mag = hamil.calc_M(&field);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(-1000, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_3d_heis_fm(0, 0, 1);
    mag = hamil.calc_M(&field);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(1000, mag[2]);

    field = gen_3d_heis_fm(0, 0, -1);
    mag = hamil.calc_M(&field);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(-1000, mag[2]);

    field = gen_3d_heis_afm(1, 0, 0);
    mag = hamil.calc_M(&field);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_3d_heis_afm(-1, 0, 0);
    mag = hamil.calc_M(&field);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_3d_heis_afm(0, 1, 0);
    mag = hamil.calc_M(&field);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_3d_heis_afm(0, -1, 0);
    mag = hamil.calc_M(&field);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_3d_heis_afm(0, 0, 1);
    mag = hamil.calc_M(&field);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_3d_heis_afm(0, 0, -1);
    mag = hamil.calc_M(&field);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);
}

TEST(FePt, 3d_submag)
{
    field_3d_h field = gen_3d_heis_fm(1, 0, 0);
    ham_FePt hamil;
    hamil.init_dim(&field);
    vector<double> mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(500, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_3d_heis_fm(-1, 0, 0);
    mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(-500, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_3d_heis_fm(0, 1, 0);
    mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(500, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_3d_heis_fm(0, -1, 0);
    mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(-500, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_3d_heis_fm(0, 0, 1);
    mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(500, mag[2]);

    field = gen_3d_heis_fm(0, 0, -1);
    mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(-500, mag[2]);

    field = gen_3d_heis_afm(1, 0, 0);
    mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(500, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_3d_heis_afm(-1, 0, 0);
    mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(-500, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_3d_heis_afm(0, 1, 0);
    mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(500, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_3d_heis_afm(0, -1, 0);
    mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(-500, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_3d_heis_afm(0, 0, 1);
    mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(500, mag[2]);

    field = gen_3d_heis_afm(0, 0, -1);
    mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(-500, mag[2]);
}

TEST(FePt, 3d_dE_consist)
{
    field_3d_h field = gen_3d_heis_fm(1, 0, 0);
    ham_FePt hamil(4);
    hamil.init_dim(&field);
    vector<int> pos(3);
    double old_E = hamil.calc_E(&field);
    for(int i = 0; i < 100; i++)
    {
        pos[0] = int(st_rand_double.gen()*10 + 1);
        pos[1] = int(st_rand_double.gen()*10 + 1);
        pos[2] = int(st_rand_double.gen()*10 + 1);
        double dE = hamil.dE(&field, pos);
        field.change_to_test(pos, &hamil);
        double new_E = hamil.calc_E(&field);
        EXPECT_NEAR(old_E + dE, new_E, abs(new_E*1e-10));
        old_E = new_E;
    }
}

#endif
