#ifndef _ISINGTEST
#define _ISINGTEST

#include "../includes/thermodynamics.hpp"

xt::xtensorf<double, xt::xshape<4>> zeroH = {0, 0, 0, 0};
xt::xtensorf<double, xt::xshape<4>> isingH = {0.1, 0, 0, 0};
particle::field::field_type isingFMField;
particle::field::field_type isingFMFieldPerio;
particle::field::field_type isingAFMField;
particle::td::functionObject exchangeOnly;

///////////////////////////////////////////////////////
// Ising model tests - 2D
///////////////////////////////////////////////////////

// TEST(Ising_model, 2d_ferromagnetic_energy_zero_field)
// {
//     field_2d_i field = gen_2d_ising_fm();
//     ham_ising hamil(0, 1);
//     hamil.init_dim(&field);
//     EXPECT_EQ(-180, hamil.calc_E(&field));
// }
//
// TEST(Ising_model, 2d_ferromagnetic_energy_ext_field)
// {
//     field_2d_i field = gen_2d_ising_fm();
//     ham_ising hamil(0.1, 1);
//     hamil.init_dim(&field);
//     EXPECT_EQ(-190, hamil.calc_E(&field));
// }
//
// TEST(Ising_model, 2d_ferromagnetic_mag)
// {
//     field_2d_i field = gen_2d_ising_fm();
//     ham_ising hamil(0, 1);
//     hamil.init_dim(&field);
//     EXPECT_EQ(100, hamil.calc_M(&field)[0]);
// }
//
// TEST(Ising_model, 2d_ferromagnetic_submag)
// {
//     field_2d_i field = gen_2d_ising_fm();
//     ham_ising hamil(0, 1);
//     hamil.init_dim(&field);
//     EXPECT_EQ(50, hamil.calc_subM(&field, 1)[0]);
// }
//
// TEST(Ising_model, 2d_ferromagnetic_dE)
// {
//     field_2d_i field = gen_2d_ising_fm();
//     ham_ising hamil(0, 1);
//     hamil.init_dim(&field);
//     std::vector<int> pos(2,5);
//     EXPECT_EQ(8, hamil.dE(&field, pos));
// }
//
// TEST(Ising_model, 2d_antiferromagnetic_energy_zero_field)
// {
//     field_2d_i field = gen_2d_ising_afm();
//     ham_ising hamil(0, 1);
//     hamil.init_dim(&field);
//     EXPECT_EQ(180, hamil.calc_E(&field));
// }
//
// TEST(Ising_model, 2d_antiferromagnetic_energy_ext_field)
// {
//     field_2d_i field = gen_2d_ising_afm();
//     ham_ising hamil(0.3, 1);
//     hamil.init_dim(&field);
//     EXPECT_EQ(180, hamil.calc_E(&field));
// }
//
// TEST(Ising_model, 2d_antiferromagnetic_mag)
// {
//     field_2d_i field = gen_2d_ising_afm();
//     ham_ising hamil(0, 1);
//     hamil.init_dim(&field);
//     EXPECT_EQ(0, hamil.calc_M(&field)[0]);
// }
//
// TEST(Ising_model, 2d_antiferromagnetic_submag)
// {
//     field_2d_i field = gen_2d_ising_afm();
//     ham_ising hamil(0, 1);
//     hamil.init_dim(&field);
//     EXPECT_EQ(50, hamil.calc_subM(&field, 1)[0]);
// }
//
// TEST(Ising_model, 2d_antiferromagnetic_dE)
// {
//     field_2d_i field = gen_2d_ising_afm();
//     ham_ising hamil(0, 1);
//     hamil.init_dim(&field);
//     std::vector<int> pos(2,5);
//     EXPECT_EQ(-8, hamil.dE(&field, pos));
// }
//
// TEST(Ising_model, 2d_dE_consist)
// {
//     field_2d_i field = gen_2d_ising_fm();
//     ham_ising hamil(0, 1);
//     hamil.init_dim(&field);
//     std::vector<int> pos(2);
//     int old_E = hamil.calc_E(&field);
//     for(int i = 0; i < 1000; i++)
//     {
//         pos[0] = int(st_rand_double.gen()*10 + 1);
//         pos[1] = int(st_rand_double.gen()*10 + 1);
//         int dE = hamil.dE(&field, pos);
//         field.change_to_test(pos, &hamil);
//         int new_E = hamil.calc_E(&field);
//         EXPECT_EQ(old_E + dE, new_E);
//         old_E = new_E;
//     }
// }
//
// ///////////////////////////////////////////////////////
// // Ising model tests - 3D
// ///////////////////////////////////////////////////////

TEST(Ising_model, 3d_ferromagnetic_energy_zero_field)
{
    EXPECT_DOUBLE_EQ(-2700, exchangeOnly.calc_E(isingFMField, zeroH));
    EXPECT_DOUBLE_EQ(-3000, exchangeOnly.calc_E(isingFMFieldPerio, zeroH));
}

TEST(Ising_model, 3d_ferromagnetic_energy_ext_field)
{
    EXPECT_DOUBLE_EQ(-2800, exchangeOnly.calc_E(isingFMField, isingH));
    EXPECT_DOUBLE_EQ(-3100, exchangeOnly.calc_E(isingFMFieldPerio, isingH));
}

TEST(Ising_model, 3d_ferromagnetic_mag)
{
    EXPECT_DOUBLE_EQ(1000, exchangeOnly.calc_M(isingFMField)[0]);
    EXPECT_DOUBLE_EQ(1000, exchangeOnly.calc_M(isingFMFieldPerio)[0]);
}

TEST(Ising_model, 3d_ferromagnetic_submag)
{
    EXPECT_DOUBLE_EQ(500, exchangeOnly.calc_subM(isingFMField, 0)[0]);
    EXPECT_DOUBLE_EQ(500, exchangeOnly.calc_subM(isingFMFieldPerio, 0)[0]);
}

TEST(Ising_model, 3d_antiferromagnetic_energy_zero_field)
{
    EXPECT_DOUBLE_EQ(2700, exchangeOnly.calc_E(isingAFMField, zeroH));
}

TEST(Ising_model, 3d_antiferromagnetic_energy_ext_field)
{
    EXPECT_DOUBLE_EQ(2700, exchangeOnly.calc_E(isingAFMField, isingH));
}

TEST(Ising_model, 3d_antiferromagnetic_mag)
{
    EXPECT_DOUBLE_EQ(0, exchangeOnly.calc_M(isingAFMField)[0]);
}

TEST(Ising_model, 3d_antiferromagnetic_submag)
{
    EXPECT_DOUBLE_EQ(500, exchangeOnly.calc_subM(isingAFMField, 0)[0]);
}

TEST(Ising_model, 3d_dE_consist)
{
    particle::field::field_type fieldCopy = isingFMField;
    int pos=0;
    double old_E = exchangeOnly.calc_E(fieldCopy, isingH);
    for(int i = 0; i < 1000; i++)
    {
        pos = int(st_rand_double.gen()*1000);
        fieldCopy.gen_rand();
        double dE = exchangeOnly.calc_dE(fieldCopy, pos, isingH);
        fieldCopy.set_rand(pos);
        double new_E = exchangeOnly.calc_E(fieldCopy, isingH);
        EXPECT_NEAR(old_E + dE, new_E, 1e-10);
        old_E = new_E;
    }
}

#endif
