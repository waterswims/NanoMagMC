#ifndef _ISINGTEST
#define _ISINGTEST
#include "field_type_new.hpp"
#include "hamiltonian_new.hpp"

#ifndef _TESTRANDGEN
#define _TESTRANDGEN
mkl_irand st_rand_int(1e5, 1);
mkl_drand st_rand_double(1e5, 2);
mkl_lnrand rand_ln(0, 0.25, 1e5, 3);
#endif

using namespace std;

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
    field.fill_ghost();
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
    field.fill_ghost();
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
    field.fill_ghost();
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
    field.fill_ghost();
    return field;
}

///////////////////////////////////////////////////////
// Ising model tests - 2D
///////////////////////////////////////////////////////

TEST(Ising_model, 2d_ferromagnetic_energy_zero_field)
{
    field_2d_i field = gen_2d_ising_fm();
    ham_ising hamil(0, 1);
    hamil.init_dim(&field);
    EXPECT_EQ(-180, hamil.calc_E(&field));
}

TEST(Ising_model, 2d_ferromagnetic_energy_ext_field)
{
    field_2d_i field = gen_2d_ising_fm();
    ham_ising hamil(0.1, 1);
    hamil.init_dim(&field);
    EXPECT_EQ(-190, hamil.calc_E(&field));
}

TEST(Ising_model, 2d_ferromagnetic_mag)
{
    field_2d_i field = gen_2d_ising_fm();
    ham_ising hamil(0, 1);
    hamil.init_dim(&field);
    EXPECT_EQ(100, hamil.calc_M(&field)[0]);
}

TEST(Ising_model, 2d_ferromagnetic_submag)
{
    field_2d_i field = gen_2d_ising_fm();
    ham_ising hamil(0, 1);
    hamil.init_dim(&field);
    EXPECT_EQ(50, hamil.calc_subM(&field, 1)[0]);
}

TEST(Ising_model, 2d_ferromagnetic_dE)
{
    field_2d_i field = gen_2d_ising_fm();
    ham_ising hamil(0, 1);
    hamil.init_dim(&field);
    vector<int> pos(2,5);
    EXPECT_EQ(8, hamil.dE(&field, pos));
}

TEST(Ising_model, 2d_antiferromagnetic_energy_zero_field)
{
    field_2d_i field = gen_2d_ising_afm();
    ham_ising hamil(0, 1);
    hamil.init_dim(&field);
    EXPECT_EQ(180, hamil.calc_E(&field));
}

TEST(Ising_model, 2d_antiferromagnetic_energy_ext_field)
{
    field_2d_i field = gen_2d_ising_afm();
    ham_ising hamil(0.3, 1);
    hamil.init_dim(&field);
    EXPECT_EQ(180, hamil.calc_E(&field));
}

TEST(Ising_model, 2d_antiferromagnetic_mag)
{
    field_2d_i field = gen_2d_ising_afm();
    ham_ising hamil(0, 1);
    hamil.init_dim(&field);
    EXPECT_EQ(0, hamil.calc_M(&field)[0]);
}

TEST(Ising_model, 2d_antiferromagnetic_submag)
{
    field_2d_i field = gen_2d_ising_afm();
    ham_ising hamil(0, 1);
    hamil.init_dim(&field);
    EXPECT_EQ(50, hamil.calc_subM(&field, 1)[0]);
}

TEST(Ising_model, 2d_antiferromagnetic_dE)
{
    field_2d_i field = gen_2d_ising_afm();
    ham_ising hamil(0, 1);
    hamil.init_dim(&field);
    vector<int> pos(2,5);
    EXPECT_EQ(-8, hamil.dE(&field, pos));
}

TEST(Ising_model, 2d_dE_consist)
{
    field_2d_i field = gen_2d_ising_fm();
    ham_ising hamil(0, 1);
    hamil.init_dim(&field);
    vector<int> pos(2);
    int old_E = hamil.calc_E(&field);
    for(int i = 0; i < 1000; i++)
    {
        pos[0] = int(st_rand_double.gen()*10 + 1);
        pos[1] = int(st_rand_double.gen()*10 + 1);
        int dE = hamil.dE(&field, pos);
        field.change_to_test(pos, &hamil);
        int new_E = hamil.calc_E(&field);
        EXPECT_EQ(old_E + dE, new_E);
        old_E = new_E;
    }
}

///////////////////////////////////////////////////////
// Ising model tests - 3D
///////////////////////////////////////////////////////

TEST(Ising_model, 3d_ferromagnetic_energy_zero_field)
{
    field_3d_i field = gen_3d_ising_fm();
    ham_ising hamil(0, 1);
    hamil.init_dim(&field);
    EXPECT_EQ(-2700, hamil.calc_E(&field));
}

TEST(Ising_model, 3d_ferromagnetic_energy_ext_field)
{
    field_3d_i field = gen_3d_ising_fm();
    ham_ising hamil(0.1, 1);
    hamil.init_dim(&field);
    EXPECT_EQ(-2800, hamil.calc_E(&field));
}

TEST(Ising_model, 3d_ferromagnetic_mag)
{
    field_3d_i field = gen_3d_ising_fm();
    ham_ising hamil(0, 1);
    hamil.init_dim(&field);
    EXPECT_EQ(1000, hamil.calc_M(&field)[0]);
}

TEST(Ising_model, 3d_ferromagnetic_submag)
{
    field_3d_i field = gen_3d_ising_fm();
    ham_ising hamil(0, 1);
    hamil.init_dim(&field);
    EXPECT_EQ(500, hamil.calc_subM(&field, 1)[0]);
}

TEST(Ising_model, 3d_ferromagnetic_dE)
{
    field_3d_i field = gen_3d_ising_fm();
    ham_ising hamil(0, 1);
    hamil.init_dim(&field);
    vector<int> pos(3,5);
    EXPECT_EQ(12, hamil.dE(&field, pos));
}

TEST(Ising_model, 3d_antiferromagnetic_energy_zero_field)
{
    field_3d_i field = gen_3d_ising_afm();
    ham_ising hamil(0, 1);
    hamil.init_dim(&field);
    EXPECT_EQ(2700, hamil.calc_E(&field));
}

TEST(Ising_model, 3d_antiferromagnetic_energy_ext_field)
{
    field_3d_i field = gen_3d_ising_afm();
    ham_ising hamil(0.1, 1);
    hamil.init_dim(&field);
    EXPECT_EQ(2700, hamil.calc_E(&field));
}

TEST(Ising_model, 3d_antiferromagnetic_mag)
{
    field_3d_i field = gen_3d_ising_afm();
    ham_ising hamil(0, 1);
    hamil.init_dim(&field);
    EXPECT_EQ(0, hamil.calc_M(&field)[0]);
}

TEST(Ising_model, 3d_antiferromagnetic_submag)
{
    field_3d_i field = gen_3d_ising_fm();
    ham_ising hamil(0, 1);
    hamil.init_dim(&field);
    EXPECT_EQ(500, hamil.calc_subM(&field, 1)[0]);
}

TEST(Ising_model, 3d_antiferromagnetic_dE)
{
    field_3d_i field = gen_3d_ising_fm();
    ham_ising hamil(0, 1);
    hamil.init_dim(&field);
    vector<int> pos(3,5);
    EXPECT_EQ(12, hamil.dE(&field, pos));
}

TEST(Ising_model, 3d_dE_consist)
{
    field_3d_i field = gen_3d_ising_fm();
    ham_ising hamil(0, 1);
    hamil.init_dim(&field);
    vector<int> pos(3);
    int old_E = hamil.calc_E(&field);
    for(int i = 0; i < 1000; i++)
    {
        pos[0] = int(st_rand_double.gen()*10 + 1);
        pos[1] = int(st_rand_double.gen()*10 + 1);
        pos[3] = int(st_rand_double.gen()*10 + 1);
        int dE = hamil.dE(&field, pos);
        field.change_to_test(pos, &hamil);
        int new_E = hamil.calc_E(&field);
        EXPECT_EQ(old_E + dE, new_E);
        old_E = new_E;
    }
}

#endif
