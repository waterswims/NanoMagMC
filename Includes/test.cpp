#include "mklrand.h"
#include "cpoints.hpp"
#include "field_type_new.hpp"
#include "hamiltonian_new.hpp"
#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <stdio.h>

using namespace std;

mkl_irand st_rand_int(1e5, 1);
mkl_drand st_rand_double(1e5, 2);
mkl_lnrand rand_ln(0, 0.25, 1e5, 3);

const double pi = 3.141592653589793;

///////////////////////////////////////////////////////
// Non-test functions
///////////////////////////////////////////////////////

double chi2(vector<int>& count, vector<double>& expect)
{
    double ans = 0;
    int skip = 0;
    for(int i=0; i < count.size(); i++)
    {
        if(count[i] == 0)
        {
            skip++;
            continue;
        }
        ans += pow((count[i] - expect[i]), 2) / expect[i];
    }
    return ans / (count.size() - 1 - skip);
}

double beta(int a, int b)
{
    return tgamma(a)*tgamma(b)/(tgamma(a+b));
}

double binomial(int n, int k)
{
    return 1 / (k * beta(k, n-k+1));
}

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

///////////////////////////////////////////////////////
// Random Number Tests
///////////////////////////////////////////////////////

TEST(Random_Numbers, Integer_Test)
{
    int N_atts = 1e6, N_toss = 100;
    vector<int> bins(N_toss/2, 0);
    for(int i=0; i < N_atts; i++)
    {
        int N_zero = 0;
        for (int j=0; j < N_toss; j++)
        {
            if(st_rand_int.gen() == 0){N_zero++;}
        }
        bins[abs(N_zero - N_toss/2)]++;
    }
    vector<double> expect(N_toss/2);
    for(int i = 0; i < N_toss/2; i++)
    {
        expect[i] = binomial(N_toss, (N_toss/2)-i) * N_atts / float(pow(2, (N_toss - 1)));
    }
    expect[0] = expect[0] / 2;
    double chi2_test = chi2(bins, expect);
    EXPECT_GT(chi2_test, 0.9);
    EXPECT_LT(chi2_test, 1.3);
}

TEST(Random_Numbers, Double_Test)
{
    int N_bins = 100, N_atts = 1e6;
    vector<int> bins(N_bins, 0);
    for(int i = 0; i < N_atts; i++)
    {
        bins[int(st_rand_double.gen()*N_bins)]++;
    }
    vector<double> expect(N_bins, (N_atts/N_bins));
    double chi2_test = chi2(bins, expect);
    EXPECT_GT(chi2_test, 0.9);
    EXPECT_LT(chi2_test, 1.3);
}

TEST(Random_Numbers, Ln_Test)
{
    int N_bins = 100, N_atts = 1e6;
    vector<int> bins(N_bins, 0);
    for(int i = 0; i < N_atts; i++)
    {
        double b_num = int(rand_ln.gen()*(N_bins/2.0));
        if (b_num >= N_bins) {i--; continue;}
        bins[b_num]++;
    }
    vector<double> expect(N_bins);
    for(int i=0; i < N_bins; i++)
    {
        double x = (i+0.5)*2.0/float(N_bins);
        double temp = pow(log(x), 2) / 0.125;
        expect[i] = N_atts * exp(-temp) / (50.08490379760093 * x * 0.25 * pow(2*pi, 0.5));
    }
    double chi2_test = chi2(bins, expect);
    EXPECT_GT(chi2_test, 0.9);
    EXPECT_LT(chi2_test, 1.3);
}

TEST(Random_Numbers, Read_Checkpoint)
{
    st_rand_double.load("Includes/Test_Files/double_load_test.cp");
    st_rand_double.fill();
    double loaded;
    ifstream in;
    in.open("Includes/Test_Files/double_load_test.result");
    for (int i = 0; i < 100; i++)
    {
        in >> loaded;
        EXPECT_NEAR(st_rand_double.gen(), loaded, 1e-6);
    }
    in.close();
}

TEST(Random_Numbers, Save_Checkpoint)
{
    vector<double> first(100);
    st_rand_double.save("Includes/Test_Files/double_save_test.cp");
    st_rand_double.fill();
    for(int i = 0; i < 100; i++)
    {
        first[i] = st_rand_double.gen();
    }
    st_rand_double.load("Includes/Test_Files/double_save_test.cp");
    st_rand_double.fill();
    for(int i = 0; i < 100; i++)
    {
        EXPECT_EQ(first[i], st_rand_double.gen());
    }
}

///////////////////////////////////////////////////////
// Checkpointing Tests
///////////////////////////////////////////////////////

TEST(Parameter_Checkpointing, Cpoint_Name)
{
    EXPECT_EQ(cpointname("pre", 1, 2, 3, 'a', 'b', 4), "Checkpoints/pre_32ab4_1.cp");
}

TEST(Parameter_Checkpointing, Read_vals)
{
    string filename = "Includes/Test_Files/val_load_test.cp";
    fstream str;
    int testvals[6];
    int corrvals[6] = {4, 2, 3, 6, 7, 8};
    read_cval(str, filename, testvals);
    for(int i = 0; i < 6; i++)
    {
        EXPECT_EQ(testvals[i], corrvals[i]);
    }
}

TEST(Parameter_Checkpointing, Read_arrs)
{
    string filename = "Includes/Test_Files/arr_load_test.cp";
    fstream str;
    double testvals[2][100];
    double corrvals[2][3];
    corrvals[0][0] = 2.4;
    corrvals[0][1] = 4.6;
    corrvals[0][2] = 6.7;
    corrvals[1][0] = 3.2;
    corrvals[1][1] = 42.6;
    corrvals[1][2] = -11.2;
    read_clist(str, filename, testvals);
    for(int i = 0; i < 2; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            EXPECT_EQ(testvals[i][j], corrvals[i][j]);
        }
    }
}

TEST(Parameter_Checkpointing, Write_vals)
{
    string filename = "Includes/Test_Files/val_save_test.cp";
    remove(filename.c_str());
    fstream stream;
    for(int i = 0; i < 100; i++)
    {
        print_cval(stream, filename, i);
    }
    int invals[100];
    EXPECT_EQ(100, read_cval(stream, filename, invals));
    for(int i = 0; i < 100; i++)
    {
        EXPECT_EQ(i, invals[i]);
    }
}

TEST(Parameter_Checkpointing, Write_arr)
{
    string filename = "Includes/Test_Files/arr_save_test.cp";
    remove(filename.c_str());
    fstream stream;
    double outvals[100];
    for(int i = 0; i < 5; i++)
    {
        for(int j = 0; j < 100; j++)
        {
            outvals[j] = i*j;
        }
        print_clist(stream, filename, outvals, 100);
    }
    double invals[5][100];
    EXPECT_EQ(5, read_clist(stream, filename, invals));
    for(int i = 0; i < 5; i++)
    {
        for(int j = 0; j < 100; j++)
        {
            EXPECT_EQ(i*j, invals[i][j]);
        }
    }
}

///////////////////////////////////////////////////////
// Ising model tests
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

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
