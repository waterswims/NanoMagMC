#ifndef _CPOINTTEST
#define _CPOINTTEST

#include "../includes/cpoints.hpp"
#include <gtest/gtest.h>
#include <fstream>

///////////////////////////////////////////////////////
// Checkpointing Tests
///////////////////////////////////////////////////////

TEST(Parameter_Checkpointing, Cpoint_Name)
{
    EXPECT_EQ(cpointname("pre", 1, 2, 3, 'a', 'b', 4), "Checkpoints/pre_32ab4_1.cp");
}

TEST(Parameter_Checkpointing, Read_vals)
{
    string filename = "tests/Test_Files/val_load_test.cp";
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
    string filename = "tests/Test_Files/arr_load_test.cp";
    fstream str;
    double** testvals = alloc_2darr<double>(2, 3);
    double** corrvals = alloc_2darr<double>(2, 3);
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
    dealloc_2darr<double>(2, testvals);
    dealloc_2darr<double>(2, corrvals);
}

TEST(Parameter_Checkpointing, Write_vals)
{
    string filename = "tests/Test_Files/val_save_test.cp";
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
    string filename = "tests/Test_Files/arr_save_test.cp";
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
    double** invals = alloc_2darr<double>(5, 100);
    EXPECT_EQ(5, read_clist(stream, filename, invals));
    for(int i = 0; i < 5; i++)
    {
        for(int j = 0; j < 100; j++)
        {
            EXPECT_EQ(i*j, invals[i][j]);
        }
    }
    dealloc_2darr<double>(5, invals);
}

#endif
