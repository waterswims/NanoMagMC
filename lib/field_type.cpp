#include "../includes/field_type.hpp"
#include "../includes/array_alloc.hpp"
#include "../includes/mklrand.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <hdf5.h>

///////////////////////
// Global Variables
///////////////////////

extern mkl_irand st_rand_int;
extern mkl_drand st_rand_double;
const double pi = 3.141592653589793;

///////////////////////
// Random Functions
///////////////////////

int boundmovedown(int test, int limit)
{
    return ((test%limit)+limit)%limit;
}

int boundmoveup(int test, int limit)
{
    return test%limit;
}

void rand_spin_h(double &x, double &y, double &z)
{
    double phi = st_rand_double.gen()*2*pi;
    double cthet = 2*st_rand_double.gen()-1;
    double sthet = pow(1 - pow(cthet, 2), 0.5);
    x = cos(phi)*sthet;
    y = sin(phi)*sthet;
    z = cthet;
}

void rand_spin_i(int &x)
{
    x = st_rand_int.gen();
}

///////////////////////
// All models
///////////////////////

void field_type::allzero()
{
    int start = (totsize - insize) / 2;
    vector<int> pos(dim, start);
    bool finished = false;

    while(!finished)
    {
        this->fill_zero(pos);
        this->next(finished, pos);
    }
}

///////////////////////
// 1d Heis-model
///////////////////////

field_cluster_h::field_cluster_h()
{
    dim = 1;
    periodic = false;
    ft = 1;
}

field_cluster_h::field_cluster_h(string filename)
{
    dim = 1;
    periodic = false;
    ft = 1;

    ifstream file;
    file.open(filename.c_str());
    if(!file.is_open())
    {
        cout << "Input file not opened" << endl;
        exit(105);
    }
    insize = 0;
    string line;
    while(getline(file, line))
    {
        insize++;
    }
    file.close();

    double temp_d;
    totsize=this->insize;
    spinx = alloc_1darr<double>(insize);
    spiny = alloc_1darr<double>(insize);
    spinz = alloc_1darr<double>(insize);

    for(int i=0; i < insize; i++)
    {
        rand_spin_h(spinx[i], spiny[i], spinz[i]);
    }
}

field_cluster_h::field_cluster_h(field_type& other)
{
    dim = 1;
    periodic = 1;
    ft = 1;
    if(other.get_ft() != 1)
    {
        cout << "Cannot copy from other field type" << endl;
        exit(104);
    }
    insize = other.get_insize();
    totsize = other.get_totsize();
    double* xoth;
    double* yoth;
    double* zoth;
    other.get_1dfield_h(xoth, yoth, zoth);
    spinx = deep_copy_1darr<double>(insize, xoth);
    spiny = deep_copy_1darr<double>(insize, yoth);
    spinz = deep_copy_1darr<double>(insize, zoth);
}

field_cluster_h::field_cluster_h(const field_cluster_h& other)
{
    dim = 1;
    periodic = 1;
    ft = 1;
    insize = other.get_insize();
    totsize = other.get_totsize();
    double* xoth;
    double* yoth;
    double* zoth;
    other.get_1dfield_h(xoth, yoth, zoth);
    spinx = deep_copy_1darr<double>(insize, xoth);
    spiny = deep_copy_1darr<double>(insize, yoth);
    spinz = deep_copy_1darr<double>(insize, zoth);
}

field_cluster_h& field_cluster_h::operator=(const field_cluster_h& other)
{
    dim = 1;
    periodic = 1;
    ft = 1;
    insize = other.get_insize();
    totsize = other.get_totsize();
    double* xoth;
    double* yoth;
    double* zoth;
    other.get_1dfield_h(xoth, yoth, zoth);
    spinx = deep_copy_1darr<double>(insize, xoth);
    spiny = deep_copy_1darr<double>(insize, yoth);
    spinz = deep_copy_1darr<double>(insize, zoth);
    return *this;
}

field_cluster_h::~field_cluster_h()
{
    dealloc_1darr<double>(spinx);
    dealloc_1darr<double>(spiny);
    dealloc_1darr<double>(spinz);
}

void field_cluster_h::h_access(vector<int>& position, vector<double>& out)
{
    out[0] = spinx[position[0]];
    out[1] = spiny[position[0]];
    out[2] = spinz[position[0]];
}

void field_cluster_h::h_next(bool &finish, vector<int> &pos, vector<double> &out)
{
    this->h_access(pos, out);
    pos[0]++;
    if(pos[0] == insize)
    {
        pos[0] = 0;
        finish = true;
    }
}

void field_cluster_h::get_1dfield_h(double* &x, double* &y, double* &z) const
{
    x = spinx;
    y = spiny;
    z = spinz;
}

void field_cluster_h::change_to_test(vector<int>& position, ham_type* hamil)
{
    hamil->get_test(spinx[position[0]],spiny[position[0]],spinz[position[0]]);
}

///////////////////////
// 2d general
///////////////////////

void field_2d::next(bool &finish, vector<int> &pos)
{
    int start = (totsize - insize) / 2;
    int end = start + insize;

    pos[1]++;
    if(pos[1] == end)
    {
        pos[1] = start;
        pos[0]++;
        if(pos[0] == end)
        {
            pos[0] = start;
            finish = true;
        }
    }
}

int field_2d::findnum()
{
    int c = 0;
    for(int i = 0; i<totsize; i++)
    {
        for(int j = 0; j<totsize; j++)
        {
            if(!(iszero[i][j]))
            {
                c++;
            }
        }
    }
    return c;
}

///////////////////////
// 2d Heis-model
///////////////////////

field_2d_h::field_2d_h()
{
    dim = 2;
    periodic = true;
    ft = 2;
}

field_2d_h::field_2d_h(int size, bool isperio)
{
    dim = 2;
    periodic = isperio;
    insize = size;
    ft = 2;
    if(periodic)
    {
        totsize = insize;
    }
    else
    {
        totsize = insize + 2;
    }
    spinx = alloc_2darr<double>(totsize, totsize);
    spiny = alloc_2darr<double>(totsize, totsize);
    spinz = alloc_2darr<double>(totsize, totsize);
    iszero = alloc_2darr<bool>(totsize, totsize);
}

field_2d_h::field_2d_h(int size, bool isperio, int p_pad)
{
    dim = 2;
    periodic = isperio;
    insize = size;
    ft = 2;
    if(periodic)
    {
        totsize = insize;
    }
    else
    {
        totsize = insize + 2 * p_pad;
    }
    spinx = alloc_2darr<double>(totsize, totsize);
    spiny = alloc_2darr<double>(totsize, totsize);
    spinz = alloc_2darr<double>(totsize, totsize);
    iszero = alloc_2darr<bool>(totsize, totsize);
}

field_2d_h::field_2d_h(field_type& other)
{
    ft = 2;
    if(other.get_ft() != 2)
    {
        cout << "Cannot copy from other field type" << endl;
        exit(104);
    }
    dim = 2;
    insize = other.get_insize();
    totsize = other.get_totsize();
    periodic = other.get_perio();
    double** xoth;
    double** yoth;
    double** zoth;
    bool** zero_oth;
    other.get_2dfield_h(xoth, yoth, zoth);
    other.get_2dzero(zero_oth);
    spinx = deep_copy_2darr<double>(totsize, totsize, xoth);
    spiny = deep_copy_2darr<double>(totsize, totsize, yoth);
    spinz = deep_copy_2darr<double>(totsize, totsize, zoth);
    iszero = deep_copy_2darr<bool>(totsize, totsize, zero_oth);
}

field_2d_h::field_2d_h(const field_2d_h& other)
{
    ft = 2;
    dim = 2;
    insize = other.get_insize();
    totsize = other.get_totsize();
    periodic = other.get_perio();
    double** xoth;
    double** yoth;
    double** zoth;
    bool** zero_oth;
    other.get_2dfield_h(xoth, yoth, zoth);
    other.get_2dzero(zero_oth);
    spinx = deep_copy_2darr<double>(totsize, totsize, xoth);
    spiny = deep_copy_2darr<double>(totsize, totsize, yoth);
    spinz = deep_copy_2darr<double>(totsize, totsize, zoth);
    iszero = deep_copy_2darr<bool>(totsize, totsize, zero_oth);
}

field_2d_h& field_2d_h::operator=(const field_2d_h& other)
{
    ft = 2;
    dim = 2;
    insize = other.get_insize();
    totsize = other.get_totsize();
    periodic = other.get_perio();
    double** xoth;
    double** yoth;
    double** zoth;
    bool** zero_oth;
    other.get_2dfield_h(xoth, yoth, zoth);
    other.get_2dzero(zero_oth);
    spinx = deep_copy_2darr<double>(totsize, totsize, xoth);
    spiny = deep_copy_2darr<double>(totsize, totsize, yoth);
    spinz = deep_copy_2darr<double>(totsize, totsize, zoth);
    iszero = deep_copy_2darr<bool>(totsize, totsize, zero_oth);
    return *this;
}

field_2d_h::~field_2d_h()
{
    dealloc_2darr<double>(totsize, spinx);
    dealloc_2darr<double>(totsize, spiny);
    dealloc_2darr<double>(totsize, spinz);
    dealloc_2darr<bool>(totsize, iszero);
}

void field_2d_h::h_access(vector<int>& position, vector<double>& out)
{
    out[0] = spinx[position[0]][position[1]];
    out[1] = spiny[position[0]][position[1]];
    out[2] = spinz[position[0]][position[1]];
}

void field_2d_h::h_next(bool &finish, vector<int> &pos, vector<double> &out)
{
    this->h_access(pos, out);
    this->next(finish, pos);
}

void field_2d_h::fill_ghost(int num_rows)
{
    for(int r = 0; r < num_rows; r++)
    {
        for(int i=0; i < totsize; i++)
        {
            spinx[i][r] = 0;
            spiny[i][r] = 0;
            spinz[i][r] = 0;
            iszero[i][r] = true;
            spinx[r][i] = 0;
            spiny[r][i] = 0;
            spinz[r][i] = 0;
            iszero[r][i] = true;
            spinx[i][totsize-1-r] = 0;
            spiny[i][totsize-1-r] = 0;
            spinz[i][totsize-1-r] = 0;
            iszero[i][totsize-1-r] = true;
            spinx[totsize-1-r][i] = 0;
            spiny[totsize-1-r][i] = 0;
            spinz[totsize-1-r][i] = 0;
            iszero[totsize-1-r][i] = true;
        }
    }
}

void field_2d_h::get_2dfield_h(double** &x, double** &y, double** &z) const
{
    x = spinx;
    y = spiny;
    z = spinz;
}

void field_2d_h::h_adjacent(vector<int>& position, double** out)
{
    dirsx[0] = boundmovedown(position[0] - 1, totsize);
    dirsx[1] = boundmoveup(position[0] + 1, totsize);
    dirsx[2] = position[0];
    dirsx[3] = position[0];

    dirsy[0] = position[1];
    dirsy[1] = position[1];
    dirsy[2] = boundmovedown(position[1] - 1, totsize);
    dirsy[3] = boundmoveup(position[1] + 1, totsize);

    for (int i = 0; i < 4; i++)
    {
        out[0][i] = spinx[dirsx[i]][dirsy[i]];
        out[1][i] = spiny[dirsx[i]][dirsy[i]];
        out[2][i] = spinz[dirsx[i]][dirsy[i]];
    }
}

void field_2d_h::fill_rand(vector<int>& position)
{
    rand_spin_h(spinx[position[0]][position[1]],
                spiny[position[0]][position[1]],
                spinz[position[0]][position[1]]);
    iszero[position[0]][position[1]] = false;
}

void field_2d_h::fill_one(vector<int>&position)
{
    spinx[position[0]][position[1]] = 0;
    spiny[position[0]][position[1]] = 0;
    spinz[position[0]][position[1]] = 1;
    iszero[position[0]][position[1]] = false;
}

void field_2d_h::fill_zero(vector<int>& position)
{
    spinx[position[0]][position[1]] = 0;
    spiny[position[0]][position[1]] = 0;
    spinz[position[0]][position[1]] = 0;
    iszero[position[0]][position[1]] = true;
}

void field_2d_h::change_to_test(vector<int>& position, ham_type* hamil)
{
    hamil->get_test(spinx[position[0]][position[1]],
                    spiny[position[0]][position[1]],
                    spinz[position[0]][position[1]]);
}

void field_2d_h::fill_val_h(vector<int>& position, double x, double y, double z)
{
    spinx[position[0]][position[1]] = x;
    spiny[position[0]][position[1]] = y;
    spinz[position[0]][position[1]] = z;
}

void field_2d_h::add_val_h(vector<int>& position, vector<double> &in)
{
    spinx[position[0]][position[1]] += in[0];
    spiny[position[0]][position[1]] += in[1];
    spinz[position[0]][position[1]] += in[2];
}

///////////////////////
// 2d Ising-model
///////////////////////

field_2d_i::field_2d_i()
{
    dim = 2;
    periodic = true;
    ft = 21;
}

field_2d_i::field_2d_i(int size, bool isperio)
{
    dim = 2;
    periodic = isperio;
    insize = size;
    ft = 21;
    if(periodic)
    {
        totsize = insize;
    }
    else
    {
        totsize = insize + 2;
    }
    spin = alloc_2darr<int>(totsize, totsize);
    iszero = alloc_2darr<bool>(totsize, totsize);
}

field_2d_i::field_2d_i(field_type& other)
{
    ft = 21;
    if(other.get_ft() != 21)
    {
        cout << "Cannot copy from other field type" << endl;
        exit(104);
    }
    dim = 2;
    insize = other.get_insize();
    totsize = other.get_totsize();
    periodic = other.get_perio();
    int** othspin;
    bool** zero_oth;
    other.get_2dfield_i(othspin);
    other.get_2dzero(zero_oth);
    spin = deep_copy_2darr<int>(totsize, totsize, othspin);
    iszero = deep_copy_2darr<bool>(totsize, totsize, zero_oth);
}

field_2d_i::field_2d_i(const field_2d_i& other)
{
    ft = 21;
    dim = 2;
    insize = other.get_insize();
    totsize = other.get_totsize();
    periodic = other.get_perio();
    int** othspin;
    bool** zero_oth;
    other.get_2dfield_i(othspin);
    other.get_2dzero(zero_oth);
    spin = deep_copy_2darr<int>(totsize, totsize, othspin);
    iszero = deep_copy_2darr<bool>(totsize, totsize, zero_oth);
}

field_2d_i& field_2d_i::operator=(const field_2d_i& other)
{
    ft = 21;
    dim = 2;
    insize = other.get_insize();
    totsize = other.get_totsize();
    periodic = other.get_perio();
    int** othspin;
    bool** zero_oth;
    other.get_2dfield_i(othspin);
    other.get_2dzero(zero_oth);
    spin = deep_copy_2darr<int>(totsize, totsize, othspin);
    iszero = deep_copy_2darr<bool>(totsize, totsize, zero_oth);
    return *this;
}

field_2d_i::~field_2d_i()
{
    dealloc_2darr<int>(totsize, spin);
    dealloc_2darr<bool>(totsize, iszero);
}

void field_2d_i::i_access(vector<int>& position, int &out)
{
    out = spin[position[0]][position[1]];
}

void field_2d_i::i_next(bool &finish, vector<int> &pos, int &out)
{
    this->i_access(pos, out);
    this->next(finish, pos);
}

void field_2d_i::fill_ghost(int num_rows)
{
    for(int r = 0; r < num_rows; r++)
    {
        for(int i=0; i < totsize; i++)
        {
            spin[i][r] = 0;
            iszero[i][r] = true;
            spin[r][i] = 0;
            iszero[r][i] = true;
            spin[i][totsize-1-r] = 0;
            iszero[i][totsize-1-r] = true;
            spin[totsize-1-r][i] = 0;
            iszero[totsize-1-r][i] = true;
        }
    }
}

void field_2d_i::get_2dfield_i(int** &x) const
{
    x = spin;
}

void field_2d_i::i_adjacent(vector<int>& position, int* out)
{
    int dirsx[4], dirsy[4];
    dirsx[0] = boundmovedown(position[0] - 1, totsize);
    dirsx[1] = boundmoveup(position[0] + 1, totsize);
    dirsx[2] = position[0];
    dirsx[3] = position[0];

    dirsy[0] = position[1];
    dirsy[1] = position[1];
    dirsy[2] = boundmovedown(position[1] - 1, totsize);
    dirsy[3] = boundmoveup(position[1] + 1, totsize);

    for (int i = 0; i < 4; i++)
    {
        out[i] = spin[dirsx[i]][dirsy[i]];
    }
}

void field_2d_i::fill_rand(vector<int>& position)
{
    rand_spin_i(spin[position[0]][position[1]]);
    iszero[position[0]][position[1]] = false;
}

void field_2d_i::fill_one(vector<int>&position)
{
    spin[position[0]][position[1]] = 1;
    iszero[position[0]][position[1]] = false;
}

void field_2d_i::fill_zero(vector<int>& position)
{
    spin[position[0]][position[1]] = 0;
    iszero[position[0]][position[1]] = true;
}

void field_2d_i::fill_val_i(vector<int>& position, int val)
{
    spin[position[0]][position[1]] = val;
}

void field_2d_i::change_to_test(vector<int>& position, ham_type* hamil)
{
    spin[position[0]][position[1]] = -spin[position[0]][position[1]];
}

///////////////////////
// 3d general
///////////////////////

void field_3d::next(bool &finish, vector<int> &pos)
{
    int start = (totsize - insize) / 2;
    int end = start + insize;

    pos[2]++;
    if(pos[2] == end)
    {
        pos[2] = start;
        pos[1]++;
        if(pos[1] == end)
        {
            pos[1] = start;
            pos[0]++;
            if(pos[0] == end)
            {
                pos[0] = start;
                finish = true;
            }
        }
    }
}

int field_3d::findnum()
{
    int c = 0;
    for(int i = 0; i<totsize; i++)
    {
        for(int j = 0; j<totsize; j++)
        {
            for(int k = 0; k<totsize; k++)
            {
                if(!(iszero[i][j][k]))
                {
                    c++;
                }
            }
        }
    }
    return c;
}

///////////////////////
// 3d Heis-model
///////////////////////

field_3d_h::field_3d_h()
{
    dim = 3;
    periodic = true;
    ft = 3;
}

field_3d_h::field_3d_h(int size, bool isperio)
{
    dim = 3;
    periodic = isperio;
    insize = size;
    ft = 3;
    if(periodic)
    {
        totsize = insize;
    }
    else
    {
        totsize = insize + 2;
    }
    spinx = alloc_3darr<double>(totsize, totsize, totsize);
    spiny = alloc_3darr<double>(totsize, totsize, totsize);
    spinz = alloc_3darr<double>(totsize, totsize, totsize);
    iszero = alloc_3darr<bool>(totsize, totsize, totsize);
    postemp = alloc_2darr<int>(3, 1358);
}

field_3d_h::field_3d_h(int size, bool isperio, int p_pad)
{
    dim = 3;
    periodic = isperio;
    insize = size;
    ft = 3;
    if(periodic)
    {
        totsize = insize;
    }
    else
    {
        totsize = insize + 2*p_pad;
    }
    spinx = alloc_3darr<double>(totsize, totsize, totsize);
    spiny = alloc_3darr<double>(totsize, totsize, totsize);
    spinz = alloc_3darr<double>(totsize, totsize, totsize);
    iszero = alloc_3darr<bool>(totsize, totsize, totsize);
    postemp = alloc_2darr<int>(3, 1358);
}

field_3d_h::field_3d_h(field_type& other)
{
    ft = 3;
    if(other.get_ft() != 3)
    {
        cout << "Cannot copy from other field type" << endl;
        exit(104);
    }
    dim = 3;
    insize = other.get_insize();
    totsize = other.get_totsize();
    periodic = other.get_perio();
    double*** xoth;
    double*** yoth;
    double*** zoth;
    bool*** zero_oth;
    other.get_3dfield_h(xoth, yoth, zoth);
    other.get_3dzero(zero_oth);
    spinx = deep_copy_3darr<double>(totsize, totsize, totsize, xoth);
    spiny = deep_copy_3darr<double>(totsize, totsize, totsize, yoth);
    spinz = deep_copy_3darr<double>(totsize, totsize, totsize, zoth);
    iszero = deep_copy_3darr<bool>(totsize, totsize, totsize, zero_oth);
    postemp = alloc_2darr<int>(3, 1358);
}

field_3d_h::field_3d_h(const field_3d_h& other)
{
    ft = 3;
    dim = 3;
    insize = other.get_insize();
    totsize = other.get_totsize();
    periodic = other.get_perio();
    double*** xoth;
    double*** yoth;
    double*** zoth;
    bool*** zero_oth;
    other.get_3dfield_h(xoth, yoth, zoth);
    other.get_3dzero(zero_oth);
    spinx = deep_copy_3darr<double>(totsize, totsize, totsize, xoth);
    spiny = deep_copy_3darr<double>(totsize, totsize, totsize, yoth);
    spinz = deep_copy_3darr<double>(totsize, totsize, totsize, zoth);
    iszero = deep_copy_3darr<bool>(totsize, totsize, totsize, zero_oth);
    postemp = alloc_2darr<int>(3, 1358);
}

field_3d_h& field_3d_h::operator=(const field_3d_h& other)
{
    ft = 3;
    dim = 3;
    insize = other.get_insize();
    totsize = other.get_totsize();
    periodic = other.get_perio();
    double*** xoth;
    double*** yoth;
    double*** zoth;
    bool*** zero_oth;
    other.get_3dfield_h(xoth, yoth, zoth);
    other.get_3dzero(zero_oth);
    spinx = deep_copy_3darr<double>(totsize, totsize, totsize, xoth);
    spiny = deep_copy_3darr<double>(totsize, totsize, totsize, yoth);
    spinz = deep_copy_3darr<double>(totsize, totsize, totsize, zoth);
    iszero = deep_copy_3darr<bool>(totsize, totsize, totsize, zero_oth);
    postemp = alloc_2darr<int>(3, 1358);
    return *this;
}

field_3d_h::~field_3d_h()
{
    dealloc_3darr<double>(totsize, totsize, spinx);
    dealloc_3darr<double>(totsize, totsize, spiny);
    dealloc_3darr<double>(totsize, totsize, spinz);
    dealloc_3darr<bool>(totsize, totsize, iszero);
    dealloc_2darr<int>(3, postemp);
}

void field_3d_h::h_access(vector<int>& position, vector<double>& out)
{
    out[0] = spinx[position[0]][position[1]][position[2]];
    out[1] = spiny[position[0]][position[1]][position[2]];
    out[2] = spinz[position[0]][position[1]][position[2]];
}

void field_3d_h::h_next(bool &finish, vector<int> &pos, vector<double> &out)
{
    this->h_access(pos, out);
    this->next(finish, pos);
}

void field_3d_h::fill_ghost(int num_rows)
{
    for(int r = 0; r < num_rows; r++)
    {
        for(int i=0; i < totsize; i++)
        {
            for(int j=0; j < totsize; j++)
            {
                spinx[i][j][r] = 0;
                spiny[i][j][r] = 0;
                spinz[i][j][r] = 0;
                iszero[i][j][r] = true;
                spinx[r][i][j] = 0;
                spiny[r][i][j] = 0;
                spinz[r][i][j] = 0;
                iszero[r][i][j] = true;
                spinx[j][r][i] = 0;
                spiny[j][r][i] = 0;
                spinz[j][r][i] = 0;
                iszero[j][r][i] = true;

                spinx[i][j][totsize-1-r] = 0;
                spiny[i][j][totsize-1-r] = 0;
                spinz[i][j][totsize-1-r] = 0;
                iszero[i][j][totsize-1-r] = true;
                spinx[totsize-1-r][i][j] = 0;
                spiny[totsize-1-r][i][j] = 0;
                spinz[totsize-1-r][i][j] = 0;
                iszero[totsize-1-r][i][j] = true;
                spinx[j][totsize-1-r][i] = 0;
                spiny[j][totsize-1-r][i] = 0;
                spinz[j][totsize-1-r][i] = 0;
                iszero[j][totsize-1-r][i] = true;
            }
        }
    }
}

void field_3d_h::get_3dfield_h(double*** &x, double*** &y, double*** &z) const
{
    x = spinx;
    y = spiny;
    z = spinz;
}

void field_3d_h::h_adjacent(vector<int>& position, double** out)
{
    dirsx[0] = boundmovedown(position[0] - 1, totsize);
    dirsx[1] = boundmoveup(position[0] + 1, totsize);
    dirsx[2] = position[0];
    dirsx[3] = position[0];
    dirsx[4] = position[0];
    dirsx[5] = position[0];

    dirsy[0] = position[1];
    dirsy[1] = position[1];
    dirsy[2] = boundmovedown(position[1] - 1, totsize);
    dirsy[3] = boundmoveup(position[1] + 1, totsize);
    dirsy[4] = position[1];
    dirsy[5] = position[1];

    dirsz[0] = position[2];
    dirsz[1] = position[2];
    dirsz[2] = position[2];
    dirsz[3] = position[2];
    dirsz[4] = boundmovedown(position[2] - 1, totsize);
    dirsz[5] = boundmoveup(position[2] + 1, totsize);

    for (int i = 0; i < 6; i++)
    {
        out[0][i] = spinx[dirsx[i]][dirsy[i]][dirsz[i]];
        out[1][i] = spiny[dirsx[i]][dirsy[i]][dirsz[i]];
        out[2][i] = spinz[dirsx[i]][dirsy[i]][dirsz[i]];
    }
}

void field_3d_h::h_2adjacent(vector<int>& position, double** out)
{
    int dirsx[6], dirsy[6], dirsz[6];
    dirsx[0] = boundmovedown(position[0] - 2, totsize);
    dirsx[1] = boundmoveup(position[0] + 2, totsize);
    dirsx[2] = position[0];
    dirsx[3] = position[0];
    dirsx[4] = position[0];
    dirsx[5] = position[0];

    dirsy[0] = position[1];
    dirsy[1] = position[1];
    dirsy[2] = boundmovedown(position[1] - 2, totsize);
    dirsy[3] = boundmoveup(position[1] + 2, totsize);
    dirsy[4] = position[1];
    dirsy[5] = position[1];

    dirsz[0] = position[2];
    dirsz[1] = position[2];
    dirsz[2] = position[2];
    dirsz[3] = position[2];
    dirsz[4] = boundmovedown(position[2] - 2, totsize);
    dirsz[5] = boundmoveup(position[2] + 2, totsize);

    for (int i = 0; i < 6; i++)
    {
        out[0][i] = spinx[dirsx[i]][dirsy[i]][dirsz[i]];
        out[1][i] = spiny[dirsx[i]][dirsy[i]][dirsz[i]];
        out[2][i] = spinz[dirsx[i]][dirsy[i]][dirsz[i]];
    }
}

void field_3d_h::h_arb_adj(vector<int>& position, vector<int>& dxs, vector<int>& dys, vector<int>& dzs, double** out, int num)
{
    #pragma simd
    for(int i = 0; i < num; i++)
    {
        postemp[0][i] = min(max(0, (position[0] + dxs[i])), totsize-1);
    }
    #pragma simd
    for(int i = 0; i < num; i++)
    {
        postemp[1][i] = min(max(0, (position[1] + dys[i])), totsize-1);
    }
    #pragma simd
    for(int i = 0; i < num; i++)
    {
        postemp[2][i] = min(max(0, (position[2] + dzs[i])), totsize-1);
    }

    // #pragma simd
    for(int i = 0; i < num; i++)
    {
        out[0][i] = spinx[postemp[0][i]][postemp[1][i]][postemp[2][i]];
        out[1][i] = spiny[postemp[0][i]][postemp[1][i]][postemp[2][i]];
        out[2][i] = spinz[postemp[0][i]][postemp[1][i]][postemp[2][i]];
    }
}

void field_3d_h::fill_rand(vector<int>& position)
{
    rand_spin_h(spinx[position[0]][position[1]][position[2]],
                spiny[position[0]][position[1]][position[2]],
                spinz[position[0]][position[1]][position[2]]);
    iszero[position[0]][position[1]][position[2]] = false;
}

void field_3d_h::fill_one(vector<int>&position)
{
    spinx[position[0]][position[1]][position[2]] = 0;
    spiny[position[0]][position[1]][position[2]] = 0;
    spinz[position[0]][position[1]][position[2]] = 1;
    iszero[position[0]][position[1]][position[2]] = false;
}

void field_3d_h::fill_zero(vector<int>& position)
{
    spinx[position[0]][position[1]][position[2]] = 0;
    spiny[position[0]][position[1]][position[2]] = 0;
    spinz[position[0]][position[1]][position[2]] = 0;
    iszero[position[0]][position[1]][position[2]] = true;
}

void field_3d_h::fill_val_h(vector<int>& position, double x, double y, double z)
{
    spinx[position[0]][position[1]][position[2]] = x;
    spiny[position[0]][position[1]][position[2]] = y;
    spinz[position[0]][position[1]][position[2]] = z;
}

void field_3d_h::change_to_test(vector<int>& position, ham_type* hamil)
{
    hamil->get_test(spinx[position[0]][position[1]][position[2]],
                    spiny[position[0]][position[1]][position[2]],
                    spinz[position[0]][position[1]][position[2]]);
}

void field_3d_h::add_val_h(vector<int>& position, vector<double> &in)
{
    spinx[position[0]][position[1]][position[2]] += in[0];
    spiny[position[0]][position[1]][position[2]] += in[1];
    spinz[position[0]][position[1]][position[2]] += in[2];
}

void field_3d_h::print(string filename, string arrname)
{
    // Open existing file
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    hid_t f_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, plist_id);
    H5Pclose(plist_id);

    // Open dataset
    hid_t dset_id = H5Dopen1(f_id, arrname.c_str());

    // Get slab and slice space
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    hsize_t count[4] = {totsize, totsize, totsize, 1};
    hsize_t offset[4] = {0, 0, 0, 0};
    hid_t slice_space_id = H5Screate_simple(4, count, NULL);
    // for x
    hid_t slab_id = H5Dget_space(dset_id);
    H5Sselect_hyperslab(slab_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, slice_space_id, slab_id, plist_id,
        spinx[0][0]);
    H5Sclose(slab_id);
    // for y
    offset[3] = 1;
    slab_id = H5Dget_space(dset_id);
    H5Sselect_hyperslab(slab_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, slice_space_id, slab_id, plist_id,
        spiny[0][0]);
    H5Sclose(slab_id);
    // for z
    offset[3] = 2;
    slab_id = H5Dget_space(dset_id);
    H5Sselect_hyperslab(slab_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, slice_space_id, slab_id, plist_id,
        spinz[0][0]);
    H5Sclose(slab_id);

    // close
    H5Pclose(plist_id);
    H5Dclose(dset_id);
    H5Fclose(f_id);
}

void field_3d_h::print_setup(const string filename, const int Tmax,
    const int Hmax)
{
    // Open existing file
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    hid_t f_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, plist_id);
    H5Pclose(plist_id);

    // create datasets
    hsize_t full_dims[4] = {totsize, totsize, totsize, 3};
    hid_t dspace_id = H5Screate_simple(4, full_dims, NULL);
    hid_t dset_id;
    stringstream nstream;
    string name;
    for(int i=0; i < Tmax; i++)
    {
        for(int j=0; j < Hmax; j++)
        {
            nstream << "/Latt_Print/T_" << i << "-H_" << j;
            nstream >> name;
            nstream.clear();
            dset_id = H5Dcreate(f_id, name.c_str(), H5T_NATIVE_DOUBLE,
                dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dclose(dset_id);
        }
    }

    // close
    H5Fclose(f_id);
}

///////////////////////
// 3d Ising-model
///////////////////////

field_3d_i::field_3d_i()
{
    dim = 3;
    periodic = true;
    ft = 31;
}

field_3d_i::field_3d_i(int size, bool isperio)
{
    dim = 3;
    periodic = isperio;
    insize = size;
    ft = 31;
    if(periodic)
    {
        totsize = insize;
    }
    else
    {
        totsize = insize + 2;
    }
    spin = alloc_3darr<int>(totsize, totsize, totsize);
    iszero = alloc_3darr<bool>(totsize, totsize, totsize);
}

field_3d_i::field_3d_i(field_type& other)
{
    ft = 31;
    if(other.get_ft() != 31)
    {
        cout << "Cannot copy from other field type" << endl;
        exit(104);
    }
    dim = 3;
    insize = other.get_insize();
    totsize = other.get_totsize();
    periodic = other.get_perio();
    int*** othspin;
    bool*** zero_oth;
    other.get_3dfield_i(othspin);
    other.get_3dzero(zero_oth);
    spin = deep_copy_3darr<int>(totsize, totsize, totsize, othspin);
    iszero = deep_copy_3darr<bool>(totsize, totsize, totsize, zero_oth);
}

field_3d_i::field_3d_i(const field_3d_i& other)
{
    ft = 31;
    dim = 3;
    insize = other.get_insize();
    totsize = other.get_totsize();
    periodic = other.get_perio();
    int*** othspin;
    bool*** zero_oth;
    other.get_3dfield_i(othspin);
    other.get_3dzero(zero_oth);
    spin = deep_copy_3darr<int>(totsize, totsize, totsize, othspin);
    iszero = deep_copy_3darr<bool>(totsize, totsize, totsize, zero_oth);
}

field_3d_i& field_3d_i::operator=(const field_3d_i& other)
{
    ft = 31;
    dim = 3;
    insize = other.get_insize();
    totsize = other.get_totsize();
    periodic = other.get_perio();
    int*** othspin;
    bool*** zero_oth;
    other.get_3dfield_i(othspin);
    other.get_3dzero(zero_oth);
    spin = deep_copy_3darr<int>(totsize, totsize, totsize, othspin);
    iszero = deep_copy_3darr<bool>(totsize, totsize, totsize, zero_oth);
    return *this;
}

field_3d_i::~field_3d_i()
{
    dealloc_3darr<int>(totsize, totsize, spin);
    dealloc_3darr<bool>(totsize, totsize, iszero);
}

void field_3d_i::i_access(vector<int>& position, int &out)
{
    out = spin[position[0]][position[1]][position[2]];
}

void field_3d_i::i_next(bool &finish, vector<int> &pos, int &out)
{
    this->i_access(pos, out);
    this->next(finish, pos);
}

void field_3d_i::fill_ghost(int num_rows)
{
    for(int r = 0; r < num_rows; r++)
    {
        for(int i=0; i < totsize; i++)
        {
            for(int j=0; j < totsize; j++)
            {
                spin[i][j][r] = 0;
                iszero[i][j][r] = true;
                spin[r][i][j] = 0;
                iszero[r][i][j] = true;
                spin[j][r][i] = 0;
                iszero[j][r][i] = true;

                spin[i][j][totsize-1-r] = 0;
                iszero[i][j][totsize-1-r] = true;
                spin[totsize-1-r][i][j] = 0;
                iszero[totsize-1-r][i][j] = true;
                spin[j][totsize-1-r][i] = 0;
                iszero[j][totsize-1-r][i] = true;
            }
        }
    }
}

void field_3d_i::get_3dfield_i(int*** &x) const
{
    x = spin;
}

void field_3d_i::i_adjacent(vector<int>& position, int* out)
{
    int dirsx[6], dirsy[6], dirsz[6];
    dirsx[0] = boundmovedown(position[0] - 1, totsize);
    dirsx[1] = boundmoveup(position[0] + 1, totsize);
    dirsx[2] = position[0];
    dirsx[3] = position[0];
    dirsx[4] = position[0];
    dirsx[5] = position[0];

    dirsy[0] = position[1];
    dirsy[1] = position[1];
    dirsy[2] = boundmovedown(position[1] - 1, totsize);
    dirsy[3] = boundmoveup(position[1] + 1, totsize);
    dirsy[4] = position[1];
    dirsy[5] = position[1];

    dirsz[0] = position[2];
    dirsz[1] = position[2];
    dirsz[2] = position[2];
    dirsz[3] = position[2];
    dirsz[4] = boundmovedown(position[2] - 1, totsize);
    dirsz[5] = boundmoveup(position[2] + 1, totsize);

    for (int i = 0; i < 6; i++)
    {
        out[i] = spin[dirsx[i]][dirsy[i]][dirsz[i]];
    }
}

void field_3d_i::fill_rand(vector<int>& position)
{
    rand_spin_i(spin[position[0]][position[1]][position[2]]);
    iszero[position[0]][position[1]][position[2]] = false;
}

void field_3d_i::fill_one(vector<int>&position)
{
    spin[position[0]][position[1]][position[2]] = 1;
    iszero[position[0]][position[1]][position[2]] = false;
}

void field_3d_i::fill_zero(vector<int>& position)
{
    spin[position[0]][position[1]][position[2]] = 0;
    iszero[position[0]][position[1]][position[2]] = true;
}

void field_3d_i::fill_val_i(vector<int>& position, int val)
{
    spin[position[0]][position[1]][position[2]] = val;
}

void field_3d_i::change_to_test(vector<int>& position, ham_type* hamil)
{
    spin[position[0]][position[1]][position[2]] = -spin[position[0]][position[1]][position[2]];
}