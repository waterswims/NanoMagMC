#include "field_type_new.hpp"
#include <fstream>
#include <sstream>

extern mkl_irand st_rand_int;
extern mkl_drand st_rand_double;

int boundmovedown(int test, int limit)
{
    return ((test%limit)+limit)%limit;
}

int boundmoveup(int test, int limit)
{
    return test%limit;
}

void field_type::rand_spin_h(double &x, double &y, double &z)
{
    double phi = st_rand_double.gen()*2*pi;
    double cthet = 2*st_rand_double.gen()-1;
    double sthet = pow(1 - pow(cthet, 2), 0.5);
    x = cos(phi)*sthet;
    y = sin(phi)*sthet;
    z = cthet;
}

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
    while(getline(file, line))
    {
        insize++;
    }
    file.close();

    string line;
    double temp_d;
    totsize=this->insize;
    spinx = alloc_1darr(insize);
    spiny = alloc_1darr(insize);
    spinz = alloc_1darr(insize);

    for(int i=0; i < insize; i++)
    {
        this->rand_spin_h(spinx[i], spiny[i], spinz[i]);
    }
}

field_cluster_h::field_cluster_h(field_type& other)
{
    dim = 1;
    periodic = 1;
    ft = 1;
    if(other.ft != 1)
    {
        cout << "Cannot copy from other field type" << endl;
        exit(104);
    }
    insize = other.get_insize();
    totsize = other.get_totsize();
    double* xoth, yoth, zoth;
    other.get_1dfield_h(xoth, yoth, zoth);
    spinx = deep_copy_1darr(insize, xoth);
    spiny = deep_copy_1darr(insize, yoth);
    spinz = deep_copy_1darr(insize, zoth);
}

field_cluster_h::field_cluster_h(const field_cluster_h& other)
{
    dim = 1;
    periodic = 1;
    ft = 1;
    insize = other.get_insize();
    totsize = other.get_totsize();
    double* xoth, yoth, zoth;
    other.get_1dfield_h(xoth, yoth, zoth);
    spinx = deep_copy_1darr(insize, xoth);
    spiny = deep_copy_1darr(insize, yoth);
    spinz = deep_copy_1darr(insize, zoth);
}

field_cluster& field_cluster_h::operator=(field_cluster& other)
{
    dim = 1;
    periodic = 1;
    ft = 1;
    insize = other.get_insize();
    totsize = other.get_totsize();
    double* xoth, yoth, zoth;
    other.get_1dfield_h(xoth, yoth, zoth);
    spinx = deep_copy_1darr(insize, xoth);
    spiny = deep_copy_1darr(insize, yoth);
    spinz = deep_copy_1darr(insize, zoth);
    return *this;
}

field_cluster_h::~field_cluster()
{
    dealloc_1darr(insize, spinx);
    dealloc_1darr(insize, spiny);
    dealloc_1darr(insize, spinz);
}

void field_cluster_h::h_access(vector<int>& postion, vector<double>& out)
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

void field_cluster_h::get_1dfield_h(double* &x, double* &y, double* &z)
{
    x = spinx;
    y = spiny;
    z = spinz;
}

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
    spinx = alloc_2darr(totsize, totsize);
    spiny = alloc_2darr(totsize, totsize);
    spinz = alloc_2darr(totsize, totsize);
}

field_2d_h::field_2d_h(field_type& other)
{
    ft = 2;
    if(other.ft != 2)
    {
        cout << "Cannot copy from other field type" << endl;
        exit(104);
    }
    dim = 2;
    insize = other.get_insize();
    totsize = other.get_totsize();
    periodic = other.get_perio();
    double** xoth, yoth, zoth;
    other.get_2dfield_h(xoth, yoth, zoth);
    spinx = deep_copy_2darr(insize, xoth);
    spiny = deep_copy_2darr(insize, yoth);
    spinz = deep_copy_2darr(insize, zoth);
}

field_2d_h::field_2d_h(const field_2d_h& other)
{
    ft = 2;
    dim = 2;
    insize = other.get_insize();
    totsize = other.get_totsize();
    periodic = other.get_perio();
    double** xoth, yoth, zoth;
    other.get_2dfield_h(xoth, yoth, zoth);
    spinx = deep_copy_2darr(insize, xoth);
    spiny = deep_copy_2darr(insize, yoth);
    spinz = deep_copy_2darr(insize, zoth);
}

field_2d_h& field_2d_h::operator=(field_2d_h& other)
{
    ft = 2;
    dim = 2;
    insize = other.get_insize();
    totsize = other.get_totsize();
    periodic = other.get_perio();
    double** xoth, yoth, zoth;
    other.get_2dfield_h(xoth, yoth, zoth);
    spinx = deep_copy_2darr(insize, xoth);
    spiny = deep_copy_2darr(insize, yoth);
    spinz = deep_copy_2darr(insize, zoth);
    return *this;
}

field_2d_h::~field_2d_h()
{
    dealloc_2darr(totsize, totsize, spinx);
    dealloc_2darr(totsize, totsize, spiny);
    dealloc_2darr(totsize, totsize, spinz);
}

void field_2d_h::h_access(vector<int>& position, vector<double>& out)
{
    out[0] = spinx[position[0]][position[1]];
    out[1] = spiny[position[0]][position[1]];
    out[2] = spinz[position[0]][position[1]];
}

void field_2d_h::h_next(bool &finish, vector<int> &pos, vector<double> &out)
{
    int start = 0;
    int end = totsize;
    if(!periodic)
    {
        start++;
        end--;
    }
    this->h_access(pos, out);
    pos[1]++;
    if(pos[1] == end)
    {
        pos[1] = start;
        pos[0]++;
        if(pos[0] == end;)
        {
            pos[0] = start;
            finish = true;
        }
    }
    return out;
}

void field_2d_h::fill_ghost()
{
    if(!periodic)
    {
        for(int i=0; i < totsize; i++)
        {
            spinx[i][0] = 0;
            spiny[i][0] = 0;
            spinz[i][0] = 0;
            spinx[0][i] = 0;
            spiny[0][i] = 0;
            spinz[0][i] = 0;
            spinx[i][totsize-1] = 0;
            spiny[i][totsize-1] = 0;
            spinz[i][totsize-1] = 0;
            spinx[totsize-1][i] = 0;
            spiny[totsize-1][i] = 0;
            spinz[totsize-1][i] = 0;
        }
    }
}

int field_2d_h::findnum()
{
    int c = 0;
    for(int i = 0; i<totsize; i++)
    {
        for(int j = 0; j<totsize; j++)
        {
            if(spinx[i][j]!=0 || spiny[i][j]!=0 || spinz[i][j]!=0)
            {
                c++;
            }
        }
    }
    return c;
}

void field_2d_h::get_2dfield_h(double** &x, double** &y, double** &z)
{
    x = spinx;
    y = spiny;
    z = spinz;
}

void field_2d_h::h_adjacent(vector<int>& position, double** &out)
{
    out = alloc_2darr(4, 3);
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
        out[i][0] = spinx[dirsx[i]][dirsy[i]];
        out[i][1] = spiny[dirsx[i]][dirsy[i]];
        out[i][2] = spinz[dirsx[i]][dirsy[i]];
    }

    return out;
}
