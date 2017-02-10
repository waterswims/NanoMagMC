#include "field_type_new.hpp"
#include <fstream>
#include <sstream>

extern mkl_irand st_rand_int;
extern mkl_drand st_rand_double;

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
    get_1dfield_h(xoth, yoth, zoth);
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
    get_1dfield_h(xoth, yoth, zoth);
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
    get_1dfield_h(xoth, yoth, zoth);
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
    out[1] = spiny[position[1]];
    out[2] = spinz[position[2]];
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
