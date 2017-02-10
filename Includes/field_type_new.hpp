#ifndef _FIELD
#define _FIELD

#include <vector>
#include <string>

using namespace std;

class field_type
{
protected:
    int ft;
    int dim, insize, totsize;
    bool periodic;
public:
    field_type(){}
    ~field_type(){}
    virtual void i_access(vector<int>& postion, int &out){}
    virtual void i_adjacent(vector<int>& position, vector<int>& out){}
    virtual void h_access(vector<int>& postion, vector<double>& out){}
    virtual void h_adjacent(vector<int>& position, vector<vector<double> >& out){}
    virtual void i_next(bool &finish, vector<int> &pos, int &out){}
    virtual void h_next(bool &finish, vector<int> &pos, vector<double> &out){}
    int get_insize(){return insize;}
    int get_totsize(){return totsize;}
    bool get_perio(){return periodic;}
    int get_dim(){return dim;}
    virtual void fill_ghost(int num){}
    virtual void new_mem(){}
    virtual int findnum(){return 0;}
    virtual void get_2dfield_i(int** &x){}
    virtual void get_3dfield_i(int*** &x){}
    virtual void get_1dfield_h(double* &x, double* &y, double* &z){}
    virtual void get_2dfield_h(double** &x, double** &y, double** &z){}
    virtual void get_3dfield_h(double*** &x, double*** &y, double*** &z){}
    virtual void print(){}
    void rand_spin_h(double &x, double &y, double &z);
};

class field_cluster_h: public field_type
{
protected:
    double* spinx;
    double* spiny;
    double* spinz;
public:
    field_cluster_h();
    field_cluster_h(string filename);
    field_cluster_h(field_type& other);
    field_cluster_h(const field_cluster_h& other);
    ~field_cluster();
    void h_access(vector<int>& postion, vector<double>& out);
    void h_next(bool &finish, vector<int> &pos, vector<double> &out);
    field_cluster& operator=(field_cluster& other);
    int findnum(){return insize;}
    void get_1dfield_h(double* &x, double* &y, double* &z);
};

#endif
