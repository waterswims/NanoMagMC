#ifndef _FIELD
#define _FIELD

#define BOOST_DISABLE_ASSERTS

#include "boost/multi_array.hpp"
#include <vector>
#include <string>

using namespace std;

class field_type_new
{
protected:
    int dim, insize, totsize;
    bool periodic;
public:
    field_type_new(){}
    ~field_type_new(){}
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
}

template <class T> class field_type
{
protected:
    int dim, insize, totsize;
    bool periodic;
public:
    field_type(){}
    ~field_type(){}
    virtual T& access(vector<int>& position){}
    virtual void adjacent(vector<int>& position, vector<T*>& out){}
    virtual T& next(bool &finish, vector<int> &pos){}
    int get_insize(){return insize;}
    int get_totsize(){return totsize;}
    bool get_perio(){return periodic;}
    int get_num(){return pow(insize, dim);}
    int get_dim(){return dim;}
    virtual void fill_ghost(int &num){}
    virtual void new_mem(){}
    virtual int findnum(){return 0;}
    virtual vector<T>* get_1dfield(){return NULL;}
    virtual boost::multi_array<T, 2>* get_2dfield(){return NULL;}
    virtual boost::multi_array<T, 3>* get_3dfield(){return NULL;}
    virtual void print(){}
};

template <class T> class field_2d: public field_type<T>
{
protected:
    boost::multi_array<T, 2>* field;
public:
    field_2d();
    field_2d(int size, bool isperio);
    field_2d(field_type<T>& other);
    field_2d(const field_2d<T>& other);
    ~field_2d(){delete field;}
    T& access(vector<int>& position);
    virtual void adjacent(vector<int>& position, vector<T*>& out);
    T& next(bool &finish, vector<int> &pos);
    void fill_ghost(int &num);
    void new_mem();
    field_2d<T>& operator=(field_2d<T>& other);
    int findnum();
    boost::multi_array<T, 2>* get_2dfield(){return field;}
    void print();
};

template <class T> class field_3d: public field_type<T>
{
protected:
    boost::multi_array<T, 3>* field;
public:
    field_3d();
    field_3d(int size, bool isperio);
    field_3d(int size, bool isperio, int pad);
    field_3d(field_type<T>& other);
    field_3d(const field_3d<T>& other);
    ~field_3d(){delete field;}
    T& access(vector<int>& position);
    virtual void adjacent(vector<int>& position, vector<T*>& out);
    T& next(bool &finish, vector<int> &pos);
    void fill_ghost(int &num);
    void new_mem();
    field_3d<T>& operator=(field_3d<T>& other);
    int findnum();
    boost::multi_array<T, 3>* get_3dfield(){return field;}
};

template <class T> class hex_2d: public field_2d<T>
{
public:
    hex_2d():field_2d<T>(){}
    hex_2d(int size, bool isperio):field_2d<T>(size, isperio){}
    hex_2d(field_type<T>& other):field_2d<T>(other){}
    hex_2d(const hex_2d<T>& other);
    ~hex_2d(){delete this->field;}
    void adjacent(vector<int> &position, vector<T*> &out);
    hex_2d<T>& operator=(hex_2d<T>& other);
    boost::multi_array<T,2>* get_2dfield(){return this->field;}
};

template <class T> class hex_3d: public field_3d<T>
{
public:
    hex_3d():field_3d<T>(){}
    hex_3d(int size, bool isperio):field_3d<T>(size, isperio){}
    hex_3d(field_type<T>& other):field_3d<T>(other){}
    hex_3d(const hex_3d<T>& other);
    ~hex_3d(){delete this->field;}
    void adjacent(vector<int> &position, vector<T*> &out);
    hex_3d<T>& operator=(hex_3d<T>& other);
    boost::multi_array<T,3>* get_3dfield(){return this->field;}
};

template <class T> class field_cluster: public field_type<T>
{
protected:
    vector<double>* xs;
    vector<double>* temp;
    vector<T>* field;
public:
    field_cluster();
    field_cluster(string filename);
    field_cluster(field_type<T>& other);
    field_cluster(const field_cluster<T>& other);
    ~field_cluster(){delete xs; delete temp; delete field;}
    T& access(vector<int>& position){return (*field)[position[0]];}
    T& next(bool &finish, vector<int> &pos);
    field_cluster<T>& operator=(field_cluster<T>& other);
    int findnum(){return (*field).size();}
    vector<T>* get_1dfield(){return field;}
};

#endif
