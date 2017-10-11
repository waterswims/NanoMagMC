#ifndef _FIELD
#define _FIELD

#include "hamiltonian.hpp"
#include "array_alloc.hpp"
#include <vector>
#include <string>

class ham_type;

using namespace std;

void rand_spin_h(double &x, double &y, double &z);

class field_type
{
protected:
    int ft;
    int dim, insize, totsize;
    bool periodic;
public:
    field_type(){}
    ~field_type(){}
    // virtual void alloc_pos(int size){}
    virtual void i_access(vector<int>& postion, int &out){}
    virtual void i_adjacent(vector<int>& position, int* out){}
    virtual void h_access(vector<int>& position, vector<double>& out){}
    virtual void h_adjacent(vector<int>& position, double** out){}
    virtual void h_2adjacent(vector<int>& position, double** out){}
    virtual void h_arb_adj(vector<int>& position, vector<int>& dxs, vector<int>& dys, vector<int>& dzs, double** out, int num){}
    virtual void next(bool &finish, vector<int> &pos){}
    virtual void i_next(bool &finish, vector<int> &pos, int &out){}
    virtual void h_next(bool &finish, vector<int> &pos, vector<double> &out){}
    int get_insize() const {return insize;}
    int get_totsize() const {return totsize;}
    bool get_perio() const {return periodic;}
    int get_dim() const {return dim;}
    virtual void fill_ghost(int num_rows){}
    virtual void fill_rand(vector<int>& position){}
    virtual void fill_zero(vector<int>& position){}
    virtual void fill_one(vector<int>&position){}
    virtual void fill_val_i(vector<int>& position, int val){}
    virtual void fill_val_h(vector<int>& position, double x, double y, double z){}
    virtual void add_val_h(vector<int>& position, vector<double> &in){}
    void allzero();
    virtual bool check_zero(vector<int>& position){return true;}
    virtual void change_to_test(vector<int>& position, ham_type* hamil){}
    virtual void new_mem(){}
    virtual int findnum(){return 0;}
    virtual void get_2dfield_i(int** &x) const{}
    virtual void get_3dfield_i(int*** &x) const{}
    virtual void get_1dfield_h(double* &x, double* &y, double* &z) const{}
    virtual void get_2dfield_h(double** &x, double** &y, double** &z) const{}
    virtual void get_3dfield_h(double*** &x, double*** &y, double*** &z) const{}
    virtual void get_2dzero(bool** &x) const{}
    virtual void get_3dzero(bool*** &x) const{}
    virtual void print(string filename){}
    int get_ft() const {return ft;}
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
    ~field_cluster_h();
    void h_access(vector<int>& position, vector<double>& out);
    void h_next(bool &finish, vector<int> &pos, vector<double> &out);
    field_cluster_h& operator=(const field_cluster_h& other);
    int findnum(){return insize;}
    void get_1dfield_h(double* &x, double* &y, double* &z) const;
    void fill_rand(vector<int>& position){}
    bool check_zero(vector<int>& position){return false;}
    void change_to_test(vector<int>& position, ham_type* hamil);
};

class field_2d: public field_type
{
protected:
    bool** iszero;
public:
    field_2d() {}
    ~field_2d() {}
    void next(bool &finish, vector<int> &pos);
    bool check_zero(vector<int>& position)
        {return iszero[position[0]][position[1]];}
    void get_2dzero(bool** &x) const {x = iszero;}
    int findnum();
};

class field_2d_h: public field_2d
{
protected:
    double** spinx;
    double** spiny;
    double** spinz;
public:
    field_2d_h();
    field_2d_h(int size, bool isperio);
    field_2d_h(int size, bool isperio, int p_pad);
    field_2d_h(field_type& other);
    field_2d_h(const field_2d_h& other);
    ~field_2d_h();
    void h_access(vector<int>& position, vector<double>& out);
    void h_next(bool &finish, vector<int> &pos, vector<double> &out);
    void fill_ghost(int num_rows);
    field_2d_h& operator=(const field_2d_h& other);
    void get_2dfield_h(double** &x, double** &y, double** &z) const;
    void h_adjacent(vector<int>& position, double** out);
    void fill_rand(vector<int>& position);
    void fill_one(vector<int>&position);
    void fill_zero(vector<int>& position);
    void fill_val_h(vector<int>& position, double x, double y, double z);
    void add_val_h(vector<int>& position, vector<double> &in);
    void change_to_test(vector<int>& position, ham_type* hamil);
};

class field_2d_i: public field_2d
{
protected:
    int** spin;
public:
    field_2d_i();
    field_2d_i(int size, bool isperio);
    field_2d_i(field_type& other);
    field_2d_i(const field_2d_i& other);
    ~field_2d_i();
    void i_access(vector<int>& position, int &out);
    void i_next(bool &finish, vector<int> &pos, int &out);
    void fill_ghost(int num_rows);
    field_2d_i& operator=(const field_2d_i& other);
    void get_2dfield_i(int** &x) const;
    void i_adjacent(vector<int>& position, int* out);
    void fill_rand(vector<int>& position);
    void fill_one(vector<int>&position);
    void fill_zero(vector<int>& position);
    void change_to_test(vector<int>& position, ham_type* hamil);
    void fill_val_i(vector<int>& position, int val);
};

class field_3d: public field_type
{
protected:
    bool*** iszero;
public:
    field_3d() {}
    ~field_3d() {}
    void next(bool &finish, vector<int> &pos);
    bool check_zero(vector<int>& position)
        {return iszero[position[0]][position[1]][position[2]];}
    void get_3dzero(bool*** &x) const {x = iszero;}
    int findnum();
};

class field_3d_h: public field_3d
{
protected:
    double*** spinx;
    double*** spiny;
    double*** spinz;
    int** postemp;
public:
    field_3d_h();
    field_3d_h(int size, bool isperio);
    field_3d_h(int size, bool isperio, int p_pad);
    field_3d_h(field_type& other);
    field_3d_h(const field_3d_h& other);
    ~field_3d_h();
    // void alloc_pos(int size) {postemp = alloc_2darr<int>(3, size);}
    void h_access(vector<int>& position, vector<double>& out);
    void h_next(bool &finish, vector<int> &pos, vector<double> &out);
    void fill_ghost(int num_rows);
    field_3d_h& operator=(const field_3d_h& other);
    void get_3dfield_h(double*** &x, double*** &y, double*** &z) const;
    void h_adjacent(vector<int>& position, double** out);
    void h_2adjacent(vector<int>& position, double** out);
    void fill_rand(vector<int>& position);
    void fill_one(vector<int>&position);
    void fill_zero(vector<int>& position);
    void fill_val_h(vector<int>& position, double x, double y, double z);
    void add_val_h(vector<int>& position, vector<double> &in);
    void change_to_test(vector<int>& position, ham_type* hamil);
    void h_arb_adj(vector<int>& position, vector<int>& dxs, vector<int>& dys, vector<int>& dzs, double** out, int num);
    void print(string filename);
};

class field_3d_i: public field_3d
{
protected:
    int*** spin;
public:
    field_3d_i();
    field_3d_i(int size, bool isperio);
    field_3d_i(field_type& other);
    field_3d_i(const field_3d_i& other);
    ~field_3d_i();
    void i_access(vector<int>& position, int &out);
    void i_next(bool &finish, vector<int> &pos, int &out);
    void fill_ghost(int num_rows);
    field_3d_i& operator=(const field_3d_i& other);
    void get_3dfield_i(int*** &x) const;
    void i_adjacent(vector<int>& position, int* out);
    void fill_rand(vector<int>& position);
    void fill_one(vector<int>&position);
    void fill_zero(vector<int>& position);
    void change_to_test(vector<int>& position, ham_type* hamil);
    void fill_val_i(vector<int>& position, int val);
};

#endif
