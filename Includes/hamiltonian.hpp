#ifndef _HAMIL
#define _HAMIL

#include "field_type.hpp"
#include "array_alloc.hpp"
#include <vector>

class field_type;

using namespace std;

class ham_type
{
protected:
    int dim;
public:
    virtual double calc_E(field_type* lattice) {return 0;}
    virtual double dE(field_type* lattice, vector<int>& position) {return 0;}
    virtual vector<double> calc_M(field_type* lattice) {}
    virtual vector<double> calc_subM(field_type* lattice, int subnumber) {}
    virtual double get_J() const {return 0;}
    virtual double get_H() const {return 0;}
    virtual vector<double> get_Js() const {}
    virtual vector<double> get_Hs() const {}
    virtual double get_K() const {return 0;}
    virtual void get_test(double &x, double &y, double &z){}
    virtual void init_dim(field_type* field) {}
    // virtual double get_Dx(){return 0;}
    // virtual double get_Dy(){return 0;}
    // virtual vector<double>* get_xs(){return NULL;}
    // virtual vector<double>* get_ys(){return NULL;}
    // virtual vector<double>* get_zs(){return NULL;}
    // virtual vector<double>* get_hs(){return NULL;}
    // virtual vector<double>* get_Rms(){return NULL;}
    // virtual vector<double>* get_Rns(){return NULL;}
    // virtual vector<heis_spin>* get_ks(){return NULL;}
    // virtual double get_ani(){return 0;}
    // virtual double get_Ms(){return 0;}
};

class ham_ising: public ham_type
{
private:
    double H, J;
    int* adj;
public:
    ham_ising(){}
    ham_ising(double Hin, double Jin){H = Hin; J = Jin;}
    ham_ising(ham_type& other);
    ~ham_ising(){dealloc_1darr<int>(adj);}
    double calc_E(field_type* lattice);
    vector<double> calc_M(field_type* lattice);
    vector<double> calc_subM(field_type* lattice, int subnumber);
    double dE(field_type* lattice, vector<int>& position);
    double get_J() const {return J;}
    double get_H() const {return H;}
    ham_ising& operator=(ham_type& other);
    void init_dim(field_type* field);
};

class ham_heis: public ham_type
{
protected:
    vector<double> H, J;
    double** adj;
    vector<double> vsum, curr, H_sum, J_sum, test;
    vector<int> pos;
public:
    ham_heis(){}
    ham_heis(double Hin, double Jin);
    ham_heis(ham_type& other);
    ~ham_heis(){dealloc_2darr<double>(4, adj);}
    virtual double calc_E(field_type* lattice);
    vector<double> calc_M(field_type* lattice);
    vector<double> calc_subM(field_type* lattice, int subnumber);
    virtual double dE(field_type* lattice, vector<int>& position);
    vector<double> get_Js() const {return J;}
    vector<double> get_Hs() const {return H;}
    ham_heis& operator=(ham_type& other);
    virtual void init_dim(field_type* field);
    void get_test(double &x, double &y, double &z)
        {x = test[0]; y = test[1]; z = test[2];}
};

class ham_FePt: public ham_heis
{
private:
    vector<int> pos2, dxs, dys, dzs;
    vector<double> Js, adj_curr, d_ijs;
    double d0;
public:
    ham_FePt();
    ham_FePt(double Hin);
    ham_FePt(ham_type& other);
    ~ham_FePt() {}
    double calc_E(field_type* lattice);
    double dE(field_type* lattice, vector<int>& position);
    ham_FePt& operator=(ham_type& other);
    void read_Js();
    void init_dim(field_type* field);
};

class ham_skyrm: public ham_heis
{
private:
    double K;
    int dirs[6];
    int mod[6];
    vector<double> cmp;
public:
    ham_skyrm(){}
    ham_skyrm(double Hin, double Jin, double Kin);
    ham_skyrm(const ham_type& other);
    ~ham_skyrm(){}
    virtual double calc_E(field_type* lattice);
    virtual double dE(field_type* lattice, vector<int>& position);
    ham_skyrm& operator=(const ham_type& other);
    double get_K() const {return K;}
    void set_dirs();
};

#endif
