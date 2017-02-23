#ifndef _HAMIL
#define _HAMIL

#include "field_type_new.hpp"
#include <vector>

class field_type;

using namespace std;

class ham_type
{
public:
    virtual double calc_E(field_type* lattice) {return 0;}
    virtual double dE(field_type* lattice, vector<int>& position) {return 0;}
    virtual vector<double> calc_M(field_type* lattice) {}
    virtual vector<double> calc_subM(field_type* lattice, int subnumber) {}
    virtual double get_J(){return 0;}
    virtual double get_H(){return 0;}
    virtual vector<double> get_Js(){}
    virtual vector<double> get_Hs(){}
    virtual void get_test(double &x, double &y, double &z){}
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
    ~ham_ising(){}
    double calc_E(field_type* lattice);
    vector<double> calc_M(field_type* lattice);
    vector<double> calc_subM(field_type* lattice, int subnumber);
    double dE(field_type* lattice, vector<int>& position);
    double get_J(){return J;}
    double get_H(){return H;}
    ham_ising& operator=(ham_type& other);
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
    ~ham_heis(){}
    virtual double calc_E(field_type* lattice);
    vector<double> calc_M(field_type* lattice);
    vector<double> calc_subM(field_type* lattice, int subnumber);
    virtual double dE(field_type* lattice, vector<int>& position);
    vector<double> get_Js(){return J;}
    vector<double> get_Hs(){return H;}
    ham_heis& operator=(ham_type& other);
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
    ham_FePt(ham_type& other);
    ~ham_FePt() {}
    double calc_E(field_type* lattice);
    double dE(field_type* lattice, vector<int>& position);
    ham_FePt& operator=(ham_type& other);
    void read_Js();
    bool check_pos(int start, int fin);
};

#endif
