#ifndef _HAMIL
#define _HAMIL

#include "field_type_new.hpp"

using namespace std;

class ham_type
{
public:
    virtual double calc_E(field_type* lattice) {return 0;}
    virtual vector<double> calc_M(field_type* lattice) {}
    virtual vector<double> calc_subM(field_type* lattice, int subnumber) {}
    double dE(field_type* lattice, vector<int>& position);
    virtual double get_J(){return 0;}
    virtual double get_H(){return 0;}
    virtual vector<double> get_Js(){}
    virtual vector<double> get_Hs(){}
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

#endif
