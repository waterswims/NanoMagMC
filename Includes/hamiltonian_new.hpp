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

#endif
