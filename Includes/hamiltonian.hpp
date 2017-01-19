#ifndef _HAMIL
#define _HAMIL

#include "field_type.hpp"
#include "FCC.hpp"
#include "spin_type.hpp"
#include "boost/assign.hpp"
#include <vector>
#include <iostream>
#include <string>

using namespace std;


template <class T> class ham_type
{
public:
    ham_type(){}
    ~ham_type(){}
    virtual double calc_E(field_type<T>* lattice){return 0;}
    virtual vector<double> calc_M(field_type<T>* lattice){}
    virtual vector<double> calc_subM(field_type<T>* lattice, int subnumber){}
    virtual double dE(field_type<T>* lattice, vector<int>& position){return 0;}
    virtual double get_J(){return 0;}
    virtual double get_H(){return 0;}
    virtual vector<double> get_Js(){}
    virtual vector<double> get_Hs(){}
    virtual T* get_test(){}
    virtual double get_Dx(){return 0;}
    virtual double get_Dy(){return 0;}
    virtual vector<double>* get_xs(){return NULL;}
    virtual vector<double>* get_ys(){return NULL;}
    virtual vector<double>* get_zs(){return NULL;}
    virtual vector<double>* get_Rms(){return NULL;}
    virtual vector<double>* get_Rns(){return NULL;}
    virtual vector<T>* get_ks(){return NULL;}
};

template <class T> class ham_ising: public ham_type<T>
{
private:
    double H, J;
    vector<ising_spin*> adj;
    T test;
public:
    ham_ising(){}
    ham_ising(double Hin, double Jin){H = Hin; J = Jin;}
    ham_ising(ham_type<ising_spin>& other);
    ham_ising(ham_type<heis_spin>& other){cout << "Error" << endl; exit(201);}
    ~ham_ising(){}
    double calc_E(field_type<ising_spin>* lattice);
    vector<double> calc_M(field_type<ising_spin>* lattice);
    vector<double> calc_subM(field_type<ising_spin>* lattice, int subnumber);
    double dE(field_type<ising_spin>* lattice, vector<int>& position);
    double get_J(){return J;}
    double get_H(){return H;}
    ham_ising<T>& operator=(ham_type<ising_spin>& other);
    T* get_test(){return &test;}
};

template <class T> class ham_heis: public ham_type<T>
{
protected:
    vector<double> H, J;
    vector<heis_spin*> adj;
    T test;
    vector<double> vsum, curr, adj_curr, H_sum, J_sum, potential;
    vector<int> pos;
public:
    ham_heis(){}
    ham_heis(double Hin, double Jin);
    ham_heis(ham_type<ising_spin>& other){cout << "Error" << endl; exit(201);}
    ham_heis(ham_type<heis_spin>& other);
    ~ham_heis(){}
    virtual double calc_E(field_type<heis_spin>* lattice);
    vector<double> calc_M(field_type<heis_spin>* lattice);
    vector<double> calc_subM(field_type<heis_spin>* lattice, int subnumber);
    virtual double dE(field_type<heis_spin>* lattice, vector<int>& position);
    vector<double> get_Js(){return J;}
    vector<double> get_Hs(){return H;}
    T* get_test(){return &test;}
    ham_heis<T>& operator=(ham_type<heis_spin>& other);
};

template <class T> class ham_FePt: public ham_type<T>
{
private:
    T test;
    double Jx_sum, Jy_sum, Jz_sum, d2_sum, d0_sum, d0;
    vector<double> vsum, curr, adj_curr, potential, Js, d_ijs;
    vector<int> pos, dxs, dys, dzs, pos2;
    vector<vector<int> > posvec;
    vector<vector<double> > adjvec;
public:
    ham_FePt();
    ham_FePt(ham_type<ising_spin>& other){cout << "Error" << endl; exit(201);}
    ham_FePt(ham_type<heis_spin>& other);
    ~ham_FePt(){}
    double calc_E(field_type<heis_spin>* lattice);
    vector<double> calc_M(field_type<heis_spin>* lattice);
    double dE(field_type<heis_spin>* lattice, vector<int>& position);
    T* get_test(){return &test;}
    ham_FePt<T>& operator=(ham_type<heis_spin>& other);
    vector<double> calc_subM(field_type<heis_spin>* lattice, int subnumber);
    void read_Js();
    bool check_pos(int start, int fin);
};

template <class T> class ham_cluster: public ham_type<T>
{
private:
    T test;
    vector<double> xs, ys, zs, Rms, Rns, Vms, Vs, h;
    vector<double> temp, potential, diff;
    vector<double> curr, curr2, curr_ani;
    vector<T> ks;
    vector<int> pos, pos2;
    double ani_const, Ms, mu0;
public:
    ham_cluster(){}
    ham_cluster(string filename);
    ham_cluster(ham_type<ising_spin>& other){cout << "Error" << endl; exit(201);}
    ham_cluster(ham_type<heis_spin>& other);
    ~ham_cluster(){}
    double calc_E(field_type<heis_spin>* lattice);
    vector<double> calc_M(field_type<heis_spin>* lattice);
    double dE(field_type<heis_spin>* lattice, vector<int>& position);
    T* get_test(){return &test;}
    ham_cluster<T>& operator=(ham_type<heis_spin>& other);
    vector<double>* get_xs(){return &xs;}
    vector<double>* get_ys(){return &ys;}
    vector<double>* get_zs(){return &zs;}
    vector<double>* get_Rms(){return &Rms;}
    vector<double>* get_Rns(){return &Rns;}
    vector<T>* get_ks(){return &ks;}
};

// template <class T> class ham_skyrm: public ham_heis<T>
// {
// private:
//     double Dx, Dy;
// public:
//     ham_skyrm(){}
//     ham_skyrm(double Hin, Jin, Dxin, Dyin);
//     ham_skyrm(ham_type<ising_spin>& other){cout << "Error" << endl; exit(201);}
//     ham_skyrm(ham_type<heis_spin>& other);
//     ~ham_skyrm(){}
//     double dE(field_type<heis_spin>* lattice, vector<int>& position);
//     double calc_E(field_type<heis_spin>* lattice);
//     double get_Dx(){return Dx;}
//     double get_Dy(){return Dy;}
//     ham_skyrm<T>& operator=(ham_type<heis_spin>& other);
// };

#endif
