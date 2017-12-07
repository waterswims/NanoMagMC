#ifndef _HAMIL
#define _HAMIL

#include "field_type.hpp"
#include "array_alloc.hpp"
#include <vector>

#ifndef _PI
#define _PI 3.141592653589793
#endif

class field_type;

class ham_type
{
protected:
    int dim;
public:
    virtual double calc_E(field_type* lattice) {return 0;}
    virtual double dE(field_type* lattice, std::vector<int>& position) {return 0;}
    virtual std::vector<double> calc_M(field_type* lattice) {}
    virtual std::vector<double> calc_subM(field_type* lattice, int subnumber) {}
    virtual double get_J() const {return 0;}
    virtual double get_H() const {return 0;}
    virtual std::vector<double> get_Js() const {}
    virtual std::vector<double> get_Hs() const {}
    virtual double get_K() const {return 0;}
    virtual void get_test(double &x, double &y, double &z){}
    virtual void init_dim(field_type* field) {}
    virtual void set_H(double Hin) {}
    virtual std::vector<double> calc_top_charge(field_type* lattice) {}
    // virtual double get_Dx(){return 0;}
    // virtual double get_Dy(){return 0;}
    // virtual std::vector<double>* get_xs(){return NULL;}
    // virtual std::vector<double>* get_ys(){return NULL;}
    // virtual std::vector<double>* get_zs(){return NULL;}
    // virtual std::vector<double>* get_hs(){return NULL;}
    // virtual std::vector<double>* get_Rms(){return NULL;}
    // virtual std::vector<double>* get_Rns(){return NULL;}
    // virtual std::vector<heis_spin>* get_ks(){return NULL;}
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
    std::vector<double> calc_M(field_type* lattice);
    std::vector<double> calc_subM(field_type* lattice, int subnumber);
    double dE(field_type* lattice, std::vector<int>& position);
    double get_J() const {return J;}
    double get_H() const {return H;}
    ham_ising& operator=(ham_type& other);
    void init_dim(field_type* field);
    void set_H(double Hin) {H = Hin;}
};

class ham_heis: public ham_type
{
protected:
    std::vector<double> H, J;
    double** adj;
    std::vector<double> vsum, curr, H_sum, J_sum, test;
    std::vector<double> tchar, s2, s3, sbuff;
    std::vector<int> pos;
    bool is3d;
    int tchar_size;
public:
    ham_heis(){}
    ham_heis(double Hin, double Jin);
    ham_heis(ham_type& other);
    ~ham_heis(){dealloc_2darr<double>(4, adj);}
    virtual double calc_E(field_type* lattice);
    std::vector<double> calc_M(field_type* lattice);
    std::vector<double> calc_subM(field_type* lattice, int subnumber);
    virtual double dE(field_type* lattice, std::vector<int>& position);
    std::vector<double> get_Js() const {return J;}
    std::vector<double> get_Hs() const {return H;}
    ham_heis& operator=(ham_type& other);
    virtual void init_dim(field_type* field);
    void get_test(double &x, double &y, double &z)
        {x = test[0]; y = test[1]; z = test[2];}
    void set_H(double Hin) {H[2] = Hin;}
    std::vector<double> calc_top_charge(field_type* lattice);
};

class ham_FePt: public ham_heis
{
private:
    std::vector<int> pos2, dxs, dys, dzs;
    std::vector<double> Js, adj_curr, d_ijs;
    double d0;
public:
    ham_FePt();
    ham_FePt(double Hin);
    ham_FePt(ham_type& other);
    ~ham_FePt() {}
    double calc_E(field_type* lattice);
    double dE(field_type* lattice, std::vector<int>& position);
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
    std::vector<double> cmp;
    std::vector<double> J2_sum;
public:
    ham_skyrm(){}
    ham_skyrm(double Hin, double Jin, double Kin);
    ham_skyrm(const ham_type& other);
    ~ham_skyrm(){}
    virtual double calc_E(field_type* lattice);
    virtual double dE(field_type* lattice, std::vector<int>& position);
    ham_skyrm& operator=(const ham_type& other);
    double get_K() const {return K;}
    void set_dirs();
};

#endif
