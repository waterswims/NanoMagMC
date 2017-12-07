#ifndef _SHAPE
#define _SHAPE

#include <vector>

class shape_type
{
public:
    shape_type(){}
    ~shape_type(){}
    virtual bool check(std::vector<int> Is, int l_size){return false;}
    virtual double get_r0(){return 0;}
    virtual double get_beta(){return 0;}
    virtual double get_a(){return 0;}
    virtual double get_b(){return 0;}
    virtual double get_c(){return 0;}
};

class shape_2d: public shape_type
{
public:
    shape_2d(){}
    ~shape_2d(){}
    virtual bool check(std::vector<int> Is, int l_size){return true;}
};

class shape_3d: public shape_type
{
public:
    shape_3d(){}
    ~shape_3d(){}
    virtual bool check(std::vector<int> Is, int l_size){return true;}
};

class weibull: public shape_type
{
private:
    double r0;
    double beta;
    double a[3];
public:
    weibull(){r0 = 0; beta = 0; a[0] = 1; a[1] = 1; a[2] = 1;}
    weibull(shape_type& other);
    weibull(double rin, double bin);
    weibull(double betain, double ain, double bin, double cin);
    ~weibull(){}
    bool check(std::vector<int> Is, int l_size);
    weibull& operator=(shape_type& other);
    double get_r0(){return r0;}
    double get_beta(){return beta;}
    double get_a(){return a[0];}
    double get_b(){return a[1];}
    double get_c(){return a[2];}
};

class square: public shape_2d
{
public:
    square(){}
    ~square(){}
};

class cube: public shape_3d
{
public:
    cube(){}
    ~cube(){}
};

class sh_cluster: public shape_3d
{
public:
    sh_cluster(){}
    ~sh_cluster(){}
};

#endif
