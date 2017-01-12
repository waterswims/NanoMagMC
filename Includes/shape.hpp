#ifndef _SHAPE
#define _SHAPE

#include <vector>

using namespace std;

class shape_type
{
public:
    shape_type(){}
    ~shape_type(){}
    virtual bool check(vector<int> Is, int l_size){return false;}
    virtual double get_r0(){return 0;}
    virtual double get_beta(){return 0;}
};

class shape_2d: public shape_type
{
public:
    shape_2d(){}
    ~shape_2d(){}
    virtual bool check(vector<int> Is, int l_size){return true;}
};

class shape_3d: public shape_type
{
public:
    shape_3d(){}
    ~shape_3d(){}
    virtual bool check(vector<int> Is, int l_size){return true;}
};

class weibull: public shape_type
{
private:
    double r0;
    double beta;
public:
    weibull(){r0 = 0; beta = 0;}
    weibull(shape_type& other){r0 = other.get_r0(); beta = other.get_beta();}
    weibull(double rin, double bin);
    ~weibull(){}
    bool check(vector<int> Is, int l_size);
    weibull& operator=(shape_type& other);
    double get_r0(){return r0;}
    double get_beta(){return beta;}
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

#endif
