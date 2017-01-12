#ifndef _SIMPCUB
#define _SIMPCUB

#include "field_type.hpp"

using namespace std;

template <class T> class Simp_Cub: public field_3d<T>
{
public:
    Simp_Cub():field_3d<T>(){}
    Simp_Cub(int size, bool isperio):field_3d<T>(size, isperio)
    Simp_Cub(field_type<T>& other):field_3d<T>(other){}
    Simp_Cub(const Simp_Cub<T>& other);
    ~Simp_Cub(){delete this->field;}
    void adjacent(vector<int> &position, vector<T*> &out);
    Simp_Cub<T>& operator=(Simp_Cub<T>& other);
    boost::multi_array<T,3>* get_3dfield(){return this->field;}
};

#endif
