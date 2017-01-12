#include "simplecubic.hpp"

template <class T> Simp_Cub<T>::Simp_Cub(const Simp_Cub<T>& other)
{
    this->dim = 3;
    this->insize = other.insize;
    this->totsize = other.totsize;
    this->periodic = other.periodic;
    this->field = new boost::multi_array<T,3>(*(other.field));
}

template <class T> Simp_Cub<T>& Simp_Cub<T>::operator=(Simp_Cub<T>& other)
{
    this->dim = other.get_dim();
    this->insize = other.get_insize();
    this->totsize = other.get_totsize();
    this->periodic = other.get_perio();
    this->field = new boost::multi_array<T, 3>(*(other.get_3dfield));
    return *this;
}

template <class T> void Simp_Cub<T>::adjacent(vector<int> &position, vector<T*> &out)
{
    out.resize(8);
    
}
