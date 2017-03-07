#include "spin_type.hpp"
#include "FCC.hpp"
#include "functions.h"
#include <iostream>

template <class T> FC_Tetr<T>::FC_Tetr(const FC_Tetr<T>& other)
{
    this->dim = 3;
    this->insize = other.insize;
    this->totsize = other.totsize;
    this->periodic = other.periodic;
    this->field = new boost::multi_array<T, 3>(*(other.field));
}

template <class T> FC_Tetr<T>& FC_Tetr<T>::operator=(FC_Tetr<T>& other)
{
    this->dim = other.get_dim();
    this->insize = other.get_insize();
    this->totsize = other.get_totsize();
    this->periodic = other.get_perio();
    this->field = new boost::multi_array<T, 3>(*(other.get_3dfield()));
    return *this;
}

template <class T> void FC_Tetr<T>::sym_plane_adj(vector<int> &position, vector<T*> &out)
{
    out.resize(4);
    if (position[0] != 0)
    {
        lshift = position[0] - 1;
    }
    else
    {
        lshift = this->totsize - 1;
    }
    if (position[0] != this->totsize-1)
    {
        rshift = position[0] + 1;
    }
    else
    {
        rshift = 0;
    }
    if (position[1] != 0)
    {
        ushift = position[1] - 1;
    }
    else
    {
        ushift = this->totsize - 1;
    }
    if (position[1] != this->totsize-1)
    {
        dshift = position[1] + 1;
    }
    else
    {
        dshift = 0;
    }
    // left up
    out[0] = &((*this->field)[lshift][ushift][position[2]]);
    // right up
    out[1] = &((*this->field)[rshift][ushift][position[2]]);
    // left down
    out[2] = &((*this->field)[lshift][dshift][position[2]]);
    // right down
    out[3] = &((*this->field)[rshift][dshift][position[2]]);
}

template <class T> void FC_Tetr<T>::sym_plane_diag(vector<int> &position, vector<T*> &out)
{
    out.resize(4);
    // left
    if (position[0] != 0)
    {
        out[0] = &((*this->field)[position[0] - 1][position[1]][position[2]]);
    }
    else
    {
        out[0] = &((*this->field)[this->totsize - 1][position[1]][position[2]]);
    }
    // right
    if (position[0] != this->totsize-1)
    {
        out[1] = &((*this->field)[position[0] + 1][position[1]][position[2]]);
    }
    else
    {
        out[1] = &((*this->field)[0][position[1]][position[2]]);
    }
    // up
    if (position[1] != 0)
    {
        out[2] = &((*this->field)[position[0]][position[1] - 1][position[2]]);
    }
    else
    {
        out[2] = &((*this->field)[position[0]][this->totsize - 1][position[2]]);
    }
    // down
    if (position[1] != this->totsize - 1)
    {
        out[3] = &((*this->field)[position[0]][position[1] + 1][position[2]]);
    }
    else
    {
        out[3] = &((*this->field)[position[0]][0][position[2]]);
    }
}

template <class T> void FC_Tetr<T>::face_diag(vector<int> &position, vector<T*> &out)
{
    out.resize(8);
    if (position[0] != this->totsize-1)
    {
        rshift = position[0] + 1;
    }
    else
    {
        rshift = 0;
    }
    if (position[1] != this->totsize-1)
    {
        dshift = position[1] + 1;
    }
    else
    {
        dshift = 0;
    }
    if (position[2] != 0)
    {
        fshift = position[2] - 1;
    }
    else
    {
        fshift = this->totsize-1;
    }
    if (position[2] != this->totsize-1)
    {
        bshift = position[2] + 1;
    }
    else
    {
        bshift = 0;
    }

    // left up f
    out[0] = &((*this->field)[position[0]][position[1]][fshift]);
    // right up f
    out[1] = &((*this->field)[rshift][position[1]][fshift]);
    // left down f
    out[2] = &((*this->field)[position[0]][dshift][fshift]);
    // right down f
    out[3] = &((*this->field)[rshift][dshift][fshift]);
    // left up b
    out[4] = &((*this->field)[position[0]][position[1]][bshift]);
    // right up b
    out[5] = &((*this->field)[rshift][position[1]][bshift]);
    // left down b
    out[6] = &((*this->field)[position[0]][dshift][bshift]);
    // right down b
    out[7] = &((*this->field)[rshift][dshift][bshift]);
}

template <class T> void FC_Tetr<T>::long_plane_adj(vector<int> &position, vector<T*> &out)
{
    out.resize(2);
    fshift = mod((position[2] - 2), this->totsize);
    bshift = (position[2] + 2)%(this->totsize);

    // forward
    out[0] = &((*this->field)[position[0]][position[1]][fshift]);
    // backward
    out[1] = &((*this->field)[position[0]][position[1]][bshift]);
}

template <class T> void FC_Tetr<T>::long_plane_diag(vector<int> &position, vector<T*> &out)
{
    out.resize(8);
    lshift = mod((position[0] - 1), this->totsize);
    rshift = (position[0] + 1)%(this->totsize);
    ushift = mod((position[1] - 1), this->totsize);
    dshift = (position[1] + 1)%(this->totsize);
    fshift = mod((position[2] - 2), this->totsize);
    bshift = (position[2] + 2)%(this->totsize);

    out[0] = &((*this->field)[lshift][ushift][fshift]);
    out[1] = &((*this->field)[lshift][dshift][fshift]);
    out[2] = &((*this->field)[rshift][ushift][fshift]);
    out[3] = &((*this->field)[rshift][dshift][fshift]);
    out[4] = &((*this->field)[lshift][ushift][bshift]);
    out[5] = &((*this->field)[lshift][dshift][bshift]);
    out[6] = &((*this->field)[rshift][ushift][bshift]);
    out[7] = &((*this->field)[rshift][dshift][bshift]);
}

template <class T> L10<T>::L10(const L10<T>& other)
{
    this->dim = 3;
    this->insize = other.insize;
    this->totsize = other.totsize;
    this->periodic = other.periodic;
    this->field = new boost::multi_array<T, 3>(*(other.field));
}

template <class T> L10<T>& L10<T>::operator=(L10<T>& other)
{
    this->dim = other.get_dim();
    this->insize = other.get_insize();
    this->totsize = other.get_totsize();
    this->periodic = other.get_perio();
    this->field = new boost::multi_array<T, 3>(*(other.get_3dfield()));
    return *this;
}

template <class T> void L10<T>::sym_plane_adj(vector<int> &position, vector<T*> &out)
{
    out.resize(4);
    if (position[0] != 0)
    {
        lshift = position[0] - 1;
    }
    else
    {
        lshift = this->totsize - 1;
    }
    if (position[0] != this->totsize-1)
    {
        rshift = position[0] + 1;
    }
    else
    {
        rshift = 0;
    }
    if (position[1] != 0)
    {
        ushift = position[1] - 1;
    }
    else
    {
        ushift = this->totsize - 1;
    }
    if (position[1] != this->totsize-1)
    {
        dshift = position[1] + 1;
    }
    else
    {
        dshift = 0;
    }
    // left up
    out[0] = &((*this->field)[lshift][ushift][position[2]]);
    // right up
    out[1] = &((*this->field)[rshift][ushift][position[2]]);
    // left down
    out[2] = &((*this->field)[lshift][dshift][position[2]]);
    // right down
    out[3] = &((*this->field)[rshift][dshift][position[2]]);
}

template <class T> void L10<T>::sym_plane_diag(vector<int> &position, vector<T*> &out)
{
    out.resize(4);
    // left
    if (position[0] != 0)
    {
        out[0] = &((*this->field)[position[0] - 1][position[1]][position[2]]);
    }
    else
    {
        out[0] = &((*this->field)[this->totsize - 1][position[1]][position[2]]);
    }
    // right
    if (position[0] != this->totsize-1)
    {
        out[1] = &((*this->field)[position[0] + 1][position[1]][position[2]]);
    }
    else
    {
        out[1] = &((*this->field)[0][position[1]][position[2]]);
    }
    // up
    if (position[1] != 0)
    {
        out[2] = &((*this->field)[position[0]][position[1] - 1][position[2]]);
    }
    else
    {
        out[2] = &((*this->field)[position[0]][this->totsize - 1][position[2]]);
    }
    // down
    if (position[1] != this->totsize - 1)
    {
        out[3] = &((*this->field)[position[0]][position[1] + 1][position[2]]);
    }
    else
    {
        out[3] = &((*this->field)[position[0]][0][position[2]]);
    }
}

template <class T> void L10<T>::long_plane_adj(vector<int> &position, vector<T*> &out)
{
    out.resize(2);
    if (position[2] != 0)
    {
        fshift = position[2] - 1;
    }
    else
    {
        fshift = this->totsize-1;
    }
    if (position[2] != this->totsize-1)
    {
        bshift = position[2] + 1;
    }
    else
    {
        bshift = 0;
    }

    // forward
    out[0] = &((*this->field)[position[0]][position[1]][fshift]);
    // backward
    out[1] = &((*this->field)[position[0]][position[1]][bshift]);
}

template <class T> void L10<T>::long_plane_diag(vector<int> &position, vector<T*> &out)
{
    out.resize(8);
    lshift = mod((position[0] - 1), this->totsize);
    rshift = (position[0] + 1)%(this->totsize);
    ushift = mod((position[1] - 1), this->totsize);
    dshift = (position[1] + 1)%(this->totsize);
    fshift = mod((position[2] - 1), this->totsize);
    bshift = (position[2] + 1)%(this->totsize);

    out[0] = &((*this->field)[lshift][ushift][fshift]);
    out[1] = &((*this->field)[lshift][dshift][fshift]);
    out[2] = &((*this->field)[rshift][ushift][fshift]);
    out[3] = &((*this->field)[rshift][dshift][fshift]);
    out[4] = &((*this->field)[lshift][ushift][bshift]);
    out[5] = &((*this->field)[lshift][dshift][bshift]);
    out[6] = &((*this->field)[rshift][ushift][bshift]);
    out[7] = &((*this->field)[rshift][dshift][bshift]);
}

template <class T> void L10<T>::common_neigh(vector<int> &position, vector<T*> &out, vector<int> &num_neigh)
{
    out.resize(26);
    num_neigh.resize(26);
    lshift = mod((position[0] - 1), this->totsize);
    rshift = (position[0] + 1)%(this->totsize);
    ushift = mod((position[1] - 1), this->totsize);
    dshift = (position[1] + 1)%(this->totsize);
    fshift = mod((position[2] - 1), this->totsize);
    bshift = (position[2] + 1)%(this->totsize);

    // centre axis
    num_neigh[0] = 4;
    num_neigh[1] = 4;
    num_neigh[2] = 4;
    num_neigh[3] = 4;
    num_neigh[4] = 2;
    num_neigh[5] = 2;
    num_neigh[6] = 2;
    num_neigh[7] = 2;

    out[0] = &((*this->field)[lshift][position[1]][position[2]]);
    out[1] = &((*this->field)[rshift][position[1]][position[2]]);
    out[2] = &((*this->field)[position[0]][ushift][position[2]]);
    out[3] = &((*this->field)[position[0]][dshift][position[2]]);
    out[4] = &((*this->field)[lshift][ushift][position[2]]);
    out[5] = &((*this->field)[rshift][ushift][position[2]]);
    out[6] = &((*this->field)[lshift][dshift][position[2]]);
    out[7] = &((*this->field)[rshift][dshift][position[2]]);

    // lower axis
    num_neigh[8] = 2;
    num_neigh[9] = 2;
    num_neigh[10] = 2;
    num_neigh[11] = 2;
    num_neigh[12] = 1;
    num_neigh[13] = 1;
    num_neigh[14] = 1;
    num_neigh[15] = 1;
    num_neigh[16] = 4;

    out[8] = &((*this->field)[lshift][position[1]][fshift]);
    out[9] = &((*this->field)[rshift][position[1]][fshift]);
    out[10] = &((*this->field)[position[0]][ushift][fshift]);
    out[11] = &((*this->field)[position[0]][dshift][fshift]);
    out[12] = &((*this->field)[lshift][ushift][fshift]);
    out[13] = &((*this->field)[rshift][ushift][fshift]);
    out[14] = &((*this->field)[lshift][dshift][fshift]);
    out[15] = &((*this->field)[rshift][dshift][fshift]);
    out[16] = &((*this->field)[position[0]][position[1]][fshift]);

    // upper axis
    num_neigh[17] = 2;
    num_neigh[18] = 2;
    num_neigh[19] = 2;
    num_neigh[20] = 2;
    num_neigh[21] = 1;
    num_neigh[22] = 1;
    num_neigh[23] = 1;
    num_neigh[24] = 1;
    num_neigh[25] = 4;

    out[17] = &((*this->field)[lshift][position[1]][bshift]);
    out[18] = &((*this->field)[rshift][position[1]][bshift]);
    out[19] = &((*this->field)[position[0]][ushift][bshift]);
    out[20] = &((*this->field)[position[0]][dshift][bshift]);
    out[21] = &((*this->field)[lshift][ushift][bshift]);
    out[22] = &((*this->field)[rshift][ushift][bshift]);
    out[23] = &((*this->field)[lshift][dshift][bshift]);
    out[24] = &((*this->field)[rshift][dshift][bshift]);
    out[25] = &((*this->field)[position[0]][position[1]][bshift]);
}

template class FC_Tetr<ising_spin>;
template class FC_Tetr<heis_spin>;
template class L10<ising_spin>;
template class L10<heis_spin>;
