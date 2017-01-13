#include "spin_type.hpp"
#include "field_type.hpp"
#include <fstream>
#include <sstream>
#include <iostream>

template <class T> field_2d<T>::field_2d()
{
    this->dim=2;
    this->periodic=true;
}

template <class T> field_3d<T>::field_3d()
{
    this->dim=3;
    this->periodic=true;
}

template <class T> field_cluster<T>::field_cluster()
{
    this->dim=1;
    this->periodic=false;
}

template <class T> field_2d<T>::field_2d(int size, bool isperio)
{
    this->dim=2;
    this->periodic=isperio;
    this->insize = size;
    if (this->periodic){
        this->totsize = this->insize;
    }
    else{
        this->totsize = this->insize + 2;
    }
    field = new boost::multi_array<T, 2>(boost::extents[this->totsize][this->totsize]);
}

template <class T> field_3d<T>::field_3d(int size, bool isperio)
{
    this->dim=3;
    this->periodic=isperio;
    this->insize = size;
    if(this->periodic){
        this->totsize = this->insize;
    }
    else{
        this->totsize = this->insize + 2;
    }
    field = new boost::multi_array<T, 3>(boost::extents[this->totsize][this->totsize][this->totsize]);
}

template <class T> field_cluster<T>::field_cluster(string filename)
{
    this->dim = 1;
    this->periodic=false;

    ifstream file;
    file.open(filename.c_str());
    if(!file.is_open())
    {
        cout << "Input file not opened" << endl;
        exit(105);
    }
    string line;
    double temp_d;
    while(getline(file, line))
    {
        istringstream stream(line);
        stream >> temp_d;
        (*xs).push_back(temp_d);
        stream >> temp_d;
        (*temp).push_back(temp_d);
        stream >> temp_d;
        (*temp).push_back(temp_d);
        stream >> temp_d;
        (*temp).push_back(temp_d);
        stream >> temp_d;
        (*temp).push_back(temp_d);
    }
    file.close();
    this->insize=(*xs).size();
    this->totsize=this->insize;
    field = new vector<T>(this->insize);
    for(int i=0; i < this->insize; i++)
    {
        (*field)[i].rand_spin();
    }
}

template <class T> field_3d<T>::field_3d(int size, bool isperio, int pad)
{
    this->dim=3;
    this->periodic=isperio;
    this->insize = size;
    if(this->periodic){
        this->totsize = this->insize + 2*pad;
    }
    else{
        this->totsize = this->insize + 2 + 2 * pad;
    }
    field = new boost::multi_array<T, 3>(boost::extents[this->totsize][this->totsize][this->totsize]);
}

template <class T> field_2d<T>::field_2d(field_type<T>& other)
{
    this->dim = 2;
    this->insize = other.get_insize();
    this->totsize = other.get_totsize();
    this->periodic = other.get_perio();
    field = new boost::multi_array<T, 2>(*(other.get_2dfield()));
}

template <class T> field_3d<T>::field_3d(field_type<T>& other)
{
    this->dim = 3;
    this->insize = other.get_insize();
    this->totsize = other.get_totsize();
    this->periodic = other.get_perio();
    field = new boost::multi_array<T, 3>(*(other.get_3dfield()));
}

template <class T> field_cluster<T>::field_cluster(field_type<T>& other)
{
    this->dim = 1;
    this->insize = other.get_insize();
    this->totsize = other.get_totsize();
    this->periodic = other.get_perio();
    xs = new vector<double>(*(other.get_xs()));
    ys = new vector<double>(*(other.get_ys()));
    zs = new vector<double>(*(other.get_zs()));
    Rms = new vector<double>(*(other.get_Rms()));
    Rns = new vector<double>(*(other.get_Rns()));
    ks = new vector<T>(*(other.get_ks()));
    field = new vector<T>(*(other.get_1dfield()));
}

template <class T> field_2d<T>::field_2d(const field_2d<T>& other)
{
    this->dim = 2;
    this->insize = other.insize;
    this->totsize = other.totsize;
    this->periodic = other.periodic;
    field = new boost::multi_array<T, 2>(*(other.field));
}

template <class T> field_3d<T>::field_3d(const field_3d<T>& other)
{
    this->dim = 3;
    this->insize = other.insize;
    this->totsize = other.totsize;
    this->periodic = other.periodic;
    field = new boost::multi_array<T, 3>(*(other.field));
}

template <class T> field_cluster<T>::field_cluster(const field_cluster<T>& other)
{
    this->dim = 1;
    this->insize = other.insize;
    this->totsize = other.totsize;
    this->periodic = other.periodic;
    xs = new vector<double>(*other.xs);
    ys = new vector<double>(*other.ys);
    zs = new vector<double>(*other.zs);
    Rms = new vector<double>(*other.Rms);
    Rns = new vector<double>(*other.Rns);
    ks = new vector<T>(*other.ks);
    field = new vector<T>(*other.field);
}

template <class T> T& field_2d<T>::access(vector<int>& position)
{
    return (*field)[position[0]][position[1]];
}

template <class T> T& field_3d<T>::access(vector<int>& position)
{
    return (*field)[position[0]][position[1]][position[2]];
}

template <class T> vector<double>& field_3d<T>::spin_access(vector<int>& position)
{
    return ((*field)[position[0]][position[1]][position[2]]).spin_access();
}

template <class T> void field_2d<T>::adjacent(vector<int>& position, vector<T*>& out)
{
    out.resize(4);
    // left
    if (position[0] != 0)
    {
        out[0] = &((*field)[position[0] - 1][position[1]]);
    }
    else
    {
        out[0] = &((*field)[this->totsize - 1][position[1]]);
    }
    // right
    if (position[0] != this->totsize-1)
    {
        out[1] = &((*field)[position[0] + 1][position[1]]);
    }
    else
    {
        out[1] = &((*field)[0][position[1]]);
    }
    // up
    if (position[1] != 0)
    {
        out[2] = &((*field)[position[0]][position[1]-1]);
    }
    else
    {
        out[2] = &((*field)[position[0]][this->totsize - 1]);
    }
    // down
    if (position[1] != this->totsize - 1)
    {
        out[3] = &((*field)[position[0]][position[1] + 1]);
    }
    else
    {
        out[3] = &((*field)[position[0]][0]);
    }
}

template <class T> void field_3d<T>::adjacent(vector<int>& position, vector<T*>& out)
{
    out.resize(6);
    // left
    if (position[0] != 0)
    {
        out[0] = &((*field)[position[0] - 1][position[1]][position[2]]);
    }
    else
    {
        out[0] = &((*field)[this->totsize - 1][position[1]][position[2]]);
    }
    // right
    if (position[0] != this->totsize-1)
    {
        out[1] = &((*field)[position[0] + 1][position[1]][position[2]]);
    }
    else
    {
        out[1] = &((*field)[0][position[1]][position[2]]);
    }
    // up
    if (position[1] != 0)
    {
        out[2] = &((*field)[position[0]][position[1] - 1][position[2]]);
    }
    else
    {
        out[2] = &((*field)[position[0]][this->totsize - 1][position[2]]);
    }
    // down
    if (position[1] != this->totsize - 1)
    {
        out[3] = &((*field)[position[0]][position[1] + 1][position[2]]);
    }
    else
    {
        out[3] = &((*field)[position[0]][0][position[2]]);
    }
    // forward
    if (position[2] != 0)
    {
        out[4] = &((*field)[position[0]][position[1]][position[2]-1]);
    }
    else
    {
        out[4] = &((*field)[position[0]][position[1]][this->totsize - 1]);
    }
    // backward
    if (position[2] != this->totsize - 1)
    {
        out[5] = &((*field)[position[0]][position[1]][position[2]+1]);
    }
    else
    {
        out[5] = &((*field)[position[0]][position[1]][0]);
    }
}

template <class T> T& field_2d<T>::next(bool &finish, vector<int> &pos)
{
    int start = 0;
    int end = this->totsize;
    if(!this->periodic)
    {
        start++;
        end--;
    }
    T* out = &(this->access(pos));
    pos[1]++;
    if(pos[1] == end)
    {
        pos[1] = start;
        pos[0]++;
        if(pos[0] == end)
        {
            pos[0] = start;
            finish = true;
        }
    }
    return *out;
}

template <class T> T& field_3d<T>::next(bool &finish, vector<int> &pos)
{
    int start = 0;
    int end = this->totsize;
    if(!this->periodic)
    {
        start++;
        end--;
    }
    T* out = &(this->access(pos));
    pos[2]++;
    if(pos[2] == end)
    {
        pos[2] = start;
        pos[1]++;
        if(pos[1] == end)
        {
            pos[1] = start;
            pos[0]++;
            if(pos[0] == end)
            {
                pos[0] = start;
                finish = true;
            }
        }
    }
    return *out;
}

template <class T> T& field_cluster<T>::next(bool &finish, vector<int> &pos)
{
    T* out = &(this->access(pos));
    pos[0]++;
    if(pos[0] == this->insize)
    {
        pos[0] = 0;
        finish = true;
    }
    return *out;
}

template <class T> void field_2d<T>::fill_ghost(int &num)
{
    if(!(this->periodic))
    {
        for (int i=0; i < this->totsize; i++)
        {
            // if((*field)[0][i].is_zero()) num--;
            // if((*field)[i][0].is_zero()) num--;
            // if((*field)[this->totsize-1][i].is_zero()) num--;
            // if((*field)[i][this->totsize-1].is_zero()) num--;

            (*field)[0][i].zero_spin();
            (*field)[i][0].zero_spin();
            (*field)[this->totsize-1][i].zero_spin();
            (*field)[i][this->totsize-1].zero_spin();
        }
    }
}

template <class T> void field_3d<T>::fill_ghost(int &num)
{
    if(!(this->periodic))
    {
        for (int i=0; i < this->totsize; i++)
        {
            for (int j=0; j < this->totsize; j++)
            {
                // if(!(*field)[0][i][j].is_zero()) num--;
                // if(!(*field)[i][j][0].is_zero()) num--;
                // if(!(*field)[j][0][i].is_zero()) num--;
                // if(!(*field)[this->totsize-1][i][j].is_zero()) num--;
                // if(!(*field)[i][j][this->totsize-1].is_zero()) num--;
                // if(!(*field)[j][this->totsize-1][i].is_zero()) num--;

                (*field)[0][i][j].zero_spin();
                (*field)[i][j][0].zero_spin();
                (*field)[j][0][i].zero_spin();
                (*field)[this->totsize-1][i][j].zero_spin();
                (*field)[i][j][this->totsize-1].zero_spin();
                (*field)[j][this->totsize-1][i].zero_spin();
            }
        }
    }
}

template <class T> void field_2d<T>::new_mem()
{
    field = new boost::multi_array<T, 2>(*field);
}

template <class T> void field_3d<T>::new_mem()
{
    field = new boost::multi_array<T, 3>(*field);
}

template <class T> field_2d<T>& field_2d<T>::operator=(field_2d<T>& other)
{
    this->dim = other.get_dim();
    this->insize = other.get_insize();
    this->totsize = other.get_totsize();
    this->periodic = other.get_perio();
    field = new boost::multi_array<T, 2>(*(other.get_2dfield()));
    return *this;
}

template <class T> field_3d<T>& field_3d<T>::operator=(field_3d<T>& other)
{
    this->dim = other.get_dim();
    this->insize = other.get_insize();
    this->totsize = other.get_totsize();
    this->periodic = other.get_perio();
    field = new boost::multi_array<T, 3>(*(other.get_3dfield()));
    return *this;
}

template <class T> field_cluster<T>& field_cluster<T>::operator=(field_cluster<T>& other)
{
    this->dim = 1;
    this->insize = other.get_insize();
    this->totsize = other.get_totsize();
    this->periodic = other.get_perio();
    xs = new vector<double>(*(other.get_xs()));
    ys = new vector<double>(*(other.get_ys()));
    zs = new vector<double>(*(other.get_zs()));
    Rms = new vector<double>(*(other.get_Rms()));
    Rns = new vector<double>(*(other.get_Rns()));
    ks = new vector<T>(*(other.get_ks()));
    field = new vector<T>(*(other.get_1dfield()));
    return *this;
}

template <class T> int field_2d<T>::findnum()
{
    int count = 0;
    for(int i = 0; i < this->totsize; i++)
    {
        for(int j=0; j < this->totsize; j++)
        {
            if(!((*field)[i][j].is_zero()))
            {
                count++;
            }
        }
    }
    return count;
}

template <class T> int field_3d<T>::findnum()
{
    int count = 0;
    for(int i = 0; i < this->totsize; i++)
    {
        for(int j=0; j < this->totsize; j++)
        {
            for(int k=0; k < this->totsize; k++)
            {
                if(!((*field)[i][j][k].is_zero()))
                {
                    count++;
                }
            }
        }
    }
    return count;
}

template <class T> void field_2d<T>::print()
{
    for(int i = 0; i < this->totsize; i++)
    {
        for(int j=0; j < this->totsize; j++)
        {
            (*field)[i][j].print();
            cout << " ";
        }
        cout << endl;
    }
}

template <class T> hex_2d<T>::hex_2d(const hex_2d<T>& other)
{
    this->dim = 2;
    this->insize = other.insize;
    this->totsize = other.totsize;
    this->periodic = other.periodic;
    this->field = new boost::multi_array<T, 2>(*(other.field));
}

template <class T> hex_2d<T>& hex_2d<T>::operator=(hex_2d<T>& other)
{
    this->dim = other.get_dim();
    this->insize = other.get_insize();
    this->totsize = other.get_totsize();
    this->periodic = other.get_perio();
    this->field = new boost::multi_array<T, 2>(*(other.get_2dfield()));
    return *this;
}

template <class T> void hex_2d<T>::adjacent(vector<int> &position, vector<T*> &out)
{
    out.resize(3);
    // left
    if((position[0]+position[1])%2 == 0)
    {
        if (position[0] != 0)
        {
            out[0] = &((*(this->field))[position[0] - 1][position[1]]);
        }
        else
        {
            out[0] = &((*(this->field))[this->totsize - 1][position[1]]);
        }
    }
    // right
    else
    {
        if (position[0] != this->totsize-1)
        {
            out[0] = &((*(this->field))[position[0] + 1][position[1]]);
        }
        else
        {
            out[0] = &((*(this->field))[0][position[1]]);
        }
    }
    // up
    if (position[1] != 0)
    {
        out[1] = &((*(this->field))[position[0]][position[1]-1]);
    }
    else
    {
        out[1] = &((*(this->field))[position[0]][this->totsize - 1]);
    }
    // down
    if (position[1] != this->totsize - 1)
    {
        out[2] = &((*(this->field))[position[0]][position[1] + 1]);
    }
    else
    {
        out[2] = &((*(this->field))[position[0]][0]);
    }
}

template <class T> hex_3d<T>::hex_3d(const hex_3d<T>& other)
{
    this->dim = 3;
    this->insize = other.insize;
    this->totsize = other.totsize;
    this->periodic = other.periodic;
    this->field = new boost::multi_array<T, 3>(*(other.field));
}

template <class T> hex_3d<T>& hex_3d<T>::operator=(hex_3d<T>& other)
{
    this->dim = other.get_dim();
    this->insize = other.get_insize();
    this->totsize = other.get_totsize();
    this->periodic = other.get_perio();
    this->field = new boost::multi_array<T, 3>(*(other.get_3dfield()));
    return *this;
}

template <class T> void hex_3d<T>::adjacent(vector<int>& position, vector<T*>& out)
{
    out.resize(5);
    // left
    if((position[0]+position[1])%2 == 0)
    {
        if (position[0] != 0)
        {
            out[0] = &((*(this->field))[position[0] - 1][position[1]][position[2]]);
        }
        else
        {
            out[0] = &((*(this->field))[this->totsize - 1][position[1]][position[2]]);
        }
    }
    // right
    else
    {
        if (position[0] != this->totsize-1)
        {
            out[0] = &((*(this->field))[position[0] + 1][position[1]][position[2]]);
        }
        else
        {
            out[0] = &((*(this->field))[0][position[1]][position[2]]);
        }
    }
    // up
    if (position[1] != 0)
    {
        out[1] = &((*(this->field))[position[0]][position[1] - 1][position[2]]);
    }
    else
    {
        out[1] = &((*(this->field))[position[0]][this->totsize - 1][position[2]]);
    }
    // down
    if (position[1] != this->totsize - 1)
    {
        out[2] = &((*(this->field))[position[0]][position[1] + 1][position[2]]);
    }
    else
    {
        out[2] = &((*(this->field))[position[0]][0][position[2]]);
    }
    // forward
    if (position[2] != 0)
    {
        out[3] = &((*(this->field))[position[0]][position[1]][position[2]-1]);
    }
    else
    {
        out[3] = &((*(this->field))[position[0]][position[1]][this->totsize - 1]);
    }
    // backward
    if (position[2] != this->totsize - 1)
    {
        out[4] = &((*(this->field))[position[0]][position[1]][position[2]+1]);
    }
    else
    {
        out[4] = &((*(this->field))[position[0]][position[1]][0]);
    }
}

template class field_2d<ising_spin>;
template class field_2d<heis_spin>;
template class field_3d<ising_spin>;
template class field_3d<heis_spin>;
template class field_cluster<ising_spin>;
template class field_cluster<heis_spin>;
template class hex_2d<heis_spin>;
template class hex_2d<ising_spin>;
template class hex_3d<heis_spin>;
template class hex_3d<ising_spin>;
