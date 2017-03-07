#ifndef _FCC
#define _FCC

#include "field_type.hpp"

template <class T> class FC_Tetr: public field_3d<T>
{
private:
    int lshift, rshift, ushift, dshift, fshift, bshift;
public:
    FC_Tetr():field_3d<T>(){}
    FC_Tetr(int size, bool isperio):field_3d<T>(size, isperio){}
    FC_Tetr(field_type<T>& other):field_3d<T>(other){}
    FC_Tetr(const FC_Tetr<T>& other);
    ~FC_Tetr(){delete this->field;}
    void sym_plane_adj(vector<int> &position, vector<T*> &out);
    void sym_plane_diag(vector<int> &position, vector<T*> &out);
    void face_diag(vector<int> &position, vector<T*> &out);
    void long_plane_adj(vector<int> &position, vector<T*> &out);
    void long_plane_diag(vector<int> &position, vector<T*> &out);
    FC_Tetr<T>& operator=(FC_Tetr<T>& other);
    boost::multi_array<T,3>* get_3dfield(){return this->field;}
};

template <class T> class L10: public field_3d<T>
{
private:
    int lshift, rshift, ushift, dshift, fshift, bshift;
public:
    L10():field_3d<T>(){}
    L10(int size, bool isperio):field_3d<T>(size, isperio){}
    L10(field_type<T>& other):field_3d<T>(other){}
    L10(const L10<T>& other);
    ~L10(){delete this->field;}
    L10<T>& operator=(L10<T>& other);
    boost::multi_array<T,3>* get_3dfield(){return this->field;}
    void sym_plane_adj(vector<int> &position, vector<T*> &out);
    void sym_plane_diag(vector<int> &position, vector<T*> &out);
    void long_plane_adj(vector<int> &position, vector<T*> &out);
    void long_plane_diag(vector<int> &position, vector<T*> &out);
    void common_neigh(vector<int> &position, vector<T*> &out, vector<int> &num_neigh);
};

#endif
