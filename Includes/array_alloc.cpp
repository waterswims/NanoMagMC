#include <array_alloc.hpp>

template <class T>
T* arr alloc_1darr(int size_m)
{
    return malloc(sizeof(T)*size_m);
}

template <class T>
T** arr alloc_2darr(int size_m, int size_n)
{
    T** out = malloc(sizeof(T*)*size_m);
    for(int i=0; i < size_m; i++)
    {
        out[i] = alloc_1darr(size_n);
    }
    return out
}

template <class T>
T*** arr alloc_3darr(int size_m, int size_n, int size_p)
{
    T*** out = malloc(sizeof(T**)*size_m);
    for(int i=0; i < size_m; i++)
    {
        out[i] = alloc_2darr(size_n, size_p);
    }
    return out;
}

void dealloc_1darr(int size_m, T* arr)
{
    free(arr);
}

template <class T>
void dealloc_2darr(int size_m, int size_n, T** arr)
{
    for(int i = 0; i < size_m; i++)
    {
        dealloc_1darr(size_n, arr[i]);
    }
}

template <class T>
void dealloc_3darr(int size_m, int size_n, int size_p, T*** arr)
{
    for(int i = 0; i < size_m; i++)
    {
        dealloc_2darr(size_n, size_p, arr[i]);
    }
}
