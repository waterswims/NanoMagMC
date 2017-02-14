#include <array_alloc.hpp>
#include <cstdlib>

template <class T>
T* arr_alloc_1darr(int size_m)
{
    return malloc(sizeof(T)*size_m);
}

template <class T>
T** arr_alloc_2darr(int size_m, int size_n)
{
    T** out = malloc(sizeof(T*)*size_m);
    for(int i=0; i < size_m; i++)
    {
        out[i] = alloc_1darr(size_n);
    }
    return out
}

template <class T>
T*** arr_alloc_3darr(int size_m, int size_n, int size_p)
{
    T*** out = malloc(sizeof(T**)*size_m);
    for(int i=0; i < size_m; i++)
    {
        out[i] = alloc_2darr(size_n, size_p);
    }
    return out;
}

template <class T>
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

template <class T>
T* deep_copy_1darr(int size_m, T* arr)
{
    T* out = alloc_1darr(size_m);
    for(int i = 0; i < size_m; i++)
    {
        out[i] = arr[i];
    }
    return out;
}

template <class T>
T** deep_copy_2darr(int size_m, int size_n, T** arr)
{
    T** out = alloc_1darr(size_m);
    for(int i = 0; i < size_m; i++)
    {
        for(int j = 0; j < size_n; j++)
        {
            out[i][j] = arr[i][j];
        }
    }
    return out;
}

template <class T>
T*** deep_copy_3darr(int size_m, int size_n, int size_p, T*** arr)
{
    T*** out = alloc_1darr(size_m);
    for(int i = 0; i < size_m; i++)
    {
        for(int j = 0; j < size_n; j++)
        {
            for(int k = 0; k < size_p; k++)
            {
                out[i][j][k] = arr[i][j][k];
            }
        }
    }
    return out;
}
