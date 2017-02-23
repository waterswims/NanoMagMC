#include <array_alloc.hpp>
#include <cstdlib>

template <class T>
T* alloc_1darr(int size_m)
{
    return (T*)malloc(sizeof(T)*size_m);
}

template <class T>
T** alloc_2darr(int size_m, int size_n)
{
    T** out = (T**)malloc(sizeof(T*)*size_m);
    for(int i=0; i < size_m; i++)
    {
        out[i] = alloc_1darr<T>(size_n);
    }
    return out;
}

template <class T>
T*** alloc_3darr(int size_m, int size_n, int size_p)
{
    T*** out = (T***)malloc(sizeof(T**)*size_m);
    for(int i=0; i < size_m; i++)
    {
        out[i] = alloc_2darr<T>(size_n, size_p);
    }
    return out;
}

template <class T>
void dealloc_1darr(T* arr)
{
    free(arr);
}

template <class T>
void dealloc_2darr(int size_m, T** arr)
{
    for(int i = 0; i < size_m; i++)
    {
        dealloc_1darr(arr[i]);
    }
    free(arr);
}

template <class T>
void dealloc_3darr(int size_m, int size_n, T*** arr)
{
    for(int i = 0; i < size_m; i++)
    {
        dealloc_2darr(size_n, arr[i]);
    }
    free(arr);
}

template <class T>
T* deep_copy_1darr(int size_m, T* arr)
{
    T* out = alloc_1darr<T>(size_m);
    for(int i = 0; i < size_m; i++)
    {
        out[i] = arr[i];
    }
    return out;
}

template <class T>
T** deep_copy_2darr(int size_m, int size_n, T** arr)
{
    T** out = alloc_2darr<T>(size_m, size_n);
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
    T*** out = alloc_3darr<T>(size_m, size_n, size_p);
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

template double* alloc_1darr<double>(int size_m);
template double** alloc_2darr<double>(int size_m, int size_n);
template double*** alloc_3darr<double>(int size_m, int size_n, int size_p);
template void dealloc_1darr<double>(double* arr);
template void dealloc_2darr<double>(int size_m, double** arr);
template void dealloc_3darr<double>(int size_m, int size_n, double*** arr);
template double* deep_copy_1darr<double>(int size_m, double* arr);
template double** deep_copy_2darr<double>(int size_m, int size_n, double** arr);
template double*** deep_copy_3darr<double>(int size_m, int size_n, int size_p, double*** arr);

template int* alloc_1darr<int>(int size_m);
template int** alloc_2darr<int>(int size_m, int size_n);
template int*** alloc_3darr<int>(int size_m, int size_n, int size_p);
template void dealloc_1darr<int>(int* arr);
template void dealloc_2darr<int>(int size_m, int** arr);
template void dealloc_3darr<int>(int size_m, int size_n, int*** arr);
template int* deep_copy_1darr<int>(int size_m, int* arr);
template int** deep_copy_2darr<int>(int size_m, int size_n, int** arr);
template int*** deep_copy_3darr<int>(int size_m, int size_n, int size_p, int*** arr);
