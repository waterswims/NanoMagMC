#ifndef _ARR_ALLOC_
#define _ARR_ALLOC_

template <class T> T* alloc_1darr(int size_m, T* arr);

template <class T> T** alloc_2darr(int size_m, int size_n, T** arr);

template <class T> T*** alloc_3darr(int size_m, int size_n, int size_p, T*** arr);

template <class T> void dealloc_1darr(int size_m, T* arr);

template <class T> void dealloc_2darr(int size_m, int size_n, T** arr);

template <class T> void dealloc_3darr(int size_m, int size_n, int size_p, T*** arr);

template <class T> deep_copy_1darr(int size_m, T* arr);

template <class T> deep_copy_2darr(int size_m, int size_n, T** arr);

template <class T> deep_copy_3darr(int size_m, int size_n, T*** arr);

#endif
