#include "mpifuncs.hpp"
#include "array_alloc.hpp"

#include <iostream>

void send1darr(int sizex, double* arr, int target)
{
    MPI_Ssend(arr, sizex, MPI_DOUBLE, target, 0, MPI_COMM_WORLD);
}

void send2darr(int sizex, int sizey, double** arr, int target)
{
    for(int i = 0; i < sizex; i++)
    {
        MPI_Ssend(arr[i], sizey, MPI_DOUBLE, target, 0, MPI_COMM_WORLD);
    }
}

void send3darr(int sizex, int sizey, int sizez, double*** arr, int target)
{
    for(int i = 0; i < sizex; i++)
    {
        for(int j = 0; j < sizey; j++)
        {
            MPI_Ssend(arr[i][j], sizez, MPI_DOUBLE, target, 0, MPI_COMM_WORLD);
        }
    }
}

double* recv1darr(int sizex, int target)
{
    MPI_Status status;
    double* out = alloc_1darr<double>(sizex);
    MPI_Recv(out, sizex, MPI_DOUBLE, target, 0, MPI_COMM_WORLD, &status);
    return out;
}

double** recv2darr(int sizex, int sizey, int target)
{
    MPI_Status status;
    double** out = alloc_2darr<double>(sizex, sizey);
    for(int i = 0; i < sizex; i++)
    {
        MPI_Recv(out[i], sizey, MPI_DOUBLE, target, 0, MPI_COMM_WORLD, &status);
    }
    return out;
}

double*** recv3darr(int sizex, int sizey, int sizez, int target)
{
    MPI_Status status;
    double*** out = alloc_3darr<double>(sizex, sizey, sizez);
    for(int i = 0; i < sizex; i++)
    {
        for(int j = 0; j < sizey; j++)
        {
            MPI_Recv(out[i][j], sizez, MPI_DOUBLE, target, 0, MPI_COMM_WORLD, &status);
        }
    }
    return out;
}

double* redc1darr(int sizex, double* arr, int root)
{
    double* out = alloc_1darr<double>(sizex);
    MPI_Reduce(arr, out, sizex, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
    return out;
}

double** redc2darr(int sizex, int sizey, double** arr, int root)
{
    double** out = alloc_2darr<double>(sizex, sizey);
    for(int i = 0; i < sizex; i++)
    {
        MPI_Reduce(arr[i], out[i], sizey, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
    }
    return out;
}

double*** redc3darr(int sizex, int sizey, int sizez, double*** arr, int root)
{
    double*** out = alloc_3darr<double>(sizex, sizey, sizez);
    for(int i = 0; i < sizex; i++)
    {
        for(int j = 0; j < sizey; j++)
        {
            MPI_Reduce(arr[i][j], out[i][j], sizez, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
        }
    }
    return out;
}

vector<int> gather1darr(int sizex, int* arr, int root, int num_ranks)
{
    vector<int> out(sizex*num_ranks);
    int* outarr = alloc_1darr<int>(sizex*num_ranks);
    MPI_Allgather(arr, sizex, MPI_INT, outarr, sizex, MPI_INT, MPI_COMM_WORLD);
    for(int i=0; i < sizex*num_ranks; i++)
    {
        out[i] = outarr[i];
    }
    dealloc_1darr<int>(outarr);
    return out;
}

vector<vector<double> > gather2darr(int sizex, int sizey, double** arr, int root, int num_ranks)
{
    vector<vector<double> > out;
    vector<double> temp(sizey*num_ranks);
    double** outarr = alloc_2darr<double>(sizex, sizey*num_ranks);
    for(int i=0; i < sizex; i++)
    {
        MPI_Allgather(arr[i], sizey, MPI_DOUBLE, outarr[i], sizey, MPI_DOUBLE, MPI_COMM_WORLD);
    }
    for(int i=0; i < sizex; i++)
    {
        for(int j=0; j < sizey*num_ranks; j++)
        {
            temp[j] = outarr[i][j];
        }
        out.push_back(temp);
    }
    dealloc_2darr<double>(sizex, outarr);
    return out;
}
