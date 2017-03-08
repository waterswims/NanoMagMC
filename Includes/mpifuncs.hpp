#ifndef _MPIFUNCS
#define _MPIFUNCS

#include <mpi.h>
#include <vector>

using namespace std;

void send1darr(int sizex, double* arr, int target);

void send2darr(int sizex, int sizey, double** arr, int target);

void send3darr(int sizex, int sizey, int sizez, double*** arr, int target);

double* recv1darr(int sizex, int target);

double** recv2darr(int sizex, int sizey, int target);

double*** recv3darr(int sizex, int sizey, int sizez, int target);

double* redc1darr(int sizex, double* arr, int root);

double** redc2darr(int sizex, int sizey, double** arr, int root);

double*** redc3darr(int sizex, int sizey, int sizez, double*** arr, int root);

vector<int> gather1darr(int sizex, int* arr, int root, int num_ranks);

vector<vector<double> > gather2darr(int sizex, int sizey, double** arr, int root, int num_ranks);

#endif
