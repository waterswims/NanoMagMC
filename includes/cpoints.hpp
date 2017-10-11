#ifndef _CHECK
#define _CHECK

#include <vector>
#include <fstream>
#include <cstring>

using namespace std;

void print_cval(fstream &stream, string fname, double cval);

void print_cval(fstream &stream, string fname, int cval);

void print_clist(fstream &stream, string fname, double* clist, int l);

int read_cval(fstream &stream, string fname, int* cvals);

int read_clist(fstream &stream, string fname, double** clist);

string cpointname(string prefix, int rank, double size, int dist, char latt, char ham, double field);

bool exists(string &name);

#endif
