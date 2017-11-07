#ifndef _FUNC
#define _FUNC

#include <vector>
#include <fstream>
#include <cstring>

using namespace std;

double mean(vector<double> &oY);

double std_dev(vector<double> &x);

double norm(vector<double> vals);

double sum(vector<double> &oY);

int sum(vector<int> &oY);

void AtoLn(double amean, double asd, double &lmean, double &lsd);

int mod(int a, int b);

#endif
