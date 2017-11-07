#include "../includes/functions.hpp"
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <sstream>

using namespace std;

double mean(vector<double> &oY){
    double s(0);
    for (int i=0; i < oY.size(); i++){
        s += oY[i];
    }
    return s / oY.size();
}

double std_dev(vector<double> &x){
    double av_x = mean(x);
    int size = x.size();
    double sum = 0;
    for(int i(0); i < size; i++){
        sum += pow(x[i] - av_x, 2);
    }
    double var = sum/(size - 1.);
    return pow(var, 0.5);
}

double norm(vector<double> vals)
{
    double s = 0;
    for(vector<double>::iterator it = vals.begin(); it != vals.end(); it++)
    {
        s += pow(*it, 2);
    }
    return pow(s, 0.5);
}

double sum(vector<double> &oY)
{
            double s = 0;
            for(vector<double>::iterator it = oY.begin(); it != oY.end(); it++)
            {
                        s += *it;
            }
            return s;
}

int sum(vector<int> &oY)
{
    int s = 0;
    for(vector<int>::iterator it = oY.begin(); it != oY.end(); it++)
    {
                s += *it;
    }
    return s;
}

void AtoLn(double amean, double asd, double &lmean, double &lsd)
{
    double av = pow(asd, 2);
    double am2 = pow(amean, 2);
    lmean = log(am2 / pow(av+am2, 0.5));
    lsd = pow(log(1 + av / am2), 0.5);
}

int mod(int a, int b)
{
    int r = a % b;
    return r < 0 ? r + b : r;
}
