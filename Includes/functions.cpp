#include "functions.h"
// #include "boost/assign.hpp"
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

double mean(double oY[]){
    int N = arr_len(oY);
    //cout << N << endl;
    double s(0);
    for (int i=0; i < N; i++){
        s += oY[i];
    }
    //cout << s << endl;
    return s / N;
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

double mag_sus(vector<double> mag, double T, vector<int> nums){
    int N = mag.size();
    double s_mag = sum(mag);
    int s_num = sum(nums);
    double s_m2ovn = 0;
    for (int i(0); i < N; i++)
    {
        s_m2ovn += pow(mag[i], 2) / nums[i];
    }
    double ms = (s_m2ovn / s_num - pow(s_mag / s_num, 2)) / T;
    return ms;
}

double mag_sus(vector<double> magx, vector<double> magy, vector<double> magz, double T, vector<int> nums)
{
    int N = magx.size();
    double s_magx = sum(magx);
    double s_magy = sum(magy);
    double s_magz = sum(magz);
    int s_num = sum(nums);
    double s_m2ovn = 0;
    vector<double> mvec(3);
    for (int i(0); i < N; i++)
    {
        mvec[0] = magx[i];
        mvec[1] = magy[i];
        mvec[2] = magz[i];
        // mvec = boost::assign::list_of(magx[i])(magy[i])(magz[i]).convert_to_container<vector<double> >();
        s_m2ovn += pow(norm(mvec), 2) / nums[i];
    }
    mvec[0] = s_magx;
    mvec[1] = s_magy;
    mvec[2] = s_magz;
    // mvec = boost::assign::list_of(s_magx)(s_magy)(s_magz).convert_to_container<vector<double> >();
    double ms = (s_m2ovn / s_num - pow(norm(mvec) / s_num, 2)) / T;
    return ms;
}

double heat_cap(vector<double> ener, double T, vector<int> nums){
    int N = ener.size();
    double s_ener = sum(ener);
    int s_num = sum(nums);
    double s_e2ovn = 0;
    for (int i(0); i < N; i++)
    {
        s_e2ovn += pow(ener[i], 2) / nums[i];
    }
    double hc = (s_e2ovn / s_num - pow(s_ener / s_num, 2)) / pow(T,2);
    return hc;
}

double binders(vector<double> mag, int size)
{
	int N = mag.size();
	vector<double> mag_sq;
	vector<double> mag_qu;
	for (int i(0); i < N; i++)
	{
		mag[i] = mag[i] / size;
		mag_sq.push_back(pow(mag[i], 2));
		mag_qu.push_back(pow(mag_sq[i], 2));
	}
	double binder = 1 - mean(mag_qu) / (3*pow(mean(mag_sq), 2));
	return binder;
}

int arr_len(double oY[])
{
            return sizeof(oY) / sizeof(oY[0]);
}

double norm(double val, int norm_typ)
{
            if(norm_typ==2)
            {
                        return abs(val);
            }
            cout << "Norm Type not defined, returning zero" << endl;
            return 0.0;
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

void printv(vector<int> &oX)
{
    cout << "[";
    for(vector<int>::iterator it = oX.begin(); it != oX.end(); it++)
    {
        cout << *it << " ";
    }
    cout << "]" << endl;
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
