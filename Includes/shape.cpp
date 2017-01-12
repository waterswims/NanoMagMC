#include "shape.hpp"
#include "mklrand.h"
#include <cmath>
#include <iostream>

extern mkl_drand st_rand_double;

weibull::weibull(double rin, double bin)
{
    beta = bin;
    r0 = rin / tgamma(1. + 1. / beta);
}

bool weibull::check(vector<int> Is, int l_size)
{
    double centre = double(l_size - 1) / 2.;
    double dist2 = 0;
    for(int i = 0; i < Is.size(); i++)
    {
        dist2 += pow((Is[i]-centre), 2);
    }
	double dist = pow(dist2, 0.5);
	double test = exp(-pow((dist/r0), beta));
	if(st_rand_double.gen() < test)
	{
		return true;
	}
	else
	{
		return false;
	}
}

weibull& weibull::operator=(shape_type& other)
{
    r0 = other.get_r0();
    beta = other.get_beta();
    return *this;
}
