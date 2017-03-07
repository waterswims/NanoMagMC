#include "spin_type.hpp"
#include "mklrand.h"

#include <iostream>
#include <cmath>

extern mkl_irand st_rand_int;
extern mkl_drand st_rand_double;

const double pi = 3.141592653589793;

ising_spin::ising_spin(){spin=0; zero_flag=true;}

ising_spin::ising_spin(const ising_spin& other)
{
    spin = other.spin;
    zero_flag = other.zero_flag;
}

void ising_spin::rand_spin()
{
    spin = st_rand_int.gen() * 2 - 1;
    zero_flag = false;
}

void ising_spin::zero_spin(){spin=0; zero_flag=true;}

int ising_spin::i_access(){return spin;}

void ising_spin::set_spin(int in){spin=in;}

bool ising_spin::is_zero()
{
    // if (spin == 0)
    // {
    //     return true;
    // }
    // else
    // {
    //     return false;
    // }
    return zero_flag;
}

ising_spin& ising_spin::operator=(const ising_spin& other)
{
    spin = other.spin;
    zero_flag = other.zero_flag;
}

void ising_spin::print()
{
    cout << spin;
}

heis_spin::heis_spin()
{
    spinx = 0;
    spiny = 0;
    spinz = 0;
    zero_flag = true;
}

heis_spin::heis_spin(const heis_spin& other)
{
    spinx = other.spinx;
    spiny = other.spiny;
    spinz = other.spinz;
    zero_flag = other.zero_flag;
}

void heis_spin::rand_spin()
{
    zero_flag = false;

    // Marsaglia
    // bool notfin = true;
    // while(notfin)
    // {
    //     double rand1 = 2*st_rand_double.gen()-1;
    //     double rand2 = 2*st_rand_double.gen()-1;
    //     double norm2 = pow(rand1, 2) + pow(rand2, 2);
    //     if (norm2 <= 1)
    //     {
    //         notfin = false;
    //         double temp = 2 * pow(1 - norm2, 0.5);
    //         spins[0] = rand1 * temp;
    //         spins[1] = rand2 * temp;
    //         spins[3] = 1 - 2 * norm2;
    //     }
    // }

    // rescaling random direction generation
    // bool notfin = true;
    // while(notfin)
    // {
    //     double randx = 2*st_rand_double.gen()-1;
    //     double randy = 2*st_rand_double.gen()-1;
    //     double randz = 2*st_rand_double.gen()-1;
    //     double norm2 = pow(randx, 2) + pow(randy, 2) + pow(randz, 2);
    //     if (norm2 <= 1 && norm2 >= 1e-10)
    //     {
    //         notfin = false;
    //         double norm = pow(norm2, 0.5);
    //         spins[0] = randx / norm;
    //         spins[1] = randy / norm;
    //         spins[2] = randz / norm;
    //     }
    // }

    // stardard random direction generation
    double phi = st_rand_double.gen()*2*pi;
    double cthet = 2*st_rand_double.gen()-1;
    double sthet = pow(1 - pow(cthet, 2), 0.5);
    spinx = cos(phi)*sthet;
    spiny = sin(phi)*sthet;
    spinz = cthet;
}

void heis_spin::zero_spin()
{
    spinx = 0;
    spiny = 0;
    spinz = 0;
    zero_flag = true;
}

void heis_spin::change_spin(spin_type* in)
{
    in->spin_access(spinx, spiny, spinz);
    zero_flag = in->is_zero();
}

bool heis_spin::is_zero()
{
    return zero_flag;
}

heis_spin& heis_spin::operator=(const heis_spin& other)
{
    spinx = 0;
    spiny = 0;
    spinz = 0;
    zero_flag = other.zero_flag;
}

void heis_spin::print()
{
    cout << "[" << spinx << "," << spiny << "," << spinz << "]";
}
