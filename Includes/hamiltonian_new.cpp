#include "hamiltonian_new.hpp"
#include "array_alloc.hpp"
#include "functions.h"

#include <iostream>
#include <cstdlib>

///////////////////////////////////
// Ising model hamiltonian
///////////////////////////////////

ham_ising::ham_ising(ham_type& other)
{
    H = other.get_J();
    J = other.get_H();
}

ham_ising& ham_ising::operator=(ham_type& other)
{
    H = other.get_J();
    J = other.get_H();
    return *this;
}

double ham_ising::calc_E(field_type* lattice)
{
    int sum=0;
    int start = 0;
    if(lattice->get_perio())
    {
        start++;
    }
    int dim = lattice->get_dim();
    vector<int> pos(dim);
    for (vector<int>::iterator it = pos.begin(); it != pos.end(); it++)
    {
        *it = start;
    }
    bool finished = false;
    int curr;
    int H_sum(0), J_sum(0);
    adj = alloc_1darr<int>(dim*2);
    while (!finished)
    {
        lattice->i_adjacent(pos, adj);
        lattice->i_next(finished, pos, curr);
        H_sum += curr;
        for (int i = 0; i < 2*dim; i++)
        {
            J_sum += curr*adj[i];
        }
    }
    dealloc_1darr<int>(dim*2, adj);

    double E = -H * H_sum - 0.5 * J * J_sum;
    return E;
}

vector<double> ham_ising::calc_M(field_type* lattice)
{
    int sum=0;
    int start = 0;
    if(lattice->get_perio())
    {
        start++;
    }
    int dim = lattice->get_dim();
    vector<int> pos(dim);
    for (vector<int>::iterator it = pos.begin(); it != pos.end(); it++)
    {
        *it = start;
    }
    bool finished = false;
    int curr;
    while (!finished)
    {
        lattice->i_next(finished, pos, curr);
        sum += curr;
    }
    vector<double> mag(1);
    mag[0] = sum;
    return mag;
}

vector<double> ham_ising::calc_subM(field_type* lattice, int subnumber)
{
    if (subnumber > 1)
    {
        cerr << "subnumber is incorrect" << endl;
        exit(304);
    }
    int tsum=0;
    int start = 0;
    if(lattice->get_perio())
    {
        start++;
    }
    int dim = lattice->get_dim();
    vector<int> pos(dim);
    for (vector<int>::iterator it = pos.begin(); it != pos.end(); it++)
    {
        *it = start;
    }
    bool finished = false;
    int possum = 0;
    int curr;
    while (!finished)
    {
        possum = sum(pos);
        lattice->i_next(finished, pos, curr);
        if (possum%2 == subnumber)
        {
            tsum += curr;
        }
    }

    vector<double> mag(1);
    mag[0] = tsum;
    return mag;
}

double ham_ising::dE(field_type* lattice, vector<int>& position)
{
    int val;
    lattice->i_access(position, val);
    double dEH = 2 * H * val;
    int sum = 0;
    int dim = lattice->get_dim();
    adj = alloc_1darr<int>(dim*2);
    lattice->i_adjacent(position, adj);
    for(int i = 0; i < dim*2; i++)
    {
        sum += adj[i];
    }
    dealloc_1darr<int>(dim*2, adj);

    double dE = dEH + 2 * J * val * sum;
    return dE;
}
