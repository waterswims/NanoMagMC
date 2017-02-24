#include "hamiltonian_new.hpp"
#include "functions.h"

#include <iostream>
#include <cstdlib>
#include <algorithm>

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
    vector<int> pos(dim);
    for (vector<int>::iterator it = pos.begin(); it != pos.end(); it++)
    {
        *it = start;
    }
    bool finished = false;
    int curr;
    int H_sum(0), J_sum(0);
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
    lattice->i_adjacent(position, adj);
    for(int i = 0; i < dim*2; i++)
    {
        sum += adj[i];
    }

    double dE = dEH + 2 * J * val * sum;
    return dE;
}

void ham_ising::init_dim(field_type* field)
{
    dim = field->get_dim();
    adj = alloc_1darr<int>(dim*2);
}

///////////////////////////////////
// Heis model hamiltonian
///////////////////////////////////

ham_heis::ham_heis(double Hin, double Jin)
{
    H.resize(4);
    J.resize(4);
    H[0] = 0;
    H[1] = 0;
    H[2] = Hin;
    J[0] = Jin;
    J[1] = Jin;
    J[2] = Jin;
    vsum.resize(4);
    curr.resize(4);
    H_sum.resize(4);
    J_sum.resize(4);
    test.resize(3);
}

ham_heis::ham_heis(ham_type& other)
{
    H = other.get_Hs();
    J = other.get_Js();
    vsum.resize(4);
    curr.resize(4);
    H_sum.resize(4);
    J_sum.resize(4);
    test.resize(3);
}

ham_heis& ham_heis::operator=(ham_type& other)
{
    H = other.get_Hs();
    J = other.get_Js();
    vsum.resize(4);
    curr.resize(4);
    H_sum.resize(4);
    J_sum.resize(4);
    test.resize(3);

    return *this;
}

double ham_heis::calc_E(field_type* lattice)
{
    int sum=0;
    int start = 0;
    if(lattice->get_perio())
    {
        start++;
    }
    pos.resize(dim);
    for (vector<int>::iterator it = pos.begin(); it != pos.end(); it++)
    {
        *it = start;
    }
    bool finished = false;
    H_sum[0] = 0;
    H_sum[1] = 0;
    H_sum[2] = 0;
    J_sum[0] = 0;
    J_sum[1] = 0;
    J_sum[2] = 0;
    int arrsize = dim*2;
    while (!finished)
    {
        lattice->h_adjacent(pos, adj);
        lattice->h_next(finished, pos, curr);
        #pragma simd
        for (int j = 0; j < 4; j++)
        {
            H_sum[j] += curr[j];
            for (int i = 0; i < arrsize; i++)
            {
                J_sum[j] += curr[j]*adj[j][i];
            }
        }
    }
    #pragma simd
    for (int j = 0; j < 4; j++)
    {
        H_sum[j] = H_sum[j] * H[j];
        J_sum[j] = J_sum[j] * J[j];
    }
    double E = -(H_sum[0] + H_sum[1] + H_sum[2]) - 0.5*(J_sum[0] + J_sum[1] + J_sum[2]);

    return E;
}

vector<double> ham_heis::calc_M(field_type* lattice)
{
    #pragma simd
    for (int i = 0; i < 4; i++) {vsum[i] += 0;}
    int start = 0;
    if(lattice->get_perio())
    {
        start++;
    }
    pos.resize(dim);
    for (int i = 0; i < dim; i++) {pos[i] += 0;}
    bool finished = false;
    while (!finished)
    {
        lattice->h_next(finished, pos, curr);
        #pragma simd
        for (int i = 0; i < 4; i++) {vsum[i] += curr[i];}
    }
    return vsum;
}

vector<double> ham_heis::calc_subM(field_type* lattice, int subnumber)
{
    #pragma simd
    for (int i = 0; i < 4; i++) {vsum[i] += 0;}
    int start = 0;
    if(lattice->get_perio())
    {
        start++;
    }
    pos.resize(dim);
    for (int i = 0; i < dim; i++) {pos[i] += 0;}
    bool finished = false;
    int possum = 0;
    while (!finished)
    {
        possum = sum(pos);
        lattice->h_next(finished, pos, curr);
        if (possum%2 == subnumber)
        {
            #pragma simd
            for (int i = 0; i < 4; i++) {vsum[i] += curr[i];}
        }
    }
    return vsum;
}

double ham_heis::dE(field_type* lattice, vector<int>& position)
{
    rand_spin_h(test[0], test[1], test[2]);
    lattice->h_access(position, curr);
    double cmp1 = curr[0] - test[0];
    double cmp2 = curr[1] - test[1];
    double cmp3 = curr[2] - test[2];
    double dEH = H[0] * cmp1 + H[1] * cmp2 + H[2] * cmp3;

    int arrsize = dim*2;
    lattice->h_adjacent(position, adj);
    #pragma simd
    for(int j = 0; j < 4; j++)
    {
        vsum[j] = 0;
        for(int i = 0; i < arrsize; i++)
        {
            vsum[j] += adj[j][i];
        }
    }

    double dE = dEH + (J[0] * cmp1 * vsum[0] + J[1] * cmp2 * vsum[1] + J[2] * cmp3 * vsum[2]);

    return dE;
}

void ham_heis::init_dim(field_type* field)
{
    dim = field->get_dim();
    adj = alloc_2darr<double>(4, dim*2);
}

///////////////////////////////////
//  FePt hamiltonian
///////////////////////////////////

ham_FePt::ham_FePt()
{
    this->read_Js();
    vsum.resize(3);
    curr.resize(3);
    adj_curr.resize(3);
    d_ijs.resize(3);
    Js.resize(3);
    test.resize(3);
}

ham_FePt::ham_FePt(ham_type& other)
{
    this->read_Js();
    vsum.resize(3);
    curr.resize(3);
    adj_curr.resize(3);
    d_ijs.resize(3);
    Js.resize(3);
    test.resize(3);
}

ham_FePt& ham_FePt::operator=(ham_type& other)
{
    this->read_Js();
    vsum.resize(3);
    curr.resize(3);
    adj_curr.resize(3);
    d_ijs.resize(3);
    Js.resize(3);
    test.resize(3);
    return *this;
}

double ham_FePt::calc_E(field_type* lattice)
{
    double Jx_sum = 0;
    double Jy_sum = 0;
    double d2_sum = 0;
    double d0_sum = 0;
    int t_size = lattice->get_totsize();
    int i_size = lattice->get_insize();
    int start = (t_size-i_size)/2;
    int fin = start+i_size;

    pos2.resize(dim);
    for (vector<int>::iterator it = pos2.begin(); it != pos2.end(); it++)
    {
        *it = start;
    }
    pos.resize(dim);
    int nN = dxs.size();
    bool finished = false;
    while (!finished)
    {
        lattice->h_next(finished, pos2, curr);
        // All neighbour interactions
        for(int i = 0; i < nN; i++)
        {
            // if version
            // pos[0] = pos2[0] + dxs[i];
            // pos[1] = pos2[1] + dys[i];
            // pos[2] = pos2[2] + dzs[i];
            //
            // if(this->check_pos(start, fin))
            // {
            //     lattice->h_access(pos, adj_curr);
            //
            //     Jx_sum += curr[0] * adj_curr[0] * Js[i];
            //     Jy_sum += curr[1] * adj_curr[1] * Js[i];
            //
            //     d2_sum += curr[2] * adj_curr[2] * d_ijs[i];
            // }

            // ternary version
            pos[0] = max(0, (pos2[0] + dxs[i])%t_size);
            pos[1] = max(0, (pos2[1] + dys[i])%t_size);
            pos[2] = max(0, (pos2[2] + dzs[i])%t_size);
            lattice->h_access(pos, adj_curr);
            Jx_sum += curr[0] * adj_curr[0] * Js[i];
            Jy_sum += curr[1] * adj_curr[1] * Js[i];
            d2_sum += curr[2] * adj_curr[2] * d_ijs[i];
        }

        // 1 ion anisotropy
        d0_sum += curr[2]*curr[2];
    }

    // totals
    double E = -(d0 * d0_sum + 0.5 * (Jx_sum + Jy_sum + d2_sum));

    return E;
}

double ham_FePt::dE(field_type* lattice, vector<int>& position)
{
    double Jx_sum = 0;
    double Jy_sum = 0;
    double d2_sum = 0;
    int t_size = lattice->get_totsize();
    // int i_size = lattice->get_insize();
    // int start = (t_size-i_size)/2;
    // int fin = start+i_size;

    rand_spin_h(test[0], test[1], test[2]);
    lattice->h_access(position, curr);
    double cmp1 = curr[0] - test[0];
    double cmp2 = curr[1] - test[1];
    double cmp3 = curr[2] - test[2];

    int nN = dxs.size();
    // All neighbour interactions
    for(int i = 0; i < nN; i++)
    {
        // if version
        // pos[0] = position[0] + dxs[i];
        // pos[1] = position[1] + dys[i];
        // pos[2] = position[2] + dzs[i];
        //
        // if(this->check_pos(start, fin))
        // {
        //     lattice->h_access(pos, adj_curr);
        //
        //     Jx_sum += adj_curr[0] * Js[i];
        //     Jy_sum += adj_curr[1] * Js[i];
        //
        //     d2_sum += adj_curr[2] * d_ijs[i];
        // }

        // ternary version
        pos[0] = max(0, (position[0] + dxs[i])%t_size);
        pos[1] = max(0, (position[1] + dys[i])%t_size);
        pos[2] = max(0, (position[2] + dzs[i])%t_size);
        lattice->h_access(pos, adj_curr);
        Jx_sum += curr[0] * adj_curr[0] * Js[i];
        Jy_sum += curr[1] * adj_curr[1] * Js[i];
        d2_sum += curr[2] * adj_curr[2] * d_ijs[i];
    }

    // 1 ion anisotropy
    double d0_sum = curr[2]*curr[2] - test[2]*test[2];

    // totals
    double dE = (d0 * d0_sum + cmp1 * Jx_sum + cmp2 * Jy_sum + cmp3 * d2_sum);

    return dE;
}

void ham_FePt::read_Js()
{
    ifstream Jstream, d_ijstream;
    Jstream.open("Includes/Js/FePt_cut_1e-3.txt");
    d_ijstream.open("Includes/Js/FePt_2ion_cut_1e-3.txt");
    int icurr;
    double dcurr;
    while(Jstream >> icurr)
    {
        dxs.push_back(icurr);
        Jstream >> icurr;
        dys.push_back(icurr);
        Jstream >> icurr;
        dzs.push_back(icurr);
        Jstream >> dcurr;
        Js.push_back(dcurr);

        d_ijstream >> icurr;
        d_ijstream >> icurr;
        d_ijstream >> icurr;
        d_ijstream >> dcurr;
        d_ijs.push_back(dcurr);
    }
    Jstream.close();
    d_ijstream.close();

    d0 = 8.42603732e-5;
}

bool ham_FePt::check_pos(int start, int fin)
{
    if(pos[0] < start) return false;
    if(pos[1] < start) return false;
    if(pos[2] < start) return false;
    if(pos[0] >= fin) return false;
    if(pos[1] >= fin) return false;
    if(pos[2] >= fin) return false;
    return true;
}
