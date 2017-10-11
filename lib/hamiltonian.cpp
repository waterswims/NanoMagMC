#include "../includes/hamiltonian.hpp"
#include "../includes/functions.hpp"

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
    int start = (lattice->get_totsize() - lattice->get_insize()) / 2;
    vector<int> pos(dim, start);
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
    int start = (lattice->get_totsize() - lattice->get_insize()) / 2;
    vector<int> pos(dim, start);
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
    int start = (lattice->get_totsize() - lattice->get_insize()) / 2;
    vector<int> pos(dim, start);
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
    int start = (lattice->get_totsize() - lattice->get_insize()) / 2;
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
    for (int i = 0; i < 4; i++) {vsum[i] = 0;}
    int start = (lattice->get_totsize() - lattice->get_insize()) / 2;
    pos.resize(dim);
    for (vector<int>::iterator it = pos.begin(); it != pos.end(); it++)
    {
        *it = start;
    }
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
    for (int i = 0; i < 4; i++) {vsum[i] = 0;}
    int start = (lattice->get_totsize() - lattice->get_insize()) / 2;
    pos.resize(dim);
    for (vector<int>::iterator it = pos.begin(); it != pos.end(); it++)
    {
        *it = start;
    }
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
    // field->alloc_pos(1);
}

///////////////////////////////////
//  FePt hamiltonian
///////////////////////////////////

ham_FePt::ham_FePt()
{
    this->read_Js();
    vsum.resize(4);
    curr.resize(4);
    adj_curr.resize(4);
    test.resize(4);
    H.resize(4);
    H[2] = 0;
}

ham_FePt::ham_FePt(double Hin)
{
    this->read_Js();
    vsum.resize(4);
    curr.resize(4);
    adj_curr.resize(4);
    test.resize(4);
    H.resize(4);
    H[2] = Hin;
}

ham_FePt::ham_FePt(ham_type& other)
{
    this->read_Js();
    vsum.resize(4);
    curr.resize(4);
    adj_curr.resize(4);
    test.resize(4);
    H = other.get_Hs();
}

ham_FePt& ham_FePt::operator=(ham_type& other)
{
    this->read_Js();
    vsum.resize(4);
    curr.resize(4);
    adj_curr.resize(4);
    test.resize(4);
    H = other.get_Hs();
    return *this;
}

double ham_FePt::calc_E(field_type* lattice)
{
    double Jx_sum = 0;
    double Jy_sum = 0;
    double d2_sum = 0;
    double d0_sum = 0;
    double H_sum = 0;
    int t_size = lattice->get_totsize();
    int i_size = lattice->get_insize();
    int start = (t_size-i_size)/2;
    pos2.resize(dim);
    for (vector<int>::iterator it = pos2.begin(); it != pos2.end(); it++)
    {
        *it = start;
    }
    int nN = dxs.size();
    bool finished = false;
    while (!finished)
    {
        lattice->h_access(pos2, curr);
        lattice->h_arb_adj(pos2, dxs, dys, dzs, adj, nN);
        #pragma simd
        for(int i = 0; i < nN; i++)
        {
            Jx_sum += curr[0] * adj[0][i] * Js[i];
            Jy_sum += curr[1] * adj[1][i] * Js[i];
            d2_sum += curr[2] * adj[2][i] * d_ijs[i];
        }
        lattice->next(finished, pos2);
        // 1 ion anisotropy
        d0_sum += curr[2]*curr[2];
        H_sum += curr[2];
    }

    // totals
    double E = -(H[2] * H_sum + d0 * d0_sum + 0.5 * (Jx_sum + Jy_sum + d2_sum));

    return E;
}

double ham_FePt::dE(field_type* lattice, vector<int>& position)
{
    double Jx_sum = 0;
    double Jy_sum = 0;
    double d2_sum = 0;
    int t_size = lattice->get_totsize();

    rand_spin_h(test[0], test[1], test[2]);
    lattice->h_access(position, curr);
    double cmp1 = curr[0] - test[0];
    double cmp2 = curr[1] - test[1];
    double cmp3 = curr[2] - test[2];
    int nN = dxs.size();
    pos.resize(dim);
    lattice->h_arb_adj(position, dxs, dys, dzs, adj, nN);
    #pragma simd
    for(int i = 0; i < nN; i++)
    {
        Jx_sum += adj[0][i] * Js[i];
        Jy_sum += adj[1][i] * Js[i];
        d2_sum += adj[2][i] * d_ijs[i];
    }

    // 1 ion anisotropy
    double d0_sum = curr[2]*curr[2] - test[2]*test[2];

    // totals
    double dE = (H[2] * cmp3 + d0 * d0_sum + cmp1 * Jx_sum + cmp2 * Jy_sum + cmp3 * d2_sum);

    return dE;
}

void ham_FePt::read_Js()
{
    ifstream Jstream, d_ijstream;
    string Jname = "Includes/Js/FePt_nocut.txt";
    string dname = "Includes/Js/FePt_2ion_nocut.txt";
    Jstream.open(Jname.c_str());
    d_ijstream.open(dname.c_str());
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

void ham_FePt::init_dim(field_type* field)
{
    dim = field->get_dim();
    adj = alloc_2darr<double>(4, dxs.size());
    // field->alloc_pos(dxs.size());
}

///////////////////////////////////
//  Skyrmion hamiltonian
///////////////////////////////////

void ham_skyrm::set_dirs()
{
    dirs[0] = 0;
    dirs[1] = 0;
    dirs[2] = 1;
    dirs[3] = 1;
    dirs[4] = 2;
    dirs[5] = 2;
    mod[0] = -1;
    mod[1] = 1;
    mod[2] = -1;
    mod[3] = 1;
    mod[4] = -1;
    mod[5] = 1;
}

ham_skyrm::ham_skyrm(double Hin, double Jin, double Kin)
{
    H.resize(4);
    J.resize(4);
    H[0] = 0;
    H[1] = 0;
    H[2] = Hin;
    J[0] = Jin;
    J[1] = Jin;
    J[2] = Jin;
    K = Kin;
    vsum.resize(4);
    curr.resize(4);
    H_sum.resize(4);
    J_sum.resize(4);
    J2_sum.resize(4);
    cmp.resize(4);
    test.resize(3);
    this->set_dirs();
}

ham_skyrm::ham_skyrm(const ham_type& other)
{
    H = other.get_Hs();
    J = other.get_Js();
    K = other.get_K();
    vsum.resize(4);
    curr.resize(4);
    H_sum.resize(4);
    J_sum.resize(4);
    J2_sum.resize(4);
    test.resize(3);
    cmp.resize(4);
    this->set_dirs();
}

ham_skyrm& ham_skyrm::operator=(const ham_type& other)
{
    H = other.get_Hs();
    J = other.get_Js();
    K = other.get_K();
    vsum.resize(4);
    curr.resize(4);
    H_sum.resize(4);
    J_sum.resize(4);
    J2_sum.resize(4);
    cmp.resize(4);
    test.resize(3);
    this->set_dirs();

    return *this;
}

double ham_skyrm::calc_E(field_type* lattice)
{
    int sum=0;
    int start = (lattice->get_totsize() - lattice->get_insize()) / 2;
    double D_sum = 0, D2_sum = 0;
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
    J2_sum[0] = 0;
    J2_sum[1] = 0;
    J2_sum[2] = 0;
    int arrsize = dim*2;
    while (!finished)
    {
        lattice->h_adjacent(pos, adj);
        lattice->h_access(pos, curr);
        #pragma simd
        for (int j = 0; j < 4; j++)
        {
            H_sum[j] += curr[j];
            for (int i = 0; i < arrsize; i++)
            {
                J_sum[j] += curr[j]*adj[j][i];
            }
        }
        // #pragma simd
        for (int i = 0; i < arrsize; i++)
        {
            D_sum += mod[i]*(curr[(dirs[i]+1)%3]*adj[(dirs[i]+2)%3][i] -
                     curr[(dirs[i]+2)%3]*adj[(dirs[i]+1)%3][i]);
        }

        // lattice->h_2adjacent(pos, adj);
        // #pragma simd
        // for (int j = 0; j < 4; j++)
        // {
        //     for (int i = 0; i < arrsize; i++)
        //     {
        //         J2_sum[j] += curr[j]*adj[j][i];
        //     }
        // }
        // // #pragma simd
        // for (int i = 0; i < arrsize; i++)
        // {
        //     D2_sum += mod[i]*(curr[(dirs[i]+1)%3]*adj[(dirs[i]+2)%3][i] -
        //              curr[(dirs[i]+2)%3]*adj[(dirs[i]+1)%3][i]);
        // }
        //
        // lattice->next(finished, pos);
    }
    #pragma simd
    for (int j = 0; j < 4; j++)
    {
        H_sum[j] = H_sum[j] * H[j];
        J_sum[j] = J_sum[j] * J[j];
        J2_sum[j] = J2_sum[j] * J[j] / 16.;
    }

    double E = -(H_sum[0] + H_sum[1] + H_sum[2]) -
                0.5*(J_sum[0] + J_sum[1] + J_sum[2] + J2_sum[0] + J2_sum[1] +
                     J2_sum[2]) - 0.5*K*(D_sum + D2_sum/8.);

    return E;
}

double ham_skyrm::dE(field_type* lattice, vector<int>& position)
{
    double D_sum = 0, D2_sum = 0;
    rand_spin_h(test[0], test[1], test[2]);
    lattice->h_access(position, curr);
    cmp[0] = curr[0] - test[0];
    cmp[1] = curr[1] - test[1];
    cmp[2] = curr[2] - test[2];
    double dEH = H[0] * cmp[0] + H[1] * cmp[1] + H[2] * cmp[2];

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
    // #pragma simd
    for (int i = 0; i < arrsize; i++)
    {
        D_sum += mod[i]*(cmp[(dirs[i]+1)%3]*adj[(dirs[i]+2)%3][i] -
                 cmp[(dirs[i]+2)%3]*adj[(dirs[i]+1)%3][i]);
    }

    // lattice->h_2adjacent(position, adj);
    // #pragma simd
    // for(int j = 0; j < 4; j++)
    // {
    //     J2_sum[j] = 0;
    //     for(int i = 0; i < arrsize; i++)
    //     {
    //         J2_sum[j] += adj[j][i];
    //     }
    // }
    // // #pragma simd
    // for (int i = 0; i < arrsize; i++)
    // {
    //     D2_sum += mod[i]*(cmp[(dirs[i]+1)%3]*adj[(dirs[i]+2)%3][i] -
    //              cmp[(dirs[i]+2)%3]*adj[(dirs[i]+1)%3][i]);
    // }

    double dEJ = J[0] * cmp[0] * vsum[0] + J[1] * cmp[1] * vsum[1] + J[2] * cmp[2] * vsum[2];

    double dEJ2 = (J[0] * cmp[0] * J2_sum[0] + J[1] * cmp[1] * J2_sum[1] + J[2] * cmp[2] * J2_sum[2]) / 16.;

    double dE = dEH + dEJ + dEJ2 + K * (D_sum + D2_sum / 8.);

    return dE;
}
