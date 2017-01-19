#include "hamiltonian.hpp"
#include "functions.h"
#include <sstream>
#include <iostream>
#include <fstream>

const double pi = 3.141592653589793;

template <class T> ham_ising<T>::ham_ising(ham_type<ising_spin>& other)
{
    H = other.get_H();
    J = other.get_J();
}

template <class T> double ham_ising<T>::calc_E(field_type<ising_spin>* lattice)
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
    int curr, adj_curr;
    int H_sum(0), J_sum(0);
    while (!finished)
    {
        lattice->adjacent(pos, adj);
        curr = (lattice->next(finished, pos)).i_access();
        H_sum += curr;
        for (vector<ising_spin*>::iterator it = adj.begin(); it != adj.end(); it++)
        {
            adj_curr = (*it)->i_access();
            J_sum += curr*adj_curr;
        }
    }

    double E = -H * H_sum - 0.5 * J * J_sum;
}

template <class T> vector<double> ham_ising<T>::calc_M(field_type<ising_spin>* lattice)
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
    while (!finished)
    {
        sum += (lattice->next(finished, pos)).i_access();
    }

    vector<double> mag = boost::assign::list_of(sum);
    return mag;
}

template <class T> vector<double> ham_ising<T>::calc_subM(field_type<ising_spin>* lattice, int subnumber)
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
    T temp;
    while (!finished)
    {
        possum = sum(pos);
        temp = lattice->next(finished, pos);
        if (possum%2 == subnumber)
        {
            tsum += (temp).i_access();
        }
    }

    vector<double> mag = boost::assign::list_of(tsum);
    return mag;
}

template <class T> double ham_ising<T>::dE(field_type<ising_spin>* lattice, vector<int>& position)
{
    int val = (lattice->access(position)).i_access();
    double dEH = 2 * H * val;

    int sum = 0;
    lattice->adjacent(position, adj);
    //printv(position);
    //cout << adj.size() << endl;
    for(vector<ising_spin*>::iterator it = adj.begin(), end = adj.end(); it != end; it++)
    {
        sum += (*it)->i_access();
    }

    double dE = dEH + 2 * J * val * sum;
    return dE;
}

template <class T> ham_ising<T>& ham_ising<T>::operator=(ham_type<ising_spin>& other)
{
    H = other.get_H();
    J = other.get_J();
    return *this;
}

template <class T> ham_heis<T>::ham_heis(double Hin, double Jin)
{
    H = vector<double>(3);
    J = vector<double>(3);
    H[0] = 0;
    H[1] = 0;
    H[2] = Hin;
    J[0] = Jin;
    J[1] = Jin;
    J[2] = Jin;
    test = T();
    vsum.resize(3);
}

template <class T> ham_heis<T>::ham_heis(ham_type<heis_spin>& other)
{
    H = other.get_Hs();
    J = other.get_Js();
    test = *(other.get_test());
    vsum.resize(3);
}

template <class T> double ham_heis<T>::calc_E(field_type<heis_spin>* lattice)
{
    int sum=0;
    int start = 0;
    if(lattice->get_perio())
    {
        start++;
    }
    int dim = lattice->get_dim();
    pos.resize(dim);
    for (vector<int>::iterator it = pos.begin(); it != pos.end(); it++)
    {
        *it = start;
    }
    bool finished = false;
    H_sum.resize(3);
    J_sum.resize(3);
    H_sum[0] = 0;
    H_sum[1] = 0;
    H_sum[2] = 0;
    J_sum[0] = 0;
    J_sum[1] = 0;
    J_sum[2] = 0;
    while (!finished)
    {
        lattice->adjacent(pos, adj);
        curr = (lattice->next(finished, pos)).spin_access();
        H_sum[0] += curr[0];
        H_sum[1] += curr[1];
        H_sum[2] += curr[2];
        for (vector<heis_spin*>::iterator it = adj.begin(); it != adj.end(); it++)
        {
            adj_curr = (*it)->spin_access();
            J_sum[0] += curr[0]*adj_curr[0];
            J_sum[1] += curr[1]*adj_curr[1];
            J_sum[2] += curr[2]*adj_curr[2];
        }
    }
    H_sum[0] = H_sum[0] * H[0];
    H_sum[1] = H_sum[1] * H[1];
    H_sum[2] = H_sum[2] * H[2];
    J_sum[0] = J_sum[0] * J[0];
    J_sum[1] = J_sum[1] * J[1];
    J_sum[2] = J_sum[2] * J[2];
    double E = -(H_sum[0] + H_sum[1] + H_sum[2]) - 0.5*(J_sum[0] + J_sum[1] + J_sum[2]);

    return E;
}

template <class T> vector<double> ham_heis<T>::calc_M(field_type<heis_spin>* lattice)
{
    vsum.resize(3);
    vsum[0] = 0;
    vsum[1] = 0;
    vsum[2] = 0;
    int start = 0;
    if(lattice->get_perio())
    {
        start++;
    }
    int dim = lattice->get_dim();
    pos.resize(dim);
    for (vector<int>::iterator it = pos.begin(); it != pos.end(); it++)
    {
        *it = start;
    }
    bool finished = false;
    while (!finished)
    {
        curr = (lattice->next(finished, pos)).spin_access();
        vsum[0] += curr[0];
        vsum[1] += curr[1];
        vsum[2] += curr[2];
    }
    return vsum;
}

template <class T> vector<double> ham_heis<T>::calc_subM(field_type<heis_spin>* lattice, int subnumber)
{
    vsum.resize(3);
    vsum[0] = 0;
    vsum[1] = 0;
    vsum[2] = 0;
    int start = 0;
    if(lattice->get_perio())
    {
        start++;
    }
    int dim = lattice->get_dim();
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
        curr = (lattice->next(finished, pos)).spin_access();
        if (possum%2 == subnumber)
        {
            vsum[0] += curr[0];
            vsum[1] += curr[1];
            vsum[2] += curr[2];
        }
    }
    return vsum;
}

template <class T> double ham_heis<T>::dE(field_type<heis_spin>* lattice, vector<int>& position)
{
    test.rand_spin();
    // potential = test.spin_access();
    // curr = (lattice->access(position)).spin_access();
    double cmp1 = (lattice->access(position)).spin_access()[0] - test.spin_access()[0];
    double cmp2 = (lattice->access(position)).spin_access()[1] - test.spin_access()[1];
    double cmp3 = (lattice->access(position)).spin_access()[2] - test.spin_access()[2];
    double dEH = H[0] * cmp1 + H[1] * cmp2 + H[2] * cmp3;

    lattice->adjacent(position, adj);
    vsum[0] = 0;
    vsum[1] = 0;
    vsum[2] = 0;
    for(vector<heis_spin*>::iterator it = adj.begin(), end = adj.end(); it != end; it++)
    {
        vsum[0] += (*it)->spin_access()[0];
        vsum[1] += (*it)->spin_access()[1];
        vsum[2] += (*it)->spin_access()[2];
    }

    double dE = dEH + (J[0] * cmp1 * vsum[0] + J[1] * cmp2 * vsum[1] + J[2] * cmp3 * vsum[2]);

    return dE;
}

template <class T> ham_heis<T>& ham_heis<T>::operator=(ham_type<heis_spin>& other)
{
    H = other.get_Hs();
    J = other.get_Js();
    test = *(other.get_test());
    vsum.resize(3);
    return *this;
}

template class ham_ising<ising_spin>;
template class ham_ising<heis_spin>;
template class ham_heis<ising_spin>;
template class ham_heis<heis_spin>;

template <class T> double ham_FePt<T>::calc_E(field_type<heis_spin>* lattice)
{
    Jx_sum = 0;
    Jy_sum = 0;
    d2_sum = 0;

    int t_size = lattice->get_totsize();
    int i_size = lattice->get_insize();
    int start = (t_size-i_size)/2;
    int fin = start+i_size;

    int dim = lattice->get_dim();
    pos2.resize(dim);
    for (vector<int>::iterator it = pos2.begin(); it != pos2.end(); it++)
    {
        *it = start;
    }
    pos = pos2;
    int nN = dxs.size();
    bool finished = false;
    while (!finished)
    {
        curr = (lattice->next(finished, pos2)).spin_access();
        // All neighbour interactions
        for(int i = 0; i < nN; i++)
        {
            pos[0] = pos2[0] + dxs[i];
            pos[1] = pos2[1] + dys[i];
            pos[2] = pos2[2] + dzs[i];

            if(this->check_pos(start, fin))
            {
                adj_curr = (lattice->access(pos2)).spin_access();

                Jx_sum += curr[0] * adj_curr[0] * Js[i];
                Jy_sum += curr[1] * adj_curr[1] * Js[i];

                d2_sum += curr[2] * adj_curr[2] * d_ijs[i];
            }
        }
        // 1 ion anisotropy
        d0_sum += curr[2]*curr[2];
    }

    // totals
    double E = -(d0 * d0_sum + 0.5 * (Jx_sum + Jy_sum + d2_sum));

    return E;
}

template <class T> vector<double> ham_FePt<T>::calc_M(field_type<heis_spin>* lattice)
{
    vsum.resize(3);
    vsum[0] = 0;
    vsum[1] = 0;
    vsum[2] = 0;
    int start = 0;
    if(lattice->get_perio())
    {
        start++;
    }
    int dim = lattice->get_dim();
    pos.resize(dim);
    for (vector<int>::iterator it = pos.begin(); it != pos.end(); it++)
    {
        *it = start;
    }
    bool finished = false;
    while (!finished)
    {
        curr = (lattice->next(finished, pos)).spin_access();
        vsum[0] += curr[0];
        vsum[1] += curr[1];
        vsum[2] += curr[2];
    }
    return vsum;
}

template <class T> vector<double> ham_FePt<T>::calc_subM(field_type<heis_spin>* lattice, int subnumber)
{
    vsum.resize(3);
    vsum[0] = 0;
    vsum[1] = 0;
    vsum[2] = 0;
    int start = 0;
    if(lattice->get_perio())
    {
        start++;
    }
    int dim = lattice->get_dim();
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
        curr = (lattice->next(finished, pos)).spin_access();
        if (possum%2 == subnumber)
        {
            vsum[0] += curr[0];
            vsum[1] += curr[1];
            vsum[2] += curr[2];
        }
    }
    return vsum;
}

template <class T> double ham_FePt<T>::dE(field_type<heis_spin>* lattice, vector<int>& position)
{
    Jx_sum = 0;
    Jy_sum = 0;
    // Jz_sum = 0;
    d2_sum = 0;

    int t_size = lattice->get_totsize();
    int i_size = lattice->get_insize();
    int start = (t_size-i_size)/2;
    int fin = start+i_size;

    test.rand_spin();
    potential = test.spin_access();
    // curr = (lattice->access(position)).spin_access();
    curr = lattice->spin_access(position);
    double cmp1 = curr[0] - potential[0];
    double cmp2 = curr[1] - potential[1];
    double cmp3 = curr[2] - potential[2];

    int nN = dxs.size();

    // All neighbour interactions

    for(int i = 0; i < nN; i++)
    {
        pos[0] = position[0] + dxs[i];
        pos[1] = position[1] + dys[i];
        pos[2] = position[2] + dzs[i];

        if(this->check_pos(start, fin))
        {
            // adj_curr = (lattice->access(pos)).spin_access();
            adj_curr = lattice->spin_access(pos);

            Jx_sum += adj_curr[0] * Js[i];
            Jy_sum += adj_curr[1] * Js[i];
            // Jz_sum += adj_curr[2] * Js[i];

            d2_sum += adj_curr[2] * d_ijs[i];
        }
    }

    // 1 ion anisotropy
    d0_sum = (curr[2]*curr[2] - potential[2]*potential[2]) * d0;

    // totals
    double dE = d0_sum + cmp1 * Jx_sum + cmp2 * Jy_sum + cmp3 * (d2_sum);

    return dE;
}

template <class T> ham_FePt<T>::ham_FePt(ham_type<heis_spin>& other)
{
    this->read_Js();
    vsum.resize(3);
    test = *(other.get_test());

}

template <class T> ham_FePt<T>::ham_FePt()
{
    this->read_Js();
    test = T();
    vsum.resize(3);
}

template <class T> ham_FePt<T>& ham_FePt<T>::operator=(ham_type<heis_spin>& other)
{
    this->read_Js();
    vsum.resize(3);
    test = *(other.get_test());
    return *this;
}

template <class T> void ham_FePt<T>::read_Js()
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

    // posvec.resize(dxs.size());
    // adjvec.resize(dxs.size());
    // for (int i = 0; i < dxs.size(); i++) {
    //     posvec[i].resize(3);
    //     adjvec[i].resize(3);
    // }
}

template <class T> bool ham_FePt<T>::check_pos(int start, int fin)
{
    if(pos[0] < start) return false;
    if(pos[1] < start) return false;
    if(pos[2] < start) return false;
    if(pos[0] >= fin) return false;
    if(pos[1] >= fin) return false;
    if(pos[2] >= fin) return false;
    return true;
}

template class ham_FePt<ising_spin>;
template class ham_FePt<heis_spin>;

template <class T> ham_cluster<T>::ham_cluster(string filename)
{
    ifstream file;
    file.open(filename.c_str());
    if(!file.is_open())
    {
        cout << "Input file not opened" << endl;
        exit(105);
    }
    string line;
    double temp_d;
    while(getline(file, line))
    {
        istringstream stream(line);
        stream >> temp_d;
        (xs).push_back(temp_d);
        stream >> temp_d;
        (ys).push_back(temp_d);
        stream >> temp_d;
        (zs).push_back(temp_d);
        stream >> temp_d;
        (Rms).push_back(temp_d);
        stream >> temp_d;
        (Rns).push_back(temp_d);
    }
    file.close();

    ks = vector<T>((xs).size());
    for(int i=0; i < this->insize; i++)
    {
        (ks)[i].rand_spin();
        (Vms).push_back(4.0*pi*pow(Rms[i], 3.0) / 3.0);
    }
    mu0 = 4e-7 * pi;
}

template <class T> double ham_cluster<T>::calc_E(field_type<heis_spin>* lattice)
{
    double ani_sum = 0, mag_sum = 0, di_sum = 0;

    int n_parts = lattice->get_totsize();
    pos.resize(1);
    pos2.resize(1);
    temp.resize(3);

    for (pos[0] = 0; pos[0] < n_parts; pos[0]++)
    {
        curr = (lattice->access(pos)).spin_access();
        curr_ani = ks[pos[0]].spin_access();

        ani_sum += (pow((curr_ani[1]*curr[2] - curr[1]*curr_ani[2]), 2) +
                   pow((curr_ani[2]*curr[0] - curr[2]*curr_ani[0]), 2) +
                   pow((curr_ani[0]*curr[1] - curr[0]*curr_ani[1]), 2)) *
                   Vms[pos[0]];

        mag_sum -= Vms[pos[0]] * (curr[0] * (*h)[0] + curr[1] *
                                  (*h)[1] + curr[2] * (*h)[2]);

        temp[0] = 0;
        temp[1] = 0;
        temp[2] = 0;
        for(pos2[0] = 0; pos2[0] < pos[0]; pos2[0]++)
        {
            curr2 = (lattice->access(pos2)).spin_access();
            double dx = xs[pos[0]] - xs[pos2[0]];
            double dy = ys[pos[0]] - ys[pos2[0]];
            double dz = zs[pos[0]] - zs[pos2[0]];

            double modr = pow((pow(dx, 2) + pow(dy, 2) + pow(dz, 2)), 0.5);

            dx = dx / modr;
            dy = dy / modr;
            dz = dz / modr;

            double invr3 = pow(modr, -3);

            double mr = curr2[0] * dx + curr2[1] * dy + curr2[2] * dz;

            temp[0] += Vms[pos2[0]] * invr3 * (-curr2[0] + dx * mr);
            temp[1] += Vms[pos2[0]] * invr3 * (-curr2[1] + dy * mr);
            temp[2] += Vms[pos2[0]] * invr3 * (-curr2[2] + dz * mr);
        }

        di_sum -= Vms[pos[0]] * Ms * (curr[0] * temp[0] + curr[1] * temp[1] +
                                      curr[2] * temp[2]);
    }

    double E = ani_const * ani_sum + mu0 * Ms * (mag_sum + di_sum * Ms /
                                                 (4.0 * pi));

    return E;
}

template <class T> double ham_cluster<T>::dE(field_type<heis_spin>* lattice, vector<int>& position)
{
    int n_parts = lattice->get_totsize();
    pos2.resize(1);
    temp.resize(3);
    diff.resize(3);

    test.rand_spin();
    potential = test.spin_access();

    curr = (lattice->access(position)).spin_access();
    curr_ani = ks[position[0]].spin_access();

    diff[0] = potential[0] - curr[0];
    diff[1] = potential[1] - curr[1];
    diff[2] = potential[2] - curr[2];

    double diff_ani_sum = (pow((curr_ani[1]*diff[2] - diff[1]*curr_ani[2]), 2) +
                          pow((curr_ani[2]*diff[0] - diff[2]*curr_ani[0]), 2) +
                          pow((curr_ani[0]*diff[1] - diff[0]*curr_ani[1]), 2)) *
                          Vms[position[0]];

    double diff_mag_sum = -Vms[position[0]] * (diff[0] * (*h)[0] + diff[1] *
                                          (*h)[1] + diff[2] * (*h)[2]);

    temp[0] = 0;
    temp[1] = 0;
    temp[2] = 0;
    for(pos2[0] = 0; pos2[0] < n_parts; pos2[0]++)
    {
        if (pos2[0] == position[0]) {continue;}

        curr2 = (lattice->access(pos2)).spin_access();
        double dx = xs[position[0]] - xs[pos2[0]];
        double dy = ys[position[0]] - ys[pos2[0]];
        double dz = zs[position[0]] - zs[pos2[0]];

        double modr = pow((pow(dx, 2) + pow(dy, 2) + pow(dz, 2)), 0.5);

        dx = dx / modr;
        dy = dy / modr;
        dz = dz / modr;

        double invr3 = pow(modr, -3);

        double mr = curr2[0] * dx + curr2[1] * dy + curr2[2] * dz;

        temp[0] += Vms[pos2[0]] * invr3 * (-curr2[0] + dx * mr);
        temp[1] += Vms[pos2[0]] * invr3 * (-curr2[1] + dy * mr);
        temp[2] += Vms[pos2[0]] * invr3 * (-curr2[2] + dz * mr);
    }

    double diff_di_sum = -Vms[position[0]] * Ms * (diff[0] * temp[0] +
                                                  diff[1] * temp[1] +
                                                  diff[2] * temp[2]);

    double dE = ani_const * diff_ani_sum + mu0 * Ms * (diff_mag_sum +
                                                       diff_di_sum * Ms /
                                                       (4.0 * pi));

    return dE;
}

template <class T> vector<double> ham_cluster<T>::calc_M(field_type<heis_spin>* lattice)
{
    temp.resize(3);
    vector<T>* field_point = lattice->get_1dfield();

    temp[0] = 0;
    temp[1] = 0;
    temp[2] = 0;

    for(vector<T>::iterator it = (*field_point).begin();
        it != (*field_point).end(); i++)
    {
        curr = (*it).spin_access();
        temp[0] += curr[0];
        temp[1] += curr[1];
        temp[2] += curr[2];
    }

    return temp;
}

// template <class T> ham_skyrm<T>::ham_skyrm(double Hin, double Jin, double Dxin, double Dyin)
// {
//     H = vector<double>(3);
//     J = vector<double>(3);
//     H[0] = 0;
//     H[1] = 0;
//     H[2] = Hin;
//     J[0] = Jin;
//     J[1] = Jin;
//     J[2] = Jin;
//     test = T();
//     vsum.resize(3);
//     Dx = Dxin;
//     Dy = Dyin;
// }
//
// template <class T> ham_skyrm<T>::ham_skyrm(ham_type<heis_spin>& other)
// {
//     H = other.get_Hs();
//     J = other.get_Js();
//     test = *(other.get_test());
//     vsum.resize(3);
//     Dx = other.get_Dx();
//     Dy = other.get_Dy();
// }
//
// template <class T> ham_skyrm<T>& ham_skyrm<T>::operator=(ham_type<heis_spin>& other)
// {
//     H = other.get_Hs();
//     J = other.get_Js();
//     test = *(other.get_test());
//     vsum.resize(3);
//     Dx = other.get_Dx();
//     Dy = other.get_Dy();
//
//     return *this;
// }
//
// template <class T> double ham_skyrm<T>::dE(field_type<heis_spin>* lattice, vector<int>& position)
// {
//     test.rand_spin();
//     double cmp1 = (lattice->access(position)).spin_access()[0] - test.spin_access()[0];
//     double cmp2 = (lattice->access(position)).spin_access()[1] - test.spin_access()[1];
//     double cmp3 = (lattice->access(position)).spin_access()[2] - test.spin_access()[2];
//     double dEH = H[0] * cmp1 + H[1] * cmp2 + H[2] * cmp3;
//
//     lattice->adjacent(position, adj);
//     vsum[0] = 0;
//     vsum[1] = 0;
//     vsum[2] = 0;
//     for(vector<heis_spin*>::iterator it = adj.begin(), end = adj.end(); it != end; it++)
//     {
//         vsum[0] += (*it)->spin_access()[0];
//         vsum[1] += (*it)->spin_access()[1];
//         vsum[2] += (*it)->spin_access()[2];
//     }
//
//     double dE = dEH + (J[0] * cmp1 * vsum[0] + J[1] * cmp2 * vsum[1] + J[2] * cmp3 * vsum[2]);
//
//     return dE;
// }
//
// template <class T> double ham_skyrm<T>::calc_E(field_type<heis_spin>* lattice)
// {
//     int sum=0;
//     int start = 0;
//     if(lattice->get_perio())
//     {
//         start++;
//     }
//     int dim = lattice->get_dim();
//     pos.resize(dim);
//     for (vector<int>::iterator it = pos.begin(); it != pos.end(); it++)
//     {
//         *it = start;
//     }
//     bool finished = false;
//     H_sum.resize(3);
//     J_sum.resize(3);
//     H_sum[0] = 0;
//     H_sum[1] = 0;
//     H_sum[2] = 0;
//     J_sum[0] = 0;
//     J_sum[1] = 0;
//     J_sum[2] = 0;
//     while (!finished)
//     {
//         lattice->adjacent(pos, adj);
//         curr = (lattice->next(finished, pos)).spin_access();
//         H_sum[0] += curr[0];
//         H_sum[1] += curr[1];
//         H_sum[2] += curr[2];
//         for (vector<heis_spin*>::iterator it = adj.begin(); it != adj.end(); it++)
//         {
//             adj_curr = (*it)->spin_access();
//             J_sum[0] += curr[0]*adj_curr[0];
//             J_sum[1] += curr[1]*adj_curr[1];
//             J_sum[2] += curr[2]*adj_curr[2];
//         }
//     }
//     H_sum[0] = H_sum[0] * H[0];
//     H_sum[1] = H_sum[1] * H[1];
//     H_sum[2] = H_sum[2] * H[2];
//     J_sum[0] = J_sum[0] * J[0];
//     J_sum[1] = J_sum[1] * J[1];
//     J_sum[2] = J_sum[2] * J[2];
//     double E = -(H_sum[0] + H_sum[1] + H_sum[2]) - 0.5*(J_sum[0] + J_sum[1] + J_sum[2]);
//
//     return E;
// }
