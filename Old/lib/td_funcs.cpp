#include "../includes/td_funcs.hpp"
#include "../includes/functions.hpp"

#include <cmath>
#include <iostream>

#include <xtensor/xeval.hpp>
#include <xtensor/xview.hpp>

xt::xtensorf<double, xt::xshape<4>> J_sum, H_sum, D_sum, diff, temp, temp2,
 loc_diff;

double particle::funcs::calc_E(field::field_type& lattice,
    xt::xtensorf<double, xt::xshape<4>>& H)
{
    int Nspins = lattice.get_size();

    J_sum = {0, 0, 0, 0};
    temp2 = {0, 0, 0, 0};

    for(int i = 0; i < Nspins; i++)
    {
        temp = {0, 0, 0, 0};
        for(int j = 0; j < lattice.get_neigh(i).size(); j++)
        {
            if(lattice.use_J())
            {
            temp += lattice.get_J(i, j) *
                lattice.access(lattice.get_neigh(i)[j]);
            }
            if(lattice.use_D())
            {
                c_prod(lattice.get_D_vec(i, j),
                    lattice.access(lattice.get_neigh(i)[j]), temp2);

                temp += temp2;
            }
        }

        J_sum += temp * lattice.access(i);
    }

    H_sum = calc_M(lattice);

    double E = xt::eval(xt::sum(-H * H_sum - 0.5 * J_sum))[0];

    return E;
}

double particle::funcs::calc_dE(particle::field::field_type& lattice,
    int position, xt::xtensorf<double, xt::xshape<4>>& H)
{
    lattice.gen_rand();
    diff = lattice.access(position) - lattice.get_rand();
    J_sum = {0, 0, 0, 0};
    D_sum = {0, 0, 0, 0};

    for(int j = 0; j < lattice.get_neigh(position).size(); j++)
    {
        J_sum += lattice.get_J(position, j) *
            lattice.access(lattice.get_neigh(position)[j]);
        c_prod(lattice.access(lattice.get_neigh(position)[j]),
            lattice.get_D_vec(position, j), temp);
        D_sum += temp;
    }

    double dE = xt::eval(xt::sum((H + J_sum + D_sum) * diff))[0];

    return dE;
}

xt::xtensorf<double, xt::xshape<4>> particle::funcs::calc_M(
    particle::field::field_type& lattice)
{
    int Nspins = lattice.get_size();
    temp = {0, 0, 0, 0};

    for(int i = 0; i < Nspins; i++)
    {
        temp += lattice.access(i);
    }

    return temp;
}

xt::xtensorf<double, xt::xshape<4>> particle::funcs::calc_subM(
    particle::field::field_type& lattice,
    int subnumber)
{
    int Nspins = lattice.get_size();
    temp = {0, 0, 0, 0};

    for(int i = 0; i < Nspins; i++)
    {
        int possum = xt::eval(xt::sum(lattice.get_loc(i)))[0];
        if (possum%2 == subnumber)
        {
            temp += lattice.access(i);
        }
    }

    return temp;
}

std::vector<double> particle::funcs::calc_TC(
    particle::field::field_type& lattice)
{
    int Nspins = lattice.get_size();
    std::vector<double> out;
    int left, right, forward, backward;

    for(int i = 0; i < Nspins; i++)
    {
        if(out.size() < lattice.get_loc(i)[2]+1)
        {
            out.resize(out.size());
        }
        for(int j = 0; j < Nspins; j++)
        {
            loc_diff = lattice.get_loc(j) - lattice.get_loc(i);
            if(loc_diff[2] == 0)
            {
                if (loc_diff[0] == 1 && loc_diff[1] == 0)
                {
                    left = j;
                }
                else if (loc_diff[0] == -1 && loc_diff[1] == 0)
                {
                    right = j;
                }
                else if (loc_diff[0] == 0 && loc_diff[1] == 1)
                {
                    forward = j;
                }
                else if (loc_diff[0] == 0 && loc_diff[1] == -1)
                {
                    backward = j;
                }
            }
        }
        out[lattice.get_loc(i)[2]] += solid_angle(lattice.access(i),
            lattice.access(right), lattice.access(forward));
        out[lattice.get_loc(i)[2]] += solid_angle(lattice.access(i),
            lattice.access(left), lattice.access(backward));
    }
    for(int i = 0; i < out.size(); i++)
    {
        out[i] /= 2 * M_PI;
    }

    return out;
}

double particle::funcs::solid_angle(const xt::xtensorf<double, xt::xshape<4>>& s1,
                const xt::xtensorf<double, xt::xshape<4>>& s2,
                const xt::xtensorf<double, xt::xshape<4>>& s3)
{
    double s1s2 = xt::eval(xt::sum(s1 * s2))[0];
    double s1s3 = xt::eval(xt::sum(s1 * s3))[0];
    double s2s3 = xt::eval(xt::sum(s2 * s3))[0];

    double rho = 2 * (1 + s1s2) * (1 + s1s3) * (1 + s2s3);
    double dotsum = 1 + s1s2 + s1s3 + s2s3;

    c_prod(s1, s3, temp);
    double crosssum = xt::eval(xt::sum(temp * s2))[0];

    double ang = atan2(crosssum / rho, dotsum / rho);

    return ang;
}
