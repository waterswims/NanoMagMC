#include "array_alloc.hpp"
#include "state.hpp"
#include "mklrand.h"
#include "functions.h"
#include "output.hpp"
#include "param_read.hpp"
#include "mpifuncs.hpp"

#include <iostream>
#include <vector>
#include <fstream>
#include <cstring>

mkl_irand st_rand_int(1e7, 1);
mkl_drand st_rand_double(1e8, 1);

int main(int argc, char **argv)
{
    //Start MPI off
	MPI_Init(&argc, &argv);
	int rank, comm_size;
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status status;
  	if(rank==0)
	{
		cout << "Total number of processes is " << comm_size << endl;
	}

	// setup arguments
	int padding, size, N_av, Nsingle;
	double J, H, K, k, beta, amean, asd, lmean, lsd;
	bool periodic, distributed;
	char shape, hamil;
	string temp_name;
	read_all_vars(argv[1], size, J, H, k, periodic, shape, hamil, N_av,
		Nsingle, padding, beta, distributed, amean, asd, temp_name, K);
	double args[] = {beta};
	AtoLn(amean, asd, lmean, lsd);

	// load temperatures
	double* Ts;
	int num_Ts;
	load_temps(temp_name, Ts, num_Ts);
	double Tmin(Ts[num_Ts-1]), Tmax(Ts[0]);

	// Check correct inputs
	if(num_Ts != 1)
	{
		cout << "Skyrmion can only be printed for a single temperature" <<endl;
		exit(1000);
	}
	if(distributed)
	{
		cout << "Skyrmion can only be printed for non-distributed sizes" <<endl;
		exit(1000);
	}

	//Set the random seeds of the generators
	st_rand_int.change_seed(rank);
	st_rand_double.change_seed(600+rank);

	// New Generator for Lognormal
    mkl_lnrand rand_ln(lmean, lsd, N_av, 1200+rank);

    if(rank==0)
    {
		cout << num_Ts << " temperatures" << endl;
        cout << N_av << " lattices per process" << endl;
        cout << Nsingle << " equillibriation steps per spin per temperature step" << endl;
    }

	// Loop through temperatures
	if (rank == 0)
	{
		cout << "Running MC..." << endl;
	}

	field_3d_h totfield(size, periodic);
	totfield.allzero();

	// main loop
	for (int j = 0; j < N_av; j++)
	{
        state curr_state(size, periodic, shape, hamil, J, H, k, Tmax, K, args);
		if (rank == 0)
		{
			cout << "Generated State..." << endl;
		}
		curr_state.equil(Nsingle*curr_state.num_spins());
		if (rank == 0)
		{
			cout << "Completed Equillibriation..." << endl;
		}

		curr_state.add_to_av(&totfield);

        if (rank==0)
        {
            cout << "Lattice " << j << " of " << N_av << " completed" << endl;
        }
	}

	double*** xsout;
	double*** ysout;
	double*** zsout;
	totfield.get_3dfield_h(xsout, ysout, zsout);
	int tsize = totfield.get_totsize();
	double*** gxspin = redc3darr(tsize, tsize, tsize, xsout, 0);
	double*** gyspin = redc3darr(tsize, tsize, tsize, ysout, 0);
	double*** gzspin = redc3darr(tsize, tsize, tsize, zsout, 0);
	vector<int> pos(3);
	for(int i = 0; i < tsize; i++)
	{
		pos[0] = i;
		for(int j = 0; j < tsize; j++)
		{
			pos[1] = j;
			for(int k = 0; k < tsize; k++)
			{
				pos[2] = k;
				totfield.fill_val_h(pos, gxspin[i][j][k], gyspin[i][j][k], gzspin[i][j][k]);
			}
		}
	}
	dealloc_3darr<double>(tsize, tsize, gxspin);
	dealloc_3darr<double>(tsize, tsize, gyspin);
	dealloc_3darr<double>(tsize, tsize, gzspin);

	if(rank == 0)
	{
		totfield.print("Output/Skyrm_Print/latt");
	}

    // Finish program
	MPI_Finalize();

    return 0;
}
