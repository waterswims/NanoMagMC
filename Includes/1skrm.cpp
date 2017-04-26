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
#include <sstream>

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
	int padding, N_av, Nsingle;
	double J, H, K, k, beta, amean, asd, lmean, lsd, size;
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

	field_3d_h totfield(size, periodic, 2);
	state curr_state(size, periodic, shape, hamil, J, 0, k, Tmax, K, args);

	for(int i=0; i < num_Ts; i++)
	{
		double T = Ts[i];
		curr_state.change_temp(T);
		totfield.allzero();
		curr_state.equil(3000*curr_state.num_spins());
		state field_state(curr_state);
		for(int j = 0; j < 20; j++)
		{
			double thisH = (j+1)*H/20;
			field_state.change_field(thisH);
			field_state.equil(100*curr_state.num_spins());
		}

		// main loop
		for (int j = 0; j < N_av; j++)
		{
			field_state.equil(Nsingle*curr_state.num_spins());
			field_state.add_to_av(&totfield);
	        cout << "Lattice " << j + 1 << " of " << N_av << " completed" << endl;
		}
		stringstream sstream;
		sstream << "mkdir -p Output/Skyrm_Print/" << T << endl;
		string runstring;
		getline(sstream, runstring);
		system(runstring.c_str());
		cout << runstring << endl;
		sstream << "Output/Skyrm_Print/" << T << "/latt" << endl;
		string folderstring;
		sstream >> folderstring;
		cout << folderstring << endl;
		totfield.print(folderstring);
		cout << "Temp " << i + 1 << " of " << num_Ts << " completed" << endl;
	}

    // Finish program
	MPI_Finalize();

    return 0;
}
