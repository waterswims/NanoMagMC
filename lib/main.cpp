#include "../includes/array_alloc.hpp"
#include "../includes/state.hpp"
#include "../includes/mklrand.hpp"
#include "../includes/functions.hpp"
#include "../includes/cpoints.hpp"
#include "../includes/output.hpp"
#include "../includes/param_read.hpp"
#include "../includes/mpifuncs.hpp"
#include "../includes/protocol.hpp"

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
	int padding, Samp_steps, Eq_steps, N_samp, protocol;
	double J=1, K=0, k=1, beta=50, amean, asd, lmean, lsd, size;
	bool periodic, distributed, print_latt;
	char shape, hamil;
	string temp_name, field_name;
	read_all_vars(argv[1], size, J, k, periodic, shape, hamil, Samp_steps,
		N_samp, Eq_steps, padding, beta, distributed, amean, asd, temp_name,
		field_name, protocol, K, print_latt);
	double args[] = {beta};
	AtoLn(amean, asd, lmean, lsd);

	// load temperatures and fields
	double* Ts;
	int num_Ts;
	double* Hs;
	int num_Hs;
	load_Hs_Ts(temp_name, Ts, num_Ts, field_name, Hs, num_Hs);
	double Tmin(Ts[num_Ts-1]), Tmax(Ts[0]);
	double Hmin(Hs[num_Hs-1]), Hmax(Hs[0]);

	//Set the random seeds of the generators
	st_rand_int.change_seed(rank);
	st_rand_double.change_seed(comm_size+rank);

	// New Generator for Lognormal
    mkl_lnrand rand_ln(lmean, lsd, N_samp, 2*comm_size+rank);

    if(rank==0)
    {
		cout << num_Ts << " temperatures" << endl;
		cout << num_Hs << " field strengths" << endl;
        cout << N_samp << " samples per configuration" << endl;
    }

	// Storage Variables
	int nums, s_nums;
	double* mag1 = alloc_1darr<double>(N_samp);
	double* ener1 = alloc_1darr<double>(N_samp);
	double* magx1 = alloc_1darr<double>(N_samp);
	double* magy1 = alloc_1darr<double>(N_samp);
	double* magz1 = alloc_1darr<double>(N_samp);
	double* smag1 = alloc_1darr<double>(N_samp);
	double* smagx1 = alloc_1darr<double>(N_samp);
	double* smagy1 = alloc_1darr<double>(N_samp);
	double* smagz1 = alloc_1darr<double>(N_samp);

	vector<double> mtemp;

	// set protocol
	double* var1_list;
	double* var2_list;
	int v1_size, v2_size;
	set_protocol(protocol, var1_list, var2_list, v1_size, v2_size,Hs, Ts,
		num_Hs, num_Ts);

	// Start MC
	if (rank == 0)
	{
		cout << "Running MC..." << endl;
	}

	// Check the Checkpoints
	if (rank == 0)
	{
		cout << "Passed Checkpoint Reading..." << endl;
	}

	// Generate the state
	state base_state(size, periodic, shape, hamil, J, Hmax, k, Tmax, K, args);
	state curr_state(base_state);
	nums = base_state.num_spins();
	s_nums = base_state.sub_num(0);
	if (rank == 0)
	{
		cout << "Generated State..." << endl;
	}

	// Create Output Folders
	create_folders(J, size, k, K, shape, hamil);

	// main loop
	for (int i = 0; i < v1_size; i++)
	{
		// Change Temp/field
		base_state.change_v1(protocol, var1_list[i]);
		// Equillibriation
		base_state.equil(Eq_steps*nums);
		// Copy
		curr_state = base_state;
        for (int j = 0; j < v2_size; j++)
        {
			// Change Temp/field
			curr_state.change_v2(protocol, var2_list[j]);
			// Equillibriation
			curr_state.equil(Eq_steps*nums);

			// Get Samples
			for (int ns = 0; ns < N_samp; ns++)
			{
				curr_state.equil(Samp_steps*nums);
				mtemp = curr_state.magnetisation();
				if (hamil != 'i' && hamil != 'I')
				{
					magx1[ns] = mtemp[0];
					magy1[ns] = mtemp[1];
					magz1[ns] = mtemp[2];
				}
				else
				{
					magz1[ns] = mtemp[0];
				}

				mag1[ns] = norm(mtemp);
				ener1[ns] = curr_state.energy();

				mtemp = curr_state.submag(0);
				if (hamil != 'i' && hamil != 'I')
				{
					smagx1[ns] = mtemp[0];
					smagy1[ns] = mtemp[1];
					smagz1[ns] = mtemp[2];
				}
				else
				{
					smagz1[ns] = mtemp[0];
				}
				smag1[ns] = norm(mtemp);
			}
			print_sngl_HT(magx1, magy1, magz1, mag1, ener1, smagx1, smagy1,
				smagz1, smag1, N_samp, protocol, var1_list[i],
				var2_list[j], J, size, k, K, shape, hamil);
			// print the lattice
			if(print_latt)
			{
				continue;
			}
        }
		// Checkpoint data
	}
	cout << endl;

	// Destroy Checkpoints

	// Destroy arrays
	dealloc_1darr<double>(Ts);
	dealloc_1darr<double>(Hs);
	dealloc_1darr<double>(mag1);
	dealloc_1darr<double>(ener1);
	dealloc_1darr<double>(magx1);
	dealloc_1darr<double>(magy1);
	dealloc_1darr<double>(magz1);
	dealloc_1darr<double>(smag1);
	dealloc_1darr<double>(smagx1);
	dealloc_1darr<double>(smagy1);
	dealloc_1darr<double>(smagz1);

    // Finish program
	MPI_Finalize();

    return 0;
}
