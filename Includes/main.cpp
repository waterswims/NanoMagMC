#include "array_alloc.hpp"
#include "state.hpp"
#include "mklrand.h"
#include "functions.h"
#include "cpoints.hpp"
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
	int padding, N_av, Nsingle;
	double J=1, H=0, K=0, k=1, beta=50, amean, asd, lmean, lsd, size;
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

	// Storage Variables
	int* nums = alloc_1darr<int>(N_av);
	int* s_nums = alloc_1darr<int>(N_av);
	double** mag1 = alloc_2darr<double>(N_av, num_Ts);
	double** ener1 = alloc_2darr<double>(N_av, num_Ts);
	double** magx1 = alloc_2darr<double>(N_av, num_Ts);
	double** magy1 = alloc_2darr<double>(N_av, num_Ts);
	double** magz1 = alloc_2darr<double>(N_av, num_Ts);
	double** smag1 = alloc_2darr<double>(N_av, num_Ts);
	double** smagx1 = alloc_2darr<double>(N_av, num_Ts);
	double** smagy1 = alloc_2darr<double>(N_av, num_Ts);
	double** smagz1 = alloc_2darr<double>(N_av, num_Ts);

	vector<vector<double> > allmag, allener, allmagx, allmagy, allmagz, allsmag,
		allsmagx, allsmagy, allsmagz;
	vector<int> allnums, allsnums;
	vector<double> mtemp;

	// Loop through temperatures
	if (rank == 0)
	{
		cout << "Running MC..." << endl;
	}

	// checkpoint names
	double csize = 0;
	if (distributed)
	{
		csize = asd;
	}
	else
	{
		csize = size;
	}
	fstream cpoint;
    string magpoint = cpointname("mag", rank, csize, distributed, shape, hamil, H);
	string magxpoint = cpointname("magx", rank, csize, distributed, shape, hamil, H);
	string magypoint = cpointname("magy", rank, csize, distributed, shape, hamil, H);
	string magzpoint = cpointname("magz", rank, csize, distributed, shape, hamil, H);
	string enpoint = cpointname("en", rank, csize, distributed, shape, hamil, H);
	string seedipoint = cpointname("seedi", rank, csize, distributed, shape, hamil, H);
	string seeddpoint = cpointname("seedd", rank, csize, distributed, shape, hamil, H);
	string seedlnpoint = cpointname("seedln", rank, csize, distributed, shape, hamil, H);
	string numpoint = cpointname("num", rank, csize, distributed, shape, hamil, H);
	string smagpoint = cpointname("smag", rank, csize, distributed, shape, hamil, H);
	string smagxpoint = cpointname("smagx", rank, csize, distributed, shape, hamil, H);
	string smagypoint = cpointname("smagy", rank, csize, distributed, shape, hamil, H);
	string smagzpoint = cpointname("smagz", rank, csize, distributed, shape, hamil, H);
	string snumpoint = cpointname("snum", rank, csize, distributed, shape, hamil, H);

	// read checkpoint
	int j=0;
	if(exists(magpoint))
	{
		read_clist(cpoint, magpoint, mag1);
		read_clist(cpoint, smagpoint, smag1);
		if (hamil != 'i' && hamil != 'I')
		{
			read_clist(cpoint, magxpoint, magx1);
			read_clist(cpoint, magypoint, magy1);
			read_clist(cpoint, magzpoint, magz1);
			read_clist(cpoint, smagxpoint, smagx1);
			read_clist(cpoint, smagypoint, smagy1);
			read_clist(cpoint, smagzpoint, smagz1);
		}
		read_clist(cpoint, enpoint, ener1);
		j = read_cval(cpoint, numpoint, nums);
		read_cval(cpoint, snumpoint, s_nums);

		st_rand_int.load(seedipoint.c_str());
		st_rand_double.load(seeddpoint.c_str());
		rand_ln.load(seedlnpoint.c_str());
		st_rand_int.fill();
		st_rand_double.fill();
		rand_ln.fill();
	}

	if (rank == 0)
	{
		cout << "Passed Checkpoint Reading..." << endl;
	}

	// main loop
	for (; j < N_av; j++)
	{
		double s_size = 0;
		if(distributed)
		{
			s_size = rand_ln.gen();
		}
		else
		{
			s_size = size;
		}
        state curr_state(s_size, periodic, shape, hamil, J, H, k, Tmax, K, args);
		if (rank == 0)
		{
			cout << "Generated State..." << endl;
		}
		nums[j] = curr_state.num_spins();
		s_nums[j] = curr_state.sub_num(0);
		curr_state.equil(5*Nsingle*nums[j]);
		// cout << nums[j] << endl;
		// curr_state.change_temp(300);
		// for(int i = 0; i < 1000; i++)
		// {
		// 	curr_state.equil(10*nums[j]);
		// 	cerr << curr_state.magnetisation()[0] << " " << curr_state.energy() << endl;
		// }
		// return 0;
		if (rank == 0)
		{
			cout << "Completed Equillibriation..." << endl;
		}
        for (int i = 0; i < num_Ts; i++)
        {
            double T = Ts[i];
            curr_state.change_temp(T);

            curr_state.equil(Nsingle*nums[j]);
            mtemp = curr_state.magnetisation();
			if (hamil != 'i' && hamil != 'I')
			{
				magx1[j][i] = mtemp[0];
				magy1[j][i] = mtemp[1];
				magz1[j][i] = mtemp[2];
			}

			mag1[j][i] = norm(mtemp);
			ener1[j][i] = curr_state.energy();

			mtemp = curr_state.submag(0);
			if (hamil != 'i' && hamil != 'I')
			{
				smagx1[j][i] = mtemp[0];
				smagy1[j][i] = mtemp[1];
				smagz1[j][i] = mtemp[2];
			}
			smag1[j][i] = norm(mtemp);

			if (rank==0)
	        {
	            cout << "\rTemp " << i + 1 << " of " << num_Ts << " completed" << flush;
	        }
        }
		if (rank==0)
		{
			cout << endl;
		}

		// Checkpoint data
		print_clist(cpoint, magpoint, mag1[j], num_Ts);
		print_clist(cpoint, smagpoint, smag1[j], num_Ts);
		if (hamil != 'i' && hamil != 'I')
		{
			print_clist(cpoint, magxpoint, magx1[j], num_Ts);
			print_clist(cpoint, magypoint, magy1[j], num_Ts);
			print_clist(cpoint, magzpoint, magz1[j], num_Ts);
			print_clist(cpoint, smagxpoint, smagx1[j], num_Ts);
			print_clist(cpoint, smagypoint, smagy1[j], num_Ts);
			print_clist(cpoint, smagzpoint, smagz1[j], num_Ts);
		}
		print_clist(cpoint, enpoint, ener1[j], num_Ts);
		print_cval(cpoint, numpoint, nums[j]);
		print_cval(cpoint, snumpoint, s_nums[j]);
		st_rand_int.save(seedipoint.c_str());
		st_rand_double.save(seeddpoint.c_str());
		rand_ln.save(seedlnpoint.c_str());

        if (rank==0)
        {
            cout << "Lattice " << j + 1 << " of " << N_av << " completed" << endl;
        }
	}

	// Prep for printing and stuff
	int tmax = num_Ts;

	double** mag = alloc_2darr<double>(num_Ts, N_av);
	double** ener = alloc_2darr<double>(num_Ts, N_av);
	double** magx = alloc_2darr<double>(num_Ts, N_av);
	double** magy = alloc_2darr<double>(num_Ts, N_av);
	double** magz = alloc_2darr<double>(num_Ts, N_av);
	double** smag = alloc_2darr<double>(num_Ts, N_av);
	double** smagx = alloc_2darr<double>(num_Ts, N_av);
	double** smagy = alloc_2darr<double>(num_Ts, N_av);
	double** smagz = alloc_2darr<double>(num_Ts, N_av);

	for (int i = 0; i < N_av; i++)
	{
		for (int j = 0; j < num_Ts; j++)
		{
			mag[j][i] = mag1[i][j];
			magx[j][i] = magx1[i][j];
			magy[j][i] = magy1[i][j];
			magz[j][i] = magz1[i][j];
			ener[j][i] = ener1[i][j];
			smag[j][i] = smag1[i][j];
			smagx[j][i] = smagx1[i][j];
			smagy[j][i] = smagy1[i][j];
			smagz[j][i] = smagz1[i][j];
		}
	}

	if(rank==0)
	{
		cout << "Reducing Data..." << endl;
	}

	allmag = gather2darr(num_Ts, N_av, mag, 0, comm_size);
	allener = gather2darr(num_Ts, N_av, ener, 0, comm_size);
	allmagx = gather2darr(num_Ts, N_av, magx, 0, comm_size);
	allmagy = gather2darr(num_Ts, N_av, magy, 0, comm_size);
	allmagz = gather2darr(num_Ts, N_av, magz, 0, comm_size);
	allsmag = gather2darr(num_Ts, N_av, smag, 0, comm_size);
	allsmagx = gather2darr(num_Ts, N_av, smagz, 0, comm_size);
	allsmagy = gather2darr(num_Ts, N_av, smagy, 0, comm_size);
	allsmagz = gather2darr(num_Ts, N_av, smagz, 0, comm_size);
	allnums = gather1darr(N_av, nums, 0, comm_size);
	allsnums = gather1darr(N_av, s_nums, 0, comm_size);

	if(rank==0)
	{
		cout << "Finding profiles and printing to file..." << endl;
	}

	// print to file
	if(rank==0)
	{
		// Total lattice Size
		int g_lattsize = sum(allnums);
		int g_slattsize = sum(allsnums);

		// ready for printing
		string avname;
		if (distributed)
		{
			avname = main_name_d(H, J, amean, asd, k, shape, hamil);
			init_avs(avname, H, J, k, amean, Tmin, Tmax);
			create_folder_d(H, J, amean, asd, k, shape, hamil);
		}
		else
		{
			avname = main_name(H, J, size, k, shape, hamil);
			init_avs(avname, H, J, k, size, Tmin, Tmax);
			create_folder(H, J, size, k, shape, hamil);
		}

		int t = 0;
		for(int j = 0; j < num_Ts; j++)
		{
			double T = Ts[j];

			string datname;
			if(distributed)
			{datname = full_name_d(H, J, amean, asd, k, shape, hamil, T);}
			else
			{datname = full_name(H, J, size, k, shape, hamil, T);}

			print_full(datname, T, allener[t], allmag[t], allmagx[t], allmagy[t],
				allmagz[t], allsmag[t], allsmagx[t], allsmagy[t], allsmagz[t],
				allnums);

			if (t == tmax) {break;}

			print_avs(avname, allener[t], allmag[t], allmagx[t], allmagy[t],
				allmagz[t], allsmag[t], allsmagx[t], allsmagy[t], allsmagz[t],
				allnums, H, hamil, g_lattsize, g_slattsize, T);

			t++;
		}
	}

	// Destroy Checkpoints
	remove(magpoint.c_str());
	remove(magxpoint.c_str());
	remove(magypoint.c_str());
	remove(magzpoint.c_str());
	remove(enpoint.c_str());
	remove(seedipoint.c_str());
	remove(seeddpoint.c_str());
	remove(seedlnpoint.c_str());
	remove(numpoint.c_str());
	remove(smagpoint.c_str());
	remove(smagxpoint.c_str());
	remove(smagypoint.c_str());
	remove(smagzpoint.c_str());
	remove(snumpoint.c_str());

	dealloc_1darr<double>(Ts);
	dealloc_1darr<int>(nums);
	dealloc_1darr<int>(s_nums);
	dealloc_2darr<double>(num_Ts, mag);
	dealloc_2darr<double>(num_Ts, ener);
	dealloc_2darr<double>(num_Ts, magx);
	dealloc_2darr<double>(num_Ts, magy);
	dealloc_2darr<double>(num_Ts, magz);
	dealloc_2darr<double>(num_Ts, smag);
	dealloc_2darr<double>(num_Ts, smagx);
	dealloc_2darr<double>(num_Ts, smagy);
	dealloc_2darr<double>(num_Ts, smagz);
	dealloc_2darr<double>(N_av, mag1);
	dealloc_2darr<double>(N_av, ener1);
	dealloc_2darr<double>(N_av, magx1);
	dealloc_2darr<double>(N_av, magy1);
	dealloc_2darr<double>(N_av, magz1);
	dealloc_2darr<double>(N_av, smag1);
	dealloc_2darr<double>(N_av, smagx1);
	dealloc_2darr<double>(N_av, smagy1);
	dealloc_2darr<double>(N_av, smagz1);

    // Finish program
	MPI_Finalize();

    return 0;
}
