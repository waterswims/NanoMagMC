
#include "state_new.hpp"
#include "mklrand.h"
#include "functions.h"
#include "cpoints.hpp"
#include "output.hpp"
#include "param_read.hpp"

#include <iostream>
#include <vector>
#include <mpi.h>
#include <fstream>
#include <cstring>

mkl_irand st_rand_int(1e7, 1);
mkl_drand st_rand_double(1e8, 1);

const double Ry = 13.606;

// define the type of spin
typedef heis_spin sp_typ;

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
	double J, H, k, beta, amean, asd, lmean, lsd;
	bool periodic, distributed;
	char shape, hamil;
	string temp_name;
	read_all_vars(argv[1], size, J, H, k, periodic, shape, hamil, N_av,
		Nsingle, padding, beta, distributed, amean, asd, temp_name);
	double args[] = {beta, padding};
	AtoLn(amean, asd, lmean, lsd);

	// load temperatures
	double Ts[100];
	load_temps(temp_name, Ts);
	int num_Ts = 100;
	double Tmin(Ts[num_Ts-1]), Tmax(Ts[0]);

	//Set the random seeds of the generators
	st_rand_int.change_seed(rank);
	st_rand_double.change_seed(600+rank);

	// New Generator for Lognormal
    mkl_lnrand rand_ln(lmean, lsd, N_av, 1200+rank);

    if(rank==0)
    {
        cout << N_av << " lattices per process" << endl;
        cout << Nsingle << " equillibriation steps per spin per temperature step" << endl;
    }

	// Storage Variables
	int nums[N_av], s_nums[N_av];
	double mag1[N_av][num_Ts], ener1[N_av][num_Ts],
		magx1[N_av][num_Ts], magy1[N_av][num_Ts], magz1[N_av][num_Ts],
		smag1[N_av][num_Ts], smagx1[N_av][num_Ts], smagy1[N_av][num_Ts],
		smagz1[N_av][num_Ts];
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
        state<sp_typ> curr_state(s_size, periodic, shape, hamil, J, H, k, Tmax, args);
        nums[j] = curr_state.num_spins();
		s_nums[j] = curr_state.sub_num(0);
		curr_state.equil(5*Nsingle*nums[j]);
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
	            cout << "Temp " << i << " of " << num_Ts << " completed" << endl;
	        }
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
            cout << "Lattice " << j << " of " << N_av << " completed" << endl;
        }
	}

	// Prep for printing and stuff
	int tmax = num_Ts;
	double mag[num_Ts][N_av], ener[num_Ts][N_av],
		magx[num_Ts][N_av], magy[num_Ts][N_av], magz[num_Ts][N_av],
		smag[num_Ts][N_av], smagx[num_Ts][N_av], smagy[num_Ts][N_av],
		smagz[num_Ts][N_av];
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
		// all_states[0].ptf("latt2.txt");
		cout << "Reducing Data..." << endl;
	}

	vector<double> temp;
	double temparr[N_av];
	int temparr2[N_av];

	// Reduce data
	for(int t=0; t < num_Ts; t++)
	{
		if(rank==0)
		{
			allmag.push_back(temp);
			allener.push_back(temp);
			allmagx.push_back(temp);
			allmagy.push_back(temp);
			allmagz.push_back(temp);
			allsmag.push_back(temp);
			allsmagx.push_back(temp);
			allsmagy.push_back(temp);
			allsmagz.push_back(temp);
			for(int i=0; i < N_av; i++)
			{
				allmag[t].push_back(mag[t][i]);
				allener[t].push_back(ener[t][i]);
				allmagx[t].push_back(magx[t][i]);
				allmagy[t].push_back(magy[t][i]);
				allmagz[t].push_back(magz[t][i]);

				allsmag[t].push_back(smag[t][i]);
				allsmagx[t].push_back(smagx[t][i]);
				allsmagy[t].push_back(smagy[t][i]);
				allsmagz[t].push_back(smagz[t][i]);
				if(t==0)
				{
					allnums.push_back(nums[i]);
					allsnums.push_back(s_nums[i]);
				}
			}
			for(int i=1; i < comm_size; i++)
			{
				MPI_Recv(temparr, N_av, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
				for(int j=0; j < N_av; j++)
				{
					allmag[t].push_back(temparr[j]);
					//cout << "Got here" << endl;
				}
				MPI_Recv(temparr, N_av, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
				for(int j=0; j < N_av; j++)
				{
					allener[t].push_back(temparr[j]);
				}
				MPI_Recv(temparr, N_av, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
				for(int j=0; j < N_av; j++)
				{
					allmagx[t].push_back(temparr[j]);
					//cout << "Got here" << endl;
				}
				MPI_Recv(temparr, N_av, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
				for(int j=0; j < N_av; j++)
				{
					allmagy[t].push_back(temparr[j]);
					//cout << "Got here" << endl;
				}
				MPI_Recv(temparr, N_av, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
				for(int j=0; j < N_av; j++)
				{
					allmagz[t].push_back(temparr[j]);
					//cout << "Got here" << endl;
				}

				MPI_Recv(temparr, N_av, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
				for(int j=0; j < N_av; j++)
				{
					allsmag[t].push_back(temparr[j]);
					//cout << "Got here" << endl;
				}
				MPI_Recv(temparr, N_av, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
				for(int j=0; j < N_av; j++)
				{
					allsmagx[t].push_back(temparr[j]);
					//cout << "Got here" << endl;
				}
				MPI_Recv(temparr, N_av, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
				for(int j=0; j < N_av; j++)
				{
					allsmagy[t].push_back(temparr[j]);
					//cout << "Got here" << endl;
				}
				MPI_Recv(temparr, N_av, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
				for(int j=0; j < N_av; j++)
				{
					allsmagz[t].push_back(temparr[j]);
					//cout << "Got here" << endl;
				}
				//cout << "Two recieved" << endl;
				if(t==0)
				{
					MPI_Recv(temparr2, N_av, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
					for(int j=0; j < N_av; j++)
					{
						allnums.push_back(temparr2[j]);
					}

					MPI_Recv(temparr2, N_av, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
					for(int j=0; j < N_av; j++)
					{
						allsnums.push_back(temparr2[j]);
					}
				}
			}
		}
		else
		{
			MPI_Ssend(&mag[t][0], N_av, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			MPI_Ssend(&ener[t][0], N_av, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			MPI_Ssend(&magx[t][0], N_av, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			MPI_Ssend(&magy[t][0], N_av, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			MPI_Ssend(&magz[t][0], N_av, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

			MPI_Ssend(&smag[t][0], N_av, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			MPI_Ssend(&smagx[t][0], N_av, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			MPI_Ssend(&smagy[t][0], N_av, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			MPI_Ssend(&smagz[t][0], N_av, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			if(t==0)
			{
				MPI_Ssend(nums, N_av, MPI_INT, 0, 0, MPI_COMM_WORLD);
				MPI_Ssend(s_nums, N_av, MPI_INT, 0, 0, MPI_COMM_WORLD);
			}
		}
	}

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

    // Finish program
	MPI_Finalize();

    return 0;
}
