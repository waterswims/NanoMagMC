#!/bin/bash

# set default resource requirements for job
# - these can be overridden on the qsub command line (this is for a 4 hour job)
#PBS -l walltime=60:00:00
#PBS -l nodes=32:ppn=16
#PBS -j oe

# Change to the directory from which the job was submitted.
# (The actual name is held in the PBS environment variable $PBS_O_WORKDIR)
cd $PBS_O_WORKDIR
module load intel/2016
module load intel/mkl/2016
module load intel/mpi/2016
nprocs=`wc -l $PBS_NODEFILE | awk '{ print $1 }'`

# Run my_prog taking input from file input_file and send output to output_file
mpirun -bootstrap rsh -f $PBS_NODEFILE -n $nprocs ./run D_INPUT_5