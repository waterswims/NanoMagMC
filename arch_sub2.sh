#!/bin/bash

# set default resource requirements for job
# - these can be overridden on the qsub command line (this is for a 4 hour job)
#PBS -l walltime=24:00:00
#PBS -l select=42
#PBS -A e05-nandef-kra
#PBS -j oe

# Change to the directory from which the job was submitted.
# (The actual name is held in the PBS environment variable $PBS_O_WORKDIR)
cd $PBS_O_WORKDIR
module swap Prg-Env-cray PrgEnv-intel

# Run my_prog taking input from file input_file and send output to output_file
aprun -n 1008 ./run ARC_INPUT2.dat
