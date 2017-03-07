# MC_HAMR: Monte Carlo Simulation Software for Atomistic Models of Magnetic Materials for Heat Assisted Magnetic Recording

## Requirements

* At the moment the code has only been tested with the Intel compilers due to the use of the Intel MKL library for random number generation.

## Compiling the Code

There is a makefile in the root folder of the code for compiling which calls a makefile in the "Includes" folder. Typing "make run" will correctly compile the code.

## Running the Code and the Input File

After compilation an executable file will be placed in the root folder named "run". This can be executed by typing "./run INPUT_FILE" where "INPUT_FILE" is the name of the input file which takes a number of arguments. An example input file is provided called "INPUT.dat"

### Compulsory Settings

* BOLTZMANN: The value of the Boltzmann constant which is being used.

* LATTSHAPE: The shape of the atomic lattice which is being used. This can take a number of inputs.

    * s: Square 2D lattice

    * w: Weibull (swiss cheese) circle (Requires WEIBULLFACT)

    * h: Hexagonal lattice in a rough square shape

    * c: Cubic lattice

    * x: Weibull (swiss cheese) sphere (Requires WEIBULLFACT)

    * o: Hexagonal columns in a rough cube

    * f: Face Centered Tetragonal in a cube

* ISPERIO: Boolean value to determine if the boundary conditions of the lattice are periodic or not. 0 for non=periodic, 1 for periodic. (Note, for some shapes this makes little or no difference. i.e. A Weibull shape.)

* LATTPADDING: Adds a 'pad' of empty lattice sites around your lattice. Takes a number signifying the depth of this pad.

* MCSWEEPS: The number of Monte Carlo sweeps that should be performed on each lattice at each temperature step. A sweep is defined as 1 Monte Carlo step for every filled lattice site.

* LATTSPERPROC: This program is implemented to MPI communication, meaning that multiple independent processes are run in tandom, each one running distinct realisations of the lattice. This value indicates the number of lattices that will be simulated on each process.

* HAMILTONIAN: The Hamiltonian which is being used to calculate the total energy and energy change per Monte Carlo step. This can take a number of inputs:

    * i: Ising Model

    * h: Heisenburg Model

    * f: FePt parameterised Hamiltonian

* ISDISTRIB: A boolean value which defines whether or not the sizes of the realisations of each lattice are fixed or distributed. The only distribution of sizes which has so far been implemented is a log-normal distribution.

* TEMPNAME: In the folder "Includes/Temps/" there are a number of example text files containing the temperatures which the simulation will run. This option is used to specify which file is to be used without the ".txt".

### Optional/Situational Settings

* SIZE: If ISDISTRIB is set to false then the fixed size of the lattice being run must be given.

* MEANSIZE: If ISDISTRIB is set to true then the arithmetic mean of the sizes of all lattices must be given.

* SIZEDEV: If ISDISTRIB is set to true then the standard deviation of the sizes of all lattices must be given.

* WEIBULLFACT: If either the circular or spherical Weibull distributed grain is chosen then WEIBULLFACT denotes the weibull factor used in the grain generation.

* EXCHANGE: The inter-atomic exchange coupling constant J between coupled atoms in the Ising and Heisenberg models.

* MAGFIELD: The external magnetic field that is present when using the Ising, Heisenberg and FePt models.

* DMISTREN: The strength of the DMI interaction for a skyrmion hamiltonian. 
