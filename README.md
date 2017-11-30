# MC_HAMR: Monte Carlo Simulation Software for Atomistic Models of Magnetic Materials for Heat Assisted Magnetic Recording

## Requirements

* Requires a parallel build of HDF5 for data storage.
* At the moment the code has only been tested with the Intel compilers due to the use of the Intel MKL library for random number generation. Future versions of the code will include non-Intel RNG.

## Compiling the Code

There is a makefile in the root folder. Typing "make" will correctly compile the code. CC can be set within this makefile or with environment variables. There is also some support in the makefile for Mac and ARCHER builds.

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

* ISPERIO: Boolean value to determine if the boundary conditions of the lattice are periodic or not. 0 for non=periodic, 1 for periodic. (Note, for some shapes this makes little or no difference. i.e. A Weibull shape.)

* EQSWEEPS: The number of Monte Carlo sweeps that should be performed on each lattice before samples are taken. A sweep is defined as 1 Monte Carlo step for every filled lattice site.

* SAMPSWEEPS: The number of Monte Carlo sweeps that should be performed between samples. A sweep is defined as 1 Monte Carlo step for every filled lattice site.

* NSAMPS: The number of Monte Carlo samples to be taken per lattice.

* HAMILTONIAN: The Hamiltonian which is being used to calculate the total energy and energy change per Monte Carlo step. This can take a number of inputs:

    * i: Ising Model

    * h: Heisenburg Model

    * f: FePt parameterised Hamiltonian

    * s: Skyrmion Hamiltonian

* ISDISTRIB: A boolean value which defines whether or not the sizes of the realisations of each lattice are fixed or distributed. The only distribution of sizes which has so far been implemented is a log-normal distribution.

* TEMPNAME: In the folder "Temps" there are a number of example text files containing the temperatures at which the simulation will run. This option is used to specify which file is to be used.

* FIELDNAME: In the folder "Fields" there are a number of example text files containing the magnetic fields at which the simulation will run. This option is used to specify which file is to be used.

* PROTOCOL: The protocol which defines the path through the phase diagram that the simulation will take. The options are:

    * 1: The lattices move through the magnetic fields initially then the temperatures.

    * 2: The lattices move through the temperatures followed by the magnetic fields.

    * 3: The lattices move through the magnetic fields in reverse order then the temperatures in forward order. (Not yet implemented.)

    * 4: The lattices move through the temperatures in reverse order followed by the magnetic fields in forward order.

* PRINT_LATT: Boolean option signifying whether the average lattice should be stored. This may significantly increase file storage requirements. (This option has only been implemented for 3D continuous spin systems so far. Using it in other situations will result in large, blank datasets within the output.)

### Optional/Situational Settings

* SIZE: If ISDISTRIB is set to false then the fixed size of the lattice being run must be given.

* MEANSIZE: If ISDISTRIB is set to true then the arithmetic mean of the sizes of all lattices must be given.

* SIZEDEV: If ISDISTRIB is set to true then the standard deviation of the sizes of all lattices must be given.

* NLATTS: If ISDISTRIB is set to true then the number of different lattices to use is required. NSAMPS samples will be taken for each of these lattices.

* WEIBULLFACT: If either the circular or spherical Weibull distributed grain is chosen then WEIBULLFACT denotes the weibull factor used in the grain generation.

* EXCHANGE: The inter-atomic exchange coupling constant J between coupled atoms in the Ising and Heisenberg models.

* MAGFIELD: The external magnetic field that is present when using the Ising, Heisenberg and FePt models.

* DMISTREN: The strength of the DMI interaction for a skyrmion hamiltonian.
