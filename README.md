# Homework 3, Due Monday, April 11

## Submission Instructions

* Before class, all code and related materials should be pushed to your private GitHub sc3260s16 repository within a subdirectory called ```hwk3```.

* Additionally, you should also hand in a summary of benchmarks and other requested output to me at the beginning of class on the due date.

* You can review the late submission policy at the course syllabus: https://github.com/sc3260s16/syllabus

  * Late submissions must be brought to me in my office or submitted in subsequent classes.

## Compiling Instructions

Use Intel's MPI C compiler (*mpicc -cc=icc*) to build all programs for this homework. If you haven't done so already, stick 
this line in your ~/.bashrc:

```
setpkgs -a intel_cluster_studio_compiler
```

Apply the following command-line flags:

```
mpicc -cc=icc -Iinc -Wall -xHost -O3 -vec-report2
```

Also ensure that you link against Intel's imf library within your Makefile:

```
LDFLAGS=-limf
```

-xHost will use the highest level of vectorization support for the processor where you are
building the program. So if you build on a system that supports AVX but then try to run on a system that only support SSE then
you'll get errors:

```
$ bash simulate.sh
running: bin/run_md -N 125

Please verify that both the operating system and the processor support Intel(R) AVX1 instructions.
.
.
.
```

A copy of the Makefile is included in this repo.

## Benchmark Instructions

It's important to run all benchmarks from the same host, otherwise you will not get an apples-to-apples comparison. Specifically,
you should run your benchmarks from one of the 7 systems below, which are assigned based on your birth month. 

* vmps11 (Jan - Feb)
* vmps12 (Mar - Apr)
* vmps13 (May - June)
* vmp902 (July - Aug)
* vmp903 (Sept - Oct)
* vmp904 (Nov)
* vmp905 (Dec)

Once you have logged into the cluster, you can reach any of these systems with ssh. For example:

```
$ ssh vmp902
```

You will be prompted for your cluster password when doing this.

### Verifying Correct Results

In parallel programs, the onus is often on the programmer to ensure the correctness of results. For the molecular dynamics code, you should check
the thermodyanmic output against the serial version of the code to ensure correctness. Note that floating-point operations are NOT associative [1], so
even slight changes to code can lead to completely different trajectories and thermodynamic output as the simulation proceeds. However, the initial and 
average properties should be the same:

* Check that the initial PE and Pressure values are the same.
* Check that the average and standard deviation of the properties vary by no more than +/- 20% or so (average pressure deviation may be even slightly larger). Run the bash script (https://github.com/sc3260s16/SimpleMD/blob/master/simulate.sh) to see these statistics at the end of your run. See this issue for more context: https://github.com/sc3260s16/announce/issues/11

[1] http://stackoverflow.com/questions/6430448/why-doesnt-gcc-optimize-aaaaaa-to-aaaaaa?lq=1

## MPI Implementation of SimpleMD

Starting with your [vectorized version of the program from part 2 of homework 2](https://github.com/sc3260s16/announce/tree/master/hwk2#part-2-vectorized-version-of-simplemd-13-of-grade), now implement multiprocess execution by adding MPI support into the program.

Specifically, divide the atoms in your simulation as evenly as possible across multiple processes to perform work in parallel. This will likely amount to assigning sequential slices of the outer loop in energy_force.c to each process, however you are free to implement the algorithm in parallel in any fashion you wish. Here are some hints:

* Using if statements, ensure that output that goes to the screen is only printed once by the master process (rank 0).
* Determine each process's sequential subset of atoms as generically as possible, so that work gets divided evenly regardless of the system size and number of processes. You can then pass the upper and lower bounds of the subset of atoms to functions that will be updating atomic information (force, velocity, position) in parallel by looping over the subset of atoms that a process is responsible for.
* Each process should have its own copy of the ```Atoms``` struct including data for all atoms in the entire simulation. Each process updates a slice of values that correspond to the atoms that process is reponsible for. MPI is used to update properties across all processes.
* It is acceptable to allow all processes to initialize atomic positions and velocities.
* You should not need any point-to-point communication in your program, only collective MPI communication (I used MPI_Allgatherv, MPI_Reduce, and MPI_Allreduce in my implementation). 
* Although not optimal from a performance standpoint, I will not count off if you make successive MPI calls rather than bundling data together into a single data structure using MPI derived types or a single large array. For instance, if you need to update the atomic coordinates across all MPI processes, it is acceptable to do so with three consecutive MPI calls (one for each dimension).
* You might find the ```MPI_IN_PLACE``` option helpful. This option allows you to use the same send and receive buffer (i.e. array), which is helpful if you want processes to send their data to all other processes while keeping their own data. Note that this option is supported in MPI_Allgatherv and MPI_Allreduce, but not in MPI_Reduce. See this thread for an example: http://stackoverflow.com/questions/34452947/same-send-and-receive-buffers-in-mpi
* Think carefully about how data gets used across the multiple functions that are called from within the driver() function in mddriver.c. Some data may need to be updated every time step, while other data may only need to be updated on timesteps where properties must be printed to the screen.
* Also be careful in compute_energy_and_force() force within energy_force.c. Remember that in the vectorized version of the code we are updating the force on atomj and not atomi, so each process will actually be updating forces on all atoms and not just those atoms its reponsible for in the outer loop. Each process should therefore be (a) doing any required initialization across all atoms and (b) multiplying forces by ```24.0 * len_jo->eps``` across all atoms, not just the subset that the process is responsible for.
* Make sure you understand how atomic positions and velocities are being updated (this occurs in multiple steps), and consider that this should be done in parallel.

After verifying correctness of the thermodynamic output, complete the following tasks:

* Run benchmarks for system sizes of 125, 1000, and 4000. Plot the speedup (relative to the serial, vectorized version of the code you completed in hwk2, part 2) as a function of the number of processes on a single plot. Use process counts of 1, 4, 8, and 12 (for 12-core hosts) OR 1, 4, 8, and 16 (for 16-core hosts). What host did you run the benchmarks on? Comment on the results.

* Also run benchmarks across multiple nodes. You can run these only for the 1000-atom system. To ensure availability of resources, submit a SLURM script to the ```mic``` partition. Make sure that you use ```srun``` rather than ```mpirun``` when executing your program. A SLURM script is included in this repo for reference. Plot the speedup as a function of the number of processes, and also the efficiency as a function of the number of processes. Comment on the results and hand in the results. Benchmark the following conditions:

	- 1 process, 1 node
	- 2 processes, 1 node
	- 4 processes, 1 node
	- 8 processes, 1 node
	- 16 processes, 1 node
	- 2 nodes, 1 process (i.e. "task" in SLURM) per node
	- 2 nodes, 2 processes per node
	- 2 nodes, 4 processes per node
	- 2 nodes, 8 processes per node
	- 4 nodes, 1 process per node
	- 4 nodes, 2 processes per node
	- 4 nodes, 4 processes per node

* Note about measuring execution time. You cannot trust the results of some of the timer information now that multiple processes are running simultaneously. However, the total wall time printed by rank 0 should be fairly reliable, so it's safe to use this value. If I were cruel I would ask you to use MPI's timer functions, but please don't worry about that. :-)

## Reference Materials

The following are included in this repo:

- Makefile to be used for building this program. Don't forget to do a ```make cleanall``` if you have logged into a new host.
- Bash scripts for automating execution interactively and through SLURM.
- Example SLURM script. You'll need one of these for each of the conditions above. Submit to SLURM using ```sbatch bmarks.slurm```.
#hwk3_2
