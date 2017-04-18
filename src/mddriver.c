#include "atoms.h"
#include "cl_parse.h"
#include "params.h"
#include "initialization.h"
#include "print_traj.h"
#include "energy_force.h"
#include "props.h"
#include "integrator.h"
#include "timer.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>

#define IF_FAIL_ERROR(expr)                                                    \
  do {                                                                         \
    if ((expr) != MPI_SUCCESS) {                                               \
      fprintf(stderr, "mpi error at %s:%d\n", __FILE__, __LINE__);             \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

// timer[0]: total simulation time
// timer[1]: energy/force computation time
// timer[2]: particle position/velocity update time
// timer[3]: thermo and trajectory print time
double timer[4]; // global to improve readability of code

//***************************************************************************
// driver() function
//   - Main body of program.
//   - Creates Atoms structure, allocates memory, initializes parameters,
//     and runs simulation.
//   - Arguments:
//       - arg_count: number of arguments passed from command line.
//       - arg_list: array of strings containing arguments from command line.
//***************************************************************************
void driver(const int arg_count, char **arg_list) {
  int nprocs, rank;
  IF_FAIL_ERROR(MPI_Comm_size(MPI_COMM_WORLD, &nprocs));
  IF_FAIL_ERROR(MPI_Comm_rank(MPI_COMM_WORLD, &rank));

  double t_0 = 0, t_f = 0;
  MPI_Barrier(MPI_COMM_WORLD); /* sync before time */
  if (0 == rank) {
    t_0 = MPI_Wtime();
  }

  // parse command line
  //      Arguments:
  //                (1) Number of particles
  //                (2) Number of timesteps
  //                (3) xyz output frequency
  //                (4) thermo output frequency
  args cl;
  MPI_Datatype args_type;
  int lengths[4] = {1, 1, 1, 1};
  MPI_Aint disps[4] = {0, sizeof(int), sizeof(int) * 2, sizeof(int) * 3};
  MPI_Datatype types[4] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT};
  MPI_Type_create_struct(4, lengths, disps, types, &args_type);
  MPI_Type_commit(&args_type);

  int counts[nprocs], displs[nprocs];
  if (0 == rank) {
    cl = parse_command_line(arg_count, arg_list);
    /* create mpi "v" vectors */
    memset(counts, 0, nprocs * sizeof(int));
    for (int cur = 0; cur < cl.N; ++cur) {
      ++counts[cur % nprocs];
    }
    displs[0] = 0;
    for (int cur = 0; cur < nprocs - 1; ++cur) {
      displs[cur + 1] = displs[cur] + counts[cur];
    }
  }
  /* broadcast params */
  MPI_Bcast(&cl, 1, args_type, 0, MPI_COMM_WORLD);
  MPI_Bcast(counts, nprocs, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(displs, nprocs, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Type_free(&args_type);

  const int bot = displs[rank];
  const int top = bot + counts[rank];

  // Allocate space to store atomic positions and velocities
  Atoms atoms; // see atoms.h for definition of Atoms data structure
  // allocate space on heap for storing atomic data
  allocate_atoms(&atoms, cl.N);

  // Specify Thermodynamic state
  float Temperature = 150.0; // temperature (K)
  float Vn = 1.0 / 0.008832; // specific volume (Ang^3/molecule)

  // Units are as follows:
  // length: Angstroms (1.0e-10 m)
  // time: fs (1.0e-15 s)
  // xmass = 1.0e-28 kg
  // energy = aJ (1.0e-18 J)
  // Temp = K

  // Parameter initialization
  lj_params lj;   // parameters for computing interaction potential
  misc_params mp; // other miscellaneous parameters needed throughout simulation
  set_params(&lj, &mp, cl.N, Vn); // load parameters

  // Print a few of the computed parameters
  if (0 == rank) {
    printf("box length: %.3e Angstrom\n", mp.side);
    printf("density: %.3e molecules/Ang^3\n", mp.density);
  }

  // Compute the long range correction in energy and force
  // These values are constant throughout the simulation and compensate for the
  // total force/energy lost due to using a finite interaction cutoff
  float ulong, vlong;
  compute_long_range_correction(&lj, &mp, &ulong, &vlong);

  // initialization of atomic positions
  initialize_positions(&atoms, mp.side, mp.sideh, counts, displs, bot, top);

  // initialization of atomic velocities
  initialize_velocities(&atoms, &mp, Temperature, counts, displs, bot, top);

  // open file for writing trajectory
  FILE *fp_out = NULL;
  if (0 == rank && cl.xyz_freq != 0) {
    fp_out = fopen("traj.xyz", "w");
  }

  // props[0]: kinetic energy
  // props[1]: potential energy
  // props[2]: total energy
  // props[3]: temperature
  // props[4]: pressure
  float props[5];

  // compute initial energy/force
  compute_energy_and_force(&atoms, &lj, &mp, counts, displs, bot, top);

  if (0 == rank) {
    printf("Beginning simulation....\n");
    print_header();
  }
  int istep;

  initialize_timer();
  timeit(0, 0); // start timer
  for (istep = 0; istep <= cl.n_timesteps; istep++) {

    if (0 == rank && istep % cl.xyz_freq == 0)
      print_xyz(fp_out, &atoms); // prints to traj.xyz in current directory

    if (0 == rank && (istep % cl.thermo_freq == 0 || istep == cl.n_timesteps)) {
      calc_props(&atoms, &mp, ulong, vlong, props);
      print_props(props, istep); // prints thermo output to screen
    }

    // update particle positions for next timestep
    update_positions(&atoms, &mp, counts, displs, bot, top);
    // impose periodic boundary conditions
    pbc(&atoms, mp.side, mp.sideh, counts, displs, bot, top);
    // compute energy/force for next timestep
    compute_energy_and_force(&atoms, &lj, &mp, counts, displs, bot, top);
    // update particle velocities for next timestep
    update_velocities(&atoms, &mp, counts, displs, bot, top);
  }
  timeit(0, 1); // end timer

  if (0 == rank) {
    printf("Simulation Complete!\n");
    print_timer();
    if (cl.xyz_freq != 0) {
      fclose(fp_out);
    }
  }

  free_atoms(&atoms);

  MPI_Barrier(MPI_COMM_WORLD); /* sync before time */
  if (0 == rank) {
    t_f = MPI_Wtime();
    printf("Time taken: %f seconds\n", t_f - t_0);
    fprintf(stderr, "%f\n", t_f - t_0);
  }
}

//***********************************************************************
// main() function
//   - Required in all C programs.
//   - Execution of entire program begins here.
//   - Arguments:
//       - argc: number of arguments passed from command line.
//       - argv: array of strings containing arguments from command line.
//***********************************************************************
int main(int argc, char **argv) {
  IF_FAIL_ERROR(MPI_Init(&argc, &argv));

  driver(argc, argv);

  assert(sizeof(unsigned int) == sizeof(float));

  IF_FAIL_ERROR(MPI_Finalize());

  return 0;
}
