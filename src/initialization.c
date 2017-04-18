#include "initialization.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdio.h>

//**********************************************************************
// initialize_positions() function
//   - Initializes atomic positions on a simple lattice.
//   - Arguments:
//       - myatoms: struct containing all atomic information.
//       - box_length: length of simulation cell.
//       - half_box_length: half the length of the simulation cell.
//**********************************************************************
void initialize_positions(Atoms *restrict myatoms, const float box_length,
                          const float half_box_length, const int *counts,
                          const int *displs, const int bot, const int top) {
  // determine number of atoms to place along each direction
  int n_atoms = myatoms->N;
  float float_N = (float)n_atoms;
  float atoms1d_float = cbrtf(float_N); // N^(1/3), from <math.h>
  int atoms1d_round =
      (int)atoms1d_float; // number of atoms along single direction
  float atoms1d_round_float = (float)atoms1d_round;
  if (atoms1d_float - atoms1d_round_float > 1.0E-14) {
    atoms1d_round++;
    atoms1d_round_float = (float)atoms1d_round;
  }

  // position atoms on simple 3d lattice
  const float atom_offset = box_length / atoms1d_round_float;
  const int round_sq = atoms1d_round * atoms1d_round;
  float *restrict xp = myatoms->xx, *restrict yp = myatoms->yy,
                  *restrict zp = myatoms->zz;
  #pragma vector aligned
  for (int atom_cnt = bot; atom_cnt < top; ++atom_cnt) {
    int iz = atom_cnt % atoms1d_round;
    int ix = atom_cnt / round_sq;
    int iy = atom_cnt / atoms1d_round - ix * atoms1d_round;
    xp[atom_cnt] = atom_offset * (float)ix - half_box_length;
    yp[atom_cnt] = atom_offset * (float)iy - half_box_length;
    zp[atom_cnt] = atom_offset * (float)iz - half_box_length;
  }
  /* see comments in energy_force.c about casts */
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_FLOAT, xp, (int *)counts, (int *)displs,
                 MPI_FLOAT, MPI_COMM_WORLD);
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_FLOAT, yp, (int *)counts, (int *)displs,
                 MPI_FLOAT, MPI_COMM_WORLD);
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_FLOAT, zp, (int *)counts, (int *)displs,
                 MPI_FLOAT, MPI_COMM_WORLD);
}

//**************************************************************************
// initialize_velocities() function
//   - Scale velocities randomly, and scale to target temperature of system.
//   - Arguments:
//       - myatoms: struct containing all atomic information.
//       - m_pars: struct containing misc. parameters.
//       - temp: temperature of system.
//**************************************************************************
void initialize_velocities(Atoms *myatoms, const misc_params *m_pars,
                           const float temp, const int *counts,
                           const int *displs, const int bot, const int top) {
  // temperature factor for velocity scaling
  float tempfac = 3.0 * m_pars->float_N * m_pars->kb * m_pars->xmassi * temp;

  // generate random seed based on the current time
  // for initializing the random number generator
  srand((unsigned int)time(NULL));

  const int natoms = myatoms->N;
  float *restrict vxp = myatoms->vx, *restrict vyp = myatoms->vy,
                  *restrict vzp = myatoms->vz;
  // random velocities from -1 to 1
  #pragma vector aligned
  /* can't be vectorized, but that's ok since this only happens once */
  for (int i = bot; i < top; i++) {
    vxp[i] = 2.0 * (float)rand() / (float)(RAND_MAX)-1.0;
    vyp[i] = 2.0 * (float)rand() / (float)(RAND_MAX)-1.0;
    vzp[i] = 2.0 * (float)rand() / (float)(RAND_MAX)-1.0;
  }

  float sumx = 0.0, sumy = 0.0, sumz = 0.0;
  #pragma vector aligned
  // enforce zero net momentum
  for (int i = bot; i < top; i++) {
    sumx += vxp[i];
    sumy += vyp[i];
    sumz += vzp[i];
  }
  MPI_Allreduce(MPI_IN_PLACE, &sumx, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &sumy, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &sumz, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

  #pragma vector aligned
  for (int i = bot; i < top; i++) {
    vxp[i] -= sumx / (float)natoms;
    vyp[i] -= sumy / (float)natoms;
    vzp[i] -= sumz / (float)natoms;
  }

  // scale velocities to set point temperature
  float sumvsq = 0.0;
  #pragma vector aligned
  for (int i = bot; i < top; i++) {
    sumvsq += vxp[i] * vxp[i] + vyp[i] * vyp[i] + vzp[i] * vzp[i];
  }
  MPI_Allreduce(MPI_IN_PLACE, &sumvsq, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  float scaling = sqrt(tempfac / sumvsq);
  #pragma vector aligned
  for (int i = bot; i < top; i++) {
    vxp[i] *= scaling;
    vyp[i] *= scaling;
    vzp[i] *= scaling;
  }
  /* see comments in energy_force.c about casts */
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_FLOAT, vxp, (int *)counts, (int *)displs,
                 MPI_FLOAT, MPI_COMM_WORLD);
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_FLOAT, vyp, (int *)counts, (int *)displs,
                 MPI_FLOAT, MPI_COMM_WORLD);
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_FLOAT, vzp, (int *)counts, (int *)displs,
                 MPI_FLOAT, MPI_COMM_WORLD);
}
