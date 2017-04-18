#include "params.h"
#include "atoms.h"
#include "integrator.h"
#include "timer.h"
#include <mpi.h>

//**********************************************************************
// update_positions() function
//   - Update particle positions by numerical integration.
//   - Arguments:
//       - myatoms: struct containing all atomic information.
//       - m_pars: struct containing misc. parameters.
//**********************************************************************
void update_positions(Atoms *myatoms, const misc_params *m_pars,
                      const int *counts, const int *displs, const int bot,
                      const int top) {

  // Update particle positions with velocity Verlet algorithm
  timeit(2, 0);
  int atomi;
  const float const_fac = 0.5 * m_pars->xmassi * m_pars->dt;
  #pragma vector aligned
  for (atomi = bot; atomi < top; atomi++) {
    float *restrict xx = myatoms->xx + atomi,
                    *restrict yy = myatoms->yy + atomi,
                    *restrict zz = myatoms->zz + atomi,
                    *restrict vx = myatoms->vx + atomi,
                    *restrict vy = myatoms->vy + atomi,
                    *restrict vz = myatoms->vz + atomi,
                    *restrict fx = myatoms->fx + atomi,
                    *restrict fy = myatoms->fy + atomi,
                    *restrict fz = myatoms->fz + atomi;
    // compute factor for subsequent updates
    float vx_halfts_fac = const_fac * (*fx);
    float vy_halfts_fac = const_fac * (*fy);
    float vz_halfts_fac = const_fac * (*fz);

    // update particle positions
    *xx += (*vx) * m_pars->dt + vx_halfts_fac * m_pars->dt;
    *yy += (*vy) * m_pars->dt + vy_halfts_fac * m_pars->dt;
    *zz += (*vz) * m_pars->dt + vz_halfts_fac * m_pars->dt;

    // compute velocity at half timestep, and store it as the
    // particle velocity
    *vx += vx_halfts_fac;
    *vy += vy_halfts_fac;
    *vz += vz_halfts_fac;
  }
  /* see comments in energy_force.c about casts */
  /* i'm almost certain this would benefit from performing these all async */
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_FLOAT, myatoms->fx, (int *)counts,
                 (int *)displs, MPI_FLOAT, MPI_COMM_WORLD);
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_FLOAT, myatoms->fy, (int *)counts,
                 (int *)displs, MPI_FLOAT, MPI_COMM_WORLD);
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_FLOAT, myatoms->fz, (int *)counts,
                 (int *)displs, MPI_FLOAT, MPI_COMM_WORLD);
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_FLOAT, myatoms->vx, (int *)counts,
                 (int *)displs, MPI_FLOAT, MPI_COMM_WORLD);
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_FLOAT, myatoms->vx, (int *)counts,
                 (int *)displs, MPI_FLOAT, MPI_COMM_WORLD);
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_FLOAT, myatoms->vx, (int *)counts,
                 (int *)displs, MPI_FLOAT, MPI_COMM_WORLD);
}

//**********************************************************************
// update_velocities() function
//   - Update velocities of atoms by numerical integration.
//   - Arguments:
//       - myatoms: struct containing all atomic information.
//       - misc_params: struct containing misc. parameters.
//**********************************************************************
void update_velocities(Atoms *myatoms, const misc_params *m_pars,
                       const int *counts, const int *displs, const int bot,
                       const int top) {

  // Update particle velocities with velocity Verlet algorithm
  timeit(2, 0);
  int atomi;
  float const_fac = 0.5 * m_pars->xmassi * m_pars->dt;
#pragma vector aligned
  for (atomi = bot; atomi < top; atomi++) {
    float *restrict vx = myatoms->vx + atomi,
                    *restrict vy = myatoms->vy + atomi,
                    *restrict vz = myatoms->vz + atomi,
                    *restrict fx = myatoms->fx + atomi,
                    *restrict fy = myatoms->fy + atomi,
                    *restrict fz = myatoms->fz + atomi;
    *vx += const_fac * (*fx);
    *vy += const_fac * (*fy);
    *vz += const_fac * (*fz);
  }
  /* see comments in energy_force.c about casts */
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_FLOAT, myatoms->vx, (int *)counts,
                 (int *)displs, MPI_FLOAT, MPI_COMM_WORLD);
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_FLOAT, myatoms->vy, (int *)counts,
                 (int *)displs, MPI_FLOAT, MPI_COMM_WORLD);
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_FLOAT, myatoms->vz, (int *)counts,
                 (int *)displs, MPI_FLOAT, MPI_COMM_WORLD);
  timeit(2, 1);
}

//**********************************************************************
// pbc() function
//   - Impose periodic boundary conditions.
//   - Arguments:
//       - myatoms: struct containing all atomic information.
//       - box_length: length of simulation cell.
//       - half_box_length: half of length of simulation cell.
//**********************************************************************
void pbc(Atoms *myatoms, const float box_length, const float half_box_length,
         const int *counts, const int *displs, const int bot, const int top) {

  int atomi;
#pragma vector aligned
  for (atomi = bot; atomi < top; atomi++) {
    float *restrict xx = myatoms->xx + atomi,
                    *restrict yy = myatoms->yy + atomi,
                    *restrict zz = myatoms->zz + atomi;
    if (*xx > half_box_length) {
      *xx -= box_length;
    }
    if (*xx < -half_box_length) {
      *xx += box_length;
    }
    if (*yy > half_box_length) {
      *yy -= box_length;
    }
    if (*yy < -half_box_length) {
      *yy += box_length;
    }
    if (*zz > half_box_length) {
      *zz -= box_length;
    }
    if (*zz < -half_box_length) {
      *zz += box_length;
    }
  }
  /* see comments in energy_force.c about casts */
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_FLOAT, myatoms->xx, (int *)counts,
                 (int *)displs, MPI_FLOAT, MPI_COMM_WORLD);
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_FLOAT, myatoms->yy, (int *)counts,
                 (int *)displs, MPI_FLOAT, MPI_COMM_WORLD);
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_FLOAT, myatoms->zz, (int *)counts,
                 (int *)displs, MPI_FLOAT, MPI_COMM_WORLD);
  timeit(2, 1);
}
