#include "energy_force.h"
#include "params.h"
#include "atoms.h"
#include "timer.h"
#include <string.h>
#include <mpi.h>
#include <stdio.h>

//************************************************************************
// compute_long_range_correction() function
//   - Calculates long range correction due to finite interaction cutoff.
//   - Arguments:
//       - len_jo: struct containing leonard jones interaction parameters.
//       - m_pars: struct containing misc. parameters.
//       - energy_long: long-range correction to energy.
//       - force_long: long-range correction to force.
//************************************************************************
void compute_long_range_correction(const lj_params *len_jo,
                                   const misc_params *m_pars,
                                   float *energy_long, float *force_long) {

  float ulongpre =
      m_pars->float_N * 8.0 * len_jo->eps * m_pars->pi * m_pars->density;
  *energy_long = ulongpre * (len_jo->sig12 / (9.0 * len_jo->rcut9) -
                             len_jo->sig6 / (6.0 * len_jo->rcut3));

  float vlongpre = 96.0 * len_jo->eps * m_pars->pi * m_pars->density;
  *force_long = -1.0 * vlongpre * (len_jo->sig12 / (9.0 * len_jo->rcut9) -
                                   len_jo->sig6 / (6.0 * len_jo->rcut3));
}

//************************************************************************
// compute_energy_and_force() function
//   - Calculates energy and force acting on each atom.
//   - Arguments:
//       - myatoms: struct containing all atomic information.
//       - len_jo: struct containing lennard jones interaction parameters.
//       - m_pars: struct containing misc. parameters.
//************************************************************************
void compute_energy_and_force(Atoms *restrict myatoms,
                              const lj_params *restrict len_jo,
                              const misc_params *restrict m_pars,
                              const int *counts, const int *displs,
                              const int bot, const int top) {
  timeit(1, 0);
  int atomi, atomj;
  float *restrict fx = myatoms->fx, *restrict fy = myatoms->fy,
                  *restrict fz = myatoms->fz;
  memset(fx, 0, sizeof(float) * myatoms->N);
  memset(fy, 0, sizeof(float) * myatoms->N);
  memset(fz, 0, sizeof(float) * myatoms->N);
  myatoms->pot_energy = 0.0;
  myatoms->virial = 0.0;

  float *restrict xx = myatoms->xx, *restrict yy = myatoms->yy,
                  *restrict zz = myatoms->zz,
                  *restrict pot = &myatoms->pot_energy,
                  *restrict virial = &myatoms->virial;
  const float *restrict sig12 = &len_jo->sig12, *restrict sig6 = &len_jo->sig6;
#ifdef __INTEL_COMPILER
#pragma vector aligned
#endif
  for (atomi = bot; atomi < top; ++atomi) {
    float *restrict fxi = fx + atomi, *restrict fyi = fy + atomi,
                    *restrict fzi = fz + atomi, *restrict xx_i = xx + atomi,
                    *restrict yy_i = yy + atomi, *restrict zz_i = zz + atomi;
#ifdef __INTEL_COMPILER
#pragma vector aligned
#endif
    for (atomj = atomi + 1; atomj < myatoms->N; ++atomj) {
      float *restrict xx_j = xx + atomj, *restrict yy_j = yy + atomj,
                      *restrict zz_j = zz + atomj;
      float xxi = minimum_image((*xx_i) - (*xx_j), m_pars->side, m_pars->sideh);
      float yyi = minimum_image((*yy_i) - (*yy_j), m_pars->side, m_pars->sideh);
      float zzi = minimum_image((*zz_i) - (*zz_j), m_pars->side, m_pars->sideh);
      const float dis2 = xxi * xxi + yyi * yyi + zzi * zzi;
      const float dis2i = 1.0 / dis2;
      const float dis6i = dis2i * dis2i * dis2i;
      const float dis12i = dis6i * dis6i;
      const float fterm = dis2i * (2.0f * (*sig12) * dis12i - (*sig6) * dis6i);
      float pot_inc, virial_dec, xft, yft, zft;
      if (dis2 <= len_jo->rcut2) {
        pot_inc = (*sig12) * dis12i - (*sig6) * dis6i;
        virial_dec = fterm * dis2;
        xft = fterm * xxi;
        yft = fterm * yyi;
        zft = fterm * zzi;
      } else {
        pot_inc = 0;
        virial_dec = 0;
        xft = 0;
        yft = 0;
        zft = 0;
      }
      *pot += pot_inc;
      *virial -= virial_dec;

      *fxi += xft;
      *fyi += yft;
      *fzi += zft;
      float *restrict fxj = fx + atomj, *restrict fyj = fy + atomj,
                      *restrict fzj = fz + atomj;
      *fxj -= xft;
      *fyj -= yft;
      *fzj -= zft;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, fx, myatoms->N, MPI_FLOAT, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, fy, myatoms->N, MPI_FLOAT, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, fz, myatoms->N, MPI_FLOAT, MPI_SUM,
                MPI_COMM_WORLD);
#ifdef __INTEL_COMPILER
#pragma simd
#pragma vector aligned
#endif
  for (atomi = bot; atomi < top; ++atomi) {
    float *restrict fxi = fx + atomi, *restrict fyi = fy + atomi,
                    *restrict fzi = fz + atomi;
    *fxi *= 24.0f * len_jo->eps;
    *fyi *= 24.0f * len_jo->eps;
    *fzi *= 24.0f * len_jo->eps;
  }
  /* the cast to non-const is required for using the "openmpi_intel14" package;
     using the openmpi_intel package makes icc fail with an unknown error, and
     the openmpi_intel14 package declares Allgatherv to only accept non-const
     counts and displs */
  /* could definitely use MPI_Iallgatherv here, but that's not available until
     mpi 1.8.8, which isn't installed for icc on these machines */
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_FLOAT, fx, (int *)counts, (int *)displs,
                 MPI_FLOAT, MPI_COMM_WORLD);
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_FLOAT, fy, (int *)counts, (int *)displs,
                 MPI_FLOAT, MPI_COMM_WORLD);
  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_FLOAT, fz, (int *)counts, (int *)displs,
                 MPI_FLOAT, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &*pot, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &*virial, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  *pot *= 4.0 * len_jo->eps;
  *virial *= 24.0 * len_jo->eps;
  timeit(1, 1);
}

//**********************************************************************
// minimum_image() function
//   - Finds the nearest images of atom i and j, and returns distance.
//   - Arguments:
//       - dist: 1d distance between atoms i and j in central sim. cell.
//       - box_length: length of simulation cell.
//       - half_box_length: half of length of simulation cell.
//**********************************************************************
float minimum_image(const float dist, const float box_length,
                    const float half_box_length) {

  float min_dist = dist;
  if (dist > half_box_length)
    min_dist = dist - box_length;
  if (dist < -half_box_length)
    min_dist = dist + box_length;
  return min_dist;
}
