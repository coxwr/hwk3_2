#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "params.h"
#include "atoms.h"

void update_positions(Atoms *,
                      const misc_params *,
                      const int *,
                      const int *,
                      const int,
                      const int);
void update_velocities(Atoms *,
                       const misc_params *,
                       const int *,
                       const int *,
                       const int,
                       const int);
void pbc(Atoms *,
         const float,
         const float,
         const int *,
         const int *,
         const int,
         const int);

#endif
