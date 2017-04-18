#ifndef INITIALIZATION_H
#define INITIALIZATION_H

#include "atoms.h"
#include "params.h"

void initialize_positions(Atoms *,
                          const float,
                          const float,
                          const int *,
                          const int *,
                          const int,
                          const int);
void initialize_velocities(Atoms *,
                           const misc_params *,
                           const float,
                           const int *,
                           const int *,
                           const int,
                           const int);

#endif
