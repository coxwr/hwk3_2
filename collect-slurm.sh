#!/bin/sh

# run this after all slurm jobs finish (when running "make test")

out=bench-slurm.out
echo "num_nodes,procs_per_node,time" > $out

for numnodes in 2 4; do
  for numprocs in 1 2 4 8; do
    echo -n "$numnodes","$numprocs", >> $out
    cat bench-slurm-"$numnodes"-"$numprocs".out >> $out
  done
done
