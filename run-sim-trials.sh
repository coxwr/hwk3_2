#!/bin/sh

setpkgs -a intel_cluster_studio_compiler
export TIMEFORMAT='%E'

outfile=sim-trials
echo "size,vectorized,time" > "$outfile"

for size in 125 500 1000 2000 4000; do
  make clean && make VECTORIZE=0
  echo -n "$size",false, >> "$outfile"
  (time bin/run_md -N "$size") 2>>"$outfile" 1>/dev/null
  make clean && make VECTORIZE=1
  echo -n "$size",true, >> "$outfile"
  (time bin/run_md -N "$size") 2>>"$outfile" 1>/dev/null
done
