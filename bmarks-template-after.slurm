# -*- mode: sh; sh-shell: bash -*-

setpkgs -a intel_cluster_studio_compiler

echo "$SLURM_NNODES" >>newfile
echo "$SLURM_JOB_NUM_NODES" >>newfile
echo "$SLURM_JOB_NODELIST" >>newfile
echo "$SLURM_NODELIST" >>newfile
echo "$SLURM_NTASKS" >>newfile
echo "$SLURM_NPROCS" >>newfile
echo "PROCS:$1" >>newfile

./simulate-slurm.sh "$1" 4000 2>&1 >>err2-slurm | \
  tee -a errfile-slurm | grep -Po "^[0-9]+\\.[0-9]+$" | head -n1
