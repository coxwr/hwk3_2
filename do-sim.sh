#!/bin/bash

out=bench.out
outserial=bench-serial.out

echo 'size,procnum,elapsed' > $out
echo 'size,elapsed' > $outserial

(cd ../hwk2/md-vectorized && make clean && make)
for size in 125 1000 4000; do
  time_serial="$(cd ../hwk2/md-vectorized && ./simulate.sh "$size" | \
                    grep -i 'Total sim' | grep -Po "[0-9]+\\.[0-9]+")"
  echo "serial: $size,$time_serial"
  echo "$size,$time_serial" >> $outserial
  for procnum in 1 2 4 8 16; do
    echo "$size,$procnum,?"
    time="$(./simulate.sh "$procnum" "$size" 2>&1 >>err2 | \
                tee -a errfile | \
                grep -Po "^[0-9]+\\.[0-9]+$" | head -n1)"
    echo "time: $time sec"
    echo "$size,$procnum,$time" >> $out
  done
done

# slurm time (yum)
for numnodes in 1 2 4; do
  for numprocs in 1 2 4 8; do
    overall_num_procs="$(echo "$numnodes" '*' "$numprocs" | bc)"
    outfile=bmarks-"$numnodes"-"$numprocs".slurm
    cat bmarks-template-before.slurm > "$outfile"
    # specifying -N,--ntasks-per-node on command line doesn't seem to work, so
    # generate file
    echo "#SBATCH --nodes=$numnodes" >> "$outfile"
    echo "#SBATCH --ntasks-per-node=$numprocs" >> "$outfile"
    echo "#SBATCH --output=bench-slurm-$numnodes-$numprocs.out" >> "$outfile"
    cat bmarks-template-after.slurm >> "$outfile"
    sbatch "$outfile" "$overall_num_procs"
  done
done
