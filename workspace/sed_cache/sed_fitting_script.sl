#!/bin/bash -l
#SBATCH --image=docker:lsstdesc/stack-jupyter:prod
#SBATCH -t 4:00:00
#SBATCH -q regular
#SBATCH -A m1727
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH -o sed_cache_v1.1.4_output.txt
#SBATCH -e sed_cache_v1.1.4_err.txt
#SBATCH -C haswell

#SBATCH -N 3

work_dir=$HOME/sims_GCRCatSimInterface_master/workspace/sed_cache/
out_dir=$SCRATCH/sed_test_190925

if [ ! -d ${out_dir} ]; then
    mkdir -p ${out_dir}
fi

srun -N 1 -n 1 -c 64 shifter ${work_dir}runshift_fitsed.sh $out_dir \
test_name_8913.h5 60 8913 &

srun -N 1 -n 1 -c 64 shifter ${work_dir}runshift_fitsed.sh $out_dir \
test_name_8794.h5 60 8794 &

srun -N 1 -n 1 -c 64 shifter ${work_dir}runshift_fitsed.sh $out_dir \
test_name_8793.h5 60 8793 &

wait
