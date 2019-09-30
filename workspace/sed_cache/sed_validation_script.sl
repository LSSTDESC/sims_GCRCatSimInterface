#!/bin/bash -l
#SBATCH --image=docker:lsstdesc/stack-jupyter:prod
#SBATCH -t 0:30:00
#SBATCH -q debug
#SBATCH -A m1727
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH -o sed_validation_output.txt
#SBATCH -e sed_validation_err.txt
#SBATCH -C haswell

#SBATCH -N 25

work_dir=$HOME/sims_GCRCatSimInterface_master/workspace/sed_cache/
out_dir=$SCRATCH/sed_test_190927

if [ ! -d ${out_dir} ]; then
    mkdir -p ${out_dir}
fi

for hpid in 9302  9428  9429  9430  9431  9556  9557  9558  9559 \
        9684  9685  9686  9687  9812  9813  9814  9815  9816 \
        9940  9941  9942  9943 10069 10070 10071
do

    srun -N 1 -n 1 -c 64 shifter ${work_dir}runshift_validatesed.sh ${hpid} &

done

wait
