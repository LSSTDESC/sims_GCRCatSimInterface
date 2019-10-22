#!/bin/bash -l
#SBATCH --image=docker:lsstdesc/stack-jupyter:prod
#SBATCH -t 2:00:00
#SBATCH -q regular
#SBATCH -A m1727
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH -o sed_cache_191015_output.txt
#SBATCH -e sed_cache_191015_err.txt
#SBATCH -C haswell

#SBATCH -N 106



work_dir=$HOME/sims_GCRCatSimInterface_master/workspace/sed_cache/
out_dir=${SCRATCH}/sed_cache_v1.1.4/

if [ ! -d ${out_dir} ]; then
    mkdir -p ${out_dir}
fi


for hpid in 8786 8787 8788 8789 8790 8791 8792 8793 8794 8913 \
 8914 8915 8916 8917 8918 8919 8920 8921 9042 9043 \
 9044 9045 9046 9047 9048 9049 9050 9169 9170 9171 \
 9172 9173 9174 9175 9176 9177 9178 9298 9299 9300 \
 9301 9303 9304 9305 9306 9425 9426 9427 \
 9432 9433 9434 9554 9555 \
 9560 9561 9562 9681 9682 9683 \
 9688 9689 9690 9810 9811 \
 9817 9818 9937 9938 9939 \
 9944 9945 9946 10066 10067 10068 \
 10072 10073 10074 10193 10194 10195 10196 10197 10198 10199 \
 10200 10201 10202 10321 10322 10323 10324 10325 10326 10327 \
 10328 10329 10444 10445 10446 10447 10448 10449 10450 10451 \
 10452
do
    srun -N 1 -n 1 -c 64 shifter ${work_dir}runshift_fitsed.sh $out_dir \
    sed_fit_${hpid}.h5 60 ${hpid} &
done
wait
