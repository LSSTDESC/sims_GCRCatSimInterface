#!/bin/bash

out_dir="${@:1:1}"
config_file="${@:2:1}"
n_parallel_jobs="${@:3:1}"
obsid_list=("${@:4}")
n_obsid=${#obsid_list[@]}
declare -i n_per_job=${n_obsid}/${n_parallel_jobs}

export OMP_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export MKL_NUM_THREADS=1

export HDF5_USE_FILE_LOCKING=FALSE

source /opt/lsst/software/stack/loadLSST.bash
setup lsst_distrib
setup lsst_sims
setup -j -r $HOME/sims_GCRCatSimInterface_master/
setup -j -r $HOME/throughputs
setup -j -r $HOME/Twinkles

declare -i i_start
i_start=0

while [ $i_start -lt $n_obsid ];
do
    id_list=${obsid_list[@]:${i_start}:${n_per_job}}
    i_start+=${n_per_job}

    python -W'ignore' \
    $SIMS_GCRCATSIMINTERFACE_DIR/bin.src/generateInstCat.py \
    --out_dir ${out_dir} \
    --config_file ${config_file} \
    --ids ${id_list[@]} --n_jobs 1 --suppress_warnings &
done
wait
