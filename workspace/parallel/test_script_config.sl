#!/bin/bash -l
#SBATCH -N 2
#SBATCH -t 00:25:00
#SBATCH -A m1727
#SBATCH -q debug
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH -C haswell

source /global/common/software/lsst/cori-haswell-gcc/Run2.0p_setup_test.bash

setup -j -r $HOME/sims_GCRCatSimInterface_desc
setup -j -r $HOME/sims_catUtils

export PYTHONPATH=$HOME/gcr-catalogs-desc/:$PYTHONPATH

export OMP_NUM_THREADS=6
export NUMEXPR_NUM_THREADS=6
export MKL_NUM_THREADS=1

out_dir=$SCRATCH/parallel_instcat_config
job_file=${out_dir}job_log.txt

if [ ! -d ${out_dir} ]; then
    mkdir -p ${out_dir}
fi

if [ -e ${job_file} ]; then
    rm ${job_file}
fi

srun -N 1 -n 1 --exclusive --mem-per-cpu 20000 \
python $SIMS_GCRCATSIMINTERFACE_DIR/bin.src/generateInstCat.py \
--config_file config_file_test.json \
--out_dir ${out_dir} \
--ids 1472 1474 --suppress_warnings \
--job_log ${job_file} &

srun -N 1 -n 1 --exclusive --mem-per-cpu 20000 \
python $SIMS_GCRCATSIMINTERFACE_DIR/bin.src/generateInstCat.py \
--config_file config_file_test.json \
--out_dir ${out_dir} \
--ids 6826 6830 --suppress_warnings \
--job_log ${job_file} &

wait
