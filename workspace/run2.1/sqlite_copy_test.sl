#!/bin/bash -l
#SBATCH -N 1
#SBATCH -o sqlite_copy_out.txt
#SBATCH -e sqlite_copy_err.txt
#SBATCH -t 0:30:00
#SBATCH -q debug
#SBATCH -A m1727
#SBATCH -C haswell
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=24

date

SFD_HOME=/global/homes/d/danielsf/
SFD_SCRATCH=/global/cscratch1/sd/danielsf/

source /global/common/software/lsst/cori-haswell-gcc/Run2.0p_setup_test.bash
setup -j -r $HOME/sims_GCRCatSimInterface_sfd
setup -j -r $SFD_HOME/sims_catUtils
setup -j -r $SFD_HOME/sims_photUtils
setup -j -r $SFD_HOME/sims_utils
setup -j -r $SFD_HOME/Twinkles
setup -j -r $SFD_HOME/sims_catalogs

export PYTHONPATH=$SFD_HOME/gcr-catalogs-desc/:$PYTHONPATH

export OMP_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export MKL_NUM_THREADS=1

export HDF5_USE_FILE_LOCKING=FALSE

config_dir=$SFD_SCRATCH/test_config_production/
config_file=test_config_production_file.json

python make_config.py ${config_dir} ${config_file}

out_dir=$SFD_SCRATCH/instcat_test_190307small/
if [ ! -d ${out_dir} ]; then
    mkdir -p ${out_dir}
fi


srun -N 1 -n 1 -c 24 --exclusive \
bash instcat_runner.sh ${out_dir} ${config_file} 68022 219928 1037341 2141478 &

wait
date
