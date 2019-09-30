#!/bin/bash

export OMP_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export MKL_NUM_THREADS=1

export HDF5_USE_FILE_LOCKING=FALSE

source /opt/lsst/software/stack/loadLSST.bash
setup lsst_distrib
setup lsst_sims
setup -j -r $HOME/sims_GCRCatSimInterface_master/
setup -j -r $HOME/throughputs

echo $OMP_NUM_THREADS
echo $NUMEXPR_NUM_THREADS
echo $MKL_NUM_THREADS
echo 'starting to process'

python validate_fitting.py --healpix $1 --seed 119 \
--in_dir $SCRATCH/sed_test_190927 --d_gal 20000 \
--nsamples 2000000 --n_threads 60

