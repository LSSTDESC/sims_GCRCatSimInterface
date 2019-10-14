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

python validate_fitting.py --healpix $1 --seed 119 \
--in_dir /global/projecta/projectdirs/lsst/groups/SSim/DC2/cosmoDC2_v1.1.4/sedLookup \
--d_gal 20000 \
--nsamples -1 --n_threads 60

