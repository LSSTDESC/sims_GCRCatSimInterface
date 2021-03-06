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

python fit_sed.py --out_dir $1 --out_name $2 --catalog cosmoDC2_v1.1.4_image \
--n_threads $3 --healpix $4
