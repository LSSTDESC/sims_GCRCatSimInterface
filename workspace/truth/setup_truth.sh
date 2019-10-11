export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
setup lsst_sims -t b4271
setup -j -r /local/lsst/danielsf/throughputs
setup -j -r /local/lsst/danielsf/sims_utils
