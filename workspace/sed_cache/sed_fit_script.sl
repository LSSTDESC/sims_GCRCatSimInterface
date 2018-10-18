#!/bin/bash -l
#SBATCH -N 2
#SBATCH -t 0:30:00
#SBATCH -q debug
#SBATCH -C haswell
#SBATCH -A m1727
#SBATCH --tasks-per-node=2
#SBATCH -o sed_fit_181017_output.txt
#SBATCH -e sed_fit_181017_err.txt

date

source /global/common/software/lsst/cori-haswell-gcc/Run2.0p_setup_test.bash
setup -j -r $HOME/sims_GCRCatSimInterface_master
setup -j -r $HOME/sims_catUtils
setup -j -r $HOME/sims_photUtils
setup -j -r $HOME/sims_utils

export OMP_NUM_THREADS=6
export NUMEXPR_NUM_THREADS=6
export MKL_NUM_THREADS=1

out_dir=${SCRATCH}/sed_cache_181017/

if [ ! -d ${out_dir} ]; then
    mkdir -p ${out_dir}
fi

date

#for hp in 9940 10068 10069 10195 10196 10197 10323 10324 10325 10447 10448;
for hp in 9940 10068 10069 10195;
do
    srun -N 1 -n 1 \
    python -m cProfile -o profile_${hp}.sav \
    fit_sed.py --healpix ${hp} \
    --out_dir ${out_dir} \
    --out_name test_${hp}.h5 \
    --n_threads 30 \
    --lim 10000 &
done

wait


#srun -N 1 -n 1 --exclusive --mem-per-cpu 20000 \
#python -m cProfile -o profile_282447_282401.sav \
#$SIMS_GCRCATSIMINTERFACE_DIR/bin.src/generateInstCat.py \
#--config_file config_file_agn.json \
#--out_dir ${out_dir} \
#--ids 282447 282401  \
#--suppress_warnings --job_log ${job_file} &


#wait
