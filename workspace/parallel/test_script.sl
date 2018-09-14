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

out_dir=$SCRATCH/parallel_instcat_3
job_file=${out_dir}job_log.txt
project_dir=/global/projecta/projectdirs/lsst/groups/SSim/DC2/

if [ ! -d ${out_dir} ]; then
    mkdir -p ${out_dir}
fi

if [ -e ${job_file} ]; then
    rm ${job_file}
fi

srun -N 1 -n 1 --exclusive --mem-per-cpu 20000 \
python $SIMS_GCRCATSIMINTERFACE_DIR/bin.src/generateInstCat.py \
--db ${project_dir}minion_1016_desc_dithered_v4.db \
--agn_db_name ${project_dir}agn_db_mbh_7.0_m_i_30.0_cosmodc2_180912.db \
--sn_db_name ${project_dir}sne_params_wfd_cosmodc2_180911.db \
--descqa_catalog cosmoDC2_v1.0_image_addon_knots \
--out_dir ${out_dir} \
--ids 1472 1474 --fov 0.02 --suppress_warnings \
--n_jobs 2 --job_log ${job_file} &

srun -N 1 -n 1 --exclusive --mem-per-cpu 20000 \
python $SIMS_GCRCATSIMINTERFACE_DIR/bin.src/generateInstCat.py \
--db ${project_dir}minion_1016_desc_dithered_v4.db \
--agn_db_name ${project_dir}agn_db_mbh_7.0_m_i_30.0_cosmodc2_180912.db \
--sn_db_name ${project_dir}sne_params_wfd_cosmodc2_180911.db \
--descqa_catalog cosmoDC2_v1.0_image_addon_knots \
--out_dir ${out_dir} \
--ids 6826 6830 --fov 0.02 --suppress_warnings \
--n_jobs 2 --job_log ${job_file}

wait
