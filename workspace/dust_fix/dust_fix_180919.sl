#!/bin/bash -l
#SBATCH -N 8
#SBATCH -t 3:00:00
#SBATCH -q premium
#SBATCH -A m1727
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH -o dust_fix_180919_output.txt
#SBATCH -e dust_fix_180919_err.txt

date

source /global/common/software/lsst/cori-haswell-gcc/Run2.0p_setup_test.bash
out_dir=/global/cscratch1/sd/desc/DC2/Run2.0i/instCat/fixed_dust_180919/
bin_dir=/global/homes/d/danielsf/sims_GCRCatSimInterface_desc/bin.src/

if [ ! -d ${out_dir} ]; then
    mkdir -p ${out_dir}
fi

date


srun -N 1 -n 1 --exclusive \
python ${bin_dir}instance_crawler.py \
data_180919/inst_cat_list_0.txt \
${out_dir} &


srun -N 1 -n 1 --exclusive \
python ${bin_dir}instance_crawler.py \
data_180919/inst_cat_list_1.txt \
${out_dir} &


srun -N 1 -n 1 --exclusive \
python ${bin_dir}instance_crawler.py \
data_180919/inst_cat_list_2.txt \
${out_dir} &


srun -N 1 -n 1 --exclusive \
python ${bin_dir}instance_crawler.py \
data_180919/inst_cat_list_3.txt \
${out_dir} &


srun -N 1 -n 1 --exclusive \
python ${bin_dir}instance_crawler.py \
data_180919/inst_cat_list_4.txt \
${out_dir} &


srun -N 1 -n 1 --exclusive \
python ${bin_dir}instance_crawler.py \
data_180919/inst_cat_list_5.txt \
${out_dir} &


srun -N 1 -n 1 --exclusive \
python ${bin_dir}instance_crawler.py \
data_180919/inst_cat_list_6.txt \
${out_dir} &


srun -N 1 -n 1 --exclusive \
python ${bin_dir}instance_crawler.py \
data_180919/inst_cat_list_7.txt \
${out_dir} &

wait
