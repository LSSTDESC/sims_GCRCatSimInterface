import GCRCatalogs
image_cat_config = GCRCatalogs.get_catalog_config('cosmoDC2_v1.1.4_image')

healpix_pixels = image_cat_config['healpix_pixels']

header = '''#!/bin/bash -l
#SBATCH -t 4:00:00
#SBATCH -q premium
#SBATCH -A m1727
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH -o sed_cache_v1.1.4_output.txt
#SBATCH -e sed_cache_v1.1.4_err.txt
#SBATCH -C haswell
'''

header += '\n#SBATCH -N %d\n\n' % len(healpix_pixels)

header += '''
python /global/common/software/lsst/common/miniconda/start-kernel-cli.py desc-stack
setup -j -r $HOME/sims_data
setup -j -r $HOME/sims_GCRCatSimInterface_master
setup -j -r $HOME/sims_catalogs
setup -j -r $HOME/sims_utils
setup -j -r $HOME/sims_catUtils
setup -j -r $HOME/sims_photUtils
setup -j -r $HOME/Twinkles
setup -j -r $HOME/throughputs

export HDF5_USE_FILE_LOCKING=FALSE

export OMP_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export MKL_NUM_THREADS=1

out_dir=${SCRATCH}/sed_cache_v1.1.4/

if [ ! -d ${out_dir} ]; then
    mkdir -p ${out_dir}
fi\n\n
'''

with open('sed_fitting_script.sl', 'w') as out_file:
    out_file.write(header)
    for hp in healpix_pixels:
        out_file.write('srun -N 1 -n 1 -c 24 \\\n')
        out_file.write('python fit_sed.py --healpix %d \\\n' % hp)
        out_file.write('--catalog cosmoDC2_v1.1.4_image \\\n')
        out_file.write('--out_dir ${out_dir} \\\n')
        out_file.write('--out_name sed_fit_%d.h5 \\\n' % hp)
        out_file.write('--n_threads 24 &\n')
        out_file.write('\n')

    out_file.write('wait\n')

exit()

#for hp in 9940 10068 10069 10195 10196 10197 10323 10324 10325 10447 10448;
#for hp in 8786 8787 8788 8913 8914 8915 8916 9042 9043 9044 9045 9169 9170 9171 9172 9298 9299 9300 9426 9427;
#do
#    srun -N 1 -n 1 \
#    python -m cProfile -o profile_${hp}_edison.sav \
#    fit_sed.py --healpix ${hp} \
#    --out_dir ${out_dir} \
#    --out_name sed_fit_${hp}.h5 \
#    --n_threads 23 &
#done

#wait


#srun -N 1 -n 1 --exclusive --mem-per-cpu 20000 \
#python -m cProfile -o profile_282447_282401.sav \
#$SIMS_GCRCATSIMINTERFACE_DIR/bin.src/generateInstCat.py \
#--config_file config_file_agn.json \
#--out_dir ${out_dir} \
#--ids 282447 282401  \
#--suppress_warnings --job_log ${job_file} &


#wait
