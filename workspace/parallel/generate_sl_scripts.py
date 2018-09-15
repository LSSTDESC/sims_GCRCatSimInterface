import json
import numpy as np

def write_script(script_name, obshistid_list):
    n_obs = len(obshistid_list)
    with open(script_name, 'w') as out_file:
        out_file.write("#!/bin/bash -l\n")
        out_file.write("#SBATCH -N %d\n" % (n_obs//2))
        out_file.write("#SBATCH -t 7:00:00\n")
        out_file.write("#SBATCH -q premium\n")
        out_file.write("#SBATCH -A m1727\n")
        out_file.write("#SBATCH --tasks-per-node=1\n")
        out_file.write("#SBATCH --cpus-per-task=2\n")
        out_file.write("#SBATCH -o %s_output.txt\n" % script_name)
        out_file.write("#SBATCH -e %s_err.txt\n" % script_name)

        out_file.write("\n\n\n")

        out_file.write("project_scratch=/global/cscratch1/sd/desc/DC2/Run2.0i/instCat\n")

        out_file.write("\n\n")

        out_file.write("source /global/common/software/lsst/cori-haswell-gcc/Run2.0p_setup_test.bash\n")
        out_file.write("setup -j -r $HOME/sims_GCRCatSimInterface_desc\n")
        out_file.write("setup -j -r $HOME/sims_catUtils\n")
        out_file.write("\n\n")
        out_file.write("export PYTHONPATH=$HOME/gcr-catalogs-desc/:$PYTHONPATH\n")
        out_file.write("\n\n")
        out_file.write("export OMP_NUM_THREADS=6\n")
        out_file.write("export NUMEXPR_NUM_THREADS=6\n")
        out_file.write("export MKL_NUM_THREADS=1\n")
        out_file.write("\n\n")
        out_file.write("out_dir=${project_scratch}180914/\n")
        out_file.write("job_file=${out_dir}%s_job_log.txt\n" % script_name)
        out_file.write("\n\n")

        out_file.write("if [ ! -d ${out_dir} ]; then\n")
        out_file.write("    mkdir -p ${out_dir}\n")
        out_file.write("fi\n")
        out_file.write("\n\n")
        out_file.write("if [ -e ${job_file} ]; then\n")
        out_file.write("    rm ${job_file}\n")
        out_file.write("fi\n")

        for ii in range(0,len(obshistid_list),2):
            out_file.write("\nsrun -N 1 -n 1 --exclusive --mem-per-cpu 20000 \ \n")
            out_file.write("python $SIMS_GCRCATSIMINTERFACE_DIR/bin.src/generateInstCat.py \ \n")
            out_file.write("--config_file config_file_edison.json \ \n")
            out_file.write("--out_dir ${out_dir} \ \n")
            out_file.write("--ids ")
            for obs in obshistid_list[ii:ii+2]:
                out_file.write("%d " % obs)
            out_file.write(" \ \n")
            out_file.write("--suppress_warnings --job_log ${job_file} & \n")

        out_file.write("wait\n")

with open('sorted_visit_list.json', 'r') as in_file:
    sorted_visit_list = json.load(in_file)
    obs_list = list([int(vv[1]) for vv in sorted_visit_list])

first_thousand = obs_list[:1000]

for ii in range(0,1000,200):
    write_script('batch_scripts/instcat_batch_%d.sl' % ii, obs_list[ii:ii+200])
