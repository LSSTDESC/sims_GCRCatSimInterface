import os
import sys
import numpy as np
import argparse

if __name__ == "__main__":

    config_file = 'config_file_2.1.wfd.json'
    out_dir_root = '/global/cscratch1/sd/desc/DC2/Run2.0i/cosmoDC2_v1.1.4/instCat/batch_'
    out_name_root = 'slurm_scripts/batch_script_'


    parser = argparse.ArgumentParser()
    parser.add_argument('--n_obs', type=int, default=1000,
                        help='Number of pointings per slurm script '
                        '(default=1000)')
    parser.add_argument('--d_obs', type=int, default=7,
                        help='Number of pointings per node (default=7)')
    args = parser.parse_args()

    obs_hist_id = []
    with open('data/wfd_obshistid_list.txt', 'r') as in_file:
        for line in in_file:
            ii=int(line.strip())
            obs_hist_id.append(ii)

    obs_hist_id = np.array(obs_hist_id)
    rng = np.random.RandomState(88123)
    rng.shuffle(obs_hist_id)

    for i_file, i_start in enumerate(range(0,len(obs_hist_id), args.n_obs)):
        batch_slice = slice(i_start, i_start+args.n_obs)
        batch = obs_hist_id[batch_slice]
        out_name = out_name_root+'%d.sl' % i_file
        out_dir = out_dir_root+'%d' % i_file
        with open(out_name, 'w') as out_file:
            n_srun = int(np.ceil(len(batch)/args.d_obs))
            out_file.write('#!/bin/bash -l\n')
            out_file.write('#SBATCH -N %d\n' % n_srun)
            out_file.write('#SBATCH -o slurm_out/batch_%d_out.txt\n' % i_file)
            out_file.write('#SBATCH -e slurm_err/batch_%d_err.txt\n' % i_file)
            with open('header.txt', 'r') as in_file:
                for line in in_file:
                    out_file.write(line)

            out_file.write('\n')
            out_file.write('out_dir=%s\n' % out_dir)
            out_file.write('config_file=%s\n' % config_file)
            out_file.write("if [ ! -d ${out_dir} ]; then\n")
            out_file.write("    mkdir -p ${out_dir}\n")
            out_file.write("fi\n")
            out_file.write('\n')

            for i_0 in range(0, len(batch), args.d_obs):
                s = slice(i_0, i_0+args.d_obs)
                these_obs = batch[s]
                out_file.write('\n')
                out_file.write('srun -N 1 -n 1 -c 24 --exclusive \\\n')
                out_file.write('bash instcat_runner.sh ${out_dir} ${config_file}')
                for ii in these_obs:
                    out_file.write(' %d' % ii)
                out_file.write(' &\n\n')

            out_file.write('\nwait\n')
