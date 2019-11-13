import os
import sys
import numpy as np
import argparse

if __name__ == "__main__":

    out_name_root = 'slurm_scripts/batch_script_'

    parser = argparse.ArgumentParser()
    parser.add_argument('--candidate_file', type=str,
                        default='data/wfd_obshistid_list.txt',
                        help='File listing the obsHistID to simulate')
    parser.add_argument('--n_obs', type=int, default=1000,
                        help='Number of obsHistID per slurm script '
                        '(default=1000)')
    parser.add_argument('--n_jobs', type=int, default=7,
                        help='number of jobs to start per node '
                        '(default 7)')
    parser.add_argument('--d_obs', type=int, default=28,
                        help='Number of obsHistID per node (default=28); '
                        'this will result in (d_obs/n_jobs) catalogs '
                        'being generated in serial by each of the n_jobs '
                        ' per node')
    parser.add_argument('--already_done', type=str, default=None,
                        help='file containing list of already '
                        'completed obsHistID')
    parser.add_argument('--max_obs', type=int, default=993348,
                        help='maximum allowed obsHistID '
                        '(default = 993348, corresponding to 4 yrs of survey)')
    parser.add_argument('--min_obs', type=int, default=-1,
                        help='minimum allowed obsHistID; '
                        'obsHistID must be < min_obs (default==-1)')
    parser.add_argument('--out_dir', type=str, default=None,
                        help='Where InstanceCatalogs will be written')
    parser.add_argument('--config_file_name', type=str, default=None,
                        help='Name of config file to use for InstanceCatalog '
                             'generation')
    args = parser.parse_args()

    if args.config_file_name is None:
        raise RuntimeError("Must specify configfile")

    if args.out_dir is None:
        raise RuntimeError("Must specify an output directory")

    obs_already_done = set()
    if args.already_done is not None:
        with open(args.already_done, 'r') as in_file:
            for line in in_file:
                obsid = int(line)
                obs_already_done.add(obsid)

    print('N already done %d' % (len(obs_already_done)))
    obs_hist_id = []
    with open(args.candidate_file, 'r') as in_file:
        for line in in_file:
            try:
                ii=int(line.strip().split(',')[0])
                if ((ii <= args.max_obs) and
                    (ii not in obs_already_done) and
                    (ii > args.min_obs)):

                    obs_hist_id.append(ii)
            except ValueError:
                pass

    assert len(obs_hist_id) > 0

    obs_hist_id = np.array(obs_hist_id)

    i_file_offset = 0
    for i_file, i_start in enumerate(range(0,len(obs_hist_id), args.n_obs)):
        batch_slice = slice(i_start, i_start+args.n_obs)
        batch = obs_hist_id[batch_slice]
        out_name = None
        while out_name is None or os.path.isfile(out_name):
            out_name = out_name_root+'%d.sl' % (i_file+i_file_offset)
            if os.path.isfile(out_name):
                i_file_offset += 1

        print('writing ',out_name)
        with open(out_name, 'w') as out_file:
            file_id = i_file+i_file_offset
            n_srun = int(np.ceil(len(batch)/float(args.d_obs)))
            n_hrs = 2+int(np.ceil(args.d_obs/float(args.n_jobs)))
            out_file.write('#!/bin/bash -l\n')
            out_file.write('#SBATCH ')
            out_file.write('--image=docker:lsstdesc/stack-sims-cat:w_2019_37-sims_w_2019_37-v1\n')
            out_file.write('#SBATCH -N %d\n' % n_srun)
            out_file.write('#SBATCH -t %d:00:00\n' % n_hrs)
            out_file.write('#SBATCH -o slurm_out/batch_%d_out.txt\n' % file_id)
            out_file.write('#SBATCH -e slurm_err/batch_%d_err.txt\n' % file_id)
            with open('header.txt', 'r') as in_file:
                for line in in_file:
                    out_file.write(line)

            out_file.write('\n')
            out_file.write('work_dir=%s\n' % os.environ['PWD'])

            out_file.write('\n')
            out_file.write('out_dir=%s\n' % args.out_dir)
            out_file.write('config_file=${work_dir}/%s\n' %
                           args.config_file_name)
            out_file.write("if [ ! -d ${out_dir} ]; then\n")
            out_file.write("    mkdir -p ${out_dir}\n")
            out_file.write("fi\n")
            out_file.write('\n')

            for i_0 in range(0, len(batch), args.d_obs):
                s = slice(i_0, i_0+args.d_obs)
                these_obs = batch[s]
                out_file.write('\n')
                out_file.write('srun -N 1 -n 1 -c 64 --exclusive \\\n')
                out_file.write('shifter ${work_dir}/runshift_instcat.sh \\\n')
                out_file.write('${out_dir} ')
                out_file.write('${config_file} %d' % args.n_jobs)
                for ii in these_obs:
                    out_file.write(' %d' % ii)
                out_file.write(' &\n\n')

            out_file.write('\nwait\n')
            out_file.write("\necho 'master all done for %s (%d)'\n" % (args.out_dir, file_id))
            out_file.write('date\n')

