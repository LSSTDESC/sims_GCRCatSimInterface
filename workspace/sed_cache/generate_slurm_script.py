import sys
import numpy as np


if __name__ == "__main__":

    n_obs = int(sys.argv[1])
    d_obs = int(sys.argv[2])
    obs_hist_id = []
    with open('obshistid_in_ddf.txt', 'r') as in_file:
        for line in in_file:
            ii=int(line.strip())
            obs_hist_id.append(ii)

    obs_hist_id = np.array(obs_hist_id)
    rng = np.random.RandomState(88123)
    obs_to_simulate = rng.choice(obs_hist_id, n_obs, replace=False)

    with open('instance_catalog_generation.sl', 'w') as out_file:
        with open('header.txt', 'r') as in_file:
            for line in in_file:
                out_file.write(line)

        out_file.write('\n')

        batches = len(obs_to_simulate)/d_obs
        print('batches %e' % batches)
        for i_start in range(0, len(obs_to_simulate), d_obs):
            s = slice(i_start, i_start+d_obs)
            these_obs = obs_to_simulate[s]
            out_file.write('\n')
            out_file.write('srun -N 1 -n 1 -c 24 --exclusive \\\n')
            out_file.write('bash instcat_runner.sh ${out_dir} ${config_file}')
            for ii in these_obs:
                out_file.write(' %d' % ii)
            out_file.write(' &\n\n')

        out_file.write('\nwait\n')            
