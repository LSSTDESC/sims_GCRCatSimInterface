import os
import numpy as np

from ic_sne import validate_sne

import time

import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--seed', type=int, default=98)

    args = parser.parse_args()

    rng = np.random.RandomState(args.seed)

    project_dir = os.path.join('/global/projecta/projectdirs',
                               'lsst/groups/SSim/DC2/cosmoDC2_v1.1.4')

    agn_db = os.path.join(project_dir,
                          'agn_db_mbh7_mi30_sf4.db')

    already_run = set()

    parent_dir = os.path.join('/global/cscratch1/sd/desc',
                              'DC2/Run2.1i/instCat')

    assert os.path.isdir(parent_dir)

    sub_dirs = os.listdir(parent_dir)


    t_start = time.time()
    with open('sne_scatter.txt', 'w') as scatter_file:
        scatter_file.write('# SALT_mag instcat_mag d_mag\n')
        with open('sne_vanished.txt', 'w') as vanishing_file:
            vanishing_file.write('# obsid snid mag c x0 x1 t0 z expmjd expmjd-t0\n')
            while len(already_run)<200:
                target_dir = rng.choice(sub_dirs, size=1)[0]
                if not os.path.isdir(os.path.join(parent_dir, target_dir)):
                    continue
                obs_list = os.listdir(os.path.join(parent_dir, target_dir))
                obs_dir = rng.choice(obs_list, size=1)[0]
                obsid = int(obs_dir)
                if obsid in already_run:
                    continue

                cat_dir = os.path.join(parent_dir, target_dir)

                validate_sne(cat_dir, obsid,
                             scatter_file=scatter_file,
                             vanishing_file=vanishing_file)

                already_run.add(obsid)
                print('%d in %e' % (len(already_run), time.time()-t_start))
