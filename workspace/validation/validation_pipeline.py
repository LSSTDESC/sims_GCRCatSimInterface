import os
import numpy as np

from ic_mags import validate_instance_catalog_magnitudes
from ic_pos import validate_instance_catalog_positions
from ic_agn import validate_agn_mags

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
    while len(already_run)<4:
        target_dir = rng.choice(sub_dirs, size=1)[0]
        obs_list = os.listdir(os.path.join(parent_dir, target_dir))
        obs_dir = rng.choice(obs_list, size=1)[0]
        obsid = int(obs_dir)
        if obsid in already_run:
            continue

        cat_dir = os.path.join(parent_dir, target_dir)
        print('\n\n%s\n\n' % cat_dir)
        #validate_instance_catalog_magnitudes(cat_dir, obsid, nrows=-1)
        validate_instance_catalog_positions(cat_dir, obsid, 2.1)
        validate_agn_mags(cat_dir, obsid, agn_db)
        already_run.add(obsid)
