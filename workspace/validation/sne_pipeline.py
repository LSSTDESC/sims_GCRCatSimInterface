import os
import numpy as np

from ic_sne import validate_sne

import time

import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--in_list', type=str, default=None)
    parser.add_argument('--out_file', type=str, default=None)

    args = parser.parse_args()

    project_dir = os.path.join('/global/projecta/projectdirs',
                               'lsst/groups/SSim/DC2/cosmoDC2_v1.1.4')

    agn_db = os.path.join(project_dir,
                          'agn_db_mbh7_mi30_sf4.db')

    ct = 0
    with open(args.out_file, 'w') as vanishing_file:
        vanishing_file.write('# obsid bandpass snid mag c x0 x1 t0 z expmjd expmjd-t0\n')
        with open(args.in_list, 'r') as in_file:
            for line in in_file:
                p = line.strip().split()
                cat_dir = p[0]
                obsid = int(p[1])
                validate_sne(cat_dir, obsid,
                             scatter_file=None,
                             vanishing_file=vanishing_file)
            ct += 1
            if ct>=10:
                break

    print('all done with %s' % args.in_list)
