import os
import numpy as np

from ic_sne import validate_sne

import multiprocessing

import time

def scan_sne(in_list, out_file):

    project_dir = os.path.join('/global/projecta/projectdirs',
                               'lsst/groups/SSim/DC2/cosmoDC2_v1.1.4')

    agn_db = os.path.join(project_dir,
                          'agn_db_mbh7_mi30_sf4.db')

    ct = 0
    with open(out_file, 'w') as vanishing_file:
        vanishing_file.write('# obsid bandpass snid mag c x0 x1 t0 z expmjd expmjd-t0\n')
        with open(in_list, 'r') as in_file:
            for line in in_file:
                p = line.strip().split()
                cat_dir = p[0]
                obsid = int(p[1])
                validate_sne(cat_dir, obsid,
                             scatter_file=None,
                             vanishing_file=vanishing_file)
                ct += 1

    print('all done with %s' % in_list)


if __name__ == "__main__":

    t_start = time.time()
    p_list = []
    for ii in range(14):
        in_name = 'instcat_for_sne_scanning_%d.txt' % ii
        out_name = 'sne_vanished_%d.txt' % ii
        p = multiprocessing.Process(target=scan_sne,
                                    args=(in_name, out_name))

        p.start()
        p_list.append(p)

    for p in p_list:
        p.join()
    print('that took %e' % (time.time()-t_start))
