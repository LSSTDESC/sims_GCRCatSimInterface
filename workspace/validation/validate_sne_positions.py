import os
import pandas as pd
import numpy as np

import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--cat_dir', type=str, default=None,
                        help='directory containing $obsHistID/ dir')
    parser.add_argument('--obs', type=int, default=None,
                        help='obsHistID of catalog')

    args = parser.parse_args()

    inst_cat_dir = os.path.join(args.cat_dir, '%.8d' % args.obs)
    if not os.path.isdir(inst_cat_dir):
        raise RuntimeError("\n%s\nis not a dir\n" % inst_cat_dir)

    colnames = ['obj', 'uniqueID', 'ra', 'dec',
                'magnorm', 'sed', 'redshift', 'g1', 'g2',
                'kappa', 'dra', 'ddec', 'src_type',
                'major_arcsec', 'minor_arcsec',
                'pa_deg', 'sindex', 'dust_rest', 'rest_av', 'rest_rv',
                'dust_obs', 'obs_av', 'obs_rv']

    col_types = {'magnorm': float, 'redshift': float,
                 'rest_av': float, 'rest_rv': float,
                 'sed': bytes, 'uniqueID': int,
                 'major_arcsec': float, 'minor_arcsec': float,
                 'pa_deg': float}

    disk_name = os.path.join(inst_cat_dir, 'disk_gal_cat_%d.txt.gz' % args.obs)
    if not os.path.isfile(disk_name):
        raise RuntimeError('\n%s\nis not a file\n' % disk_name)
    bulge_name = os.path.join(inst_cat_dir, 'bulge_gal_cat_%d.txt.gz' %
                              args.obs)
    if not os.path.isfile(bulge_name):
        raise RuntimeError('\n%s\nis not a file\n' % bulge_name)

    disk_df = pd.read_csv(disk_name, delimiter=' ',
                          compression='gzip', names=colnames,
                          dtype=col_types, nrows=None)

    disk_df['galaxy_id'] = pd.Series(disk_df['uniqueID']//1024,
                                     index=disk_df.index)

    disk_df.set_index('galaxy_id')

    bulge_df = pd.read_csv(bulge_name, delimiter=' ',
                           compression='gzip', names=colnames,
                           dtype=col_types, nrows=None)

    bulge_df['galaxy_id'] = pd.Series(bulge_df['uniqueID']//1024,
                                      index=bulge_df.index)

    bulge_df.set_index('galaxy_id')

    wanted_col = ['pa_deg', 'major_arcsec', 'minor_arcsec']

    galaxy_df = disk_df[wanted_col].join(bulge_df[wanted_col], how='outer',
                                         lsuffix='_disk', rsuffix='_bulge')


