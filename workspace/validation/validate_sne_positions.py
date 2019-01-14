import os
import sqlite3
import numpy as np
import GCRCatalogs
from GCR import GCRQuery

import argparse
import time

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--cosmoDC2', type=str,
                        default='cosmoDC2_v1.1.4_image',
                        help='extragalactic catalog to load '
                        '(default: cosmoDC2_v1.1.4_image)')
    parser.add_argument('--sn_db', type=str,
                        default=None,
                        help='SNe parameter database to check '
                        'against')

    args = parser.parse_args()

    if not os.path.isfile(args.sn_db):
        raise RuntimeError('\n%s\nis not a file\n')

    with sqlite3.connect(args.sn_db) as sn_conn:
        sn_c = sn_conn.cursor()
        query = "SELECT galaxy_id, snra_in, sndec_in FROM sne_params"
        sn_data = sn_c.execute(query).fetchall()
        sn_data = np.array(sn_data).transpose()

    cat_config = GCRCatalogs.get_catalog_config(args.cosmoDC2)
    healpix_list = np.sort(cat_config['healpix_pixels'])

    gal_cat = GCRCatalogs.load_catalog(args.cosmoDC2)

    t_start = time.time()
    for i_pix, healpix_pixel in enumerate(healpix_list):
        spatial_query = GCRQuery('healpix_pixel==%d' % healpix_pixel)

        gal_q = gal_cat.get_quantities(['galaxy_id', 'ra', 'dec',
                                        'position_angle_true',
                                        'size_disk_true',
                                        'size_minor_disk_true',
                                        'size_bulge_true',
                                        'size_minor_bulge_true'],
                                       native_filters = [spatial_query])

        sn_in_gal = np.where(np.in1d(sn_data[0], gal_q['galaxy_id']))
        sn_gid = sn_data[0][sn_in_gal]
        sn_ra = sn_data[1][sn_in_gal]
        sn_dec = sn_data[2][sn_in_gal]

        gal_in_sn = np.where(np.in1d(gal_q['galaxy_id'], sn_gid))
        for k in gal_q:
            gal_q[k] = gal_q[k][gal_in_sn]

        sorted_dex = np.argsort(gal_q['galaxy_id'])
        for k in gal_q:
            gal_q[k] = gal_q[k][sorted_dex]

        sorted_dex = np.argsort(sn_gid)
        sn_gid = sn_gid[sorted_dex]
        sn_ra = sn_ra[sorted_dex]
        sn_dec = sn_dec[sorted_dex]

        np.testing.assert_array_equal(sn_gid, gal_q['galaxy_id'])
        duration = time.time()-t_start
        per = duration/(i_pix+1)
        print('did %d in %e; %e' % (i_pix+1, duration, per))
