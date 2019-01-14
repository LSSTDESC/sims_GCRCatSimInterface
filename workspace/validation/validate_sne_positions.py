import os
import sqlite3
import numpy as np
import GCRCatalogs

import argparse

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

    gal_cat = GCRCatalogs.load_catalog(args.cosmoDC2)
    gal_q = gal_cat.get_quantities(['galaxy_id', 'ra', 'dec',
                                    'position_angle_true',
                                    'size_disk_true',
                                    'size_minor_disk_true',
                                    'size_bulge_true',
                                    'size_minor_bulge_true'])

    with sqlite3.connect(args.sn_db) as sn_conn:
        sn_c = sn_conn.cursor()
        query = "SELECT galaxy_id, snra_in, sndec_in FROM sne_params"
        sn_data = sn_c.execute(query).fetchall()
        sn_data = np.array(sn_data).transpose()
        sn_gid = sn_data[0]
        sn_ra = sn_data[1]
        sn_dec = sn_data[2]

    valid = np.where(np.in1d(gal_q['galaxy_id'], sn_gid))
    for k in gal_q:
        gal_q[k] = gal_q[k][valid]
