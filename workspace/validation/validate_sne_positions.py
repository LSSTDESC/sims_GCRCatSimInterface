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
        sn_data[0] = sn_data[0].astype(int)

    cat_config = GCRCatalogs.get_catalog_config(args.cosmoDC2)
    healpix_list = np.sort(cat_config['healpix_pixels'])

    gal_cat = GCRCatalogs.load_catalog(args.cosmoDC2)

    t_start = time.time()
    normalized_radius = []
    normalized_aa = []
    normalized_bb = []
    delta_pa = []
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
        sn_gid = sn_data[0][sn_in_gal].astype(int)
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

        for i_gal in range(len(sn_gid)):
            if gal_q['size_disk_true'][i_gal] > gal_q['size_bulge_true'][i_gal]:
                comp = 'disk'
            else:
                comp = 'bulge'

            ra = gal_q['ra'][i_gal]
            dec = gal_q['dec'][i_gal]
            pa = np.radians(gal_q['position_angle_true'][i_gal])
            # vectors will be in (RA, Dec) form
            major_axis = np.array([np.sin(pa), np.cos(pa)])
            minor_axis = np.array([np.cos(pa), -np.sin(pa)])

            sn_vec = np.array([(sn_ra[i_gal]-ra)/np.cos(np.radians(dec)),
                               sn_dec[i_gal]-dec])

            aa = np.dot(sn_vec, major_axis)*3600.0 # converting to arcsec
            bb = np.dot(sn_vec, minor_axis)*3600.0

            aa /= gal_q['size_%s_true' % comp][i_gal]
            bb /= gal_q['size_minor_%s_true' % comp][i_gal]

            sn_vec /= np.sqrt(np.dot(sn_vec, sn_vec))
            sne_pa = np.arctan2(sn_vec[0], sn_vec[1])%(2.0*np.pi)

            d_pa = np.degrees(pa-sne_pa)
            delta_pa.append(d_pa)

            rrsq = (aa**2 + bb**2)
            normalized_radius.append(rrsq)
            normalized_aa.append(aa)
            normalized_bb.append(bb)

        duration = time.time()-t_start
        per = duration/(i_pix+1)
        print('did %d in %e; %e' % (i_pix+1, duration, per))

    normalized_radius = np.sqrt(np.array(normalized_radius))
    with open('rr_out.txt', 'w') as out_file:
        for rr, aa, bb, d_pa in \
        zip(normalized_radius, normalized_aa, normalized_bb, delta_pa):
            out_file.write('%e %e %e %e\n' % (rr, aa, bb, d_pa))
