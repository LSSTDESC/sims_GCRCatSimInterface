"""
make reference catalog for DC2 run 2.1i
"""

import numpy as np
import os
import sqlite3
import argparse

import GCRCatalogs
from GCR import GCRQuery

from lsst.sims.utils import sphericalFromCartesian
from lsst.sims.utils import cartesianFromSpherical

out_dir = os.path.join(os.environ['SCRATCH'], 'dc2_run2.1i_reference')
if not os.path.isdir(out_dir):
    raise RuntimeError('\n\n%s\nis not a dir\n\n' % out_dir)

cat_name = 'cosmoDC2_v1.1.4_image'
out_name = os.path.join(out_dir, 'reference_catalog_190304.txt')

ast_err_deg = 0.0001/3600.0
ast_err_rad = np.radians(ast_err_deg)  # 0.1 milliarcsec in radians
phot_err = 0.001  # in mags

cat_config = GCRCatalogs.get_catalog_config(cat_name)
cat = GCRCatalogs.load_catalog(cat_name)

agn_db_name = '/global/projecta/projectdirs/lsst/groups/SSim/DC2'
agn_db_name = os.path.join(agn_db_name, 'cosmoDC2_v1.1.4',
                           'agn_db_mbh7_mi30_sf4.db')

if not os.path.isfile(agn_db_name):
    raise RuntimeError('\n\n%s\nis not a file\n' % agn_db_name)

with sqlite3.connect(agn_db_name) as conn:
    c = conn.cursor()
    r = c.execute('SELECT galaxy_id FROM agn_params').fetchall()
    agn_galaxy_id = np.array(r).flatten()

print('%d agn' % len(agn_galaxy_id))

healpix_list = cat_config['healpix_pixels']

if os.path.exists(out_name):
    os.unlink(out_name)

q_names = ['galaxy_id', 'ra', 'dec']
for bp in 'ugrizy':
    q_names.append('mag_%s_lsst' % bp)

rng = np.random.RandomState(4223157)

with open(out_name, 'w') as out_file:
    out_file.write('# id ra dec sigma_ra sigma_dec ra_smeared dec_smeared ')
    out_file.write('u sigma_u g sigma_g r sigma_r i sigma_i z sigma_z ')
    out_file.write('y sigma_y u_smeared g_smeared r_smeared i_smeared ')
    out_file.write('z_smeared y_smeared isresolved isagn ')
    out_file.write('properMotionRa properMotionDec parallax radialVelocity\n')

    for healpix in healpix_list:
        healpix_filter = GCRQuery('healpix_pixel==%d' % healpix)
        q = cat.get_quantities(q_names, native_filters=[healpix_filter],
                               filters=[(lambda r: r<=23.0, 'mag_r_lsst')])
        is_agn = np.in1d(q['galaxy_id'], agn_galaxy_id).astype(int)
        dec_rad = np.radians(q['dec'])
        ra_rad = np.radians(q['ra'])
        xyz = cartesianFromSpherical(ra_rad, dec_rad)

        delta_ast = rng.normal(0.0, ast_err_rad, len(ra_rad))
        delta_mag = rng.normal(0.0, phot_err, size=(6, len(ra_rad)))
        unq = q['galaxy_id']*1024 + 1

        for ii in range(len(ra_rad)):
            random_vec = rng.normal(0.0, 1.0, 3)
            parallel = np.dot(random_vec, xyz[ii])
            random_vec -= parallel*xyz[ii]
            random_vec /= np.sqrt(np.dot(random_vec, random_vec))
            cos_d = np.cos(delta_ast[ii])
            sin_d = np.sin(delta_ast[ii])
            new_xyz = cos_d*xyz[ii] + sin_d*random_vec
            new_xyz /= np.sqrt(np.dot(new_xyz, new_xyz))
            new_ra, new_dec = sphericalFromCartesian(new_xyz)
            new_ra = np.degrees(new_ra)
            new_dec = np.degrees(new_dec)
            out_file.write('%d, %.10f, %.10f, ' %
                           (unq[ii], q['ra'][ii], q['dec'][ii]))
            out_file.write('%.12f, %.12f, ' % (ast_err_deg, ast_err_deg))
            out_file.write('%.10f, %.10f, ' % (new_ra, new_dec))
            for bp in 'ugrizy':
                out_file.write('%.5f, %.3f, ' %
                               (q['mag_%s_lsst' % bp][ii], phot_err))
            for i_bp, bp in enumerate('ugrizy'):
                out_file.write('%.5f, ' %
                               (q['mag_%s_lsst' % bp][ii]+delta_mag[i_bp][ii]))
            out_file.write('0, %d, 0.0, 0.0, 0.0, 0.0\n' % is_agn[ii])
        exit()
