"""
make reference catalog for DC2 run 2.1i
"""

import numpy as np
import os
import sqlite3
import argparse
import time

import GCRCatalogs
from GCR import GCRQuery

from lsst.sims.utils import defaultSpecMap
from lsst.sims.utils import sphericalFromCartesian
from lsst.sims.utils import cartesianFromSpherical
from lsst.sims.catUtils.dust import EBVbase
import lsst.sims.photUtils as photUtils

def smear_ra_dec(ra_deg, dec_deg, delta_ast, rng):
    xyz = cartesianFromSpherical(np.radians(ra_deg), np.radians(dec_deg))
    ra_out = np.zeros(len(ra_deg), dtype=float)
    dec_out = np.zeros(len(dec_deg), dtype=float)
    for ii in range(len(ra_deg)):
        random_vec = rng.normal(0.0, 1.0, 3)
        parallel = np.dot(random_vec, xyz[ii])
        random_vec -= parallel*xyz[ii]
        random_vec /= np.sqrt(np.dot(random_vec, random_vec))
        cos_d = np.cos(delta_ast[ii])
        sin_d = np.sin(delta_ast[ii])
        new_xyz = cos_d*xyz[ii] + sin_d*random_vec
        new_xyz /= np.sqrt(np.dot(new_xyz, new_xyz))
        new_ra, new_dec = sphericalFromCartesian(new_xyz)
        ra_out[ii] = np.degrees(new_ra)
        dec_out[ii] = np.degrees(new_dec)
    return ra_out, dec_out

r_mag_limit = 23.0
star_db_name = '/global/projecta/projectdirs/lsst/groups'
star_db_name = os.path.join(star_db_name, 'SSim/DC2/dc2_stellar_db.db')
if not os.path.isfile(star_db_name):
    raise RuntimeError('\n\n%s\bis not a file\n' % star_db_name)

out_dir = os.path.join(os.environ['SCRATCH'], 'dc2_run2.1i_reference')
if not os.path.isdir(out_dir):
    raise RuntimeError('\n\n%s\nis not a dir\n\n' % out_dir)

cat_name = 'cosmoDC2_v1.1.4_image'
out_name = os.path.join(out_dir, 'reference_catalog_190304.txt')

ast_err_deg = 0.0001/3600.0
ast_err_rad = np.radians(ast_err_deg)/np.sqrt(2)  # 0.1 milliarcsec in radians (total)
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

chunk_size = 100000

sed_cache_dir = os.path.join(os.environ['SCRATCH'], 'lsst_sed_cache')
if not os.path.exists(sed_cache_dir):
    os.mkdir(sed_cache_dir)
photUtils.cache_LSST_seds(wavelen_min=0.0, wavelen_max=1600.0,
                          cache_dir=sed_cache_dir)

with open(out_name, 'w') as out_file:
    out_file.write('# id ra dec sigma_ra sigma_dec ra_smeared dec_smeared ')
    out_file.write('u sigma_u g sigma_g r sigma_r i sigma_i z sigma_z ')
    out_file.write('y sigma_y u_smeared g_smeared r_smeared i_smeared ')
    out_file.write('z_smeared y_smeared isresolved isagn ')
    out_file.write('properMotionRa properMotionDec parallax radialVelocity\n')

    ebv_gen = EBVbase()
    dust_wav = np.array([0.0])
    sed_dir = os.environ['SIMS_SED_LIBRARY_DIR']
    assert os.path.isdir(sed_dir)

    bp_dict = photUtils.BandpassDict.loadTotalBandpassesFromFiles()

    where_clause = 'WHERE ra>=47.72 AND ra<=75.98 '
    where_clause += 'AND decl>=-46.61 AND decl<=-24.594'

    with sqlite3.connect(star_db_name) as star_conn:
        t_start_stars = time.time()
        star_cursor = star_conn.cursor()

        final_ct = star_cursor.execute('SELECT count(simobjid) from stars '
                                       + where_clause).fetchall()[0][0]

        star_cmd = star_cursor.execute('SELECT simobjid, ra, decl, '
                                       + 'sedFilename, magNorm FROM stars '
                                       + where_clause)

        star_r = star_cmd.fetchmany(chunk_size)
        tot_stars = 0
        while len(star_r)>0:
            n_stars = len(star_r)
            tot_stars += n_stars
            print('n_stars %d tot %d of %d' % (n_stars, tot_stars, final_ct))
            star_r = np.array(star_r).transpose()
            unq = star_r[0].astype(int)*1024 + 2
            ra = star_r[1].astype(float)
            dec = star_r[2].astype(float)
            sed_name =star_r[3].astype(str)
            mag_norm = star_r[4].astype(float)
            ebv = ebv_gen.calculateEbv(equatorialCoordinates=np.array([np.radians(ra),
                                                                       np.radians(dec)]),
                                       interp=True)

            delta_ast = rng.normal(0.0, ast_err_rad, n_stars)
            (ra_smeared,
             dec_smeared) = smear_ra_dec(ra, dec, delta_ast, rng)

            delta_mag = rng.normal(0.0, phot_err, size=(n_stars, 6))

            mags = np.zeros((n_stars, 6), dtype=float)
            for ii in range(n_stars):
                spec = photUtils.Sed()
                file_name = os.path.join(sed_dir, defaultSpecMap[sed_name[ii]])
                spec.readSED_flambda(file_name)
                if not np.array_equal(spec.wavelen, dust_wav):
                    a_x, b_x = spec.setupCCMab()
                    dust_wav = np.copy(spec.wavelen)
                fnorm = photUtils.getImsimFluxNorm(spec, mag_norm[ii])
                spec.multiplyFluxNorm(fnorm)
                spec.addDust(a_x, b_x, ebv=ebv[ii], R_v=3.1)
                mags[ii] = bp_dict.magListForSed(spec)

            for ii in range(n_stars):
                if mags[ii][2]>r_mag_limit:
                    continue
                out_file.write('%d, %.10f, %.10f, ' %
                               (unq[ii], ra[ii], dec[ii]))
                out_file.write('%.12f, %.12f, ' % (ast_err_deg, ast_err_deg))
                out_file.write('%.10f, %.10f, ' % (ra_smeared[ii], dec_smeared[ii]))
                for ibp in range(6):
                    out_file.write('%.5f, %.3f, ' %
                                   (mags[ii][ibp], phot_err))
                for ibp in range(6):
                    out_file.write('%.5f, ' %
                                   (mags[ii][ibp]+delta_mag[ii][ibp]))
                out_file.write('0, 0, 0.0, 0.0, 0.0, 0.0\n')

            duration = time.time()-t_start_stars
            duration = duration/3600.0
            proj = final_ct*duration/tot_stars
            print('can do all stars in %e hrs (%e now)' %
                  (proj, duration))

            star_r = star_cmd.fetchmany(chunk_size)

    print('done with stars')
    for healpix in healpix_list:
        healpix_filter = GCRQuery('healpix_pixel==%d' % healpix)
        q = cat.get_quantities(q_names, native_filters=[healpix_filter],
                               filters=[(lambda r: r<=r_mag_limit, 'mag_r_lsst')])
        is_agn = np.in1d(q['galaxy_id'], agn_galaxy_id).astype(int)
        n_is_agn_true = len(np.where(is_agn==1)[0])
        print('gal %d is_agn %d' % (len(is_agn), n_is_agn_true))

        n_gal = len(q['ra'])
        delta_ast = rng.normal(0.0, ast_err_rad, n_gal)
        delta_mag = rng.normal(0.0, phot_err, size=(6, n_gal))
        unq = q['galaxy_id']*1024 + 1

        (ra_smeared,
         dec_smeared) = smear_ra_dec(q['ra'], q['dec'], delta_ast, rng)

        for ii in range(n_gal):
            out_file.write('%d, %.10f, %.10f, ' %
                           (unq[ii], q['ra'][ii], q['dec'][ii]))
            out_file.write('%.12f, %.12f, ' % (ast_err_deg, ast_err_deg))
            out_file.write('%.10f, %.10f, ' % (ra_smeared[ii], dec_smeared[ii]))
            for bp in 'ugrizy':
                out_file.write('%.5f, %.3f, ' %
                               (q['mag_%s_lsst' % bp][ii], phot_err))
            for i_bp, bp in enumerate('ugrizy'):
                out_file.write('%.5f, ' %
                               (q['mag_%s_lsst' % bp][ii]+delta_mag[i_bp][ii]))
            out_file.write('1, %d, 0.0, 0.0, 0.0, 0.0\n' % is_agn[ii])
