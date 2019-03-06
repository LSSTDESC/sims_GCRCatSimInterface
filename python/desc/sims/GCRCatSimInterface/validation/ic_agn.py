import os
import sqlite3
import numpy as np
import pandas as pd
import json

import GCRCatalogs
from GCR import GCRQuery
import healpy

from lsst.sims.catUtils.mixins import ExtraGalacticVariabilityModels

import argparse

__all__ = ["validate_agn_mags"]

def validate_agn_mags(cat_dir, obsid, agn_db):
    """
    Parameters
    ----------
    cat_dir is the parent dir of $obsid

    obsid is the obsHistID of the pointing (an int)

    agn_db is the database of AGN parameters
    """
    if not os.path.isfile(agn_db):
        raise RuntimeError('\n%s\nis not a file\n' % agn_db)

    inst_cat_dir = os.path.join(cat_dir, '%.8d' % obsid)
    if not os.path.isdir(inst_cat_dir):
        raise RuntimeError('\n%s\nis not a dir\n' % inst_cat_dir)

    agn_name = os.path.join(inst_cat_dir, 'agn_gal_cat_%d.txt.gz' % obsid)
    if not os.path.isfile(agn_name):
        raise RuntimeError('\n%s\nis not a file\n' % agn_name)

    phosim_name = os.path.join(inst_cat_dir, 'phosim_cat_%d.txt' % obsid)
    if not os.path.isfile(agn_name):
        raise RuntimeError('\n%s\nis not a file\n' % phosim_name)

    bandpass = None
    vistime = None
    with open(phosim_name, 'r') as in_file:
        for line in in_file:
            params = line.strip().split()
            if params[0] == 'filter':
                bandpass = int(params[1])
            elif params[0] == 'vistime':
                vistime = float(params[1])

            if (bandpass is not None and
                vistime is not None):

                break

    if bandpass is None:
        raise RuntimeError("Did not read bandpass")

    if vistime is None:
        raise RuntimeError("Did not read vistime")

    opsim_db = os.path.join('/global/projecta/projectdirs/lsst',
                            'groups/SSim/DC2/minion_1016_desc_dithered_v4_sfd.db')

    if not os.path.isfile(opsim_db):
        raise RuntimeError('\n%s\nis not a file' % opsim_db)

    with sqlite3.connect(opsim_db) as conn:
        c = conn.cursor()
        r = c.execute('SELECT expMJD, descDitheredRA, descDitheredDec '
                      'FROM Summary WHERE obsHistID==%d' % obsid).fetchall()
        mjd = float(r[0][0])
        pointing_ra = float(r[0][1])
        pointing_dec = float(r[0][2])

    agn_colnames = ['obj', 'uniqueID', 'ra', 'dec',
                    'magnorm', 'sed', 'redshift', 'g1', 'g2',
                    'kappa', 'dra', 'ddec', 'src_type',
                    'dust_rest', 'dust_obs', 'obs_av', 'obs_rv']

    agn_col_types = {'ra': float, 'dec':float,
                     'magnorm': float, 'redshift': float,
                     'sed': bytes, 'uniqueID': int}

    agn_df = pd.read_csv(agn_name, delimiter=' ',
                         compression='gzip', names=agn_colnames,
                         dtype=agn_col_types, nrows=None)

    agn_df['galaxy_id'] = pd.Series(agn_df['uniqueID']//1024,
                                    index=agn_df.index)

    vv = np.array([np.cos(pointing_dec)*np.cos(pointing_ra),
                   np.cos(pointing_dec)*np.sin(pointing_ra),
                   np.sin(pointing_dec)])
    hp_list = healpy.query_disc(32,vv,np.radians(2.2),nest=False,inclusive=True)

    hp_query = GCRQuery('healpix_pixel==%d' % hp_list[0])
    for hp in hp_list[1:]:
        hp_query |= GCRQuery('healpix_pixel==%d' % hp)

    with sqlite3.connect(agn_db) as agn_params_conn:
        agn_params_cursor = agn_params_conn.cursor()
        query = 'SELECT galaxy_id, magNorm, varParamStr FROM agn_params'
        agn_params = agn_params_cursor.execute(query).fetchall()
        agn_params = np.array(agn_params).transpose()
        agn_gid = agn_params[0].astype(int)
        agn_magnorm = agn_params[1].astype(float)
        agn_varParamStr = agn_params[2]
        valid_agn = np.where(np.in1d(agn_gid, agn_df['galaxy_id'].values))
        agn_gid = agn_gid[valid_agn]
        agn_magnorm = agn_magnorm[valid_agn]
        agn_varParamStr = agn_varParamStr[valid_agn]

        del agn_params

    sorted_dex = np.argsort(agn_gid)
    agn_gid = agn_gid[sorted_dex]
    agn_magnorm = agn_magnorm[sorted_dex]
    agn_varParamStr = agn_varParamStr[sorted_dex]

    instcat_gid = agn_df['galaxy_id'].values
    instcat_magnorm = agn_df['magnorm'].values
    instcat_z = agn_df['redshift'].values

    valid = np.where(instcat_gid<1.0e11)
    instcat_gid = instcat_gid[valid]
    instcat_magnorm = instcat_magnorm[valid]
    instcat_z = instcat_z[valid]
    sorted_dex = np.argsort(instcat_gid)
    instcat_gid = instcat_gid[sorted_dex]
    instcat_magnorm = instcat_magnorm[sorted_dex]
    instcat_z = instcat_z[sorted_dex]

    cat = GCRCatalogs.load_catalog('cosmoDC2_v1.1.4_image')
    cat_q = cat.get_quantities(['galaxy_id', 'redshift_true'], native_filters=[hp_query])

    valid = np.in1d(cat_q['galaxy_id'], agn_df['galaxy_id'])
    for k in cat_q:
        cat_q[k] = cat_q[k][valid]

    sorted_dex = np.argsort(cat_q['galaxy_id'])
    for k in cat_q:
        cat_q[k] = cat_q[k][sorted_dex]

    if not np.array_equal(cat_q['galaxy_id'], instcat_gid):
        msg = "GCR gid not equal to InstCat\n"
        msg += "len gcr %d\n" % len(cat_q['galaxy_id'])
        msg += "len instcat %d\n" % len(instcat_gid)
        msg += "other comparison %s\n" % str(np.array_equal(instcat_gid, agn_gid))
        raise RuntimeError(msg)

    if not np.array_equal(instcat_gid, agn_gid):
        raise RuntimeError("galaxy_id arrays are not equal")

    if len(instcat_gid) == 0:
        raise RuntimeError("no AGN to test")

    agn_params = None
    for var in agn_varParamStr:
        var_dict = json.loads(var)
        if agn_params is None:
            agn_params = {}
            for k in var_dict['p']:
                agn_params[k] = []
        for k in var_dict['p']:
            agn_params[k].append(var_dict['p'][k])

    for k in agn_params:
        agn_params[k] = np.array(agn_params[k])

    agn_simulator = ExtraGalacticVariabilityModels()
    agn_simulator._agn_threads = 3
    d_mag = agn_simulator.applyAgn([np.arange(len(agn_gid), dtype=int)],
                                   agn_params, mjd, redshift=cat_q['redshift_true'])

    d_mag_instcat = instcat_magnorm - agn_magnorm
    error = np.abs(d_mag[bandpass]-d_mag_instcat)
    max_error = error.max()
    violation = np.where(error>1.0e-5)
    for ii in violation[0]:
        print("%e -- %e %e %e" % (error[ii], d_mag[bandpass][ii],
                                  d_mag_instcat[ii],
                                  instcat_magnorm[ii]))

        for k in agn_params:
            print('    %s: %e' % (k,agn_params[k][ii]))

    valid = np.where(error<=1.0e-5)
    d_mag_valid = d_mag_instcat[valid]
    mag_valid = instcat_magnorm[valid]

    if np.max(error)>1.0e-5:
        raise RuntimeError("\n%s\nAGN validation failed: max mag error %e" %
                           (agn_name, max_error))


if __name__ == "__main__":

    project_dir = os.path.join('/global/projecta/projectdirs',
                               'lsst/groups/SSim/DC2/cosmoDC2_v1.1.4')

    default_agn_db = os.path.join(project_dir,
                                  'agn_db_mbh7_mi30_sf4.db')

    parser = argparse.ArgumentParser()
    parser.add_argument('--cat_dir', type=str, default=None,
                        help='parent directory of $obsHistID/')
    parser.add_argument('--obs', type=int, default=None,
                        help='obsHistID of pointing')
    parser.add_argument('--agn_db', type=str,
                        default=default_agn_db,
                        help='Name of agn parameters db\n'
                        '(default %s)' % default_agn_db)

    args = parser.parse_args()

    validate_agn_mags(args.cat_dir, args.obs, args.agn_db)
