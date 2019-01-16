import os
import sqlite3
import numpy as np
import pandas as pd
import json

from lsst.sims.catUtils.mixins import ExtraGalacticVariabilityModels

import argparse

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

    if not os.path.isfile(args.agn_db):
        raise RuntimeError('\n%s\nis not a file\n' % args.agn_db)

    inst_cat_dir = os.path.join(args.cat_dir, '%.8d' % args.obs)
    if not os.path.isdir(inst_cat_dir):
        raise RuntimeError('\n%s\nis not a dir\n' % inst_cat_dir)

    agn_name = os.path.join(inst_cat_dir, 'agn_gal_cat_%d.txt.gz' % args.obs)
    if not os.path.isfile(agn_name):
        raise RuntimeError('\n%s\nis not a file\n' % agn_name)

    phosim_name = os.path.join(inst_cat_dir, 'phosim_cat_%d.txt' % args.obs)
    if not os.path.isfile(agn_name):
        raise RuntimeError('\n%s\nis not a file\n' % phosim_name)

    mjd = None
    bandpass = None
    with open(phosim_name, 'r') as in_file:
        for line in in_file:
            params = line.strip().split()
            if params[0] == 'mjd':
                mjd = float(params[1])
            elif params[0] == 'filter':
                bandpass = int(params[1])

            if mjd is not None and bandpass is not None:
                break

    if mjd is None:
        raise RuntimeError("Did not read MJD")

    if bandpass is None:
        raise RuntimeError("Did not read bandpass")

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

    with sqlite3.connect(args.agn_db) as agn_params_conn:
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
    d_mag = agn_simulator.applyAgn([np.arange(len(agn_gid), dtype=int)],
                                   agn_params, mjd, redshift=instcat_z)

    d_mag_instcat = instcat_magnorm - agn_magnorm
    print('max err %e' % (np.max(np.abs(d_mag-d_mag_instcat))))
    np.testing.assert_array_equal(d_mag_instcat, d_mag[bandpass])
