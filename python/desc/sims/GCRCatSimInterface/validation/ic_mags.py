import os
import numpy as np
import pandas as pd
import healpy

import multiprocessing

import GCRCatalogs
from GCR import GCRQuery

from lsst.sims.utils import angularSeparation
from lsst.sims.photUtils import BandpassDict, Sed, Bandpass

import argparse

import time

__all__ = ["validate_instance_catalog_magnitudes"]

def get_sed(name, magnorm, redshift, av, rv):
    if not hasattr(get_sed, '_rest_dict'):
        get_sed._rest_dict = {}
        get_sed._imsim_bp = Bandpass()
        get_sed._imsim_bp.imsimBandpass()
        get_sed._sed_dir = os.environ['SIMS_SED_LIBRARY_DIR']
        get_sed._ccm_w = None

    tag = '%s_%.2f_%.2f' % (name, av, rv)
    if tag not in get_sed._rest_dict:
        ss = Sed()
        ss.readSED_flambda(os.path.join(get_sed._sed_dir, name))
        if get_sed._ccm_w is None or not np.array_equal(ss.wavelen, get_sed._ccm_w):
            get_sed._ccm_w = np.copy(ss.wavelen)
            get_sed._ax, get_sed._bx = ss.setupCCM_ab()

        mn = ss.calcMag(get_sed._imsim_bp)
        ss.addDust(get_sed._ax, get_sed._bx, A_v=av, R_v=rv)
        get_sed._rest_dict[tag] = (ss, mn)

    base_sed = get_sed._rest_dict[tag][0]
    ss = Sed(wavelen = base_sed.wavelen, flambda=base_sed.flambda)
    dmag = magnorm-get_sed._rest_dict[tag][1]
    fnorm = np.power(10.0,-0.4*dmag)
    ss.multiplyFluxNorm(fnorm)
    ss.redshiftSED(redshift, dimming=True)
    return ss


def validate_mag(row, mag_true, bandpass):

    if np.isnan(row['magnorm_disk']):
        disk_flux = 0.0
    else:
        ss = get_sed(row['sed_disk'], row['magnorm_disk'],
                     row['redshift_disk'], row['rest_av_disk'],
                     row['rest_rv_disk'])
        disk_mag = ss.calcMag(bandpass)
        disk_flux = np.power(10.0,-0.4*disk_mag)

    if np.isnan(row['magnorm_bulge']):
        bulge_flux = 0.0
    else:
        ss = get_sed(row['sed_bulge'], row['magnorm_bulge'],
                     row['redshift_bulge'], row['rest_av_bulge'],
                     row['rest_rv_bulge'])
        bulge_mag = ss.calcMag(bandpass)
        bulge_flux = np.power(10.0,-0.4*bulge_mag)

    if np.isnan(row['magnorm_knots']):
        knots_flux = 0.0
    else:
        ss = get_sed(row['sed_knots'], row['magnorm_knots'],
                     row['redshift_knots'], row['rest_av_knots'],
                     row['rest_rv_knots'])
        knots_mag = ss.calcMag(bandpass)
        knots_flux = np.power(10.0,-0.4*knots_mag)

    tot_mag = -2.5*np.log10(disk_flux+bulge_flux+knots_flux)
    d_mag = np.abs(tot_mag-mag_true)
    return d_mag

def validate_batch(mag_true_arr, galaxy_arr, bandpass, out_dict):
    pid = os.getpid()
    d_mag_max = 0.0
    for mag_true, (index, galaxy) in zip(mag_true_arr, galaxy_arr.iterrows()):
        d_mag = validate_mag(galaxy, mag_true, bandpass)
        if d_mag>d_mag_max:
            d_mag_max = d_mag

    out_dict[pid] = d_mag_max


def validate_instance_catalog_magnitudes(cat_dir, obsid, seed=99, nrows=-1):
    """
    Parameters
    ----------
    cat_dir is the parent dir of $obsid

    obsid is the obsHistID of the pointing

    seed is the seed for a random number generator

    nrows is the number of galaxies to test (if <0, test all of them)
    """
    agn_dtype = np.dtype([('galaxy_id', int), ('twinkles_id', int)])
    agn_cache = np.genfromtxt(os.path.join(os.environ['TWINKLES_DIR'], 'data',
                                           'cosmoDC2_v1.1.4_agn_cache.csv'),
                              dtype=agn_dtype,
                              delimiter=',',
                              skip_header=1)

    sne_cache = np.genfromtxt(os.path.join(os.environ['TWINKLES_DIR'], 'data',
                                           'cosmoDC2_v1.1.4_sne_cache.csv'),
                              dtype=agn_dtype,
                              delimiter=',',
                              skip_header=1)

    sprinkled_gid = np.append(agn_cache['galaxy_id'], sne_cache['galaxy_id'])



    colnames = ['obj', 'uniqueID', 'ra', 'dec',
                'magnorm', 'sed', 'redshift', 'g1', 'g2',
                'kappa', 'dra', 'ddec', 'src_type', 'major', 'minor',
                'positionAngle', 'sindex', 'dust_rest', 'rest_av', 'rest_rv',
                'dust_obs', 'obs_av', 'obs_rv']

    to_drop = ['obj', 'g1', 'g2', 'kappa', 'dra', 'ddec', 'src_type',
               'major', 'minor', 'positionAngle', 'sindex', 'dust_rest',
               'dust_obs']

    col_types = {'magnorm': float, 'redshift': float,
                 'rest_av': float, 'rest_rv': float,
                 'sed': bytes, 'uniqueID': int}

    assert os.path.isdir(cat_dir)
    data_dir = os.path.join(cat_dir,'%.8d' % obsid)
    if not os.path.isdir(data_dir):
        raise RuntimeError('\n\n%s\nis not a dir\n\n' % data_dir)

    phosim_file = os.path.join(data_dir, 'phosim_cat_%d.txt' % obsid)
    assert os.path.isfile(phosim_file)
    bandpass_name = None
    bandpass_name_list = 'ugrizy'
    with open(phosim_file, 'r') as in_file:
        for line in in_file:
            params = line.strip().split()
            if params[0] == 'filter':
                bandpass_name = bandpass_name_list[int(params[1])]

    assert bandpass_name is not None

    (tot_dict,
     hw_dict) = BandpassDict.loadBandpassesFromFiles()

    bandpass = hw_dict[bandpass_name]

    disk_file = os.path.join(data_dir, 'disk_gal_cat_%d.txt.gz' % obsid)
    if not os.path.isfile(disk_file):
        raise RuntimeError("%s is not a file" % disk_file)

    bulge_file = os.path.join(data_dir, 'bulge_gal_cat_%d.txt.gz' % obsid)
    assert os.path.isfile(bulge_file)

    knots_file = os.path.join(data_dir, 'knots_cat_%d.txt.gz' % obsid)
    assert os.path.isfile(knots_file)

    print('reading disks')
    disk_df = pd.read_csv(disk_file, delimiter=' ',
                          compression='gzip', names=colnames, dtype=col_types, nrows=None)
    disk_df.drop(labels=to_drop, axis='columns', inplace=True)
    print('read disks')

    disk_df['galaxy_id'] = pd.Series(disk_df['uniqueID']//1024,
                                     index=disk_df.index)
    disk_df = disk_df.set_index('galaxy_id')

    print('reading bulges')
    bulge_df = pd.read_csv(bulge_file, delimiter=' ',
                           compression='gzip', names=colnames,
                           dtype=col_types, nrows=None)

    bulge_df.drop(labels=to_drop, axis='columns', inplace=True)
    print('read bulges')

    bulge_df['galaxy_id'] = pd.Series(bulge_df['uniqueID']//1024,
                                      index=bulge_df.index)
    bulge_df = bulge_df.set_index('galaxy_id')

    for ii in range(len(colnames)):
        colnames[ii] = colnames[ii]+'_knots'
    for ii in range(len(to_drop)):
        to_drop[ii] = to_drop[ii]+'_knots'

    print('reading knots')
    knots_df = pd.read_csv(knots_file, delimiter=' ',
                           compression='gzip', names=colnames,
                           dtype=col_types, nrows=None)
    knots_df.drop(labels=to_drop, axis='columns', inplace=True)
    print('read knots')

    knots_df['galaxy_id'] = pd.Series(knots_df['uniqueID_knots']//1024,
                                      index=knots_df.index)
    knots_df = knots_df.set_index('galaxy_id')

    wanted_col = ['sed', 'magnorm', 'redshift',
                  'rest_av', 'rest_rv', 'ra', 'dec']

    galaxy_df = disk_df[wanted_col].join(bulge_df[wanted_col], how='outer',
                                         lsuffix='_disk', rsuffix='_bulge')

    for ii in range(len(wanted_col)):
        wanted_col[ii] = wanted_col[ii]+'_knots'
    galaxy_df = galaxy_df.join(knots_df[wanted_col],
                               how='outer',
                               rsuffix='_knots')

    valid_galaxies = np.where(np.logical_not(np.in1d(galaxy_df.index,
                                                     sprinkled_gid)))

    galaxy_df = galaxy_df.iloc[valid_galaxies]

    ra_center = np.nanmedian(galaxy_df['ra_disk'].values)
    dec_center = np.nanmedian(galaxy_df['dec_disk'].values)

    dd = angularSeparation(ra_center, dec_center,
                           galaxy_df['ra_disk'].values,
                           galaxy_df['dec_disk'].values)
    radius_deg = np.nanmax(dd)
    ra_rad = np.radians(ra_center)
    dec_rad = np.radians(dec_center)
    vv = np.array([np.cos(ra_rad)*np.cos(dec_rad),
                   np.sin(ra_rad)*np.cos(dec_rad),
                   np.sin(dec_rad)])

    healpix_list = healpy.query_disc(32, vv, np.radians(radius_deg),
                                     nest=False, inclusive=True)

    gal_id_values = galaxy_df.index.values

    cat = GCRCatalogs.load_catalog('cosmoDC2_v1.1.4_image')
    cat_qties = {}
    cat_qties['galaxy_id'] = []
    cat_qties['ra'] = []
    cat_qties['dec'] = []
    for hp in healpix_list:
        hp_query = GCRQuery('healpix_pixel==%d' % hp)
        local_qties = cat.get_quantities(['galaxy_id', 'ra', 'dec'],
                                       native_filters=[hp_query])
        valid = np.in1d(local_qties['galaxy_id'], gal_id_values)
        if valid.any():
            for k in local_qties:
                cat_qties[k].append(local_qties[k][valid])

    for k in cat_qties:
        cat_qties[k] = np.concatenate(cat_qties[k])

    cat_dexes = np.arange(len(cat_qties['galaxy_id']), dtype=int)

    if nrows>0:
        rng = np.random.RandomState(seed)
        dexes = rng.choice(galaxy_df.index.values,
                           size=nrows,
                           replace=False)
        galaxy_df = galaxy_df.loc[dexes]

    galaxy_df = galaxy_df.sort_index()
    invalid_knots = np.where(np.logical_not(
                             np.isfinite(
                             galaxy_df['magnorm_knots'].values.astype(np.float))))

    dd = angularSeparation(ra_center, dec_center,
                           cat_qties['ra'], cat_qties['dec'])

    dd_cut = np.where(dd<(radius_deg+0.05))
    gid = cat_qties['galaxy_id'][dd_cut]
    cat_dexes = cat_dexes[dd_cut]

    in1d_valid_dexes = np.where(np.in1d(gid,
                                galaxy_df.index.values,assume_unique=True))
    valid_dexes = cat_dexes[in1d_valid_dexes]
    gid = gid[in1d_valid_dexes]

    sorted_dex = np.argsort(gid)
    valid_dexes = valid_dexes[sorted_dex]

    assert len(gid) == len(galaxy_df.index.values)
    np.testing.assert_array_equal(gid[sorted_dex], galaxy_df.index.values)

    mag_name = 'mag_true_%s_lsst' % bandpass_name
    qties = {}
    qties['galaxy_id'] = []
    qties[mag_name] = []
    for hp in healpix_list:
        hp_query = GCRQuery('healpix_pixel==%d' % hp)
        local_qties = cat.get_quantities(['galaxy_id', mag_name],
                                         native_filters=[hp_query])

        valid = np.in1d(local_qties['galaxy_id'], gal_id_values)
        if valid.any():
            for k in local_qties:
                qties[k].append(local_qties[k][valid])

    for k in qties:
        qties[k] = np.concatenate(qties[k])

    np.testing.assert_array_equal(qties['galaxy_id'], cat_qties['galaxy_id'])

    mags = qties[mag_name][valid_dexes]
    gid = qties['galaxy_id'][valid_dexes]

    assert len(gid) == len(mags)
    assert len(mags) > 0
    if nrows>0:
        assert len(mags) == nrows

    t_start = time.time()
    n_proc = 3
    d_proc = len(gid)//n_proc
    mgr = multiprocessing.Manager()
    out_dict = mgr.dict()
    p_list = []
    for i_start in range(0, len(gid), d_proc):
        mag_true = mags[i_start:i_start+d_proc]
        galaxy_arr = galaxy_df.iloc[i_start:i_start+d_proc]
        p = multiprocessing.Process(target=validate_batch,
                                    args=(mag_true, galaxy_arr,
                                          bandpass, out_dict))
        p.start()
        p_list.append(p)

    for p in p_list:
        p.join()

    assert len(list(out_dict.keys())) > 0

    d_mag_max = 0.0
    for k in out_dict.keys():
        if out_dict[k] > d_mag_max:
            d_mag_max = out_dict[k]

    if d_mag_max > 1.0e-5:
        raise RuntimeError("\nobsHistID failed magnitud validation\n"
                           "d_mag_max %e" % d_mag_max)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--obs', type=int, default=None,
                        help='obsHistID')
    parser.add_argument('--nrows', type=int, default=100000,
                        help='number of (randomly chosen) rows to '
                        'validate (if negative, will validate all '
                        'rows; default=10**5)')
    parser.add_argument('--seed', type=int, default=9812,
                        help='seed for random number generator (default=9812; '
                        'only needed if nrows>0)')
    parser.add_argument('--cat_dir', type=str, default=None,
                        help='parent directory of $obsHistID/ directory')

    args = parser.parse_args()

    validate_instance_catalog_magnitudes(args.cat_dir, args.obs,
                                         seed=args.seed,
                                         nrows=args.nrows)
