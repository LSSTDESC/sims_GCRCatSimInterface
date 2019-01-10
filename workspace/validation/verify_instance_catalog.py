import os
import numpy as np
import pandas as pd
import healpy

import GCRCatalogs
from GCR import GCRQuery

from lsst.sims.utils import angularSeparation
from lsst.sims.photUtils import BandpassDict, Sed, Bandpass

import argparse

import time

d_mag_max = -1.0
ct_knots = 0

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


def validate_mag(row, mag_true):
    global ct_knots
    global d_mag_max

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
        ct_knots += 1

    tot_mag = -2.5*np.log10(disk_flux+bulge_flux+knots_flux)
    d_mag = np.abs(tot_mag-mag_true)
    if d_mag>d_mag_max:
        d_mag_max = d_mag
        print('d_mag_max %e -- InstanceCatalog %e truth %e -- %d' %
              (d_mag_max, tot_mag, mag_true,index))
        #print(row)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--obs', type=int, default=None)
    parser.add_argument('--nrows', type=int, default=100000)
    parser.add_argument('--seed', type=int, default=9812)
    parser.add_argument('--cat_dir', type=str, default=None)

    args = parser.parse_args()

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



    colnames = ['obj', 'uniqueID', 'ra', 'dec', 'magnorm', 'sed', 'redshift', 'g1', 'g2',
                'kappa', 'dra', 'ddec', 'src_type', 'major', 'minor',
                'positionAngle', 'sindex', 'dust_rest', 'rest_av', 'rest_rv',
                'dust_obs', 'obs_av', 'obs_rv']

    col_types = {'magnorm': float, 'redshift': float,
                 'rest_av': float, 'rest_rv': float,
                 'sed': bytes, 'uniqueID': int}

    assert os.path.isdir(args.cat_dir)
    data_dir = os.path.join(args.cat_dir,'%.8d' % args.obs)
    assert os.path.isdir(data_dir)

    phosim_file = os.path.join(data_dir, 'phosim_cat_%d.txt' % args.obs)
    assert os.path.isfile(phosim_file)
    bandpass_name = None
    bandpass_name_list = 'ugrizy'
    with open(phosim_file, 'r') as in_file:
        for line in in_file:
            params = line.strip().split()
            if params[0] == 'filter':
                bandpass_name = bandpass_name_list[int(params[1])]

    assert bandpass_name is not None
    print('bandpass is ',bandpass_name)

    (tot_dict,
     hw_dict) = BandpassDict.loadBandpassesFromFiles()

    bandpass = hw_dict[bandpass_name]

    disk_file = os.path.join(data_dir, 'disk_gal_cat_%d.txt.gz' % args.obs)
    if not os.path.isfile(disk_file):
        raise RuntimeError("%s is not a file" % disk_file)

    bulge_file = os.path.join(data_dir, 'bulge_gal_cat_%d.txt.gz' % args.obs)
    assert os.path.isfile(bulge_file)

    knots_file = os.path.join(data_dir, 'knots_cat_%d.txt.gz' % args.obs)
    assert os.path.isfile(knots_file)

    disk_df = pd.read_csv(disk_file, delimiter=' ',
                          compression='gzip', names=colnames, dtype=col_types, nrows=None)

    disk_df['galaxy_id'] = pd.Series(disk_df['uniqueID']//1024, index=disk_df.index)
    disk_df = disk_df.set_index('galaxy_id')
    print('read disks %e' % (len(disk_df)))

    bulge_df = pd.read_csv(bulge_file, delimiter=' ',
                           compression='gzip', names=colnames, dtype=col_types, nrows=None)
    bulge_df['galaxy_id'] = pd.Series(bulge_df['uniqueID']//1024, index=bulge_df.index)
    bulge_df = bulge_df.set_index('galaxy_id')
    print('read bulges %e' % (len(bulge_df)))

    for ii in range(len(colnames)):
        colnames[ii] = colnames[ii]+'_knots'

    knots_df = pd.read_csv(knots_file, delimiter=' ',
                           compression='gzip', names=colnames, dtype=col_types, nrows=None)
    knots_df['galaxy_id'] = pd.Series(knots_df['uniqueID_knots']//1024, index=knots_df.index)
    knots_df = knots_df.set_index('galaxy_id')

    wanted_col = ['sed', 'magnorm', 'redshift', 'rest_av', 'rest_rv', 'ra', 'dec']
    galaxy_df = disk_df[wanted_col].join(bulge_df[wanted_col], how='outer',
                                         lsuffix='_disk', rsuffix='_bulge')
    for ii in range(len(wanted_col)):
        wanted_col[ii] = wanted_col[ii]+'_knots'
    galaxy_df = galaxy_df.join(knots_df[wanted_col], how='outer', rsuffix='_knots')
    print('read knots')

    valid_galaxies = np.where(np.logical_not(np.in1d(galaxy_df.index,
                                                     sprinkled_gid)))

    galaxy_df = galaxy_df.iloc[valid_galaxies]
    print('threw out sprinkled systems; left with %e' % len(galaxy_df))

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

    print('healpix list')
    print(healpix_list)
    print(ra_center, dec_center)

    hp_query = None
    for hp in healpix_list:
        print(hp)
        if hp_query is None:
            hp_query = GCRQuery('healpix_pixel==%d' % hp)
        else:
            hp_query |= GCRQuery('healpix_pixel==%d' % hp)

    print('len(galaxy_df) ',len(galaxy_df))
    print('built final df')
    cat = GCRCatalogs.load_catalog('cosmoDC2_v1.1.4_image')
    cat_qties = cat.get_quantities(['galaxy_id', 'ra', 'dec'], native_filters=[hp_query])
    print('loaded galaxy_id %e' % len(cat_qties['galaxy_id']))
    cat_dexes = np.arange(len(cat_qties['galaxy_id']), dtype=int)

    gid_max = cat_qties['galaxy_id'].max()

    valid_galaxies = np.where(galaxy_df.index.values<=gid_max)
    galaxy_df = galaxy_df.iloc[valid_galaxies]
    print('cut on gid_max')

    rng = np.random.RandomState(args.seed)
    dexes = rng.choice(galaxy_df.index.values, size=args.nrows, replace=False)

    galaxy_df = galaxy_df.loc[dexes]

    galaxy_df = galaxy_df.sort_index()
    invalid_knots = np.where(np.logical_not(np.isfinite(galaxy_df['magnorm_knots'].values.astype(np.float))))

    print('no knots %e' % len(invalid_knots[0]))

    dd = angularSeparation(ra_center, dec_center,
                           cat_qties['ra'], cat_qties['dec'])

    dd_cut = np.where(dd<(radius_deg+0.05))
    gid = cat_qties['galaxy_id'][dd_cut]
    cat_dexes = cat_dexes[dd_cut]
    print('did spatial cut %d of %d' % (len(cat_dexes), len(cat_qties['ra'])))

    in1d_valid_dexes = np.where(np.in1d(gid, galaxy_df.index.values,assume_unique=True))
    valid_dexes = cat_dexes[in1d_valid_dexes]
    gid = gid[in1d_valid_dexes]

    print('got valid_dexes')
    print(len(valid_dexes))
    print(len(galaxy_df))

    sorted_dex = np.argsort(gid)
    valid_dexes = valid_dexes[sorted_dex]

    assert len(gid) == len(galaxy_df.index.values)
    np.testing.assert_array_equal(gid[sorted_dex], galaxy_df.index.values)

    mag_name = 'mag_true_%s_lsst' % bandpass_name
    qties = cat.get_quantities(['galaxy_id', mag_name], native_filters=[hp_query])
    mags = qties[mag_name][valid_dexes]
    gid = qties['galaxy_id'][valid_dexes]

    t_start = time.time()
    for g, mag_true, (index, row) in zip(gid, mags, galaxy_df.iterrows()):
        assert g==index

        validate_mag(row, mag_true)

    print('\nall done %d' % args.obs)
    print('knots %e' % ct_knots)
    print('took %e' % (time.time()-t_start))
