#!/usr/bin/env python
import argparse
import sqlite3
import numpy as np
import os
import gc
import time

import GCRCatalogs
from desc.sims.GCRCatSimInterface import M_i_from_L_Mass
from desc.sims.GCRCatSimInterface import log_Eddington_ratio
from desc.sims.GCRCatSimInterface import k_correction
from desc.sims.GCRCatSimInterface import tau_from_params
from desc.sims.GCRCatSimInterface import SF_from_params

from lsst.utils import getPackageDir
from lsst.sims.photUtils import Sed, BandpassDict, CosmologyObject
from lsst.sims.photUtils import Bandpass
from lsst.sims.utils import findHtmid

def create_k_corr_grid():
    """
    Returns a grid of redshifts and K corrections on the
    LSST Simulations AGN SED that can be used for K correction
    interpolation.
    """
    bp_dict = BandpassDict.loadTotalBandpassesFromFiles()
    bp_i = bp_dict['i']
    sed_dir = os.path.join(getPackageDir('sims_sed_library'),
                           'agnSED')
    sed_name = os.path.join(sed_dir, 'agn.spec.gz')
    if not os.path.exists(sed_name):
        raise RuntimeError('\n\n%s\n\nndoes not exist\n\n' % sed_name)
    base_sed = Sed()
    base_sed.readSED_flambda(sed_name)
    z_grid = np.arange(0.0, redshift.max(), 0.01)
    k_grid = np.zeros(len(z_grid),dtype=float)

    for i_z, zz in enumerate(z_grid):
        ss = Sed(flambda=base_sed.flambda, wavelen=base_sed.wavelen)
        ss.redshiftSED(zz, dimming=True)
        k = k_correction(ss, bp_i, zz)
        k_grid[i_z] = k

    return z_grid, k_grid


def get_m_i(abs_mag_i, redshift):
    """
    Take numpy arrays of absolute i-band magnitude and
    cosmological redshift.  Return a numpy array of
    observed i-band magnitudes
    """
    z_grid, k_grid = create_k_corr_grid()
    k_corr = np.interp(redshift, z_grid, k_grid)

    dc2_cosmo = CosmologyObject(H0=71.0, Om0=0.265)
    distance_modulus = dc2_cosmo.distanceModulus(redshift=redshift)
    obs_mag_i = abs_mag_i + distance_modulus + k_corr
    return obs_mag_i

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--mbh_cut', type=float, default=7.0,
                        help="log10 of the minimum black hole mass "
                        "necessary to be considered an AGN "
                        "(in solar masses).  Default=7")

    parser.add_argument('--m_i_cut', type=float, default=30.0,
                        help="Dimmest apparent i-band magnitude "
                        "to be considered an AGN.  Default = 30 "
                        "(Note: setting this to None results in "
                        "nonsensical variability parameters that "
                        "can cause some AGN to reach magnitudes "
                        "of ~ -200)")

    parser.add_argument('--out_dir', type=str, default='.',
                        help="Directory in which to write the output "
                        "database of AGN parameters.  Will be created "
                        "if it does not exist. Default = .")

    parser.add_argument('--out_file', type=str, default=None,
                        help="Name of file in which to write database "
                        "of AGN parameters.  Default = None")

    parser.add_argument('--clobber', type=str, default='False',
                        help="If 'True', will overwrite existing "
                        "out_dir/out_file.  Default='False'")

    parser.add_argument('--yaml_file', type=str,
                        default='protoDC2',
                        help="yaml file to load with GCRCatalogs. "
                        "Default = 'protoDC2'")

    parser.add_argument('--seed', type=int, default=81,
                        help="Seed value for random number generator "
                        "used to introduce scatter into AGN parameter "
                        "distributions. Default=81")

    parser.add_argument('--max_sf', type=float, default=4.0,
                        help="Maximum allowed value for the structure "
                        "function of the random walk driving AGN "
                        "variability.  Default=4 (in magnitudes)")

    args = parser.parse_args()

    if args.out_file is None:
        raise RuntimeError('Must specify an out_file')

    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)
    else:
        if not os.path.isdir(args.out_dir):
            raise RuntimeError('%s is not a dir' % args.out_dir)

    rng = np.random.RandomState(args.seed)

    out_file_name = os.path.join(args.out_dir, args.out_file)
    if os.path.exists(out_file_name):
        if not args.clobber.lower() == 'true':
            raise RuntimeError('%s already exists; clobber set to %s' %
                               (out_file_name, args.clobber))
        else:
            os.unlink(out_file_name)

    cat = GCRCatalogs.load_catalog(args.yaml_file)

    qty_list = cat.list_all_quantities(include_native=True)

    use_direct_eddington = ('blackHoleEddingtonRatio' in qty_list)

    qty_names= ['redshift_true', 'blackHoleMass', 'galaxy_id', 'ra', 'dec']
    filters = [(lambda x: x>0.0, 'blackHoleMass'),
               (lambda x: np.log10(x)>args.mbh_cut, 'blackHoleMass')]

    if use_direct_eddington:
        print('using native Eddington ratio')
        qty_names.append('blackHoleEddingtonRatio')
        filters.append((lambda x:x>0.0, 'blackHoleEddingtonRatio')

    else:
        qty_names.append('blackHoleAccretionRate')
        filters.append((lambda x: x>0.0, 'blackHoleAccretionRate'))

    cat_qties = cat.get_quantities(qty_names, filters=filters)
    if not use_direct_eddington:
        accretion_rate_full = cat_qties['blackHoleAccretionRate']

    redshift_full = cat_qties['redshift_true']
    bhm_full = cat_qties['blackHoleMass']
    galaxy_id_full = cat_qties['galaxy_id']
    ra_full = cat_qties['ra']
    dec_full = cat_qties['dec']

    if use_direct_eddington:
        log_edd_ratio_full = np.log10(cat_qties['blackHoleEddingtonRatio'])
    else:
        print('solving for Eddington ratio')
        log_edd_ratio_full = log_Eddington_ratio(bhm_full,
                                                 accretion_rate_full)

    del cat_qties
    gc.collect()

    # sort by galaxy_id so that random scatter in AGN parameters
    # is reproducible
    sorted_dex = np.argsort(galaxy_id_full)
    redshift_full = redshift_full[sorted_dex]
    bhm_full = bhm_full[sorted_dex]
    galaxy_id_full = galaxy_id_full[sorted_dex]
    ra_full = ra_full[sorted_dex]
    dec_full = dec_full[sorted_dex]
    log_edd_ratio_full = log_edd_ratio_full[sorted_dex]

    bp_dict = BandpassDict.loadTotalBandpassesFromFiles()

    sed_dir = os.path.join(getPackageDir('sims_sed_library'),
                           'agnSED')
    sed_name = os.path.join(sed_dir, 'agn.spec.gz')
    if not os.path.exists(sed_name):
        raise RuntimeError('\n\n%s\n\nndoes not exist\n\n' % sed_name)

    base_sed = Sed()
    base_sed.readSED_flambda(sed_name)

    imsimband = Bandpass()
    imsimband.imsimBandpass()

    z_grid = np.arange(0.0, redshift_full.max(), 0.01)
    m_i_grid = np.zeros(len(z_grid), dtype=float)
    mag_norm_grid = np.zeros(len(z_grid), dtype=float)
    for i_z, zz in enumerate(z_grid):
        ss = Sed(wavelen=base_sed.wavelen, flambda=base_sed.flambda)
        ss.redshiftSED(zz, dimming=True)
        m_i_grid[i_z] = ss.calcMag(bp_dict['i'])
        mag_norm_grid[i_z] = ss.calcMag(imsimband)

    htmid_level = 6
    with sqlite3.connect(out_file_name) as connection:
        cursor = connection.cursor()
        cursor.execute('''CREATE TABLE agn_params
                          (galaxy_id int, htmid_%d int, magNorm real, varParamStr text)''' % htmid_level)

        connection.commit()

        chunk_size = 100000
        full_size = len(redshift_full)
        print('starting iteration')
        t_start = time.time()
        ct_simulated = 0
        for i_start in range(0, full_size, chunk_size):
            i_end = i_start + chunk_size
            duration = (time.time()-t_start)/3600.0
            per = duration/(1+i_start)
            predicted = full_size*per
            print("    %d through %d -- took %.2e hrs; predict %.2e" %
            (i_start,i_end,duration,predicted))

            selection = slice(i_start:i_end)
            galaxy_id = galaxy_id_full[selection]
            ra = ra_full[selection]
            dec = dec_full[selection]
            redshift = redshift_full[selection]
            log_edd_ratio = log_edd_ratio_full[selection]
            bhm = bhm_full[selection]

            ct_simulated += len(galaxy_id)

            abs_mag_i = M_i_from_L_Mass(log_edd_ratio, np.log10(bhm))

            obs_mag_i = get_m_i(abs_mag_i, redshift)

            if args.m_i_cut is not None:
                valid = np.where(obs_mag_i <= args.m_i_cut)
                redshift = redshift[valid]
                bhm = bhm[valid]
                log_edd_ratio = log_edd_ratio[valid]
                abs_mag_i = abs_mag_i[valid]
                galaxy_id = galaxy_id[valid]
                ra = ra[valid]
                dec = dec[valid]
                obs_mag_i = obs_mag_i[valid]

            tau = tau_from_params(redshift,
                                  abs_mag_i,
                                  bhm, rng=rng)

            sf_dict = {}
            for bp in ('u', 'g', 'r', 'i', 'z', 'y'):
                eff_wavelen = 10.0*bp_dict[bp].calcEffWavelen()[0]
                sf_dict[bp] = SF_from_params(redshift, abs_mag_i,
                                             bhm, eff_wavelen,
                                             rng=rng)

            # Cut on structure function value.
            # Specifically, we are looping over all LSST bandpasses
            # and cutting out any AGN that has a structure function
            # value greater than args.max_sf in that bandpass.  At
            # the end of this block of code, the only AGN that should
            # remain are those for whom all structure function values
            # are less than args.max_sf.
            for bp in 'ugrizy':
                valid = np.where(sf_dict[bp]<args.max_sf)

                # update all data structures, including sf_dict
                # to reflect the current cut on structure function
                # values
                redshift = redshift[valid]
                tau = tau[valid]
                log_edd_ratio = log_edd_ratio[valid]
                abs_mag_i = abs_mag_i[valid]
                obs_mag_i = obs_mag_i[valid]
                bhm = bhm[valid]
                galaxy_id = galaxy_id[valid]
                ra = ra[valid]
                dec = dec[valid]
                for other_bp in 'ugrizy':
                    sf_dict[other_bp] = sf_dict[other_bp][valid]

            interpolated_m_i = np.interp(redshift, z_grid, m_i_grid)
            interpolated_mag_norm = np.interp(redshift, z_grid, mag_norm_grid)
            mag_norm = interpolated_mag_norm + (obs_mag_i - interpolated_m_i)

            seed_arr = rng.randint(1,high=10000000, size=len(tau))

            htmid = findHtmid(ra, dec, htmid_level)

            vals = ((int(ii), int(hh), mm, '{"m": "applyAgn", '
                          + '"p": {"seed": %d, "agn_tau": %.3e, "agn_sfu": %.3e, ' % (ss, tt, sfu)
                          + '"agn_sfg": %.3e, "agn_sfr": %.3e, "agn_sfi": %.3e, ' % (sfg, sfr, sfi)
                          + '"agn_sfz": %.3e, "agn_sfy": %.3e}}' % (sfz, sfy))
                     for ii, hh, mm, ss, tt, sfu, sfg, sfr, sfi, sfz, sfy in
                     zip(galaxy_id, htmid, mag_norm, seed_arr, tau,
                         sf_dict['u'], sf_dict['g'], sf_dict['r'], sf_dict['i'],
                         sf_dict['z'], sf_dict['y']))

            cursor.executemany('INSERT INTO agn_params VALUES(?, ?, ?, ?)', vals)
            connection.commit()

        assert ct_simulated == full_size

        print('creating index')
        cursor.execute('CREATE INDEX htmid ON agn_params (htmid_6)')
        connection.commit()

    print('all done')
