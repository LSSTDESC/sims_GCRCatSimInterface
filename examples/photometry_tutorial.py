import os
import numpy as np
import h5py

from GCR import GCRQuery
import GCRCatalogs

import lsst.sims.photUtils as sims_photUtils
from lsst.sims.catUtils.dust import EBVbase

sed_fit_dir = '/global/projecta/projectdirs/lsst/groups/SSim'
assert os.path.isdir(sed_fit_dir)
sed_fit_dir = os.path.join(sed_fit_dir,'DC2/cosmoDC2_v1.1.4/sedLookup')
assert os.path.isdir(sed_fit_dir)

hpid = 10069  # an example healpix pixel that has been fit

sed_fit_name = os.path.join(sed_fit_dir, 'sed_fit_%d.h5' % hpid)
assert os.path.isfile(sed_fit_name)

# load cosmoDC2
cat = GCRCatalogs.load_catalog('cosmoDC2_v1.1.4_image')

# get galaxy_id and redshift for crossmatching with SED fit files;
# we will also get the magnitudes that should be reproduced
# by our synthetic photometry
# (hp_query makes sure we only load the healpixel we are interested in)
hp_query = GCRQuery('healpix_pixel==%d' % hpid)
cosmoDC2_data = cat.get_quantities(['galaxy_id', 'redshift', 'ra', 'dec',
                                    'mag_true_u_lsst', 'mag_true_g_lsst',
                                    'mag_true_r_lsst', 'mag_true_i_lsst',
                                    'mag_true_z_lsst', 'mag_true_y_lsst',
                                    'mag_u_lsst', 'mag_g_lsst', 'mag_r_lsst',
                                    'mag_i_lsst', 'mag_z_lsst', 'mag_y_lsst',
                                    'shear_1', 'shear_2', 'convergence'],
                                    native_filters=[hp_query])

# make sure cosmoDC2_data is sorted by galaxy_id
sorted_dex = np.argsort(cosmoDC2_data['galaxy_id'])
for colname in cosmoDC2_data.keys():
    cosmoDC2_data[colname] = cosmoDC2_data[colname][sorted_dex]

# read in LSST bandpasses
lsst_bp_dict = sims_photUtils.BandpassDict.loadTotalBandpassesFromFiles()

# the parent directory for the SED library
sed_lib_dir = os.environ['SIMS_SED_LIBRARY_DIR']

# open the file with the SED fits
with h5py.File(sed_fit_name, 'r') as sed_fit_file:

    sed_names = sed_fit_file['sed_names'][()]
    sed_names = [s.decode() for s in sed_names]  # because stored as bytes

    # we will just consider the first 100 galaxies for this example
    n_test_gals = 100
    subset = slice(0, n_test_gals)

    galaxy_id = sed_fit_file['galaxy_id'][()][subset]

    # get the cross-match between the sed fit and cosmoDC2
    crossmatch_dex = np.searchsorted(cosmoDC2_data['galaxy_id'], galaxy_id)
    np.testing.assert_array_equal(galaxy_id,
                                  cosmoDC2_data['galaxy_id'][crossmatch_dex])

    ra = sed_fit_file['ra'][()][subset]
    dec = sed_fit_file['dec'][()][subset]
    np.testing.assert_array_equal(ra, cosmoDC2_data['ra'][crossmatch_dex])
    np.testing.assert_array_equal(dec, cosmoDC2_data['dec'][crossmatch_dex])

    # Calculate E(B-V) for dust extinction in Milky Way along relevant
    # lines of sight
    equatorial_coords = np.array([np.radians(ra), np.radians(dec)])
    ebv_model = EBVbase()
    ebv_vals = ebv_model.calculateEbv(equatorialCoordinates=equatorial_coords)

    ccm_w = None  # so that we only have to initialize internal dust once
    for i_bp, bp in enumerate('ugrizy'):
        fluxes_noMW = {}
        fluxes = {}
        for component in ['disk', 'bulge']:
            fluxes_noMW[component] = np.zeros(n_test_gals, dtype=float)
            fluxes[component] = np.zeros(n_test_gals, dtype=float)

        for component in ['disk', 'bulge']:
            sed_arr = sed_fit_file['%s_sed' % component][()][subset]
            av_arr = sed_fit_file['%s_av' % component][()][subset]
            rv_arr = sed_fit_file['%s_rv' % component][()][subset]
            mn_arr = sed_fit_file['%s_magnorm' % component][()][i_bp,:][subset]
            z_arr = cosmoDC2_data['redshift'][crossmatch_dex]
            for i_gal, (s_dex, mn, av,
                        rv, zz, ebv) in enumerate(zip(sed_arr,
                                                      mn_arr,
                                                      av_arr,
                                                      rv_arr,
                                                      z_arr,
                                                      ebv_vals)):

                # read in the SED file from the library
                sed_file_name = os.path.join(sed_lib_dir, sed_names[s_dex])
                sed = sims_photUtils.Sed()
                sed.readSED_flambda(sed_file_name)

                # find and apply normalizing flux
                fnorm = sims_photUtils.getImsimFluxNorm(sed, mn)
                sed.multiplyFluxNorm(fnorm)

                # add internal dust
                if ccm_w is None or not np.array_equal(sed.wavelen, ccm_w):
                    ccm_w = np.copy(sed.wavelen)
                    a_x, b_x = sed.setupCCM_ab()
                sed.addDust(a_x, b_x, A_v=av, R_v=rv)

                # apply redshift
                sed.redshiftSED(zz, dimming=True)

                # flux, in Janskys, without Milky Way dust extinction
                f_noMW = sed.calcFlux(lsst_bp_dict[bp])

                # apply Milky Way dust
                # (cannot reuse a_x, b_x because wavelength grid changed
                # when we called redshiftSED)
                a_x_mw, b_x_mw = sed.setupCCM_ab()
                sed.addDust(a_x_mw, b_x_mw, R_v=3.1, ebv=ebv)

                f_MW = sed.calcFlux(lsst_bp_dict[bp])

                fluxes_noMW[component][i_gal] = f_noMW
                fluxes[component][i_gal] = f_MW

        total_fluxes = fluxes_noMW['disk'] + fluxes_noMW['bulge']

        dummy_sed = sims_photUtils.Sed()
        cosmoDC2_mags = cosmoDC2_data['mag_true_%s_lsst' % bp][crossmatch_dex]
        cosmoDC2_fluxes = dummy_sed.fluxFromMag(cosmoDC2_mags)

        d_flux = np.abs(1.0-total_fluxes/cosmoDC2_fluxes)

        # add magnification due to weak lensing
        kappa = cosmoDC2_data['convergence'][crossmatch_dex]
        gamma_sq = (cosmoDC2_data['shear_1'][crossmatch_dex]**2
                    + cosmoDC2_data['shear_2'][crossmatch_dex]**2)
        magnification = 1.0/((1.0-kappa)**2-gamma_sq)
        magnified_fluxes = magnification*total_fluxes
        cosmoDC2_mags_mag = cosmoDC2_data['mag_%s_lsst' % bp][crossmatch_dex]
        cosmoDC2_fluxes_mag = dummy_sed.fluxFromMag(cosmoDC2_mags_mag)

        d_flux_mag = np.abs(1.0-magnified_fluxes/cosmoDC2_fluxes_mag)

        print('max fractional delta_flux for %s: %e, %e with lensing' %
              (bp, d_flux.max(), d_flux_mag.max()))
