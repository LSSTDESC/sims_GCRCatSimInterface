import os
import numpy as np
import h5py
import sqlite3

from GCR import GCRQuery
import GCRCatalogs

import lsst.sims.photUtils as sims_photUtils
from lsst.sims.catUtils.dust import EBVbase

sed_fit_dir = '/global/projecta/projectdirs/lsst/groups/SSim'
assert os.path.isdir(sed_fit_dir)
sed_fit_dir = os.path.join(sed_fit_dir,'DC2/cosmoDC2_v1.1.4/sedLookup')
assert os.path.isdir(sed_fit_dir)

### This should be an argument
# hpid = 10069  # an example healpix pixel that has been fit
hpid = 9430

sed_fit_name = os.path.join(sed_fit_dir, 'sed_fit_%d.h5' % hpid)
assert os.path.isfile(sed_fit_name)

# Hard-code output file path for now
#the_output = '/global/cscratch1/sd/jrbogart/truth/summary_table_10069.sqlite3'
the_output = '/global/cscratch1/sd/jrbogart/truth/summary_table_hp{}_n{}.sqlite3'

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

galaxy_ids = []
ra = []
dec = []
redshift = []
n_test_gals = 0
# open the file with the SED fits
with h5py.File(sed_fit_name, 'r') as sed_fit_file:

    sed_names = sed_fit_file['sed_names'][()]
    sed_names = [s.decode() for s in sed_names]  # because stored as bytes

    # we will just consider the first 100000 galaxies for this example
    n_test_gals = 100000
    subset = slice(0, n_test_gals)
    galaxy_ids = sed_fit_file['galaxy_id'][()][subset]
    #galaxy_ids = sed_fit_file['galaxy_id'][()]

    #n_test_gals = len(galaxy_ids)

    print("Number of fit galaxies: ", n_test_gals)

    the_output = the_output.format(hpid, n_test_gals)
    print("writing to output /n{}".format(the_output))

    # get the cross-match between the sed fit and cosmoDC2
    crossmatch_dex = np.searchsorted(cosmoDC2_data['galaxy_id'], galaxy_ids)
    np.testing.assert_array_equal(galaxy_ids,
                                  cosmoDC2_data['galaxy_id'][crossmatch_dex])

    print("Did crossmatch")
    ra = sed_fit_file['ra'][()][subset]
    dec = sed_fit_file['dec'][()][subset]
    np.testing.assert_array_equal(ra, cosmoDC2_data['ra'][crossmatch_dex])
    np.testing.assert_array_equal(dec, cosmoDC2_data['dec'][crossmatch_dex])

    flux_by_band_MW = {}
    flux_by_band_noMW = {}

    # Some columns can be written immediately:
    #     galaxy_id as id
    #     ra
    #     dec
    #     redshift
    #     host_galaxy = -1
    #     is_variable = 0
    #     is_pointsource = 0
    
    # Calculate E(B-V) for dust extinction in Milky Way along relevant
    # lines of sight
    equatorial_coords = np.array([np.radians(ra), np.radians(dec)])
    ebv_model = EBVbase()
    ebv_vals = ebv_model.calculateEbv(equatorialCoordinates=equatorial_coords)

    print("Calculated internal extinction")
    ccm_w = None  # so that we only have to initialize internal dust once
    for i_bp, bp in enumerate('ugrizy'):
        print("Processing band ", bp)
        fluxes_noMW = {}
        fluxes = {}
        for component in ['disk', 'bulge']:
            fluxes_noMW[component] = np.zeros(n_test_gals, dtype=float)
            fluxes[component] = np.zeros(n_test_gals, dtype=float)

        for component in ['disk', 'bulge']:
            print("   Processing component ", component)
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
            if (component == 'disk') and (bp == 'r'):
                redshift = z_arr
                print("picked a redshift array: component disk, band r")
        total_fluxes = fluxes_noMW['disk'] + fluxes_noMW['bulge']

        total_fluxes_MW = fluxes['disk'] + fluxes['bulge']

        print("fluxes computed, with and without extinction")
        #dummy_sed = sims_photUtils.Sed()
        #cosmoDC2_mags = cosmoDC2_data['mag_true_%s_lsst' % bp][crossmatch_dex]
        #cosmoDC2_fluxes = dummy_sed.fluxFromMag(cosmoDC2_mags)

        #d_flux = np.abs(1.0-total_fluxes/cosmoDC2_fluxes)

        # add magnification due to weak lensing
        kappa = cosmoDC2_data['convergence'][crossmatch_dex]
        gamma_sq = (cosmoDC2_data['shear_1'][crossmatch_dex]**2
                    + cosmoDC2_data['shear_2'][crossmatch_dex]**2)
        magnification = 1.0/((1.0-kappa)**2-gamma_sq)
        magnified_fluxes = magnification*total_fluxes
        magnified_fluxes_MW = magnification*total_fluxes_MW
        flux_by_band_noMW[bp] = magnified_fluxes
        flux_by_band_MW[bp] = magnified_fluxes_MW

        print("added magnification")
#  Open connection to sqlite db and write     
#  Input arguments should include directory and healpix id.
#  Default output file path will be derived from these.

with sqlite3.connect(the_output) as conn:
    cursor = conn.cursor()
    cmd = '''CREATE TABLE truth_summary
          (id BIGINT, host_galaxy BIGINT, ra DOUBLE, dec DOUBLE,
          redshift FLOAT, is_variable INT, is_pointsource INT,
          flux_u FLOAT, flux_g FLOAT, flux_r FLOAT, 
          flux_i FLOAT, flux_z FLOAT, flux_y FLOAT,
          flux_u_noMW FLOAT, flux_g_noMW FLOAT, flux_r_noMW FLOAT, 
          flux_i_noMW FLOAT, flux_z_noMW FLOAT, flux_y_noMW FLOAT)'''
    cursor.execute(cmd)
    conn.commit()
    print("Created table truth_summary")
    print("Count of galaxy_id: ", len(galaxy_ids))

    values = ((int(galaxy_ids[i_obj]),int(-1),
               ra[i_obj],dec[i_obj],
               redshift[i_obj], 0, 0,
               flux_by_band_MW['u'][i_obj], flux_by_band_MW['g'][i_obj],
               flux_by_band_MW['r'][i_obj], flux_by_band_MW['i'][i_obj],
               flux_by_band_MW['z'][i_obj], flux_by_band_MW['y'][i_obj],
               flux_by_band_noMW['u'][i_obj], flux_by_band_noMW['g'][i_obj],
               flux_by_band_noMW['r'][i_obj], flux_by_band_noMW['i'][i_obj],
               flux_by_band_noMW['z'][i_obj], flux_by_band_noMW['y'][i_obj])
              for i_obj in range(len(galaxy_ids)))
    #print("# rows to insert: ", len(values))
    cursor.executemany('''INSERT INTO truth_summary
                       VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)''',
                       values)
    conn.commit()

    print("inserted entries")
    # Add index
    cmd = '''CREATE INDEX ra_dec ON truth_summary (ra,dec)'''
    cursor.execute(cmd)
    cmd = '''CREATE INDEX gal_id ON truth_summary (id)'''
    cursor.execute(cmd)
    conn.commit
    print("created indexes")

    print('done')
