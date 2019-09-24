import os
import numpy as np
import h5py
import time
import multiprocessing

from lsst.sims.catUtils.dust import EBVbase
import lsst.sims.photUtils as sims_photUtils

import argparse

def process_component(sed_names,
                      ebv_arr,
                      redshift,
                      sed_dexes,
                      av,
                      rv,
                      magnorm,
                      my_lock):

    bp_dict = sims_photUtils.BandpassDict.loadTotalBandpassesFromFiles()
    n_obj = len(sed_dexes)
    fluxes = {}
    fluxes_noMW = {}
    for bp in 'ugrizy':
        fluxes[bp] = np.zeros(n_obj, dtype=float)
        fluxes_noMW[bp] = np.zeros(n_obj, dtype=float)

    print('starting flux calculation')
    t_start = time.time()
    ct = 0
    for sed_id in np.unique(sed_dexes):
        valid = np.where(sed_dexes==sed_id)
        sed_full_name = os.path.join(os.environ['SIMS_SED_LIBRARY_DIR'],
                                     sed_names[sed_id])

        spec_rest = sims_photUtils.Sed()
        spec_rest.readSED_flambda(sed_full_name, cache_sed=False)
        a_x_rest, b_x_rest = spec_rest.setupCCM_ab()
        fnorm_at_21 = sims_photUtils.getImsimFluxNorm(spec_rest, 21.0)
        for g_dex in valid[0]:
            spec_rest_dust = sims_photUtils.Sed(wavelen=spec_rest.wavelen,
                                                flambda=spec_rest.flambda)
            spec_rest_dust.addDust(a_x_rest, b_x_rest,
                                   A_v=av[g_dex],
                                   R_v=rv[g_dex])

            spec_rest_dust.redshiftSED(redshift[g_dex], dimming=True)
            flux_list_noMW = bp_dict.fluxListForSed(spec_rest_dust)
            a_x, b_x = spec_rest_dust.setupCCM_ab()
            spec_rest_dust.addDust(a_x, b_x, R_v=3.1, ebv=ebv_arr[g_dex])
            flux_list = bp_dict.fluxListForSed(spec_rest_dust)
            for i_bp, bp in enumerate('ugrizy'):
                factor = fnorm_at_21*np.power(10.0,
                                     -0.4*(magnorm[i_bp][g_dex]-21.0))

                fluxes[bp][g_dex] = factor*flux_list[i_bp]
                fluxes_noMW[bp][g_dex] = factor*flux_list_noMW[i_bp]
            ct += 1
            if ct%1000 == 0:
                with my_lock as context:
                    duration = (time.time()-t_start)/3600.0
                    per = duration/ct
                    pred = per*n_obj
                    print('ran ',ct,duration,per,pred)

    return fluxes, fluxes_noMW


def calculate_fluxes(in_name, out_name, healpix_id, my_lock):

    base_dir = os.path.join('/astro/store/pogo4/danielsf/desc_dc2_truth')

    t_start = time.time()
    with my_lock as context:
        bp_dict = sims_photUtils.BandpassDict.loadTotalBandpassesFromFiles()

        redshift_file_name = os.path.join(base_dir, 'redshift',
                                          'redshift_%d.h5' % healpix_id)
        with h5py.File(redshift_file_name, 'r') as in_file:
            redshift = in_file['redshift'][()]
            galaxy_id_z = in_file['galaxy_id'][()]


        with h5py.File(in_name, 'r') as in_file:
            sed_names = in_file['sed_names'][()].astype(str)

            galaxy_id_disk = in_file['galaxy_id'][()]
            sorted_dex = np.argsort(galaxy_id_disk)
            galaxy_id_disk = galaxy_id_disk[sorted_dex]

            ra = in_file['ra'][()][sorted_dex]
            dec = in_file['dec'][()][sorted_dex]
            disk_sed = in_file['disk_sed'][()][sorted_dex]
            disk_magnorm = in_file['disk_magnorm'][()]
            for ii in range(6):
                disk_magnorm[ii] = disk_magnorm[ii][sorted_dex]
            disk_av = in_file['disk_av'][()][sorted_dex]
            disk_rv = in_file['disk_rv'][()][sorted_dex]
            disk_fluxes_in = in_file['disk_fluxes'][()]
            for ii in range(6):
                disk_fluxes_in[ii] = disk_fluxes_in[ii][sorted_dex]

        print('ra ',ra.min(),ra.max())
        print('dec ',dec.min(), dec.max())
        ebv_model = EBVbase()
        ebv_arr = ebv_model.calculateEbv(interp=True,
                           equatorialCoordinates=np.array([np.radians(ra),
                                                           np.radians(dec)]))
        del ebv_model
        del ra
        del dec

    print('loaded disk data in %e' % (time.time()-t_start))

    sorted_dex_z = np.argsort(galaxy_id_z)
    galaxy_id_z = galaxy_id_z[sorted_dex_z]
    redshift = redshift[sorted_dex_z]

    if len(galaxy_id_z) != len(galaxy_id_disk):
        galaxy_id_z = galaxy_id_z[:len(galaxy_id_disk)]
        redshift = redshift[:len(galaxy_id_disk)]

    np.testing.assert_array_equal(galaxy_id_z, galaxy_id_disk)

    (disk_fluxes,
     disk_fluxes_noMW) = process_component(sed_names,
                                           ebv_arr,
                                           redshift,
                                           disk_sed,
                                           disk_av,
                                           disk_rv,
                                           disk_magnorm,
                                           my_lock)

    for i_bp, bp in enumerate('ugrizy'):
        d_flux_ratio = disk_fluxes_noMW[bp]/disk_fluxes_in[i_bp]
        print(bp,' disk flux ratio ',d_flux_ratio.max(),d_flux_ratio.min())

    del disk_sed
    del disk_av
    del disk_rv
    del disk_magnorm

    with my_lock as context:
        with h5py.File(in_name, 'r') as in_file:
            galaxy_id_bulge = in_file['galaxy_id'][()]
            sorted_dex = np.argsort(galaxy_id_bulge)
            galaxy_id_bulge = galaxy_id_bulge[sorted_dex]

            bulge_sed = in_file['bulge_sed'][()][sorted_dex]
            bulge_av = in_file['bulge_av'][()][sorted_dex]
            bulge_rv = in_file['bulge_rv'][()][sorted_dex]
            bulge_fluxes_in = in_file['bulge_fluxes'][()]
            for ii in range(6):
                bulge_fluxes_in[ii] = bulge_fluxes_in[ii][sorted_dex]
            bulge_magnorm = in_file['bulge_magnorm'][()]
            for ii in range(6):
                bulge_magnorm[ii] = bulge_magnorm[ii][sorted_dex]

    np.testing.assert_array_equal(galaxy_id_z, galaxy_id_bulge)

    (bulge_fluxes,
     bulge_fluxes_noMW) = process_component(sed_names,
                                            ebv_arr,
                                            redshift,
                                            bulge_sed,
                                            bulge_av,
                                            bulge_rv,
                                            bulge_magnorm,
                                            my_lock)

    for i_bp, bp in enumerate('ugrizy'):
        bulge_flux_ratio = bulge_fluxes_in[i_bp]/bulge_fluxes_noMW[bp]
        print(bp,' bulge flux ratio ',np.nanmax(bulge_flux_ratio),
              np.nanmin(bulge_flux_ratio))


    with my_lock as context:
        with h5py.File(out_name, 'w') as out_file:
            out_file.create_dataset('galaxy_id', data=galaxy_id_z)
            out_file.create_dataset('redshift', data=redshift)
            for bp in 'ugrizy':
                out_file.create_dataset('flux_%s' % bp,
                                 data=1.0e9*(bulge_fluxes[bp]+disk_fluxes[bp]))
                out_file.create_dataset('flux_%s_noMW' % bp,
                       data=1.0e9*(bulge_fluxes_noMW[bp]+disk_fluxes_noMW[bp]))


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--out_dir', type=str, default=None)
    parser.add_argument('--healpix', type=int, nargs='+', default=None)
    args = parser.parse_args()

    data_dir = '/astro/store/pogo4/danielsf/desc_dc2_truth/sedLookup'
    assert os.path.isdir(data_dir)

    if isinstance(args.healpix, int):
        args.healpix = [args.healpix]
    elif args.healpix is None:
        file_list = os.listdir(data_dir)
        args.healpix = []
        for name in file_list:
            p = name.split('_')
            args.healpix.append(int(p[-1].replace('.h5','')))

    print('n_healpix %d' % len(args.healpix))
    for healpix_id in args.healpix:
        file_name = os.path.join(data_dir, 'sed_fit_%d.h5' % healpix_id)
        assert os.path.isfile(file_name)

    my_lock = multiprocessing.Lock()
    in_name = os.path.join(data_dir, 'sed_fit_%d.h5' % args.healpix[0])
    out_name = os.path.join(args.out_dir,
                            'galaxy_truth_%d.h5' % args.healpix[0])
    calculate_fluxes(in_name, out_name, args.healpix[0], my_lock)
