import os
import numpy as np
import h5py
import lsst.sims.photUtils as photUtils
import GCRCatalogs
from GCR import GCRQuery

import time
import argparse

import multiprocessing

def validate_chunk(data_in,
                   in_dir, healpix,
                   read_lock, write_lock, output_dict):

    galaxy_id = data_in['galaxy_id']
    redshift = data_in['redshift']
    mag_in = {}
    for bp in 'ugrizy':
        mag_in[bp] = data_in['mag_true_%s_lsst' % bp]

    eps = 1.0e-20
    disk_to_bulge = {}
    for bp in 'ugrizy':
        disk_name = 'LSST_filters/diskLuminositiesStellar:'
        disk_name += 'LSST_%s:observed:dustAtlas' % bp

        bulge_name = 'LSST_filters/spheroidLuminositiesStellar:'
        bulge_name += 'LSST_%s:observed:dustAtlas' % bp

        denom = np.where(data_in[bulge_name]>0.0, data_in[bulge_name], eps)
        disk_to_bulge[bp] = data_in[disk_name]/denom
        assert np.isfinite(disk_to_bulge[bp]).sum()==len(disk_to_bulge[bp])

    sed_dir = os.environ['SIMS_SED_LIBRARY_DIR']

    tot_bp_dict = photUtils.BandpassDict.loadTotalBandpassesFromFiles()

    fit_file = os.path.join(in_dir, 'sed_fit_%d.h5' % healpix)
    assert os.path.isfile(fit_file)
    with read_lock:
        with h5py.File(fit_file, 'r') as in_file:
            sed_names = in_file['sed_names'][()]
            galaxy_id_fit = in_file['galaxy_id'][()]
            valid = np.in1d(galaxy_id_fit, galaxy_id)
            galaxy_id_fit = galaxy_id_fit[valid]
            np.testing.assert_array_equal(galaxy_id_fit, galaxy_id)

            disk_magnorm = in_file['disk_magnorm'][()][:,valid]
            disk_av = in_file['disk_av'][()][valid]
            disk_rv = in_file['disk_rv'][()][valid]
            disk_sed = in_file['disk_sed'][()][valid]

            bulge_magnorm = in_file['bulge_magnorm'][()][:,valid]
            bulge_av = in_file['bulge_av'][()][valid]
            bulge_rv = in_file['bulge_rv'][()][valid]
            bulge_sed = in_file['bulge_sed'][()][valid]

    sed_names = [name.decode() for name in sed_names]

    local_worst = {}
    local_worst_ratio = {}
    ccm_w = None
    dummy_sed = photUtils.Sed()
    for i_obj in range(len(disk_av)):
        for i_bp, bp in enumerate('ugrizy'):
            d_sed = photUtils.Sed()
            fname = os.path.join(sed_dir, sed_names[disk_sed[i_obj]])
            d_sed.readSED_flambda(fname)
            fnorm = photUtils.getImsimFluxNorm(d_sed,
                                               disk_magnorm[i_bp][i_obj])
            d_sed.multiplyFluxNorm(fnorm)
            if ccm_w is None or not np.array_equal(d_sed.wavelen, ccm_w):
                ccm_w = np.copy(d_sed.wavelen)
                ax, bx = d_sed.setupCCM_ab()
            d_sed.addDust(ax, bx, R_v=disk_rv[i_obj], A_v=disk_av[i_obj])
            d_sed.redshiftSED(redshift[i_obj], dimming=True)
            d_flux = d_sed.calcFlux(tot_bp_dict[bp])

            b_sed = photUtils.Sed()
            fname = os.path.join(sed_dir, sed_names[bulge_sed[i_obj]])
            b_sed.readSED_flambda(fname)
            fnorm = photUtils.getImsimFluxNorm(b_sed,
                                               bulge_magnorm[i_bp][i_obj])
            b_sed.multiplyFluxNorm(fnorm)
            if ccm_w is None or not np.array_equal(b_sed.wavelen, ccm_w):
                ccm_w = np.copy(b_sed.wavelen)
                ax, bx = b_sed.setupCCM_ab()
            b_sed.addDust(ax, bx, R_v=bulge_rv[i_obj], A_v=bulge_av[i_obj])
            b_sed.redshiftSED(redshift[i_obj], dimming=True)
            b_flux = b_sed.calcFlux(tot_bp_dict[bp])

            flux_ratio = d_flux/max(b_flux, eps)

            if flux_ratio<1.0e3 or disk_to_bulge[bp][i_obj]<1.0e3:
                delta_ratio = np.abs(flux_ratio-disk_to_bulge[bp][i_obj])
                if disk_to_bulge[bp][i_obj]>eps:
                    delta_ratio /= disk_to_bulge[bp][i_obj]
                if bp not in local_worst_ratio or delta_ratio>local_worst_ratio[bp]:
                    local_worst_ratio[bp] = delta_ratio

            tot_flux = b_flux + d_flux
            true_flux = dummy_sed.fluxFromMag(mag_in[bp][i_obj])
            delta_flux = np.abs(1.0-tot_flux/true_flux)
            if bp not in local_worst or delta_flux>local_worst[bp]:
                local_worst[bp] = delta_flux

    with write_lock:
        for bp in 'ugrizy':
            if bp not in output_dict or local_worst[bp]>output_dict[bp]:
                output_dict[bp] = local_worst[bp]
            ratio_name = '%s_ratio' % bp
            if (ratio_name not in output_dict
                or local_worst_ratio[bp]>output_dict[ratio_name]):
                output_dict[ratio_name] = local_worst_ratio[bp]

if __name__ == "__main__":

    t_start = time.time()

    parser = argparse.ArgumentParser()
    parser.add_argument('--healpix', type=int, default=None)
    parser.add_argument('--in_dir', type=str, default=None)
    parser.add_argument('--seed', type=int, default=None)
    parser.add_argument('--nsamples', type=int, default=1000000)
    parser.add_argument('--n_threads', type=int, default=10)
    parser.add_argument('--d_gal', type=int, default=10000)

    args = parser.parse_args()
    if args.healpix is None:
        raise RuntimeError("must specify healpix")
    if args.in_dir is None or not os.path.isdir(args.in_dir):
        raise RuntimeError("invalid in_dir")
    if args.seed is None:
        raise RuntimeError("must specify seed")

    cache_dir = os.path.join(os.environ['HOME'],'workspace','sed_cache')
    photUtils.cache_LSST_seds(wavelen_min=0.0, wavelen_max=1500.0,
                              cache_dir=cache_dir)

    mgr = multiprocessing.Manager()
    read_lock = mgr.Lock()
    write_lock = mgr.Lock()
    output_dict = mgr.dict()

    col_names = ['galaxy_id', 'redshift']
    for bp in 'ugrizy':
        col_names.append('mag_true_%s_lsst' % bp)
        col_names.append('LSST_filters/diskLuminositiesStellar:'
                         'LSST_%s:observed:dustAtlas' % bp)
        col_names.append('LSST_filters/spheroidLuminositiesStellar:'
                         'LSST_%s:observed:dustAtlas' % bp)

    print('loading catalog')
    cat = GCRCatalogs.load_catalog('cosmoDC2_v1.1.4_image')
    h_query = GCRQuery('healpix_pixel==%d' % args.healpix)
    data = cat.get_quantities(col_names,
                              native_filters=[h_query])
    print('got catalog')

    if args.nsamples>0:
        rng = np.random.RandomState(args.seed)
        sub_dexes = rng.choice(np.arange(len(data['galaxy_id']), dtype=int),
                              size=args.nsamples,
                              replace=False)

        for k in data.keys():
            data[k] = data[k][sub_dexes]

        sorted_dex = np.argsort(data['galaxy_id'])
        for k in data.keys():
            data[k] = data[k][sorted_dex]

    print('n galaxies %e' % len(data['galaxy_id']))

    p_list = []
    ct = 0
    for i_start in range(0,len(data['galaxy_id']), args.d_gal):
        sub_sample = slice(i_start, i_start+args.d_gal)

        sub_data = {}
        for k in data:
            sub_data[k] = data[k][sub_sample]

        p = multiprocessing.Process(target=validate_chunk,
                                    args=(sub_data,
                                          args.in_dir, args.healpix,
                                          read_lock, write_lock, output_dict))

        p.start()
        p_list.append(p)
        while len(p_list)>=args.n_threads:
            exit_list = []
            for p in p_list:
                exit_list.append(p.exitcode)
            was_popped = False
            for ii in range(len(p_list)-1,-1,-1):
                if exit_list[ii] is not None:
                    p_list.pop(ii)
                    ct += args.d_gal
                    was_popped = True
            #if was_popped:
            #    duration = (time.time()-t_start)/3600.0
            #    per = duration/ct
            #    pred = len(data['galaxy_id'])*per
            #    print('checked %d in %.2e; pred %.2e' % (ct, duration, pred))

    for p in p_list:
        p.join()

    for bp in output_dict.keys():
        print('hp %d %s %e' % (args.healpix, bp, output_dict[bp]))
    print('hp %d that took %.2e hrs' %
          (args.healpix, (time.time()-t_start)/3600.0))
