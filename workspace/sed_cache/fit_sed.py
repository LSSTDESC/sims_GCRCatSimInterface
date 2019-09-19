import os
import h5py
import numpy as np
import healpy
import multiprocessing
import time

import GCRCatalogs
from GCR import GCRQuery

from SedFitter import sed_filter_names_from_catalog
from SedFitter import sed_from_galacticus_mags
from lsst.sims.photUtils import BandpassDict, Sed, Bandpass
from lsst.sims.photUtils import cache_LSST_seds, getImsimFluxNorm
from lsst.utils import getPackageDir
#from lsst.sims.utils import htmModule as htm

import argparse


def _parallel_fitting(mag_array, redshift, redshift_true,
                      H0, Om0, wav_min, wav_width,
                      lsst_mag_array, out_dict, tag):
    pid = os.getpid()
    print('starting ',pid,len(redshift))
    (sed_names,
     mag_norms,
     av_arr,
     rv_arr) = sed_from_galacticus_mags(mag_array,
                                        redshift,
                                        redshift_true,
                                        H0, Om0,
                                        wav_min, wav_width,
                                        lsst_mag_array)

    tot_bp_dict = BandpassDict.loadTotalBandpassesFromFiles()

    sed_dir = getPackageDir('sims_sed_library')
    lsst_fit_fluxes = np.zeros((6,len(sed_names)), dtype=float)
    t_start = time.time()

    ccm_w = None
    restframe_seds = {}
    imsim_bp = Bandpass()
    imsim_bp.imsimBandpass()
    n04_ln10 = -0.4*np.log(10)

    for ii in range(len(sed_names)):
        av_val = av_arr[ii]
        rv_val = rv_arr[ii]

        sed_tag = '%s_%.3f_%.3f' % (sed_names[ii], av_val, rv_val)
        if sed_tag not in restframe_seds:
            rest_sed = Sed()
            rest_sed.readSED_flambda(os.path.join(sed_dir, sed_names[ii]))
            mag = rest_sed.calcMag(imsim_bp)
            if ccm_w is None or not np.array_equal(rest_sed.wavelen, ccm_w):
                ccm_w = np.copy(rest_sed.wavelen)
                ax, bx = rest_sed.setupCCM_ab()
            rest_sed.addDust(ax, bx, A_v=av_val, R_v=rv_val)
            restframe_seds[sed_tag] = (rest_sed, mag)

        for i_bp, bp in enumerate('ugrizy'):
            m_norm = mag_norms[i_bp][ii]
            if m_norm>0.0 and not np.isfinite(m_norm):
                continue

            spec = Sed(wavelen=restframe_seds[sed_tag][0].wavelen,
                       flambda=restframe_seds[sed_tag][0].flambda)
            fnorm = np.exp(n04_ln10*(m_norm - restframe_seds[sed_tag][1]))
            try:
                assert np.isfinite(fnorm)
                assert fnorm>0.0
            except AssertionError:
                print('\n\nmagnorm %e\n\n' % (m_norm))
                raise
            spec.multiplyFluxNorm(fnorm)
            spec.redshiftSED(redshift[ii], dimming=True)
            ff = spec.calcFlux(tot_bp_dict[bp])
            lsst_fit_fluxes[i_bp][ii] = ff

    out_dict[tag] = (sed_names, mag_norms, av_arr, rv_arr, lsst_fit_fluxes)
    print('done with ',pid)


def do_fitting(cat, component, healpix, lim, n_threads):
    """
    Fit a set of components to SEDs, Av, Rv, magNorm using sed_from_galacticus_mags

    Parameters
    ----------
    cat -- the result of GCRCatalogs.load_catalog('catalog_name')

    component -- a string; either 'disk' or 'bulge'

    healpix -- an int indicating which healpixel to fit

    lim -- an int indicating how many objects to actually fit

    Returns
    -------
    numpy arrays of:
    redshift
    galaxy_id
    sed_name
    magNorm
    Av
    Rv
    """

    filter_data = sed_filter_names_from_catalog(cat)
    filter_names = filter_data[component]['filter_name']
    lsst_filter_names = filter_data[component]['lsst_fluxes']
    wav_min = filter_data[component]['wav_min']
    wav_width = filter_data[component]['wav_width']

    H0 = cat.cosmology.H0.value
    Om0 = cat.cosmology.Om0

    healpix_query = GCRQuery('healpix_pixel==%d' % healpix)

    qties = cat.get_quantities(list(filter_names) + list(lsst_filter_names) +
                              ['redshift', 'redshift_true', 'galaxy_id'],
                               native_filters=[healpix_query])

    print("testing on %d of %d" % (lim, len(qties['galaxy_id'])))
    with np.errstate(divide='ignore', invalid='ignore'):
        mag_array = np.array([-2.5*np.log10(qties[ff][:lim])
                              for ff in filter_names])

        lsst_mag_array = np.array([-2.5*np.log10(qties[ff][:lim])
                                   for ff in lsst_filter_names])


    print('getting sed_from_galacticus_mags')
    redshift = qties['redshift'][:lim]
    redshift_true = qties['redshift_true'][:lim]
    (sed_names,
     mag_norms,
     av_arr,
     rv_arr) = sed_from_galacticus_mags(mag_array[:,:2],
                                        redshift[:2],
                                        redshift_true[:2],
                                        H0, Om0,
                                        wav_min, wav_width,
                                        lsst_mag_array[:,:2])

    print('got sed_from_galacticus_mags')

    mgr = multiprocessing.Manager()
    out_dict = mgr.dict()
    d_gal = len(redshift)//n_threads
    p_list = []
    print('should start calling parallel_fitting')
    for i_start in range(0, len(redshift), d_gal):
        s = slice(i_start, i_start+d_gal)
        p = multiprocessing.Process(target=_parallel_fitting,
                                    args=(mag_array[:,s], redshift[s],
                                          redshift_true[s],
                                          H0, Om0, wav_min, wav_width,
                                          lsst_mag_array[:,s],
                                          out_dict, i_start))
        p.start()
        p_list.append(p)

    for p in p_list:
        p.join()

    print('done running parallel fitting')
    sed_names = np.empty(len(redshift), dtype=(str,200))
    mag_norms = np.zeros((6,len(redshift)), dtype=float)
    av_arr = np.zeros(len(redshift), dtype=float)
    rv_arr = np.zeros(len(redshift), dtype=float)
    lsst_fluxes = np.zeros((6,len(redshift)), dtype=float)

    t_start_slicing = time.time()
    for i_start in out_dict.keys():
        s = slice(i_start, i_start+d_gal)
        sed_names[s] = out_dict[i_start][0]
        mag_norms[:,s] = out_dict[i_start][1]
        av_arr[s] = out_dict[i_start][2]
        rv_arr[s] = out_dict[i_start][3]
        lsst_fluxes[:,s] = out_dict[i_start][4]

    print('done slicing %e' % (time.time()-t_start_slicing))
    return (redshift, qties['galaxy_id'][:lim],
            sed_names, mag_norms, av_arr, rv_arr, lsst_fluxes)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--healpix', type=int, default=None,
                        help='The healpixel to fit')
    parser.add_argument('--out_dir', type=str, default=None,
                        help='The directory in which to write the output file')
    parser.add_argument('--lim', type=int, default=None,
                        help='The number of galaxies to fit (if you are just testing)')
    parser.add_argument('--n_threads', type=int, default=60,
                        help='The number of threads to use')
    parser.add_argument('--out_name', type=str, default=None,
                        help='The name of the output file')
    parser.add_argument('--catalog', type=str, default='cosmoDC2_v1.0_image',
                        help='The name of the extragalactic catalog '
                             '(defaults to cosmoDC2_v1.0_image)')
    parser.add_argument('--validate', type=int, default=0,
                        help='number of galaxies to randomly validate '
                        '(defaults to zero)')

    args = parser.parse_args()
    assert args.healpix is not None
    assert args.out_dir is not None
    assert args.out_name is not None
    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir)

    sed_dir = getPackageDir('sims_sed_library')

    cat = GCRCatalogs.load_catalog(args.catalog)
    h_query = GCRQuery('healpix_pixel==%d' % args.healpix)
    if args.lim is None:
        gid = cat.get_quantities('galaxy_id', native_filters=[h_query])['galaxy_id']
        args.lim = 2*len(gid)

    out_file_name = os.path.join(args.out_dir,args.out_name)

    t_start = time.time()
    ########## actually fit SED, magNorm, and dust parameters to disks and bulges

    t0 = 1539899570.0
    print('starting %d at %.2f' % (args.healpix, time.time()-t0))

    (disk_redshift, disk_id, disk_sed_name, disk_magnorm,
     disk_av, disk_rv, disk_lsst_fluxes) = do_fitting(cat, 'disk',
                                                      args.healpix, args.lim,
                                                      args.n_threads)

    print("fit disks %d at %.2f" % (args.healpix, time.time()-t0))

    (bulge_redshift, bulge_id, bulge_sed_name, bulge_magnorm,
     bulge_av, bulge_rv, bulge_lsst_fluxes) = do_fitting(cat, 'bulge',
                                                         args.healpix, args.lim,
                                                         args.n_threads)

    print("fit bulges %d at %.2f" % (args.healpix, time.time()-t0))

    np.testing.assert_array_equal(disk_id, bulge_id)
    np.testing.assert_array_equal(disk_redshift, bulge_redshift)

    tot_lsst_fluxes = bulge_lsst_fluxes + disk_lsst_fluxes

    ############ get true values of magnitudes from extragalactic catalog;
    ############ adjust magNorm to demand agreement

    q_list = ['galaxy_id', 'ra', 'dec', 'redshift']
    for bp in 'ugrizy':
        q_list.append('mag_true_%s_lsst' % bp)

    control_qties = cat.get_quantities(q_list, native_filters=[h_query])
    for kk in control_qties:
        control_qties[kk] = control_qties[kk][:args.lim]

    print("got controls %d at %.2f" % (args.healpix, time.time()-t0))

    np.testing.assert_array_equal(control_qties['galaxy_id'], disk_id)

    dummy_spec = Sed()
    for i_bp, bp in enumerate('ugrizy'):
        cosmodc2_flux = dummy_spec.fluxFromMag(control_qties['mag_true_%s_lsst' % bp])
        flux_ratio = cosmodc2_flux/tot_lsst_fluxes[i_bp]
        dmag = -2.5*np.log10(flux_ratio)
        disk_magnorm[i_bp,:] += dmag
        bulge_magnorm[i_bp,:] += dmag

    ############# save everything in an hdf5 file

    sed_dir = os.path.join(os.environ['SIMS_SED_LIBRARY_DIR'], 'galaxySED')
    list_of_seds = os.listdir(sed_dir)
    list_of_seds.sort()
    sed_names = np.empty(len(list_of_seds), dtype=(bytes, 100))
    sed_idx = np.empty(len(list_of_seds), dtype=int)
    name_to_int = {}
    for ii, name in enumerate(list_of_seds):
        full_name = os.path.join('galaxySED', name)
        name_to_int[full_name] = ii
        sed_idx[ii] = ii
        sed_names[ii] = full_name

    disk_sed_idx = np.array([name_to_int[nn] for nn in disk_sed_name])
    bulge_sed_idx = np.array([name_to_int[nn] for nn in bulge_sed_name])

    with h5py.File(out_file_name, 'w') as out_file:
        out_file.create_dataset('galaxy_id', data=control_qties['galaxy_id'])
        out_file.create_dataset('ra', data=control_qties['ra'])
        out_file.create_dataset('dec', data=control_qties['dec'])
        out_file.create_dataset('sed_names', data=sed_names)
        out_file.create_dataset('disk_sed', data=disk_sed_idx)
        out_file.create_dataset('bulge_sed', data=bulge_sed_idx)
        out_file.create_dataset('bulge_magnorm', data=bulge_magnorm)
        out_file.create_dataset('disk_magnorm', data=disk_magnorm)
        out_file.create_dataset('disk_av', data=disk_av)
        out_file.create_dataset('disk_rv', data=disk_rv)
        out_file.create_dataset('bulge_av', data=bulge_av)
        out_file.create_dataset('bulge_rv', data=bulge_rv)
        out_file.create_dataset('bulge_fluxes', data=bulge_lsst_fluxes)
        out_file.create_dataset('disk_fluxes', data=disk_lsst_fluxes)
        out_file.create_dataset('tot_fluxes', data=tot_lsst_fluxes)

    duration = (time.time()-t_start)/3600.0
    print('all done %d at %.2f duration %.4f hrs' %
          (args.healpix, time.time()-t0, duration))

    if args.validate > 0:
        print("doing some validation")
        n_to_validate = args.validate
        validate_rng = np.random.RandomState(88)
        if n_to_validate < len(control_qties['redshift']):
            dexes_to_validate = validate_rng.choice(
                       np.arange(len(control_qties['redshift']), dtype=int),
                       replace=False,
                       size=n_to_validate)
        else:
            dexes_to_validate = np.arange(len(control_qties['redshift']),
                                          dtype=int)
        max_offset = -1.0
        sed_dir = os.environ['SIMS_SED_LIBRARY_DIR']
        tot_bp_dict = BandpassDict.loadTotalBandpassesFromFiles()
        for ii in dexes_to_validate:
            disk_name = os.path.join(sed_dir,
                                     sed_names[disk_sed_idx[ii]].decode())
            bulge_name = os.path.join(sed_dir,
                                      sed_names[bulge_sed_idx[ii]].decode())
            for i_bp, bp in enumerate('ugrizy'):
                disk_s = Sed()
                disk_s.readSED_flambda(disk_name)
                fnorm = getImsimFluxNorm(disk_s, disk_magnorm[i_bp][ii])
                disk_s.multiplyFluxNorm(fnorm)
                ax, bx = disk_s.setupCCM_ab()
                disk_s.addDust(ax, bx, A_v=disk_av[ii], R_v=disk_rv[ii])
                disk_s.redshiftSED(control_qties['redshift'][ii], dimming=True)
                disk_f = disk_s.calcFlux(tot_bp_dict[bp])

                bulge_s = Sed()
                bulge_s.readSED_flambda(bulge_name)
                fnorm = getImsimFluxNorm(bulge_s, bulge_magnorm[i_bp][ii])
                bulge_s.multiplyFluxNorm(fnorm)
                ax, bx = bulge_s.setupCCM_ab()
                bulge_s.addDust(ax, bx, A_v=bulge_av[ii], R_v=bulge_rv[ii])
                bulge_s.redshiftSED(control_qties['redshift'][ii], dimming=True)
                bulge_f = bulge_s.calcFlux(tot_bp_dict[bp])

                tot_f = disk_f+bulge_f
                f_true = dummy_spec.fluxFromMag(control_qties['mag_true_%s_lsst' % bp][ii])
                offset = np.abs(1.0-(tot_f/f_true))
                if offset>max_offset:
                    print('final max_offset %e' % offset)
                    max_offset = offset
