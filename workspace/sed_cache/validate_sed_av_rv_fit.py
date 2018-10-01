import os
import h5py
import numpy as np
import healpy
import multiprocessing
import time

import GCRCatalogs
from GCR import GCRQuery

from desc.sims.GCRCatSimInterface import sed_filter_names_from_catalog
from desc.sims.GCRCatSimInterface import sed_from_galacticus_mags
from lsst.sims.photUtils import BandpassDict, Sed
from lsst.sims.photUtils import cache_LSST_seds, getImsimFluxNorm
from lsst.utils import getPackageDir

import argparse


def do_fitting(cat, component, healpix, lim):

    filter_data = sed_filter_names_from_catalog(cat)
    filter_names = filter_data[component]['filter_name']
    wav_min = filter_data[component]['wav_min']
    wav_width = filter_data[component]['wav_width']

    H0 = cat.cosmology.H0.value
    Om0 = cat.cosmology.Om0

    healpix_query = GCRQuery('healpix_pixel==%d' % healpix)

    qties = cat.get_quantities(list(filter_names) +
                              ['redshift_true', 'galaxy_id'],
                               native_filters=[healpix_query])

    print("testing on %d of %d" % (lim, len(qties['galaxy_id'])))
    with np.errstate(divide='ignore', invalid='ignore'):
        mag_array = np.array([-2.5*np.log10(qties[ff][:lim]) for ff in filter_names])

    (sed_names,
     mag_norms,
     av_arr,
     rv_arr) = sed_from_galacticus_mags(mag_array,
                                        qties['redshift_true'][:lim],
                                        H0, Om0,
                                        wav_min, wav_width)

    return (qties['redshift_true'][:lim], qties['galaxy_id'][:lim],
            sed_names, mag_norms, av_arr, rv_arr)

def calc_mags(disk_sed_list, disk_magnorm_list, disk_av_list, disk_rv_list,
              bulge_sed_list, bulge_magnorm_list, bulge_av_list, bulge_rv_list,
              out_dict, out_tag):

    bp_dict = BandpassDict.loadTotalBandpassesFromFiles()
    fit_mags = np.zeros((6,len(disk_sed_list)), dtype=float)

    ax = None
    bx = None
    ccm_w = None
    t_start = time.time()
    for ii in range(len(disk_sed_list)):
        if ii>0 and ii%1000==0:
            dur = (time.time()-t_start)/3600.0
            pred = len(disk_sed_list)*dur/ii
            print('%d of %d; %.2e hrs left' % (ii,len(disk_sed_list), pred-dur))

        disk_sed = Sed()
        disk_sed.readSED_flambda(os.path.join(sed_dir, disk_sed_list[ii]))
        fnorm = getImsimFluxNorm(disk_sed, disk_magnorm_list[ii])
        disk_sed.multiplyFluxNorm(fnorm)
        if ax is None or not np.array_equal(disk_sed.wavelen, ccm_w):
            ax, bx = disk_sed.setupCCMab()
            ccm_w = np.copy(disk_sed.wavelen)
        disk_sed.addCCMDust(ax, bx, A_v=disk_av_list[ii], R_v=disk_rv_list[ii])
        disk_fluxes = bp_dict.fluxListForSed(disk_sed)

        bulge_sed = Sed()
        bulge_sed.readSED_flambda(os.path.join(sed_dir, bulge_sed_list[ii]))
        fnorm = getImsimFluxNorm(bulge_sed, bulge_magnorm_list[ii])
        bulge_sed.multiplyFluxNorm(fnorm)
        if ax is None or not np.array_equal(bulge_sed.wavelen, ccm_w):
            ax, bx = bulge_sed.setupCCMab()
            ccm_w = np.copy(bulge_sed.wavelen)
        bulge_sed.addCCMDust(ax, bx, A_v=bulge_av_list[ii], R_v=bulge_rv_list[ii])
        bulge_fluxes = bp_dict.fluxListForSed(bulge_sed)

        fluxes = bulge_fluxes + disk_fluxes
        mags = disk_sed.magFromFlux(fluxes)
        fit_mags[:,ii] = mags

    out_dict[out_tag] = fit_mags


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--healpix', type=int, default=None)
    parser.add_argument('--out_dir', type=str, default=None)
    parser.add_argument('--lim', type=int, default=180000000)
    parser.add_argument('--out_name', type=str, default=None)
    args = parser.parse_args()
    assert args.healpix is not None
    assert args.out_dir is not None
    assert args.out_name is not None
    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir)

    sed_dir = getPackageDir('sims_sed_library')

    cat = GCRCatalogs.load_catalog('cosmoDC2_v1.0_image')

    out_file_name = os.path.join(args.out_dir,args.out_name)


    #cache_LSST_seds(wavelen_min=0.0, wavelen_max=3000.0)

    (disk_redshift, disk_id, disk_sed_name, disk_mag,
     disk_av, disk_rv) = do_fitting(cat, 'disk', args.healpix, args.lim)

    print("fit disks")

    (bulge_redshift, bulge_id, bulge_sed_name, bulge_mag,
     bulge_av, bulge_rv) = do_fitting(cat, 'bulge', args.healpix, args.lim)

    print("fit bulges")

    np.testing.assert_array_equal(disk_id, bulge_id)
    np.testing.assert_array_equal(disk_redshift, bulge_redshift)

    q_list = ['galaxy_id']
    for bp in 'ugrizy':
        q_list.append('Mag_true_%s_lsst_z0' % bp)

    h_query = GCRQuery('healpix_pixel==%d' % args.healpix)
    control_qties = cat.get_quantities(q_list, native_filters=[h_query])
    for kk in control_qties:
        control_qties[kk] = control_qties[kk][:args.lim]

    print("got controls")

    np.testing.assert_array_equal(control_qties['galaxy_id'], disk_id)

    p_list = []
    mgr = multiprocessing.Manager()
    out_dict = mgr.dict()
    fit_mags = np.zeros((6, len(disk_av)), dtype=float)
    d_gal = len(disk_av)//23
    for i_start in range(0, len(disk_av), d_gal):
        i_end = i_start + d_gal
        selection = slice(i_start, i_end)
        p = multiprocessing.Process(target=calc_mags,
                                    args=(disk_sed_name[selection],
                                          disk_mag[selection],
                                          disk_av[selection],
                                          disk_rv[selection],
                                          bulge_sed_name[selection],
                                          bulge_mag[selection],
                                          bulge_av[selection],
                                          bulge_rv[selection],
                                          out_dict,
                                          i_start))

        p.start()
        p_list.append(p)

    for p in p_list:
        p.join()

    for i_start in range(0, len(disk_av), d_gal):
        i_end = i_start+d_gal
        fit_mags[:,i_start:i_end] = out_dict[i_start]

    with h5py.File(out_file_name, 'w') as out_file:
        out_file.create_dataset('galaxy_id', data=control_qties['galaxy_id'])
        out_file.create_dataset('disk_sed', data=[s.encode('utf-8') for s in disk_sed_name])
        out_file.create_dataset('bulge_sed', data=[s.encode('utf-8') for s in bulge_sed_name])
        out_file.create_dataset('disk_av', data=disk_av)
        out_file.create_dataset('disk_rv', data=disk_rv)
        out_file.create_dataset('bulge_av', data=bulge_av)
        out_file.create_dataset('bulge_rv', data=bulge_rv)
        for i_bp, bp in enumerate('ugrizy'):
            out_file.create_dataset('fit_%s' % bp, data=fit_mags[i_bp])
            out_file.create_dataset('cosmo_%s' % bp,
                                    data=control_qties['Mag_true_%s_lsst_z0' % bp])
