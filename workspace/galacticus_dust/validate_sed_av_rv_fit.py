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

def _parallel_do_fitting(mag_array, redshift, H0, Om0, wav_min, wav_width,
                         out_dict, tag):

    (sed_names,
     mag_norms) = sed_from_galacticus_mags(mag_array,
                                        redshift,
                                        H0, Om0,
                                        wav_min, wav_width)

    out_dict[tag] = (sed_names, mag_norms)



def do_fitting(cat, component, healpix, lim):
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

    redshift = qties['redshift_true'][:lim]

    (sed_names,
     mag_norms) = sed_from_galacticus_mags(mag_array[:,:2],
                                           redshift[:2],
                                           H0, Om0,
                                           wav_min, wav_width)

    n_gal = len(redshift)
    d_gal = n_gal//23

    p_list= []
    mgr = multiprocessing.Manager()
    out_dict = mgr.dict()
    for i_start in range(0,n_gal,d_gal):
        s = slice(i_start,i_start+d_gal)
        p = multiprocessing.Process(target=_parallel_do_fitting,
                                    args=(mag_array[:,s], redshift[s],
                                          H0, Om0, wav_min, wav_width,
                                          out_dict, i_start))
        p.start()
        p_list.append(p)

    for p in p_list:
        p.join()

    sed_names = np.empty(len(redshift), dtype=(str,200))
    mag_norms = np.zeros(len(redshift), dtype=float)
    for i_start in out_dict.keys():
        s = slice(i_start,i_start+d_gal)
        sed_names[s] = out_dict[i_start][0]
        mag_norms[s] = out_dict[i_start][1]

    return (qties['redshift_true'][:lim], qties['galaxy_id'][:lim],
            sed_names, mag_norms)



def calc_mags(disk_sed_list, disk_magnorm_list,
              bulge_sed_list, bulge_magnorm_list,
              out_dict, out_tag):
    """
    Calculate the magnitudes of galaxies as fit by CatSim.
    Designed to be run on several threads at once.

    Parameters
    ----------
    disk_sed_list -- array of SED names for disks

    disk_magnorm_list -- array of magNorm for disks

    bulge_sed_list -- array of SED names for bulges

    bulge_magnorm_list -- array of magNorm for bulges

    out_dict -- a multiprocessing.Manager().dict() to store the results
    (results will be a numpy array of magnitudes of shape (6, N_galaxies))

    tag -- the key value in out_dict indicating this chunk of galaxies
    """

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

        # load the disk SED
        disk_sed = Sed()
        disk_sed.readSED_flambda(disk_sed_list[ii])

        # normalize the disk SED
        fnorm = getImsimFluxNorm(disk_sed, disk_magnorm_list[ii])
        disk_sed.multiplyFluxNorm(fnorm)

        disk_fluxes = bp_dict.fluxListForSed(disk_sed)

        # load the bluge SED
        bulge_sed = Sed()
        bulge_sed.readSED_flambda(bulge_sed_list[ii])

        # normalize the bulge SED
        fnorm = getImsimFluxNorm(bulge_sed, bulge_magnorm_list[ii])
        bulge_sed.multiplyFluxNorm(fnorm)

        bulge_fluxes = bp_dict.fluxListForSed(bulge_sed)

        # combine disk and bulge SED to get total galaxy magnitudes
        fluxes = bulge_fluxes + disk_fluxes
        mags = disk_sed.magFromFlux(fluxes)
        fit_mags[:,ii] = mags

    out_dict[out_tag] = fit_mags


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--healpix', type=int, default=None,
                        help='The healpixel to fit')
    parser.add_argument('--out_dir', type=str, default=None,
                        help='The directory in which to write the output file')
    parser.add_argument('--lim', type=int, default=180000000,
                        help='The number of galaxies to fit (if you are just testing)')
    parser.add_argument('--out_name', type=str, default=None,
                        help='The name of the output file')
    args = parser.parse_args()
    assert args.healpix is not None
    assert args.out_dir is not None
    assert args.out_name is not None
    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir)

    sed_dir = getPackageDir('sims_sed_library')

    cat = GCRCatalogs.load_catalog('cosmoDC2_v1.0_image')

    out_file_name = os.path.join(args.out_dir,args.out_name)

    ########## actually fit SED, magNorm, and dust parameters to disks and bulges

    (disk_redshift, disk_id,
     disk_sed_name, disk_mag) = do_fitting(cat, 'disk',
                                           args.healpix, args.lim)

    print("fit disks")

    (bulge_redshift, bulge_id,
     bulge_sed_name, bulge_mag) = do_fitting(cat, 'bulge',
                                             args.healpix, args.lim)

    print("fit bulges")

    np.testing.assert_array_equal(disk_id, bulge_id)
    np.testing.assert_array_equal(disk_redshift, bulge_redshift)

    ############ get true values of magnitudes from extragalactic catalog

    q_list = ['galaxy_id']
    for bp in 'ugrizy':
        q_list.append('Mag_true_%s_lsst_z0' % bp)

    h_query = GCRQuery('healpix_pixel==%d' % args.healpix)
    control_qties = cat.get_quantities(q_list, native_filters=[h_query])
    for kk in control_qties:
        control_qties[kk] = control_qties[kk][:args.lim]

    print("got controls")

    np.testing.assert_array_equal(control_qties['galaxy_id'], disk_id)

    ############# use multiprocessing to calculate the CatSim fit magnitudes of the galaxies

    p_list = []
    mgr = multiprocessing.Manager()
    out_dict = mgr.dict()
    fit_mags = np.zeros((6, len(disk_sed_name)), dtype=float)
    d_gal = len(disk_sed_name)//63
    for i_start in range(0, len(disk_sed_name), d_gal):
        i_end = i_start + d_gal
        selection = slice(i_start, i_end)
        p = multiprocessing.Process(target=calc_mags,
                                    args=(disk_sed_name[selection],
                                          disk_mag[selection],
                                          bulge_sed_name[selection],
                                          bulge_mag[selection],
                                          out_dict,
                                          i_start))

        p.start()
        p_list.append(p)

    for p in p_list:
        p.join()

    for i_start in range(0, len(disk_sed_name), d_gal):
        i_end = i_start+d_gal
        fit_mags[:,i_start:i_end] = out_dict[i_start]

    ############# save everything in an hdf5 file

    with h5py.File(out_file_name, 'w') as out_file:
        out_file.create_dataset('galaxy_id', data=control_qties['galaxy_id'])
        out_file.create_dataset('disk_sed', data=[s.encode('utf-8') for s in disk_sed_name])
        out_file.create_dataset('bulge_sed', data=[s.encode('utf-8') for s in bulge_sed_name])
        for i_bp, bp in enumerate('ugrizy'):
            out_file.create_dataset('fit_%s' % bp, data=fit_mags[i_bp])
            out_file.create_dataset('cosmo_%s' % bp,
                                    data=control_qties['Mag_true_%s_lsst_z0' % bp])
