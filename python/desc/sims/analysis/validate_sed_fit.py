import os
import h5py
import numpy as np
import healpy
import multiprocessing
import time
from os.path import expanduser

import GCRCatalogs
from GCR import GCRQuery

import sys
home = expanduser("~")
path_to_sims_GCRCatSimInterface = os.path.join(home, 'sims_GCRCatSimInterface/python')
sys.path.insert(0, path_to_sims_GCRCatSimInterface)

from desc.sims.GCRCatSimInterface import sed_filter_names_from_catalog
from desc.sims.GCRCatSimInterface import sed_from_galacticus_mags
from lsst.sims.photUtils import BandpassDict, Sed
from lsst.sims.photUtils import cache_LSST_seds, getImsimFluxNorm
from lsst.utils import getPackageDir

import argparse

alt_sed_library_dir = '/global/cscratch1/sd/danielsf/extincted_galaxy_seds'

def do_fitting(cat, component, healpix, lim, use_extincted=False, use_alt_library=False):
    """
    Fit a set of components to SEDs, magNorm using sed_from_galacticus_mags

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

    filter_data = sed_filter_names_from_catalog(cat, use_extincted=use_extincted)
    filter_names = filter_data[component]['filter_name']
    wav_min = filter_data[component]['wav_min']
    wav_width = filter_data[component]['wav_width']

    H0 = cat.cosmology.H0.value
    Om0 = cat.cosmology.Om0

    healpix_query = GCRQuery('healpix_pixel==%d' % healpix)

    qties = cat.get_quantities(list(filter_names) +
                              ['redshift_true', 'galaxy_id'],
                               native_filters=[healpix_query])

    alternate_sed_library_dir = alt_sed_library_dir if use_alt_library else None

    print("testing on %d of %d" % (lim, len(qties['galaxy_id'])))
    with np.errstate(divide='ignore', invalid='ignore'):
        mag_array = np.array([-2.5*np.log10(qties[ff][:lim]) for ff in filter_names])

    #find NN (in color space) CatSim SED to match SED colors from cosmoDC2 catalog 
    (sed_names, mag_norms, sed_dists) = sed_from_galacticus_mags(mag_array,
                                            qties['redshift_true'][:lim],
                                            H0, Om0,
                                            wav_min, wav_width, alternate_sed_library_dir=alternate_sed_library_dir)

    return (qties['redshift_true'][:lim], qties['galaxy_id'][:lim],
            sed_names, mag_norms, sed_dists)


def calc_mags(disk_sed_list, disk_magnorm_list, bulge_sed_list, bulge_magnorm_list, 
              redshift, distance_modulus, out_dict, out_tag, observer_frame=True):
    """
    Calculate the magnitudes of galaxies as fit by CatSim.
    Designed to be run on several threads at once.
    Optionally calculate magnitudes in observer frame

    Parameters
    ----------
    disk_sed_list -- array of SED names for disks

    disk_magnorm_list -- array of magNorm for disks

    bulge_sed_list -- array of SED names for bulges

    bulge_magnorm_list -- array of magNorm for bulges

    out_dict -- a multiprocessing.Manager().dict() to store the results
    (results will be numpy array(s) of magnitudes of shape (6, N_galaxies))

    tag -- the key value in out_dict indicating this chunk of galaxies
    """

    bp_dict = BandpassDict.loadTotalBandpassesFromFiles()
    fit_mags_rest = np.zeros((6,len(disk_sed_list)), dtype=float)
    fit_mags_obs = np.zeros((6,len(disk_sed_list)), dtype=float)

    ax = None
    bx = None
    ccm_w = None
    t_start = time.time()
    for ii in range(len(disk_sed_list)):
        if ii>0 and ii%10000==0:
            dur = (time.time()-t_start)/3600.0
            pred = len(disk_sed_list)*dur/ii
            print('%d of %d; %.2e hrs left' % (ii,len(disk_sed_list), pred-dur))

        # load the disk SED
        disk_sed = Sed()
        disk_sed.readSED_flambda(os.path.join(sed_dir, disk_sed_list[ii]))

        # apply redshift to SED 
        disk_sed_obs = Sed(wavelen=disk_sed.wavelen, flambda=disk_sed.flambda)
        disk_sed_obs.redshiftSED(redshift[ii], dimming=True)
        
        # normalize the disk SED
        fnorm = getImsimFluxNorm(disk_sed, disk_magnorm_list[ii])
        disk_sed.multiplyFluxNorm(fnorm)
        disk_sed_obs.multiplyFluxNorm(fnorm)

        # apply dust to the diskSED
        #if ax is None or not np.array_equal(disk_sed.wavelen, ccm_w):
        #    ax, bx = disk_sed.setupCCMab()
        #    ccm_w = np.copy(disk_sed.wavelen)
        #disk_sed.addCCMDust(ax, bx, A_v=disk_av_list[ii], R_v=disk_rv_list[ii])
        disk_fluxes = bp_dict.fluxListForSed(disk_sed)
        disk_fluxes_obs = bp_dict.fluxListForSed(disk_sed_obs)

        # load the bluge SED
        bulge_sed = Sed()
        bulge_sed.readSED_flambda(os.path.join(sed_dir, bulge_sed_list[ii]))

        # apply redshift to SED
        bulge_sed_obs = Sed(wavelen=bulge_sed.wavelen, flambda=bulge_sed.flambda)
        bulge_sed_obs.redshiftSED(redshift[ii], dimming=True)

        # normalize the bulge SED
        fnorm = getImsimFluxNorm(bulge_sed, bulge_magnorm_list[ii])
        bulge_sed.multiplyFluxNorm(fnorm)
        bulge_sed_obs.multiplyFluxNorm(fnorm)

        # apply dust to the bulge SED
        #if ax is None or not np.array_equal(bulge_sed.wavelen, ccm_w):
        #    ax, bx = bulge_sed.setupCCMab()
        #    ccm_w = np.copy(bulge_sed.wavelen)
        #bulge_sed.addCCMDust(ax, bx, A_v=bulge_av_list[ii], R_v=bulge_rv_list[ii])
        bulge_fluxes = bp_dict.fluxListForSed(bulge_sed)
        bulge_fluxes_obs = bp_dict.fluxListForSed(bulge_sed_obs)

        # combine disk and bulge SED to get total galaxy magnitudes
        fluxes = bulge_fluxes + disk_fluxes
        mags = disk_sed.magFromFlux(fluxes)
        fluxes_obs = bulge_fluxes_obs + disk_fluxes_obs
        mags_obs = disk_sed_obs.magFromFlux(fluxes_obs)
        fit_mags_rest[:,ii] = mags
        fit_mags_obs[:,ii] = mags_obs + distance_modulus[ii]

    out_dict[out_tag] = {'rest':fit_mags_rest, 'obs':fit_mags_obs}


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
    parser.add_argument('--use_extincted', help='Use extincted fluxes for SED matching', default=False, action='store_true')
    parser.add_argument('--save_unextincted', help='Save CosmoDC2 unextincted fluxes', default=False, action='store_true')
    parser.add_argument('--use_alternate_sed_library', help='Use alternate extincted sed library', default=False, action='store_true')
    args = parser.parse_args()
    assert args.healpix is not None
    assert args.out_dir is not None
    assert args.out_name is not None
    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir)
    
    sed_dir = alt_sed_library_dir if args.use_alternate_sed_library else  getPackageDir('sims_sed_library')
    print('Using sed library in {}'.format(sed_dir))

    cat = GCRCatalogs.load_catalog('cosmoDC2_v1.0_image')

    out_file_name = os.path.join(args.out_dir,args.out_name)

    ########## actually fit SED, magNorm, and dust parameters to disks and bulges

    (disk_redshift, disk_id, disk_sed_name, disk_mag,
     disk_dist) = do_fitting(cat, 'disk', args.healpix, args.lim, use_extincted=args.use_extincted, use_alt_library=args.use_alternate_sed_library)

    print("fit disks")

    (bulge_redshift, bulge_id, bulge_sed_name, bulge_mag,
     bulge_dist) = do_fitting(cat, 'bulge', args.healpix, args.lim, use_extincted=args.use_extincted, use_alt_library=args.use_alternate_sed_library)

    print("fit bulges")

    np.testing.assert_array_equal(disk_id, bulge_id)
    np.testing.assert_array_equal(disk_redshift, bulge_redshift)

    ############ get true values of magnitudes from extragalactic catalog

    q_list = ['galaxy_id', 'stellar_mass', 'stellar_mass_disk', 'stellar_mass_bulge']
    Mag_template = 'Mag_true_{}_lsst_z0_no_host_extinction' if args.save_unextincted else 'Mag_true_{}_lsst_z0'
    mag_template = 'mag_true_{}_lsst_no_host_extinction' if args.save_unextincted else 'mag_true_{}_lsst'
    for bp in 'ugrizy':
        q_list.append(Mag_template.format(bp))
        q_list.append(mag_template.format(bp))

    print('CosmoDC2 quantities: {}'.format(', '.join(q_list)))

    h_query = GCRQuery('healpix_pixel==%d' % args.healpix)
    control_qties = cat.get_quantities(q_list, native_filters=[h_query])
    for kk in control_qties:
        control_qties[kk] = control_qties[kk][:args.lim]

    print("got controls")

    distance_modulus = cat.cosmology.distmod(disk_redshift).value
    print('computed distance moduli')

    np.testing.assert_array_equal(control_qties['galaxy_id'], disk_id)

    ############# use multiprocessing to calculate the CatSim fit magnitudes of the galaxies

    p_list = []
    mgr = multiprocessing.Manager()
    out_dict = mgr.dict()
    fit_mags_rest = np.zeros((6, len(disk_mag)), dtype=float)
    fit_mags_obs = np.zeros((6, len(disk_mag)), dtype=float)
    d_gal = len(disk_mag)//23
    for i_start in range(0, len(disk_mag), d_gal):
        i_end = i_start + d_gal
        selection = slice(i_start, i_end)
        p = multiprocessing.Process(target=calc_mags,
                                    args=(disk_sed_name[selection],
                                          disk_mag[selection],
                                          bulge_sed_name[selection],
                                          bulge_mag[selection],
                                          disk_redshift[selection],
                                          distance_modulus[selection],
                                          out_dict,
                                          i_start))

        p.start()
        p_list.append(p)

    for p in p_list:
        p.join()

    for i_start in range(0, len(disk_mag), d_gal):
        i_end = i_start+d_gal
        fit_mags_rest[:,i_start:i_end] = out_dict[i_start]['rest']
        fit_mags_obs[:,i_start:i_end] = out_dict[i_start]['obs']

    ############# save everything in an hdf5 file

    with h5py.File(out_file_name, 'w') as out_file:
        out_file.create_dataset('galaxy_id', data=control_qties['galaxy_id'])
        out_file.create_dataset('stellar_mass', data=control_qties['stellar_mass'])
        out_file.create_dataset('disk_sed', data=[s.encode('utf-8') for s in disk_sed_name])
        out_file.create_dataset('bulge_sed', data=[s.encode('utf-8') for s in bulge_sed_name])
        out_file.create_dataset('stellar_mass_disk', data=control_qties['stellar_mass_disk'])
        out_file.create_dataset('disk_dist', data=disk_dist)
        out_file.create_dataset('stellar_mass_bulge', data=control_qties['stellar_mass_bulge'])
        out_file.create_dataset('bulge_dist', data=bulge_dist)
        out_file.create_dataset('redshift_true', data=disk_redshift)
        for i_bp, bp in enumerate('ugrizy'):
            out_file.create_dataset('CatSim_fit_rest_%s' % bp, data=fit_mags_rest[i_bp])
            out_file.create_dataset('CatSim_fit_obs_%s' % bp, data=fit_mags_obs[i_bp])
            out_file.create_dataset('cosmoDC2_rest_%s' % bp,
                                    data=control_qties[Mag_template.format(bp)])
            out_file.create_dataset('cosmoDC2_obs_%s' % bp,
                                    data=control_qties[mag_template.format(bp)])
