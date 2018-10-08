import os
import re
import numpy as np
import multiprocessing
import time
import GCRCatalogs
import scipy.spatial as scipy_spatial
from lsst.utils import getPackageDir
from lsst.sims.utils import defaultSpecMap
from lsst.sims.photUtils import BandpassDict, Bandpass, Sed, CosmologyObject

__all__ = ["disk_re", "bulge_re", "sed_filter_names_from_catalog", "sed_from_galacticus_mags"]

_galaxy_sed_dir = os.path.join(os.environ['SCRATCH'], 'extincted_galaxy_seds')

def sed_filter_names_from_catalog(catalog):
    """
    Takes an already-loaded GCR catalog and returns the names, wavelengths,
    and widths of the SED-defining bandpasses

    Parameters
    ----------
    catalog -- is a catalog loaded with GCR.load_catalog()

    Returns
    -------
    A dict keyed to 'bulge' and 'disk'.  The values in this dict will
    be dicts keyed to 'filter_name', 'wav_min', 'wav_width'.  The
    corresponding values are:

    filter_name -- list of the names of the columns defining the SED
    wav_min -- list of the minimum wavelengths of SED-defining bandpasses (in nm)
    wav_width -- list of the widths of the SED-defining bandpasses (in nm)

    All outputs will be returned in order of increasing wav_min
    """

    disk_names = []
    bulge_names = []
    for bp in 'ugrizy':
        disk_names.append('LSST_filters/diskLuminositiesStellar:LSST_%s:rest:dustAtlas' % bp)
        bulge_names.append('LSST_filters/bulgeLuminositiesStellar:LSST_%s:rest:dustAtlas' % bp)

    return {'disk':{'filter_name': disk_names,
                    'wav_min': None,
                    'wav_width': None},
            'bulge':{'filter_name': bulge_names,
                     'wav_min': None,
                     'wav_width': None}}

def _parallel_sed_mags(sed_file_name_list, bandpass_dict, out_dict, tag):
    imsim_bp = Bandpass()
    imsim_bp.imsimBandpass()
    sed_names = np.empty(len(sed_file_name_list), dtype=(str, 200))
    sed_mag_list = np.zeros((len(sed_file_name_list), len(bandpass_dict)), dtype=float)
    sed_mag_norm = np.zeros(len(sed_file_name_list), dtype=float)
    t_start = time.time()

    assert len(bandpass_dict)==30

    for i_sed, sed_file_name in enumerate(sed_file_name_list):
        if i_sed>0 and i_sed%100==0:
            d = (time.time()-t_start)/3600.0
            p = len(sed_file_name_list)*d/i_sed
            print('run %d of %d; %.2e hrs left' %
            (i_sed, len(sed_file_name_list), p-d))

        spec = Sed()
        spec.readSED_flambda(os.path.join(_galaxy_sed_dir, sed_file_name))
        sed_names[i_sed] = os.path.join(_galaxy_sed_dir, sed_file_name)
        sed_mag_list[i_sed] = bandpass_dict.magListForSed(spec)
        sed_mag_norm[i_sed] = spec.calcMag(imsim_bp)

    out_dict[tag] = (sed_names, sed_mag_list, sed_mag_norm)


def _create_sed_library_mags(wav_min, wav_width):
    """
    Calculate the magnitudes of the SEDs in sims_sed_library dir in the
    tophat filters specified by wav_min, wav_width

    Parameters
    ----------
    wav_min is a numpy array of the minimum wavelengths of the tophat
    filters (in nm)

    wav_width is a numpy array of the widths of the tophat filters (in nm)

    Returns
    -------
    sed_names is an array containing the names of the SED files

    sed_mag_list is MxN float array, with M = number of SED files in the library,
        and N = number of top hat filters in the catalog

    sed_mag_norm is 1d float array, with length = number of SED files in the library
    """

    (total_bp_dict,
     bandpass_dict) = BandpassDict.loadBandpassesFromFiles()

    imsim_bp = Bandpass()
    imsim_bp.imsimBandpass()

    print('creating arrays')

    list_of_sed_names = os.listdir(_galaxy_sed_dir)
    sed_names_out = np.empty(len(list_of_sed_names), dtype=(str,200))
    sed_mag_list_out = np.zeros((len(list_of_sed_names),
                             len(bandpass_dict)), dtype=float)
    sed_mag_norm_out = np.zeros(len(list_of_sed_names), dtype=float)

    t_start = time.time()
    d_sed = len(list_of_sed_names)//23
    p_list = []
    mgr = multiprocessing.Manager()
    out_dict = mgr.dict()
    for i_start in range(0, len(list_of_sed_names), d_sed):
        s = slice(i_start,i_start+d_sed)
        p = multiprocessing.Process(target=_parallel_sed_mags,
                                    args=(list_of_sed_names[s],
                                          bandpass_dict,
                                          out_dict, i_start))
        p.start()
        p_list.append(p)

    for p in p_list:
        p.join()

    for i_start in out_dict.keys():
        s=slice(i_start,i_start+d_sed)
        sed_names_out[s] = out_dict[i_start][0]
        sed_mag_list_out[s] = out_dict[i_start][1]
        sed_mag_norm_out[s] = out_dict[i_start][2]

    return sed_names_out, sed_mag_list_out, sed_mag_norm_out


def sed_from_galacticus_mags(galacticus_mags, redshift, H0, Om0,
                             wav_min, wav_width):
    """
    Fit SEDs from sims_sed_library to Galacticus galaxies based on the
    magnitudes in tophat filters.

    Parameters
    ----------

    galacticus_mags is a numpy array such that
    galacticus_mags[i][j] is the magnitude of the jth star in the ith bandpass,
    where the bandpasses are ordered in ascending order of minimum wavelength.

    redshift is an array of redshifts for the galaxies being fit

    H0 is the Hubbleparameter in units of km/s/Mpc

    Om0 is the critical density parameter for matter

    wav_min is a numpy array of the minimum wavelengths of the tophat
    filters (in nm)

    wav_grid is a numpy array of the widths of the tophat filters
    (in nm)

    Returns
    -------
    a numpy array of SED names and a numpy array of magNorms.
    """

    if not hasattr(sed_from_galacticus_mags, '_color_tree'):

        (sed_names,
         sed_mag_list,
         sed_mag_norm) = _create_sed_library_mags(wav_min, wav_width)


        sed_colors = sed_mag_list[:,1:] - sed_mag_list[:,:-1]
        sed_from_galacticus_mags._sed_names = sed_names
        sed_from_galacticus_mags._mag_norm = sed_mag_norm # N_sed
        sed_from_galacticus_mags._sed_mags = sed_mag_list # N_sed by N_mag
        print('making tree')
        sed_from_galacticus_mags._color_tree = scipy_spatial.cKDTree(sed_colors)
        print('made tree')
        sed_from_galacticus_mags._wav_min = wav_min
        sed_from_galacticus_mags._wav_width = wav_width

    if (not hasattr(sed_from_galacticus_mags, '_cosmo') or
        np.abs(sed_from_galacticus_mags._cosmo.H()-H0)>1.0e-6 or
        np.abs(sed_from_galacticus_mags._cosmo.OmegaMatter()-Om0)>1.0e-6):

        sed_from_galacticus_mags._cosmo = CosmologyObject(H0=H0, Om0=Om0)

    galacticus_mags_t = np.asarray(galacticus_mags).T # N_star by N_mag
    assert galacticus_mags_t.shape == (len(redshift), sed_from_galacticus_mags._sed_mags.shape[1])

    with np.errstate(invalid='ignore', divide='ignore'):
        galacticus_colors = galacticus_mags_t[:,1:] - galacticus_mags_t[:,:-1] # N_star by (N_mag - 1)

    (sed_dist,
     sed_idx) = sed_from_galacticus_mags._color_tree.query(galacticus_colors, k=1)

    # cKDTree returns an invalid index (==len(tree_data)) in cases
    # where the distance is not finite
    sed_idx = np.where(sed_idx<len(sed_from_galacticus_mags._sed_names),
                       sed_idx, 0)

    distance_modulus = sed_from_galacticus_mags._cosmo.distanceModulus(redshift=redshift)
    output_names = sed_from_galacticus_mags._sed_names[sed_idx]
    d_mag = (galacticus_mags_t - sed_from_galacticus_mags._sed_mags[sed_idx]).mean(axis=1)
    output_mag_norm = sed_from_galacticus_mags._mag_norm[sed_idx] + d_mag + distance_modulus

    return output_names, output_mag_norm
