import os
import re
import numpy as np
import GCRCatalogs
import scipy.spatial as scipy_spatial
from lsst.utils import getPackageDir
from lsst.sims.utils import defaultSpecMap
from lsst.sims.photUtils import BandpassDict, Bandpass, Sed, CosmologyObject

from SedFitter import disk_re, bulge_re, sed_filter_names_from_catalog

__all__ = ["sed_from_galacticus_mags_classic"]


def _create_sed_library_mags_classic(wav_min, wav_width):
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

    wav_max = max((wav0+width
                  for wav0, width in zip(wav_min, wav_width)))
    wav_grid = np.arange(wav_min.min(), wav_max, 0.1)

    bp_name_list = list()
    bp_list = list()
    for wav0, width in zip(wav_min, wav_width):
        sb_grid = ((wav_grid >= wav0) & (wav_grid <= (wav0+width))).astype(float)
        bp_list.append(Bandpass(wavelen=wav_grid, sb=sb_grid))
        bp_name_list.append('%d_%d' % (wav0, width))

    bandpass_dict = BandpassDict(bp_list, bp_name_list)

    sed_names = list()
    sed_mag_list = list()
    sed_mag_norm = list()

    imsim_bp = Bandpass()
    imsim_bp.imsimBandpass()

    for sed_file_name in os.listdir(_galaxy_sed_dir):
        spec = Sed()
        spec.readSED_flambda(os.path.join(_galaxy_sed_dir, sed_file_name))
        sed_names.append(defaultSpecMap[sed_file_name])
        sed_mag_list.append(tuple(bandpass_dict.magListForSed(spec)))
        sed_mag_norm.append(spec.calcMag(imsim_bp))

    return np.array(sed_names), np.array(sed_mag_list), np.array(sed_mag_norm)


def sed_from_galacticus_mags_classic(galacticus_mags, redshift, H0, Om0,
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

    if (not hasattr(sed_from_galacticus_mags, '_color_tree') or
        not np.allclose(wav_min, sed_from_galacticus_mags._wav_min,
                        atol=1.0e-10, rtol=0.0) or
        not np.allclose(wav_width, sed_from_galacticus_mags._wav_width,
                        atol=1.0e-10, rtol=0.0)):

        (sed_names,
         sed_mag_list,
         sed_mag_norm) = _create_sed_library_mags_classic(wav_min, wav_width)


        sed_colors = sed_mag_list[:,1:] - sed_mag_list[:,:-1]
        sed_from_galacticus_mags._sed_names = sed_names
        sed_from_galacticus_mags._mag_norm = sed_mag_norm # N_sed
        sed_from_galacticus_mags._sed_mags = sed_mag_list # N_sed by N_mag
        sed_from_galacticus_mags._color_tree = scipy_spatial.cKDTree(sed_colors)
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
