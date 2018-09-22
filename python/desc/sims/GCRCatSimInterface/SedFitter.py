import os
import re
import numpy as np
import GCRCatalogs
import time
import scipy.spatial as scipy_spatial
from lsst.utils import getPackageDir
from lsst.sims.utils import defaultSpecMap
from lsst.sims.photUtils import BandpassDict, Bandpass, Sed, CosmologyObject

__all__ = ["disk_re", "bulge_re", "sed_filter_names_from_catalog", "sed_from_galacticus_mags"]

_galaxy_sed_dir = os.path.join(getPackageDir('sims_sed_library'), 'galaxySED')

disk_re = re.compile(r'sed_(\d+)_(\d+)_disk$')
bulge_re = re.compile(r'sed_(\d+)_(\d+)_bulge$')

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

    all_quantities = catalog.list_all_quantities()

    bulge_names = []
    bulge_wav_min = []
    bulge_wav_width = []

    disk_names = []
    disk_wav_min = []
    disk_wav_width = []

    for qty_name in all_quantities:
        disk_match = disk_re.match(qty_name)
        if disk_match is not None:
            disk_names.append(qty_name)
            disk_wav_min.append(0.1*float(disk_match[1]))  # 0.1 converts to nm
            disk_wav_width.append(0.1*float(disk_match[2]))

        bulge_match = bulge_re.match(qty_name)
        if bulge_match is not None:
            bulge_names.append(qty_name)
            bulge_wav_min.append(0.1*float(bulge_match[1]))
            bulge_wav_width.append(0.1*float(bulge_match[2]))

    disk_wav_min = np.array(disk_wav_min)
    disk_wav_width = np.array(disk_wav_width)
    disk_names = np.array(disk_names)
    sorted_dex = np.argsort(disk_wav_min)
    disk_wav_width = disk_wav_width[sorted_dex]
    disk_names = disk_names[sorted_dex]
    disk_wav_min = disk_wav_min[sorted_dex]

    bulge_wav_min = np.array(bulge_wav_min)
    bulge_wav_width = np.array(bulge_wav_width)
    bulge_names = np.array(bulge_names)
    sorted_dex = np.argsort(bulge_wav_min)
    bulge_wav_width = bulge_wav_width[sorted_dex]
    bulge_names = bulge_names[sorted_dex]
    bulge_wav_min = bulge_wav_min[sorted_dex]

    return {'disk':{'filter_name': disk_names,
                    'wav_min': disk_wav_min,
                    'wav_width': disk_wav_width},
            'bulge':{'filter_name': bulge_names,
                     'wav_min': bulge_wav_min,
                     'wav_width': bulge_wav_width}}


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
    sed_names is an array containing the names of the SED files repeated over
    combinations of dust parameters (sorry; that wording is awkward)

    sed_mag_list is MxN float array, with M = number of SED file, dust parameter
    combinations in the library, and N = number of top hat filters in the catalog

    sed_mag_norm is 1d float array, with length = number of SED file, dust parameter
    combinations in the library

    av_out_list is a 1d float array of Av

    rv_out_list is a 1d float array of Rv
    """

    av_grid = np.arange(0.0, 3.0, 0.1)
    rv_grid = np.arange(1.0, 5.0, 0.1)

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
    av_out_list = list()
    rv_out_list = list()

    imsim_bp = Bandpass()
    imsim_bp.imsimBandpass()

    n_tot = len(os.listdir(_galaxy_sed_dir))*len(av_grid)*len(rv_grid)
    i_gen = 0
    t_start = time.time()
    print('\n\ncreating library')
    for sed_file_name in os.listdir(_galaxy_sed_dir):
        base_spec = Sed()
        base_spec.readSED_flambda(os.path.join(_galaxy_sed_dir, sed_file_name))
        ax, bx = base_spec.setupCCMab()
        for av in av_grid:
            for rv in rv_grid:
                if i_gen>0 and i_gen%1000==0:
                    duration = (time.time()-t_start)/3600.0
                    predicted = n_tot*duration/i_gen
                    print('%d of %d; dur %.2e pred %.2e' % (i_gen, n_tot, duration, predicted))
                spec = Sed(wavelen=base_spec.wavelen, flambda=base_spec.flambda)
                sed_names.append(defaultSpecMap[sed_file_name])
                sed_mag_norm.append(spec.calcMag(imsim_bp))
                spec.addCCMDust(ax, bx, A_v=av, R_v=rv)
                av_out_list.append(av)
                rv_out_list.append(rv)
                sed_mag_list.append(tuple(bandpass_dict.magListForSed(spec)))
                i_gen += 1

    print('made library')

    return (np.array(sed_names), np.array(sed_mag_list), np.array(sed_mag_norm),
            np.array(av_out_list), np.array(rv_out_list))


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

    if (not hasattr(sed_from_galacticus_mags, '_color_tree') or
        not np.allclose(wav_min, sed_from_galacticus_mags._wav_min,
                        atol=1.0e-10, rtol=0.0) or
        not np.allclose(wav_width, sed_from_galacticus_mags._wav_width,
                        atol=1.0e-10, rtol=0.0)):

        (sed_names,
         sed_mag_list,
         sed_mag_norm,
         av_grid, rv_grid) = _create_sed_library_mags(wav_min, wav_width)


        sed_colors = sed_mag_list[:,1:] - sed_mag_list[:,:-1]
        sed_from_galacticus_mags._sed_names = sed_names
        sed_from_galacticus_mags._mag_norm = sed_mag_norm # N_sed
        sed_from_galacticus_mags._av_grid = av_grid
        sed_from_galacticus_mags._rv_grid = rv_grid
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

    return (output_names, output_mag_norm,
            sed_from_galacticus_mags._av_grid[sed_idx],
            sed_from_galacticus_mags._rv_grid[sed_idx])
