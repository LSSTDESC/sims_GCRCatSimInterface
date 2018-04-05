import os
import re
import numpy as np
import GCRCatalogs
from lsst.utils import getPackageDir
from lsst.sims.photUtils import BandpassDict, Bandpass, Sed, CosmologyObject

__all__ = ["sed_filter_names_from_catalog", "sed_from_galacticus_mags"]

_galaxy_sed_dir = os.path.join(getPackageDir('sims_sed_library'), 'galaxySED')


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

    disk_re = re.compile(r'sed_(\d+)_(\d+)_disk$')
    bulge_re = re.compile(r'sed_(\d+)_(\d+)_bulge$')

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


def _get_sed_mags_and_cosmology(catalog_name):
    """
    Returns 5 numpy arrays: sed_names, sed_mag_list, sed_mag_norm,
                            wav_min, wav_width
    and a dictionary for cosmology
    sed_names is 1d str array, with length = number of SED files in the library
    sed_mag_list is MxN float array, with M = number of SED files in the library,
        and N = number of top hat filters in the catalog
    sed_mag_norm is 1d float array, with length = number of SED files in the library
    wav_min is a 1d array of the minimum wavelength of the catalog bandpasses in nm
    wav_grid is a 1d array of the width of the catalog bandpass in nm
    cosmology = {'H0': H0, 'Om0': Om0}
    """
    gc = GCRCatalogs.load_catalog(catalog_name, config_overwrite=dict(md5=False))
    cosmo = dict(H0=gc.cosmology.H0.value, Om0=gc.cosmology.Om0)

    sed_filter_dict = sed_filter_names_from_catalog(gc)

    disk_names = sed_filter_dict['disk']['filter_name']
    bulge_names = sed_filter_dict['bulge']['filter_name']

    try:
        for d_name in disk_names:
            b_name = d_name.replace('disk', 'bulge')
            assert b_name in bulge_names
        assert len(disk_names) == len(bulge_names)
    except AssertionError:
        for ii in range(max(len(bulge_names), len(disk_names))):
            if ii < len(disk_names):
                d_name = disk_names[ii]
            else:
                d_name = None

            if ii < len(bulge_names):
                b_name = bulge_names[ii]
            else:
                b_name = None
            print('%s %s' % (b_name, d_name))

        raise RuntimeError("disk and bulge SED filter names do not "
                           "match in catalog")

    np.testing.assert_array_almost_equal(sed_filter_dict['disk']['wav_min'],
                                         sed_filter_dict['bulge']['wav_min'],
                                         decimal=10)

    np.testing.assert_array_almost_equal(sed_filter_dict['disk']['wav_width'],
                                         sed_filter_dict['bulge']['wav_width'],
                                         decimal=10)

    wav_min_grid = sed_filter_dict['disk']['wav_min']
    wav_width_grid = sed_filter_dict['disk']['wav_width']
    wav_max = max((wav0+width
                  for wav0, width in zip(wav_min_grid, wav_width_grid)))
    wav_grid = np.arange(wav_min_grid.min(), wav_max, 0.1)

    bp_name_list = list()
    bp_list = list()
    for wav0, width in zip(wav_min_grid, wav_width_grid):
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
        sed_names.append(sed_file_name)
        sed_mag_list.append(tuple(bandpass_dict.magListForSed(spec)))
        sed_mag_norm.append(spec.calcMag(imsim_bp))

    return (np.array(sed_names), np.array(sed_mag_list),
            np.array(sed_mag_norm), wav_min_grid, wav_width_grid,
            cosmo)


def sed_from_galacticus_mags(galacticus_mags, redshift,
                             catalog_name='protoDC2',
                             wav_min=None, wav_width=None):
    """
    galacticus_mags is a numpy array such that
    galacticus_mags[i][j] is the magnitude of the jth star in the ith bandpass,
    where the bandpasses are ordered in ascending order of minimum wavelength.

    wav_min and wav_width are numpy arrays of the minimum wavelength and
    wavelength grid width (both in nanometers) of the bandpasses
    corresponding to galacticus_mags.  If passed in, this method will
    verify that it is comparing to the correct set of magnitudes

    Will return a numpy array of SED names and a numpy array of magNorms.
    """
    assert catalog_name, '`catalog_name` cannot be None or empty'

    if getattr(sed_from_galacticus_mags, '_catalog_name', None) != catalog_name:

        (sed_names, sed_mag_list, sed_mag_norm,
         wav_min_cat, wav_width_cat, cosmo) = _get_sed_mags_and_cosmology(catalog_name)

        sed_colors = sed_mag_list[:,1:] - sed_mag_list[:,:-1]
        sed_from_galacticus_mags._catalog_name = catalog_name
        sed_from_galacticus_mags._sed_names = sed_names   # N_sed
        sed_from_galacticus_mags._mag_norm = sed_mag_norm # N_sed
        sed_from_galacticus_mags._sed_mags = sed_mag_list # N_sed by N_mag
        sed_from_galacticus_mags._sed_colors = sed_colors # N_sed by (N_mag - 1)
        sed_from_galacticus_mags._cosmo = CosmologyObject(**cosmo)
        sed_from_galacticus_mags._wav_min = wav_min_cat
        sed_from_galacticus_mags._wav_width = wav_width_cat

    if wav_min is not None:
        np.testing.assert_array_almost_equal(wav_min,
                                             sed_from_galacticus_mags._wav_min,
                                             decimal=10)
        np.testing.assert_array_almost_equal(wav_width,
                                             sed_from_galacticus_mags._wav_width,
                                             decimal=10)

    galacticus_mags_t = np.asarray(galacticus_mags).T # N_star by N_mag
    assert galacticus_mags_t.shape == (len(redshift), sed_from_galacticus_mags._sed_mags.shape[1])

    def _find_closest_sed(colors_this):
        return np.argmin(np.sum((sed_from_galacticus_mags._sed_colors - colors_this)**2, axis=1))

    galacticus_colors = galacticus_mags_t[:,1:] - galacticus_mags_t[:,:-1] # N_star by (N_mag - 1)
    sed_idx = np.fromiter(
        (_find_closest_sed(colors_this) for colors_this in galacticus_colors),
        np.int,
        len(galacticus_colors),
    ) # N_star

    distance_modulus = sed_from_galacticus_mags._cosmo.distanceModulus(redshift=redshift)
    output_names = sed_from_galacticus_mags._sed_names[sed_idx]
    d_mag = (galacticus_mags_t - sed_from_galacticus_mags._sed_mags[sed_idx]).mean(axis=1)
    output_mag_norm = sed_from_galacticus_mags._mag_norm[sed_idx] + d_mag + distance_modulus

    return output_names, output_mag_norm
