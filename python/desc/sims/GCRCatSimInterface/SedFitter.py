import os
import re
import numpy as np
import GCRCatalogs
from lsst.utils import getPackageDir
from lsst.sims.photUtils import BandpassDict, Bandpass, Sed, CosmologyObject

__all__ = ["sed_from_galacticus_mags"]

_gcr_sed_re = re.compile(r'sed_(\d+)_(\d+)_(bulge|disk)$')
_galaxy_sed_dir = os.path.join(getPackageDir('sims_sed_library'), 'galaxySED')

def _get_sed_mags_and_cosmology(catalog_name):
    """
    Returns 3 numpy arrays: sed_names, sed_mag_list, sed_mag_norm
    and a dictionary for cosmology
    sed_names is 1d str array, with length = number of SED files in the library
    sed_mag_list is MxN float array, with M = number of SED files in the library,
        and N = number of top hat filters in the catalog
    sed_mag_norm is 1d float array, with length = number of SED files in the library
    cosmology = {'H0': H0, 'Om0': Om0}
    """
    gc = GCRCatalogs.load_catalog(catalog_name, config_overwrite=dict(md5=False))
    cosmo = dict(H0=gc.cosmology.H0.value, Om0=gc.cosmology.Om0)

    bp_params_raw = {'disk': set(), 'bulge': set()}
    for q in gc.list_all_quantities():
        m = _gcr_sed_re.match(q)
        if m:
            wav0, width, tag = m.groups()
            bp_params_raw[tag].add((int(wav0), int(width)))
    assert bp_params_raw['disk'] == bp_params_raw['bulge'], 'SEDs for disk and bulge do not match'
    assert bp_params_raw['disk'], 'No SED found'

    bp_params_sorted = sorted(bp_params_raw['disk'], key=lambda p: p[0])
    wav_min = bp_params_sorted[0][0] * 0.1
    wav_max = max((0.1*(wav0+width) for wav0, width in bp_params_sorted))
    wav_grid = np.arange(wav_min, wav_max, 0.1)

    bp_name_list = list()
    bp_list = list()
    for wav0, width in bp_params_sorted:
        sb_grid = ((wav_grid >= wav0*0.1) & (wav_grid <= (wav0+width)*0.1)).astype(float)
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

    return np.array(sed_names), np.array(sed_mag_list), np.array(sed_mag_norm), cosmo


def sed_from_galacticus_mags(galacticus_mags, redshift, catalog_name='protoDC2'):
    """
    galacticus_mags is a numpy array such that
    galacticus_mags[i][j] is the magnitude of the jth star in the ith bandpass,
    where the bandpasses are ordered in ascending order of minimum wavelength.

    Will return a numpy array of SED names and a numpy array of magNorms.
    """
    assert catalog_name, '`catalog_name` cannot be None or empty'

    if getattr(sed_from_galacticus_mags, '_catalog_name', None) != catalog_name:
        sed_names, sed_mag_list, sed_mag_norm, cosmo = _get_sed_mags_and_cosmology(catalog_name)
        sed_colors = sed_mag_list[:,1:] - sed_mag_list[:,:-1]
        sed_from_galacticus_mags._catalog_name = catalog_name
        sed_from_galacticus_mags._sed_names = sed_names   # N_sed
        sed_from_galacticus_mags._mag_norm = sed_mag_norm # N_sed
        sed_from_galacticus_mags._sed_mags = sed_mag_list # N_sed by N_mag
        sed_from_galacticus_mags._sed_colors = sed_colors # N_sed by (N_mag - 1)
        sed_from_galacticus_mags._cosmo = CosmologyObject(**cosmo)

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
