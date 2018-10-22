import os
import re
import numpy as np
import GCRCatalogs
import multiprocessing
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

    disk_lsst_names = []
    bulge_lsst_names = []
    for bp in 'ugrizy':
        disk_lsst_names.append('LSST_filters/diskLuminositiesStellar:LSST_%s:observed:dustAtlas' % bp)
        bulge_lsst_names.append('LSST_filters/spheroidLuminositiesStellar:LSST_%s:observed:dustAtlas' % bp)

    return {'disk':{'filter_name': disk_names,
                    'wav_min': disk_wav_min,
                    'wav_width': disk_wav_width,
                    'lsst_fluxes': np.array(disk_lsst_names)},
            'bulge':{'filter_name': bulge_names,
                     'wav_min': bulge_wav_min,
                     'wav_width': bulge_wav_width,
                     'lsst_fluxes': np.array(bulge_lsst_names)}}

def _create_library_one_sed(_galaxy_sed_dir,
                            av_grid, rv_grid,
                            i_start, i_end, bandpass_dict,
                            out_dict):

    n_dust = len(av_grid)
    sed_name_list = os.listdir(_galaxy_sed_dir)
    sed_name_list.sort()

    imsim_bp = Bandpass()
    imsim_bp.imsimBandpass()

    t_start = time.time()
    for i_sed in range(i_start, i_end, 1):
        sed_file_name = sed_name_list[i_sed]
        base_spec = Sed()
        base_spec.readSED_flambda(os.path.join(_galaxy_sed_dir, sed_file_name))
        ax, bx = base_spec.setupCCM_ab()

        mag_norm = base_spec.calcMag(imsim_bp)

        for i_dust, (av, rv) in enumerate(zip(av_grid, rv_grid)):
            spec = Sed(wavelen=base_spec.wavelen, flambda=base_spec.flambda)
            spec.addDust(ax, bx, A_v=av, R_v=rv)
            sed_mag_list = bandpass_dict.magListForSed(spec)
            i_out = i_sed*n_dust + i_dust
            out_dict['sed_names'][i_out] = i_sed
            out_dict['sed_magnorms'][i_out] = mag_norm
            out_dict['av'][i_out] = av
            out_dict['rv'][i_out] = rv
            for i_bp, bp in enumerate('ugrizy'):
                out_dict['lsst_%s' %bp][i_out] = sed_mag_list[i_bp]

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

    sed_lsst_mags is a 2D numpy array giving the magnitude of the SEDs in
    the LSST bands.  sed_lsst_mags[:,0] contains the u band magnitude of
    every SED, etc.

    av_out_list is a 1d float array of Av

    rv_out_list is a 1d float array of Rv
    """

    _av_grid = np.arange(0.0, 3.0, 0.1)
    _rv_grid = np.arange(2.0, 4.1, 0.1)

    n_dust = 0
    for av in _av_grid:
        for i_rv, rv in enumerate(_rv_grid):
            if av<0.01 and i_rv>0:
                continue
            n_dust += 1

    av_grid = np.zeros(n_dust, dtype=float)
    rv_grid = np.zeros(n_dust, dtype=float)
    i_dust = 0
    for av in av_grid:
        for i_rv, rv, in enumerate(rv_grid):
            if av<0.01 and i_rv>0:
                continue
            av_grid[i_dust] = av
            rv_grid[i_dust] = rv
            i_dust += 1

    assert i_dust == n_dust

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

    list_of_files = os.listdir(_galaxy_sed_dir)
    n_tot = len(list_of_files)*n_dust
    t_start = time.time()

    p_list = []
    n_proc = 30
    mgr = multiprocessing.Manager()
    out_dict = mgr.dict()

    out_dict['sed_names'] = mgr.Array('l', np.ones(n_tot, dtype=int))
    out_dict['sed_magnorms'] = mgr.Array('f', np.ones(n_tot, dtype=float))
    out_dict['av'] = mgr.Array('f', np.ones(n_tot, dtype=float))
    out_dict['rv'] = mgr.Array('f', np.ones(n_tot, dtype=float))
    for bp in 'ugrizy':
        out_dict['lsst_%s' % bp] = mgr.Array('f', np.ones(n_tot, dtype=float))

    i_stored = 0
    d_start = len(list_of_files)//n_proc
    i_start_list = range(0, len(list_of_files), d_start)
    for i_meta, i_start in enumerate(i_start_list):
        i_end = i_start+d_start
        if i_meta == len(i_start_list)-1:
            i_end = len(list_of_files)

        p = multiprocessing.Process(target=_create_library_one_sed,
                                    args=(_galaxy_sed_dir,
                                          av_grid, rv_grid,
                                          i_start, i_end,
                                          bandpass_dict, out_dict))

        p.start()
        p_list.append(p)

    for p in p_list:
        p.join()

    sed_names = np.array(out_dict['sed_names'])
    sed_mag_norm = np.array(out_dict['sed_magnorms'])
    av_out_list = np.array(out_dict['av'])
    rv_out_list = np.array(out_dict['rv'])

    sed_mag_list = np.zeros((6, n_tot), dtype=float)
    for i_bp, bp in enumerate('ugrizy'):
        sed_mag_list[i_bp] = np.array(out_dict['lsst_%s' % bp])

    return (sed_names, sed_mag_list, sed_mag_norm,
            av_out_list, rv_out_list)


def sed_from_galacticus_mags(galacticus_mags, redshift, H0, Om0,
                             wav_min, wav_width, obs_lsst_mags):
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

    ob_lsst_mags is a numpy array of observer frame LSST magnitudes.
    obs_lsst_mags[0] will contain the u band magnitudes of every object.

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

        assert rv_grid.min()>0.0
        assert len(np.where(np.logical_not(np.isfinite(rv_grid)))[0])==0

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

    t_start = time.time()
    (sed_dist,
     sed_idx) = sed_from_galacticus_mags._color_tree.query(galacticus_colors, k=1)

    # cKDTree returns an invalid index (==len(tree_data)) in cases
    # where the distance is not finite
    sed_idx = np.where(sed_idx<len(sed_from_galacticus_mags._sed_names),
                       sed_idx, 0)

    distance_modulus = sed_from_galacticus_mags._cosmo.distanceModulus(redshift=redshift)

    output_names = sed_from_galacticus_mags._sed_names[sed_idx]

    (tot_bp_dict,
     lsst_bp_dict) = BandpassDict.loadBandpassesFromFiles()

    output_mag_norm = np.zeros((6, len(output_names)), dtype=float)
    base_norm = sed_from_galacticus_mags._mag_norm[sed_idx]
    assert len(np.where(np.logical_not(np.isfinite(base_norm)))[0])==0
    ccm_w = None
    av_arr = sed_from_galacticus_mags._av_grid[sed_idx]
    rv_arr = sed_from_galacticus_mags._rv_grid[sed_idx]
    assert rv_arr.min()>0.0
    assert len(np.where(np.logical_not(np.isfinite(rv_arr)))[0])==0
    for i_bp in range(6):
        output_mag_norm[i_bp,:] = base_norm + distance_modulus

    sed_dir = getPackageDir('sims_sed_library')

    for i_obj in range(len(output_names)):
        spec = Sed()
        spec.readSED_flambda(os.path.join(sed_dir, output_names[i_obj]))
        if ccm_w is None or not np.array_equal(spec.wavelen, ccm_w):
            ccm_w = np.copy(spec.wavelen)
            ax, bx = spec.setupCCM_ab()
        spec.addDust(ax, bx, A_v=av_arr[i_obj], R_v=rv_arr[i_obj])
        spec.redshiftSED(redshift[i_obj], dimming=True)
        lsst_mags = lsst_bp_dict.magListForSed(spec)
        d_mag = obs_lsst_mags[:,i_obj] - lsst_mags
        output_mag_norm[:,i_obj] += d_mag

    return (output_names, output_mag_norm, av_arr, rv_arr)
