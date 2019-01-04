import os
import numpy as np
import healpy as hp
import multiprocessing as mp
import sqlite3

import time

from lsst.utils import getPackageDir
from lsst.sims.photUtils import BandpassDict
from lsst.sims.photUtils import Sed, getImsimFluxNorm
from lsst.sims.utils import defaultSpecMap


__all__ = ["write_galaxies_to_truth"]

_galaxy_query = '''SELECT b.sedFile, b.magNorm,
                   d.sedFile, d.magNorm,
                   a.sedFilepath, a.magNorm,
                   b.redshift, b.galaxy_id,
                   b.raJ2000, b.decJ2000,
                   b.is_sprinkled, a.is_agn,
                   b.shear1, b.shear2, b.kappa
                   FROM bulge as b
                   LEFT JOIN disk as d ON b.galaxy_id=d.galaxy_id
                   LEFT JOIN zpoint as a ON b.galaxy_id=a.galaxy_id
                   WHERE a.is_agn=1 OR a.galaxy_id IS NULL
                   UNION ALL
                   SELECT b.sedFile, b.magNorm,
                   d.sedFile, d.magNorm,
                   a.sedFilepath, a.magNorm,
                   d.redshift, d.galaxy_id,
                   d.raJ2000, d.decJ2000,
                   d.is_sprinkled, a.is_agn,
                   d.shear1, d.shear2, d.kappa
                   FROM disk as d
                   LEFT JOIN bulge as b ON d.galaxy_id=b.galaxy_id
                   LEFT JOIN zpoint as a on d.galaxy_id=a.galaxy_id
                   WHERE b.galaxy_id IS NULL
                   AND (a.is_agn=1 OR a.galaxy_id IS NULL)'''

_col_name_to_int = {}
_col_name_to_int['ra'] = 8
_col_name_to_int['dec'] = 9
_col_name_to_int['redshift'] = 6
_col_name_to_int['galaxy_id'] = 7


def _fluxes(sed_name, mag_norm, redshift):
    """
    Find the fluxes for a galaxy component

    Parameters
    ----------
    sed_name is an SED file name

    mag_norm is a float

    redshift is a float

    Returns
    -------
    array of fluxes in ugrizy order
    """
    if not hasattr(_fluxes, '_bp_dict'):
        bp_dir = getPackageDir('throughputs')
        bp_dir = os.path.join(bp_dir, 'imsim', 'goal')
        _fluxes._bp_dict =  BandpassDict.loadTotalBandpassesFromFiles(bandpassDir=bp_dir)

        _fluxes._sed_dir = getPackageDir('sims_sed_library')

    spec = Sed()
    full_sed_name = os.path.join(_fluxes._sed_dir, sed_name)
    if not os.path.isfile(full_sed_name):
        full_sed_name = os.path.join(_fluxes._sed_dir, defaultSpecMap[sed_name])
    spec.readSED_flambda(full_sed_name)
    fnorm = getImsimFluxNorm(spec, mag_norm)
    spec.multiplyFluxNorm(fnorm)
    spec.redshiftSED(redshift, dimming=True)
    return _fluxes._bp_dict.fluxListForSed(spec)


def write_results(conn, cursor, mag_dict, position_dict):
    """
    Write galaxy truth results to the truth table

    Parameters
    ----------
    conn is a sqlite3 connection to the database

    cursor is a sqlite3.conneciton.cursor() object

    mag_dict is a dict of mags.  It is keyed on the pid of the
    Process used to process a chunk of magnitudes.  Each value
    is a 2-D numpy array of shape (n_obj, n_bandpasses).  It is
    produced by calculate_magnitudes.

    position_dict is a dict keyed on pid of the Process used to
    process a chunk of stars.  The values are also dicts, these
    keyed on 'healpix', 'ra', 'dec', 'galaxy_id', 'redshift',
    'has_agn', 'is_sprinkled', with the values being
    arrays of those quantities for the corresponding chunk of
    stars.

    Returns
    -------
    None

    Just writes to the database
    """
    assert len(mag_dict) == len(position_dict)

    row_ct = 0
    for k in mag_dict.keys():
        mm = mag_dict[k]
        pp = position_dict[k]
        row_ct += len(pp['ra'])
        assert len(mm) == len(pp['ra'])

        values = ((int(pp['healpix'][i_obj]),
                   int(pp['galaxy_id'][i_obj]),
                   0,
                   int(pp['has_agn'][i_obj]),
                   int(pp['is_sprinkled'][i_obj]),
                   pp['ra'][i_obj], pp['dec'][i_obj],
                   pp['redshift'][i_obj],
                   mm[i_obj][0], mm[i_obj][1], mm[i_obj][2],
                   mm[i_obj][3], mm[i_obj][4], mm[i_obj][5])
                  for i_obj in range(len(pp['ra']))
                  if not np.isnan(mm[i_obj][0]))

        cursor.executemany('''INSERT INTO truth
                           VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)''', values)
        conn.commit()

    return row_ct


def calculate_mags(galaxy_list, out_dict):
    """
    Calculate the total (bulge+disk+agn) magnitudes of galaxies

    Parameters
    ----------
    galaxy_list is a list of tuples.  The tuples are of the form
    (bulge.sed_name, bulge.magNorm, disk.sed_name, disk.magNorm
     agn.sed_name, agn.magNorm)

    out_dict is a Multiprocessing.Manager.dict object that will
    store the results of this calculation
    """
    global _col_name_to_int

    i_process = mp.current_process().pid

    bulge_fluxes = np.zeros((len(galaxy_list), 6), dtype=float)
    disk_fluxes = np.zeros((len(galaxy_list), 6), dtype=float)

    magnification = np.array([1.0/((1.0-g[14])**2-g[12]**2-g[13]**2)
                              for g in galaxy_list])

    assert magnification.min() >= 0.0
    assert len(np.where(np.isfinite(magnification))[0]) == len(magnification)

    for i_gal, galaxy in enumerate(galaxy_list):
        if galaxy[0] is not None and galaxy[1] is not None:
            bulge_fluxes[i_gal] = _fluxes(galaxy[0], galaxy[1],
                                          galaxy[_col_name_to_int['redshift']])

        if galaxy[2] is not None and galaxy[3] is not None:
            disk_fluxes[i_gal] = _fluxes(galaxy[2], galaxy[3],
                                         galaxy[_col_name_to_int['redshift']])

    tot_fluxes = bulge_fluxes + disk_fluxes

    for i_filter in range(6):
        tot_fluxes[:,i_filter] *= magnification

    dummy_sed = Sed()
    valid = np.where(tot_fluxes>0.0)
    valid_mags = dummy_sed.magFromFlux(tot_fluxes[valid])
    out_mags = np.NaN*np.ones((len(galaxy_list), 6), dtype=float)
    out_mags[valid] = valid_mags
    out_dict[i_process] = out_mags


def write_galaxies_to_truth(n_side=2048, input_db=None, output=None,
                            n_procs=10, clobber=False):
    """
    Write static galaxy truth to the truth catalog

    Parameters
    ----------
    input_db is the path to the sqlite file containing extragalactic
    parameters as written by write_sprinkled_param_db

    output is the path to the output database

    n_procs is the number of Multiprocessing processes to use when
    calculating magnitudes

    n_side is the nside parameter for calculating healpix locations

    clobber is a boolean.  If True, delete any already existing databases
    with the same file name as output (default=False)

    Returns
    -------
    None

    Just writes to the database
    """
    global _galaxy_query
    global _col_name_to_int

    if input_db is None:
        raise RuntimeError("Must specify input database")

    if output is None:
        raise RuntimeError("Must specify output database")

    if output == input_db:
        raise RuntimeError("output == input_db")

    if os.path.isfile(output):
        if clobber:
            os.unlink(output)

    if not os.path.isfile(input_db):
        raise RuntimeError("%s does not exist" % input_db)

    chunk_size = 10000
    p_list = []

    mgr = mp.Manager()
    mag_dict = mgr.dict()
    position_dict = {}

    t_start = time.time()
    row_ct = 0
    iteration = 0

    is_agn_converter = {None:0, 1:1, 0:0}

    with sqlite3.connect(output) as out_conn:
        out_cursor = out_conn.cursor()

        with sqlite3.connect(input_db) as in_conn:
            in_cursor = in_conn.cursor()
            query = in_cursor.execute(_galaxy_query)

            while True:
                results = query.fetchmany(chunk_size)
                if len(results) == 0:
                    break

                proc = mp.Process(target=calculate_mags,
                                  args=(results, mag_dict))
                proc.start()
                p_list.append(proc)

                ra_arr = np.degrees(np.array([r[_col_name_to_int['ra']] for r in results]))
                dec_arr = np.degrees(np.array([r[_col_name_to_int['dec']] for r in results]))
                hp_arr = hp.ang2pix(n_side, ra_arr, dec_arr,
                                    lonlat=True,
                                    nest=True)

                local_dict = {}
                local_dict['healpix'] = hp_arr
                local_dict['ra'] = ra_arr
                local_dict['dec'] = dec_arr
                local_dict['redshift'] = np.array([r[_col_name_to_int['redshift']] for r in results])
                local_dict['galaxy_id'] = np.array([r[_col_name_to_int['galaxy_id']] for r in results])
                local_dict['is_sprinkled'] = [r[10]
                                              for r in results]
                local_dict['has_agn'] = [is_agn_converter[r[11]]
                                         for r in results]

                position_dict[proc.pid] = local_dict

                if len(p_list) >= n_procs:
                    for p in p_list:
                        p.join()
                    row_ct += write_results(out_conn, out_cursor,
                                            mag_dict, position_dict)
                    p_list = []
                    position_dict = {}
                    mag_dict = mgr.dict()
                    iteration += 1
                    duration = (time.time()-t_start)/3600.0
                    predicted = 1.0e7*duration/row_ct
                    print('output %d in %.2e hrs; 10 million in %.2e' %
                          (row_ct, duration, predicted))

            if len(p_list) > 0:
                for p in p_list:
                    p.join()
                write_results(out_conn, out_cursor,
                              mag_dict, position_dict)
