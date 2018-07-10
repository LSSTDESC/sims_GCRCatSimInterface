import os
import numpy as np
import healpy as hp
import multiprocessing as mp
import sqlite3

import time

from lsst.utils import getPackageDir
from lsst.sims.photUtils import BandpassDict
from lsst.sims.photUtils import Sed, getImsimFluxNorm


__all__ = ["write_galaxies_to_truth"]


def _fluxes(sed_name, mag_norm, redshift):
    if not hasattr(_fluxes, '_bp_dict'):
        bp_dir = getPackageDir('throughputs')
        bp_dir = os.path.join(bp_dir, 'imsim', 'goal')
        _fluxes._bp_dict =  BandpassDict.loadTotalBandpassesFromFiles(bandpassDir=bp_dir)

        _fluxes._sed_dir = getPackageDir('sims_sed_library')

    spec = Sed()
    spec.readSED_flambda(os.path.join(_fluxes._sed_dir, sed_name))
    fnorm = getImsimFluxNorm(spec, mag_norm)
    spec.multiplyFluxNorm(fnorm)
    spec.redshiftSED(redshift, dimming=True)
    return _fluxes._bp_dict.fluxListForSed(spec)


def write_results(conn, cursor, mag_dict, position_dict):

    assert len(mag_dict) == len(position_dict)

    row_ct = 0
    for k in mag_dict.keys():
        mm = mag_dict[k]
        pp = position_dict[k]
        row_ct += len(pp_arr)
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
                  for i_obj in range(len(pp['ra'])))

        cursor.executemany('''INSERT INTO truth
                           VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)''', values)
        conn.commit()

    return row_ct


def calculate_mags(galaxy_list, out_dict):
    """
    galaxy_list will be a list of tuples as returned by the sqlite3
    query.  The tuples will be:

    query = 'SELECT b.sedFile, b.magNorm, '
    query += 'd.sedFile, d.magNorm, '
    query += 'd.redshift, d.uniqueId, d.galaxy_id, '
    query += 'd.raJ2000, d.decJ2000 '
    """
    i_process = mp.current_process().pid

    bulge_fluxes = np.zeros((len(galaxy_list), 6), dtype=float)
    disk_fluxes = np.zeros((len(galaxy_list), 6), dtype=float)
    agn_fluxes = np.zeros((len(galaxy_list), 6), dtype=float)

    for i_gal, galaxy in enumerate(galaxy_list):
        if galaxy[0] is not None and galaxy[1] is not None:
            bulge_fluxes[i_gal] = _fluxes(galaxy[0], galaxy[1], galaxy[6])

        if galaxy[2] is not None and galaxy[3] is not None:
            disk_fluxes[i_gal] = _fluxes(galaxy[2], galaxy[3], galaxy[6])

        if galaxy[4] is not None and galaxy[5] is not None:
            agn_fluxes[i_gal] = _fluxes(galaxy[4], galaxy[5], galaxy[6])

    tot_fluxes = bulge_fluxes + disk_fluxes + agn_fluxes
    dummy_sed = Sed()
    valid = np.where(tot_fluxes>0.0)
    valid_mags = dummy_sed.magFromFlux(tot_fluxes[valid])
    out_mags = np.NaN*np.ones((len(galaxy_list), 6), dtype=float)
    out_mags[valid] = valid_mags
    out_dict[i_process] = out_mags


def write_galaxies_to_truth(n_side=2048, input=None, output=None
                            n_procs=10, clobber=False):

    if input is None:
        raise RuntimeError("Must specify input database")

    if output is None:
        raise RuntimeError("Must specify output database")

    if output == input:
        raise RuntimeError("output == input")

    if os.path.isfile(output):
        if clobber:
            os.unlink(output)

    if not os.path.isfile(input):
        raise RuntimeError("%s does not exist" % input)

    query = 'SELECT b.sedFile, b.magNorm, '
    query += 'd.sedFile, d.magNorm, '
    query += 'a.sedFilepath, a.magNorm, '
    query += 'b.redshift, b.galaxy_id, '
    query += 'b.raJ2000, b.decJ2000, '
    query += 'b.is_sprinkled, a.is_agn '

    query += 'FROM bulge as b '
    query += 'LEFT JOIN disk as d ON d.galaxy_id=b.galaxy_id '
    query += 'LEFT JOIN zpoint as a ON d.galaxy_id=a.galaxy_id '
    query += 'WHERE a.is_agn=1 OR a.galaxy_id IS NULL '

    query += 'UNION ALL '

    query += 'SELECT b.sedFile, b.magNorm, '
    query += 'd.sedFile, d.magNorm, '
    query += 'a.sedFilepath, a.magNorm, '
    query += 'd.redshift, d.galaxy_id, '
    query += 'd.raJ2000, d.decJ2000, '
    query += 'd.is_sprinkled, a.is_agn '

    query += 'FROM disk as d '
    query += 'LEFT JOIN bulge as b ON b.galaxy_id=d.galaxy_id '
    query += 'LEFT JOIN zpoint as a on b.galaxy_id=a.galaxy_id '
    query += 'WHERE b.galaxy_id IS NULL '
    query += 'AND (a.is_agn=1 OR a.galaxy_id IS NULL)'

    chunk_size = 10000
    p_list = []

    mgr = mp.Manager()
    mag_dict = mgr.dict()
    position_dict = {}

    t_start = time.time()
    row_ct = 0
    iteration = 0

    is_agn_converter = {None: 0, 1:1, 0:0}

    with sqlite3.connect(output) as out_conn:
        out_cursor = out_conn.cursor()

        with sqlite3.connect(input) as in_conn:
            in_cursor = in_conn.cursor()
            query = in_cursor.execute(query)

            while True:
                results = query.fetchmany(chunk_size)
                if len(results) == 0:
                    break

                proc = mp.Process(target=calculate_mags,
                                  args=(results, mag_dict))
                proc.start()
                p_list.append(proc)

                ra_arr = np.array([r[8] for r in results])
                dec_arr = np.array([r[9] for r in results])
                hp_arr = hp.ang2pix(n_side, ra_arr, dec_arr,
                                    lonlat=True,
                                    nest=True)

                local_dict = {}
                local_dict['healpix'] = hp_arr
                local_dict['ra'] = np.degrees(ra_arr)
                local_dict['dec'] = np.degrees(dec_arr)
                local_dict['redshift'] = np.array([r[6] for r in results])
                local_dict['galaxy_id'] = np.array([r[7] for r in results])
                local_dict['is_sprinkled'] = np.array([r[10] for r in results])
                local_dict['has_agn'] = np.array([r[11] for r in results])

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

        #out_cursor.execute('CREATE INDEX gid ON static_galaxies (galaxy_id)')
        #out_cursor.execute('CREATE INDEX g_hp ON static_galaxies (healpix_id)')
