import os
import numpy as np
import multiprocessing as mp
import sqlite3
import argparse

import time

from lsst.utils import getPackageDir
from lsst.sims.photUtils import BandpassDict
from lsst.sims.photUtils import Sed, getImsimFluxNorm


def write_results(conn, cursor, mag_dict, position_dict):

    assert len(mag_dict) == len(position_dict)

    row_ct = 0
    for k in mag_dict.keys():
        mm_arr = mag_dict[k]
        pp_arr = position_dict[k]
        row_ct += len(pp_arr)
        assert len(mm_arr) == len(pp_arr)
        values = ((pp[0], int(pp[1]), pp[2], pp[3],
                   mm[0], mm[1], mm[2], mm[3], mm[4], mm[5])
                  for pp, mm in zip(pp_arr, mm_arr))

        cursor.executemany('''INSERT INTO static_galaxies
                           VALUES (?,?,?,?,?,?,?,?,?,?)''', values)
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

    bp_dir = getPackageDir('throughputs')
    bp_dir = os.path.join(bp_dir, 'imsim', 'goal')
    bp_dict =  BandpassDict.loadTotalBandpassesFromFiles(bandpassDir=bp_dir)

    sed_dir = getPackageDir('sims_sed_library')

    bulge_fluxes = np.zeros((len(galaxy_list), 6), dtype=float)
    disk_fluxes = np.zeros((len(galaxy_list), 6), dtype=float)

    for i_gal, galaxy in enumerate(galaxy_list):
        if galaxy[0] is not None and galaxy[1] is not None:
            spec = Sed()
            spec.readSED_flambda(os.path.join(sed_dir, galaxy[0]))
            fnorm = getImsimFluxNorm(spec, galaxy[1])
            spec.multiplyFluxNorm(fnorm)
            spec.redshiftSED(galaxy[4], dimming=True)
            bulge_fluxes[i_gal] = bp_dict.fluxListForSed(spec)

        if galaxy[2] is not None and galaxy[3] is not None:
            spec = Sed()
            spec.readSED_flambda(os.path.join(sed_dir, galaxy[2]))
            fnorm = getImsimFluxNorm(spec, galaxy[3])
            spec.multiplyFluxNorm(fnorm)
            spec.redshiftSED(galaxy[4], dimming=True)
            disk_fluxes[i_gal] = bp_dict.fluxListForSed(spec)

    tot_fluxes = bulge_fluxes + disk_fluxes
    dummy_sed = Sed()
    valid = np.where(tot_fluxes>0.0)
    valid_mags = dummy_sed.magFromFlux(tot_fluxes[valid])
    out_mags = np.NaN*np.ones((len(galaxy_list), 6), dtype=float)
    out_mags[valid] = valid_mags
    out_dict[i_process] = out_mags


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('--input', type=str, default=None,
                        help='sqlite database containing parameters')

    parser.add_argument('--output', type=str, default=None,
                        help='sqlite database to write')

    parser.add_argument('--n_procs', type=int, default=10,
                        help='number of processes; default=10')

    parser.add_argument('--clobber', dest='clobber', action='store_true')
    parser.set_defaults(clobber=False)

    args = parser.parse_args()

    if args.input is None:
        raise RuntimeError("Must specify input database")

    if args.output is None:
        raise RuntimeError("Must specify output database")

    if args.output == args.input:
        raise RuntimeError("args.output == args.input")

    if os.path.isfile(args.output):
        if not args.clobber:
            raise RuntimeError("%s already exists" % args.output)
        os.unlink(args.output)

    if not os.path.isfile(args.input):
        raise RuntimeError("%s does not exist" % args.input)

    query = 'SELECT b.sedFile, b.magNorm, '
    query += 'd.sedFile, d.magNorm, '
    query += 'b.redshift, b.galaxy_id, '
    query += 'b.raJ2000, b.decJ2000 '

    query += 'FROM bulge as b '
    query += 'LEFT JOIN disk as d ON d.galaxy_id=b.galaxy_id '

    query += 'UNION ALL'

    query = 'SELECT b.sedFile, b.magNorm, '
    query += 'd.sedFile, d.magNorm, '
    query += 'd.redshift, d.galaxy_id, '
    query += 'd.raJ2000, d.decJ2000 '

    query += 'FROM disk as d '
    query += 'LEFT JOIN bulge as b ON d.galaxy_id=b.galaxy_id '
    query += 'WHERE b.galaxy_id IS NULL'

    chunk_size = 10000
    p_list = []

    mgr = mp.Manager()
    mag_dict = mgr.dict()
    position_dict = {}

    t_start = time.time()
    row_ct = 0
    iteration = 0

    with sqlite3.connect(args.output) as out_conn:
        out_cursor = out_conn.cursor()
        creation_cmd = 'CREATE TABLE static_galaxies '
        creation_cmd += '(redshfit float, galaxy_id int, ra float, dec float, '
        creation_cmd += 'u float, g float, r float, i float, z float, y float)'
        out_cursor.execute(creation_cmd)
        out_conn.commit()

        with sqlite3.connect(args.input) as in_conn:
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

                position_dict[proc.pid] = [(r[4], r[5], r[6], r[7])
                                           for r in results]

                if len(p_list) >= args.n_procs:
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

                    if iteration > 2:
                        exit()

            if len(p_list) > 0:
                for p in p_list:
                    p.join()
                write_results(out_conn, out_cursor,
                              mag_dict, position_dict)

        out_cursor.execute('CREATE INDEX gid ON static_galaxies (galaxy_id)')
