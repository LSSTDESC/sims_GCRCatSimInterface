import os
import numpy as np
import multiprocessing as mp
import sqlite3
import argparse

import time

from lsst.utils import getPackageDir
from lsst.sims.utils import ObservationMetaData
from lsst.sims.utils import defaultSpecMap
from lsst.sims.photUtils import BandpassDict
from lsst.sims.photUtils import Sed, getImsimFluxNorm
from lsst.sims.catUtils.baseCatalogModels import StarObj

def write_results(conn, cursor, mag_dict, position_dict):

    assert len(mag_dict) == len(position_dict)

    row_ct = 0
    for k in mag_dict.keys():
        mm_arr = mag_dict[k]
        pp_arr = position_dict[k]
        row_ct += len(pp_arr)
        if len(mm_arr) != len(pp_arr):
            raise RuntimeError('%d mm %d pp' % (len(mm_arr), len(pp_arr)))

        values = ((int(pp[0]), pp[1], pp[2],
                   mm[0], mm[1], mm[2], mm[3], mm[4], mm[5])
                  for pp, mm in zip(pp_arr, mm_arr))

        cursor.executemany('''INSERT INTO stars
                           VALUES (?,?,?,?,?,?,?,?,?)''', values)
        conn.commit()

    return row_ct


def calculate_mags(sed_name, mag_norm, out_dict):
    """
    inputs are numpy arrays 
    """
    i_process = mp.current_process().pid

    bp_dir = getPackageDir('throughputs')
    bp_dir = os.path.join(bp_dir, 'imsim', 'goal')
    bp_dict =  BandpassDict.loadTotalBandpassesFromFiles(bandpassDir=bp_dir)

    sed_dir = getPackageDir('sims_sed_library')

    out_mags = np.zeros((len(sed_name), 6), dtype=float)
    for i_star, (s_name, m_norm) in enumerate(zip(sed_name, mag_norm)):
        spec = Sed()
        spec.readSED_flambda(os.path.join(sed_dir, defaultSpecMap[s_name]))
        fnorm = getImsimFluxNorm(spec, m_norm)
        spec.multiplyFluxNorm(fnorm)
        mags = bp_dict.magListForSed(spec)
        out_mags[i_star] = mags

    out_dict[i_process] = out_mags


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('--output', type=str, default=None,
                        help='sqlite database to write')

    parser.add_argument('--n_procs', type=int, default=10,
                        help='number of processes; default=10')

    parser.add_argument('--clobber', dest='clobber', action='store_true')
    parser.set_defaults(clobber=False)

    args = parser.parse_args()

    if args.output is None:
        raise RuntimeError("Must specify output database")

    if os.path.isfile(args.output):
        if args.clobber:
            os.unlink(args.output)

    obs = ObservationMetaData(pointingRA=55.064,
                              pointingDec=-29.783,
                              mjd=59580.0,
                              boundType='circle',
                              boundLength=4.0)

    db = StarObj(database='LSSTCATSIM',
                 host='fatboy.phys.washington.edu',
                 port=1433,
                 driver='mssql+pymssql')

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
        creation_cmd = 'CREATE TABLE stars '
        creation_cmd += '(uniqueId int, ra float, dec float, '
        creation_cmd += 'u float, g float, r float, i float, z float, y float)'
        out_cursor.execute(creation_cmd)
        out_conn.commit()

        data_iter = db.query_columns(colnames=['simobjid', 'sedFilename',
                                               'magNorm', 'ra', 'decl']
                                     obs_metadata=obs)

        for star_chunk in data_iter:

            proc = mp.Process(target=calculate_mags,
                              args=(star_chunk['sedFilename'],
                                    star_chunk['magNorm'],
                                    mag_dict))
            proc.start()
            p_list.append(proc)

            position_dict[proc.pid] = [(star['simobjid'],
                                        star['ra'],
                                        star['decl'])
                                       for star in star_chunk]

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

        if len(p_list) > 0:
            for p in p_list:
                p.join()
            write_results(out_conn, out_cursor,
                          mag_dict, position_dict)

    out_cursor.execute('CREATE INDEX unqid ON stars (uniqueId)')
