import os
import healpy as hp
import numpy as np
import multiprocessing as mp
import sqlite3

import time

from lsst.utils import getPackageDir
from lsst.sims.utils import ObservationMetaData
from lsst.sims.utils import defaultSpecMap
from lsst.sims.photUtils import BandpassDict
from lsst.sims.photUtils import Sed, getImsimFluxNorm
from lsst.sims.catUtils.baseCatalogModels import StarObj


__all__ = ["write_stars_to_truth"]


def write_results(conn, cursor, mag_dict, position_dict):

    assert len(mag_dict) == len(position_dict)

    row_ct = 0
    for k in mag_dict.keys():
        mm = mag_dict[k]
        pp = position_dict[k]
        row_ct += len(pp['ra'])
        if len(mm) != len(pp['ra']):
            raise RuntimeError('%d mm %d pp' % (len(mm), len(pp['ra'])))

        values = ((int(pp['healpix'][i_obj]),
                   int(pp['id'][i_obj]), True, False, False,
                   pp['ra'][i_obj], pp['dec'][i_obj], 0.0,
                   mm[i_obj][0], mm[i_obj][1], mm[i_obj][2],
                   mm[i_obj][3], mm[i_obj][4], mm[i_obj][5])
                  for i_obj in range(len(pp['ra'])))

        cursor.executemany('''INSERT INTO truth
                           VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)''', values)
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


def write_stars_to_truth(output=None,
                         n_procs=10, n_side=2048, clobber=False):

    if output is None:
        raise RuntimeError("Must specify output database")

    if os.path.isfile(output):
        if clobber:
            os.unlink(output)

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

    with sqlite3.connect(output) as out_conn:
        out_cursor = out_conn.cursor()

        data_iter = db.query_columns(colnames=['simobjid', 'sedFilename',
                                               'magNorm', 'ra', 'decl'],
                                     obs_metadata=obs)

        for star_chunk in data_iter:

            proc = mp.Process(target=calculate_mags,
                              args=(star_chunk['sedFilename'],
                                    star_chunk['magNorm'],
                                    mag_dict))
            proc.start()
            p_list.append(proc)

            hp_arr = hp.ang2pix(n_side,
                                np.radians(star_chunk['ra']),
                                np.radians(star_chunk['decl']),
                                lonlat=True,
                                nest=True)

            position_dict[proc.pid] = {}
            position_dict[proc.pid]['healpix'] = hp_arr
            position_dict[proc.pid]['id'] = star_chunk['simobjid']
            position_dict[proc.pid]['ra'] = star_chunk['ra']
            position_dict[proc.pid]['dec'] = star_chunk['decl']

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
