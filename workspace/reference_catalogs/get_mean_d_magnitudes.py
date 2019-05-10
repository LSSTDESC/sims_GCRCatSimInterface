import os
import numpy as np
import sqlite3
import json
import multiprocessing
import time

from lsst.sims.utils import radiansFromArcsec
from lsst.sims.utils import arcsecFromRadians
from lsst.sims.photUtils import PhotometricParameters
from lsst.sims.photUtils import BandpassDict
from lsst.sims.catUtils.mixins import MLTflaringMixin
from lsst.sims.catUtils.mixins import ParametrizedLightCurveMixin
from lsst.sims.catUtils.mixins import StellarVariabilityModels
from lsst.sims.catUtils.mixins import create_variability_cache

class VariabilitySimulator(MLTflaringMixin,
                           ParametrizedLightCurveMixin,
                           StellarVariabilityModels):

    def __init__(self, ebv, parallax, uu, gg, rr, ii, zz, yy):
        """
        parallax in radians
        """
        self.photParams = PhotometricParameters(nexp=1, exptime=30)
        self.lsstBandpassDict = BandpassDict.loadTotalBandpassesFromFiles()

        self._ebv = ebv
        px_arcsec = arcsecFromRadians(parallax)
        print('px %e %e' % (1000.0*px_arcsec.min(),1000.0*px_arcsec.max()))
        self._parallax = parallax
        self._q_mag = {}
        self._q_mag['u'] = uu
        self._q_mag['g'] = gg
        self._q_mag['r'] = rr
        self._q_mag['i'] = ii
        self._q_mag['z'] = zz
        self._q_mag['y'] = yy
        self._actually_calculated_columns = list(['lsst_%s' % bp
                                                  for bp in 'ugrizy'])

    def column_by_name(self, name):
        if name == 'ebv':
            return self._ebv
        elif name == 'parallax':
            return self._parallax
        elif name.startswith('quiescent_lsst'):
            bp = name[-1]
            return self._q_mag[bp]
        raise RuntimeError("Cannot return column %s" % name)

def process_results(results, mjd_grid, variability_cache, out_name, lock):

        results = np.array(results)
        simulator = VariabilitySimulator(results[:,1].astype(float),
                                         results[:,2].astype(float),
                                         results[:,3].astype(float),
                                         results[:,4].astype(float),
                                         results[:,5].astype(float),
                                         results[:,6].astype(float),
                                         results[:,7].astype(float),
                                         results[:,8].astype(float))

        t_start = time.time()
        dmag = simulator.applyVariability(results[:,9], expmjd=mjd_grid,
                                          variability_cache=variability_cache)

        d_mean = np.zeros(dmag.shape[:2], dtype=float)
        d_rms = np.zeros(dmag.shape[:2], dtype=float)
        for ii in range(6):
            for jj in range(len(results)):
                d_mean[ii,jj] = np.mean(dmag[ii,jj,:])
                d_rms[ii,jj] = np.sqrt(np.mean(dmag[ii,jj,:]**2))

        lock.acquire()
        with open(out_name, 'a') as out_file:
            for ii in range(len(results)):
                out_file.write('%s ' % results[ii,0])
                for jj in range(6):
                    out_file.write('%e %e ' % (d_mean[jj,ii], d_rms[jj,ii]))
                out_file.write('\n')

        duration = time.time()-t_start
        print(os.getpid(),dmag.shape,duration,
              '%.2e %.2e %.2e' %
               (np.abs(dmag).max(),np.abs(d_mean).max(),d_rms.max()))

        lock.release()

if __name__ == "__main__":

    db_name = '/astro/store/pogo4/danielsf/dc2_stellar_db.db'
    assert os.path.isfile(db_name)

    out_name = '/astro/store/pogo4/danielsf/dc2_stellar_dmag_stats.txt'
    if os.path.exists(out_name):
        os.unlink(out_name)

    rng = np.random.RandomState(5612)
    mjd_grid = np.sort(59580.0+rng.random_sample(2000)*3652.5)

    ra_min=47.91
    ra_max=75.80
    dec_min=-45.32
    dec_max=-26.24

    variability_cache = create_variability_cache()
    plc = ParametrizedLightCurveMixin()
    plc.load_parametrized_light_curves(variability_cache=variability_cache)

    query = "SELECT "
    query += "simobjid, ebv, parallax, "
    query += "umag, gmag, rmag, imag, zmag, ymag, "
    query += "varParamStr "
    query += "FROM stars "
    query += "WHERE ra>=%e " % ra_min
    query += "AND ra<=%e " % ra_max
    query += "AND decl>=%e " % dec_min
    query += "AND decl<=%e " % dec_max

    chunk_size = 10000
    lock = multiprocessing.Lock()
    p_list = []
    n_procs = 15
    with sqlite3.connect(db_name) as connection:
        cursor = connection.cursor()
        query = cursor.execute(query)
        results = query.fetchmany(chunk_size)
        while len(results)>0:
            p = multiprocessing.Process(target=process_results,
                                        args=(results, mjd_grid,
                                              variability_cache,
                                              out_name, lock))

            p.start()
            p_list.append(p)
            if len(p_list)>=n_procs:
                for p in p_list:
                    p.join()
                p_list = []
            results = query.fetchmany(chunk_size)

        for p in p_list:
            p.join()
