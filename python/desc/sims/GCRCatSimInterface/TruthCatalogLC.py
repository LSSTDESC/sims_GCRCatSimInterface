import os
import numpy as np
import shutil
import sqlite3
import json
import multiprocessing
import time
from lsst.sims.catalogs.db import DBObject
from lsst.sims.catUtils.mixins import ExtraGalacticVariabilityModels
from lsst.sims.photUtils import BandpassDict, Sed, getImsimFluxNorm
from lsst.sims.catUtils.supernovae import SNObject
from desc.twinkles import TimeDelayVariability
from . import get_pointing_htmid

__all__ = ["write_sprinkled_lc"]


class SneSimulator(object):

    def __init__(self, bp_dict):
        """
        bp_dict is a BandpassDict
        """
        self._bp_dict = bp_dict

    def _calculate_batch_of_sn(self, sn_truth_params, mjd_arr, filter_arr,
                               out_dict, tag):
        """
        For processing SNe in batches using multiprocessing

        -----
        sn_truth_params is a numpy array of json-ized
        dicts defining the state of an SNObject

        mjd_arr is a numpy array of the TAI times of
        observations

        filter_arr is a numpy array of ints indicating
        the filter being observed at each time

        out_dict is a multprocessing.Manager().dict() to store
        the results

        tag is an integer denoting which batch of SNe this is.
        """

        int_to_filter = 'ugrizy'

        mags = np.NaN*np.ones((len(sn_truth_params), len(mjd_arr)), dtype=float)

        for i_obj, sn_par in enumerate(sn_truth_params):
            sn_obj = SNObject.fromSNState(json.loads(sn_par))
            for i_time in range(len(mjd_arr)):
                mjd = mjd_arr[i_time]
                filter_name = int_to_filter[filter_arr[i_time]]
                if mjd < sn_obj.mintime() or mjd > sn_obj.maxtime():
                    continue
                bp = self._bp_dict[filter_name]
                sn_sed = sn_obj.SNObjectSED(mjd, bandpass=bp,
                                            applyExtinction=False)
                ff = sn_sed.calcFlux(bp)
                if ff>1.0e-300:
                    mm = sn_sed.magFromFlux(ff)
                    mags[i_obj][i_time] = mm

        out_dict[tag] = mags

    def calculate_sn_magnitudes(self, sn_truth_params, mjd_arr, filter_arr):
        """
        sn_truth_params is a numpy array of json-ized
        dicts defining the state of an SNObject

        mjd_arr is a numpy array of the TAI times of
        observations

        filter_arr is a numpy array of ints indicating
        the filter being observed at each time
        """
        n_threads = 12
        d_sne = max(12, len(sn_truth_params)//(n_threads-1))
        p_list = []
        mgr = multiprocessing.Manager()
        out_dict = mgr.dict()
        for i_start in range(0, len(sn_truth_params), d_sne):
            i_end = i_start+d_sne
            selection = slice(i_start, i_end)
            p = multiprocessing.Process(target=self._calculate_batch_of_sn,
                                        args=(sn_truth_params[selection],
                                              mjd_arr, filter_arr,
                                              out_dict, i_start))

            p.start()
            p_list.append(p)

        for p in p_list:
            p.join()

        mags = np.NaN*np.ones((len(sn_truth_params), len(mjd_arr)), dtype=float)
        for i_start in range(0, len(sn_truth_params), d_sne):
            i_end = i_start + d_sne
            selection = slice(i_start, i_end)
            mags[selection] = out_dict[i_start]

        return mags


class AgnSimulator(ExtraGalacticVariabilityModels, TimeDelayVariability):

    def __init__(self, redshift_arr):
        self._redshift_arr = redshift_arr
        self._agn_threads = 12

    def column_by_name(self, key):
        if key!= 'redshift':
            raise RuntimeError("Do not know how to get column %s" % key)
        return self._redshift_arr


def create_sprinkled_sql_file(sql_name):
    with sqlite3.connect(sql_name) as conn:
        cursor = conn.cursor()
        cmd = '''CREATE TABLE agn_lc '''
        cmd += '''(uniqueId int,
                   obshistid int, mag float)'''
        cursor.execute(cmd)

        cmd = '''CREATE TABLE sne_lc '''
        cmd += '''(uniqueId int,
                   obshistid int, mag float)'''
        cursor.execute(cmd)

        cmd = '''CREATE TABLE obs_metadata '''
        cmd += '''(obshistid int, mjd float, filter int)'''
        cursor.execute(cmd)

        cmd = '''CREATE TABLE variables_and_transients '''
        cmd += '''(uniqueId int, galaxy_id int,
                   ra float, dec float, sprinkled int)'''
        cursor.execute(cmd)

        conn.commit()

def write_sprinkled_lc(out_file_name, total_obs_md,
                       pointing_dir, opsim_db_name,
                       ra_colname='descDitheredRA',
                       dec_colname='descDitheredDec',
                       sql_file_name=None,
                       bp_dict=None):

    """
    Create database of light curves

    Note: this is still under development.  It has not yet been
    used for a production-level truth catalog

    Parameters
    ----------
    out_file_name is the name of the sqlite file to be written

    total_obs_md is an ObservationMetaData covering the whole
    survey area

    pointing_dir contains a series of files that are two columns: obshistid, mjd.
    The files must each have 'visits' in their name.  These specify the pointings
    for which we are assembling data.  See:
        https://github.com/LSSTDESC/DC2_Repo/tree/master/data/Run1.1
    for an example.

    opsim_db_name is the name of the OpSim database to be queried for pointings

    ra_colname is the column used for RA of the pointing (default:
    descDitheredRA)

    dec_colname is the column used for the Dec of the pointing (default:
    descDitheredDec)

    sql_file_name is the name of the parameter database produced by
    write_sprinkled_param_db to be used

    bp_dict is a BandpassDict of the telescope filters to be used

    Returns
    -------
    None

    Writes out a database to out_file_name.  The tables of this database and
    their columns are:

    agn_lc/sn_lc:
        - uniqueId -- an int unique to all objects
        - obshistid -- an int unique to all pointings
        - mag -- the magnitude observed for this object at that pointing

    obs_metadata:
        - obshistid -- an int unique to all pointings
        - mjd -- the date of the pointing
        - filter -- an int corresponding to the telescope filter (0==u, 1==g..)

    variables_and_transients:
        - uniqueId -- an int unique to all objects
        - galaxy_id -- an int indicating the host galaxy
        - ra -- in degrees
        - dec -- in degrees
    """

    t0_master = time.time()

    if not os.path.isfile(sql_file_name):
        raise RuntimeError('%s does not exist' % sql_file_name)


    sn_simulator = SneSimulator(bp_dict)
    sed_dir = os.environ['SIMS_SED_LIBRARY_DIR']

    create_sprinkled_sql_file(out_file_name)

    t_start = time.time()
    (htmid_dict,
     mjd_dict,
     filter_dict) = get_pointing_htmid(pointing_dir, opsim_db_name,
                                       ra_colname=ra_colname,
                                       dec_colname=dec_colname)
    t_htmid_dict = time.time()-t_start

    bp_to_int = {'u':0, 'g':1, 'r':2, 'i':3, 'z':4 'y':5}

    with sqlite3.connect(out_file_name) as conn:
        cursor = conn.cursor()
        values = ((int(obs), mjd_dict[obs], bp_to_int[filter_dict[obs]])
                  for obs in mjd_dict)
        cursor.executemany('''INSERT INTO obs_metadata VALUES (?,?,?)''', values)
        cursor.execute('''CREATE INDEX obs_filter
                       ON obs_metadata (obshistid, filter)''')
        conn.commit()

    print('\ngot htmid_dict -- %d in %e seconds' % (len(htmid_dict), t_htmid_dict))

    db = DBObject(sql_file_name, driver='sqlite')

    query = 'SELECT DISTINCT htmid FROM zpoint WHERE is_agn=1 OR is_sn=1'
    dtype = np.dtype([('htmid', int)])

    results = db.execute_arbitrary(query, dtype=dtype)

    object_htmid = results['htmid']

    agn_dtype = np.dtype([('uniqueId', int), ('galaxy_id', int),
                          ('ra', float), ('dec', float),
                          ('redshift', float), ('sed', str, 500),
                          ('magnorm', float), ('varParamStr', str, 500),
                          ('is_sprinkled', int)])

    agn_base_query = 'SELECT uniqueId, galaxy_id, '
    agn_base_query += 'raJ2000, decJ2000, '
    agn_base_query += 'redshift, sedFilepath, '
    agn_base_query += 'magNorm, varParamStr, is_sprinkled '
    agn_base_query += 'FROM zpoint WHERE is_agn=1 '

    sn_dtype = np.dtype([('uniqueId', int), ('galaxy_id', int),
                         ('ra', float), ('dec', float),
                         ('redshift', float), ('sn_truth_params', str, 500),
                         ('is_sprinkled', int)])

    sn_base_query = 'SELECT uniqueId, galaxy_id, '
    sn_base_query += 'raJ2000, decJ2000, '
    sn_base_query += 'redshift, sn_truth_params, is_sprinkled '
    sn_base_query += 'FROM zpoint WHERE is_sn=1 '

    filter_to_int = {'u':0, 'g':1, 'r':2, 'i':3, 'z':4, 'y':5}

    n_floats = 0
    with sqlite3.connect(out_file_name) as conn:
        cursor = conn.cursor()
        t_before_htmid = time.time()
        for htmid_dex, htmid in enumerate(object_htmid):
            if htmid_dex>0:
                htmid_duration = (time.time()-t_before_htmid)/3600.0
                htmid_prediction = len(object_htmid)*htmid_duration/htmid_dex
                print('%d htmid out of %d in %e hours; predict %e hours remaining' %
                (htmid_dex, len(object_htmid), htmid_duration,htmid_prediction-htmid_duration))
            mjd_arr = []
            obs_arr = []
            filter_arr = []

            for obshistid in htmid_dict:
                is_contained = False
                for bounds in htmid_dict[obshistid]:
                    if htmid<=bounds[1] and htmid>=bounds[0]:
                        is_contained = True
                        break
                if is_contained:
                    mjd_arr.append(mjd_dict[obshistid])
                    obs_arr.append(obshistid)
                    filter_arr.append(filter_to_int[filter_dict[obshistid]])
            if len(mjd_arr) == 0:
                continue
            mjd_arr = np.array(mjd_arr)
            obs_arr = np.array(obs_arr)
            filter_arr = np.array(filter_arr)
            sorted_dex = np.argsort(mjd_arr)
            mjd_arr = mjd_arr[sorted_dex]
            obs_arr = obs_arr[sorted_dex]
            filter_arr = filter_arr[sorted_dex]

            agn_query = agn_base_query + 'AND htmid=%d' % htmid

            agn_iter = db.get_arbitrary_chunk_iterator(agn_query,
                                                       dtype=agn_dtype,
                                                       chunk_size=10000)

            for i_chunk, agn_results in enumerate(agn_iter):
                values = ((int(agn_results['uniqueId'][i_obj]),
                           int(agn_results['galaxy_id'][i_obj]),
                           np.degrees(agn_results['ra'][i_obj]),
                           np.degrees(agn_results['dec'][i_obj]),
                           int(agn_results['is_sprinkled'][i_obj]))
                          for i_obj in range(len(agn_results)))

                cursor.executemany('''INSERT INTO variables_and_transients VALUES
                                      (?,?,?,?,?)''', values)

                agn_simulator = AgnSimulator(agn_results['redshift'])

                quiescent_mag = np.zeros((len(agn_results), 6), dtype=float)
                for i_obj, (sed_name, zz, mm) in enumerate(zip(agn_results['sed'],
                                                               agn_results['redshift'],
                                                               agn_results['magnorm'])):
                    spec = Sed()
                    spec.readSED_flambda(os.path.join(sed_dir, sed_name))
                    fnorm = getImsimFluxNorm(spec, mm)
                    spec.multiplyFluxNorm(fnorm)
                    spec.redshiftSED(zz, dimming=True)
                    mag_list = bp_dict.magListForSed(spec)
                    quiescent_mag[i_obj] = mag_list

                dmag = agn_simulator.applyVariability(agn_results['varParamStr'],
                                                      expmjd=mjd_arr)


                for i_time in range(len(obs_arr)):
                    values = ((int(agn_results['uniqueId'][i_obj]),
                               int(obs_arr[i_time]),
                               quiescent_mag[i_obj][filter_arr[i_time]]+
                               dmag[filter_arr[i_time]][i_obj][i_time])
                              for i_obj in range(len(agn_results)))
                    cursor.executemany('''INSERT INTO agn_lc VALUES
                                       (?,?,?)''', values)

                conn.commit()
                n_floats += len(dmag.flatten())

            sn_query = sn_base_query + 'AND htmid=%d' % htmid

            sn_iter = db.get_arbitrary_chunk_iterator(sn_query,
                                                      dtype=sn_dtype,
                                                      chunk_size=10000)

            for sn_results in sn_iter:
                t0_sne = time.time()

                values = ((int(sn_results['uniqueId'][i_obj]),
                           int(sn_results['galaxy_id'][i_obj]),
                           np.degrees(sn_results['ra'][i_obj]),
                           np.degrees(sn_results['dec'][i_obj]),
                           int(sn_results['is_sprinkled'][i_obj]))
                          for i_obj in range(len(sn_results)))

                cursor.executemany('''INSERT INTO variables_and_transients VALUES
                                      (?,?,?,?,?)''', values)
                conn.commit()

                sn_mags = sn_simulator.calculate_sn_magnitudes(sn_results['sn_truth_params'],
                                                               mjd_arr, filter_arr)
                print('    did %d sne in %e seconds' % (len(sn_results), time.time()-t0_sne))

                for i_time in range(len(mjd_arr)):
                    valid_obj = np.where(np.isfinite(sn_mags[:,i_time]))
                    if len(valid_obj[0]) == 0:
                        continue

                    values = ((int(sn_results['uniqueId'][i_obj]),
                               int(obs_arr[i_time]),
                               sn_mags[i_obj][i_time])
                              for i_obj in valid_obj[0])

                    cursor.executemany('''INSERT INTO sne_lc VALUES (?,?,?)''', values)
                    conn.commit()
                    n_floats += len(valid_obj[0])

        cursor.execute('CREATE INDEX unq_obs_agn ON agn_lc (uniqueId, obshistid)')
        conn.commit()
        cursor.execute('CREATE INDEX unq_obs_sne ON sne_lc (uniqueId, obshistid)')
        conn.commit()

    print('n_floats %d' % n_floats)
    print('in %e seconds' % (time.time()-t0_master))
