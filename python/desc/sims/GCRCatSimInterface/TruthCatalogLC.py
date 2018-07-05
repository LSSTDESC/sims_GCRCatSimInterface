import tempfile
import os
import numpy as np
import shutil
import sqlite3
import json
import time
from lsst.sims.catalogs.db import DBObject
from lsst.sims.catUtils.mixins import ExtraGalacticVariabilityModels
from lsst.sims.photUtils import BandpassDict, Sed, getImsimFluxNorm
from lsst.sims.catUtils.supernovae import SNObject
from desc.twinkles import TimeDelayVariability
from . import write_sprinkled_truth_db
from . import get_pointing_htmid

__all__ = ["write_sprinkled_lc"]


class SneSimulator(object):

    def __init__(self, bp_dict):
        """
        bp_dict is a BandpassDict
        """
        self._bp_dict = bp_dict

    def calculate_sn_magnitudes(self, sn_truth_params, mjd_arr, filter_arr):
        """
        sn_truth_params is a numpy array of json-ized
        dicts defining the state of an SNObject

        mjd_arr is a numpy array of the TAI times of
        observations

        filter_arr is a numpy array of ints indicating
        the filter being observed at each time
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
                mm = sn_obj.catsimBandMag(self._bp_dict[filter_name], mjd)
                mags[i_obj][i_time] = mm

        return mags


class AgnSimulator(ExtraGalacticVariabilityModels, TimeDelayVariability):

    def __init__(self, redshift_arr):
        self._redshift_arr = redshift_arr

    def column_by_name(self, key):
        if key!= 'redshift':
            raise RuntimeError("Do not know how to get column %s" % key)
        return self._redshift_arr


def create_sprinkled_sql_file(sql_name):
    with sqlite3.connect(sql_name) as conn:
        cursor = conn.cursor()
        cmd = '''CREATE TABLE sprinkled_agn '''
        cmd += '''(uniqueId int,
                   obshistid int, mag float)'''
        cursor.execute(cmd)

        cmd = '''CREATE TABLE sprinkled_sne '''
        cmd += '''(uniqueId int,
                   obshistid int, mag float)'''
        cursor.execute(cmd)

        cmd = '''CREATE TABLE obs_metadata '''
        cmd += '''(obshistid int, mjd float, filter int)'''
        cursor.execute(cmd)

        cmd = '''CREATE TABLE sprinkled_objects '''
        cmd += '''(uniqueId int, galaxy_id int,
                   ra float, dec float)'''
        cursor.execute(cmd)

        conn.commit()

def write_sprinkled_lc(out_file_name, total_obs_md,
                       pointing_dir, opsim_db_name,
                       field_ra=55.064, field_dec=-29.783,
                       agn_db=None, yaml_file='proto-dc2_v4.6.1',
                       ra_colname='descDitheredRA',
                       dec_colname='descDitheredDec'):
    t0_master = time.time()

    bp_dict = BandpassDict.loadTotalBandpassesFromFiles()
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

    with sqlite3.connect(out_file_name) as conn:
        cursor = conn.cursor()
        values = ((obs, mjd_dict[obs], filter_dict[obs])
                  for obs in mjd_dict)
        cursor.executemany('''INSERT INTO obs_metadata VALUES (?,?,?)''', values)
        cursor.execute('''CREATE INDEX obs_filter
                       ON obs_metadata (obshistid, filter)''')
        conn.commit()

    print('\ngot htmid_dict -- %d in %e seconds' % (len(htmid_dict), t_htmid_dict))

    sql_dir = tempfile.mkdtemp(dir=os.environ['SCRATCH'],
                               prefix='sprinkled_sql_')

    (file_name,
     table_list) = write_sprinkled_truth_db(total_obs_md,
                                            field_ra=field_ra,
                                            field_dec=field_dec,
                                            agn_db=agn_db,
                                            yaml_file=yaml_file,
                                            out_dir=sql_dir)

    print('\nwrote db %s' % file_name)

    db = DBObject(file_name, driver='sqlite')

    query = 'SELECT DISTINCT htmid FROM zpoint WHERE is_agn=1'
    dtype = np.dtype([('htmid', int)])

    results = db.execute_arbitrary(query, dtype=dtype)

    object_htmid = results['htmid']

    agn_dtype = np.dtype([('uniqueId', int), ('galaxy_id', int),
                          ('ra', float), ('dec', float),
                          ('redshift', float), ('sed', str, 500),
                          ('magnorm', float), ('varParamStr', str, 500)])

    agn_base_query = 'SELECT uniqueId, galaxy_id, '
    agn_base_query += 'raJ2000, decJ2000, '
    agn_base_query += 'redshift, sedFilepath, '
    agn_base_query += 'magNorm, varParamStr '
    agn_base_query += 'FROM zpoint WHERE is_agn=1 '

    sn_dtype = np.dtype([('uniqueId', int), ('galaxy_id', int),
                         ('ra', float), ('dec', float),
                         ('redshift', float), ('sn_truth_params', str, 500)])

    sn_base_query = 'SELECT uniqueId, galaxy_id, '
    sn_base_query += 'raJ2000, decJ2000, '
    sn_base_query += 'redshift, sn_truth_params '
    sn_base_query += 'FROM zpoint WHERE is_sn=1 '

    filter_to_int = {'u':0, 'g':1, 'r':2, 'i':3, 'z':4, 'y':5}

    n_floats = 0
    with sqlite3.connect(out_file_name) as conn:
        cursor = conn.cursor()
        for htmid_dex, htmid in enumerate(object_htmid):
            mjd_arr = []
            obs_arr = []
            filter_arr = []
            t_start = time.time()
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
            duration = time.time()-t_start
            print('made mjd_arr in %e seconds' % duration)

            agn_query = agn_base_query + 'AND htmid=%d' % htmid

            agn_iter = db.get_arbitrary_chunk_iterator(agn_query,
                                                       dtype=agn_dtype,
                                                       chunk_size=10000)

            for i_chunk, agn_results in enumerate(agn_iter):
                values = ((int(agn_results['uniqueId'][i_obj]),
                           int(agn_results['galaxy_id'][i_obj]),
                           np.degrees(agn_results['ra'][i_obj]),
                           np.degrees(agn_results['dec'][i_obj]))
                          for i_obj in range(len(agn_results)))

                cursor.executemany('''INSERT INTO sprinkled_objects VALUES
                                      (?,?,?,?)''', values)

                agn_simulator = AgnSimulator(agn_results['redshift'])

                quiescent_mag = np.zeros((len(agn_results), 6), dtype=float)
                for i_obj, (sed_name, zz, mm) in enumerate(zip(agn_results['sed'],
                                                               agn_results['redshift'],
                                                               agn_results['magnorm'])):
                    spec = Sed()
                    spec.readSED_flambda(os.path.join(sed_dir, sed_name))
                    fnorm = getImsimFluxNorm(spec, mm)
                    spec.multiplyFluxNorm(fnorm)
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
                    cursor.executemany('''INSERT INTO sprinkled_agn VALUES
                                       (?,?,?)''', values)

                conn.commit()
                n_floats += len(dmag.flatten())

            sn_query = sn_base_query + 'AND htmid=%d' % htmid

            sn_iter = db.get_arbitrary_chunk_iterator(sn_query,
                                                      dtype=sn_dtype,
                                                      chunk_size=10000)

            for sn_results in sn_iter:
                t0_sne = time.time()
                sn_mags = sn_simulator.calculate_sn_magnitudes(sn_results['sn_truth_params'],
                                                               mjd_arr, filter_arr)
                print('    did %d sne in %e seconds' % (len(sn_results), time.time()-t0_sne))

                for i_time in range(len(mjd_arr)):
                    valid_obj = np.where(np.isfinite(sn_mags[:,i_time]))
                    if len(valid_obj[0]) == 0:
                        continue
                    print('valid sn %d' % len(valid_obj[0]))
                    values = ((int(sn_results['uniqueId'][i_obj]),
                               int(obs_arr[i_time]),
                               sn_mags[i_obj][i_time])
                              for i_obj in valid_obj[0])

                    cursor.executemany('''INSERT INTO sprinkled_sne VALUES (?,?,?)''', values)
                    conn.commit()
                    n_floats += len(valid_obj[0])

    print('n_floats %d' % n_floats)
    print('in %e seconds' % (time.time()-t0_master))

    del db
    for file_name in os.listdir(sql_dir):
        os.unlink(os.path.join(sql_dir, file_name))
    shutil.rmtree(sql_dir)
