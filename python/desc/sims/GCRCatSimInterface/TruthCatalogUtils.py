import numpy as np
import os
import sqlite3
import time

from lsst.sims.catalogs.decorators import cached, compound
from lsst.sims.catalogs.definitions import InstanceCatalog
from lsst.sims.catalogs.db import DBObject
from lsst.sims.utils import findHtmid, halfSpaceFromRaDec
from lsst.sims.utils import defaultSpecMap
from lsst.sims.photUtils import Sed, getImsimFluxNorm
from lsst.sims.catUtils.baseCatalogModels import StarObj
from . import sprinklerCompound_DC2_truth
from . import TwinklesCompoundInstanceCatalog_DC2
from . import SQLSubCatalogMixin
from . import diskDESCQAObject_protoDC2 as diskDESCQAObject
from . import bulgeDESCQAObject_protoDC2 as bulgeDESCQAObject
from . import agnDESCQAObject_protoDC2 as agnDESCQAObject
from . import PhoSimDESCQA
from . import TwinklesCatalogZPoint_DC2
from .TwinklesClasses import twinkles_spec_map

__all__ = ["write_sprinkled_param_db",
           "get_pointing_htmid",
           "write_star_truth_db"]


_truth_trixel_level = 7  # resolution of trixels in which to separate sources


class _SprinkledTruth(object):

    def get_htmid(self):
        ra = np.degrees(self.column_by_name('raJ2000'))
        dec = np.degrees(self.column_by_name('decJ2000'))
        return findHtmid(ra, dec, max_level=_truth_trixel_level)

    def write_header(self, file_handle):
        InstanceCatalog.write_header(self, file_handle)


class _SersicTruth(_SprinkledTruth):
    column_outputs = ['uniqueId', 'galaxy_id',
                      'raJ2000', 'decJ2000', 'redshift',
                      'magNorm', 'sedFile', 'internalAv',
                      'internalRv', 'galacticAv', 'galacticRv',
                      'is_sprinkled']

    _t_start = None
    _sersic_ct = 0

    @cached
    def get_sedFile(self):
        self.column_by_name('is_sprinkled')
        self.column_by_name('majorAxis')
        self.column_by_name('minorAxis')
        self.column_by_name('positionAngle')
        if self._t_start is None:
            self._t_start = time.time()

        out = self.column_by_name('sedFilepath')
        self._sersic_ct += len(out)

        if self._sersic_ct > 0:
            duration = (time.time()-self._t_start)/3600.0
            prediction = 8.0e6*duration/self._sersic_ct
            print('%d in %.2e hrs; 8 million in %.2e' %
            (self._sersic_ct, duration,prediction))

        return out


class _ZPointTruth(_SprinkledTruth):
    column_outputs = ['uniqueId', 'galaxy_id', 'htmid',
                      'raJ2000', 'decJ2000', 'redshift',
                      'sedFilepath', 'magNorm',
                      'varParamStr', 'sn_truth_params',
                      'sn_t0',
                      'is_sprinkled', 'is_sn', 'is_agn']

    override_formats = {'varParamStr': '%s', 'sn_truth_params': '%s'}

    @cached
    def get_is_sn(self):
        sn = self.column_by_name('sn_truth_params').astype(str)
        return np.where(np.char.find(sn, 'None')==0, 0, 1)

    @cached
    def get_is_agn(self):
        var = self.column_by_name('varParamStr').astype(str)
        return np.where(np.char.find(var, 'None')==0, 0, 1)


class BulgeTruth(_SersicTruth, SQLSubCatalogMixin, PhoSimDESCQA):
    cannot_be_null = ['hasBulge', 'sedFilepath']
    _file_name = 'sprinkled_objects.sqlite'
    _table_name = 'bulge'


class DiskTruth(_SersicTruth, SQLSubCatalogMixin, PhoSimDESCQA):
    cannot_be_null = ['hasDisk', 'sedFilepath']
    _file_name = 'sprinkled_objects.sqlite'
    _table_name = 'disk'


class AgnTruth(_ZPointTruth, SQLSubCatalogMixin, TwinklesCatalogZPoint_DC2):
    cannot_be_null = ['has_params']
    _file_name = 'sprinkled_objects.sqlite'
    _table_name = 'zpoint'

    def get_has_params(self):
        sn = self.column_by_name('is_sn')
        agn = self.column_by_name('is_agn')
        magnorm = self.column_by_name('magNorm')

        with np.errstate(divide='ignore', invalid='ignore'):
            return np.where(np.logical_or(sn==1,
                            np.logical_and(agn==1,
                            np.logical_and(np.isfinite(magnorm),
                                           magnorm<100.0))),
                            True,
                            None)


def write_sprinkled_param_db(obs, field_ra=55.064, field_dec=-29.783,
                             agn_db=None, yaml_file='proto-dc2_v4.6.1',
                             out_dir=None):
    """
    This method writes out a sqlite database that contains truth information
    on all of the sprinkled sources.  It will return the name of the database
    and a list of the tables in that database.

    Parameters
    ----------
    obs is an ObservationMetaData specifying the entire survey area

    field_ra, field_dec are the centers of the survey region in degrees

    agn_db is the path to the database of AGN parameters

    yaml_file is the yaml file denoting the extragalactic catalog
    to be loaded (defaulte proto-dc2_v4.6.1)

    out_dir is the directory where the output database will be written

    Returns
    -------
    The name of the written sqlite database and a list of the tables
    it contains.
    """
    if not os.path.isfile(agn_db):
        raise RuntimeError("%s is not a valid file" % agn_db)

    if not os.path.isdir(out_dir):
        raise RuntimeError("%s is not a valid dir" % out_dir)

    # needed for the supernovae
    twinkles_spec_map.subdir_map['(^specFileGLSN)'] = 'Dynamic'

    cat_class_list = [BulgeTruth, DiskTruth, AgnTruth]
    db_class_list = [bulgeDESCQAObject,
                     diskDESCQAObject,
                     agnDESCQAObject]

    table_name_list = []
    for cat_class in cat_class_list:
        assert cat_class._file_name == cat_class_list[0]._file_name
        table_name_list.append(cat_class._table_name)

    for db_class in db_class_list:
        db_class.yaml_file_name = yaml_file

    cat = TwinklesCompoundInstanceCatalog_DC2(cat_class_list,
                                              db_class_list,
                                              obs_metadata=obs,
                                              field_ra=field_ra,
                                              field_dec=field_dec,
                                              agn_params_db=agn_db,
                                              compoundDBclass=sprinklerCompound_DC2_truth)

    cat.sed_dir = None

    cat.write_catalog(os.path.join(out_dir,'params.txt'), chunk_size=10000)

    # remove empty text file
    txt_name = os.path.join(out_dir, 'params.txt')
    if os.path.exists(txt_name):
        os.unlink(txt_name)

    zpoint_file_name = os.path.join(out_dir, AgnTruth._file_name)
    with sqlite3.connect(zpoint_file_name) as connection:
        cursor = connection.cursor()
        index_cmd = 'CREATE INDEX htmid_indx ON %s (htmid, is_sn, is_agn)' % (AgnTruth._table_name)
        cursor.execute(index_cmd)
        connection.commit()

        cursor.execute('CREATE INDEX agn_gal_id ON zpoint (galaxy_id, is_sn, is_agn)')
        connection.commit()
        cursor.execute('CREATE INDEX gal_id_bulge ON bulge (galaxy_id)')
        connection.commit()
        cursor.execute('CREATE INDEX gal_id_disk ON disk (galaxy_id)')
        connection.commit()

    full_file_name = os.path.join(out_dir, cat_class_list[0]._file_name)
    if not os.path.exists(full_file_name):
        raise RuntimeError("After generating db; %s does not exist" % full_file_name)

    return full_file_name, table_name_list


def get_pointing_htmid(pointing_dir, opsim_db_name,
                       ra_colname = 'descDitheredRA',
                       dec_colname = 'descDitheredDec'):
    """
    For a list of OpSim pointings, find dicts mapping those pointings to:
    - The trixels filling the pointings
    - The MJDs of the pointings
    - The telescope filters of the pointings

    Parameters
    ----------
    pointing_dir contains a series of files that are two columns: obshistid, mjd.
    The files must each have 'visits' in their name.  These specify the pointings
    for which we are assembling data.  See:
        https://github.com/LSSTDESC/DC2_Repo/tree/master/data/Run1.1
    for an example.

    opsim_db_name is the path to the OpSim database from which to take those pointings

    ra_colname is the column used for RA of the pointing (default:
    descDitheredRA)

    dec_colname is the column used for the Dec of the pointing (default:
    descDitheredDec)

    Returns
    -------
    htmid_bound_dict -- a dict keyed on ObsHistID.  Values are the list of trixels filling
    the OpSim pointing, as returned by lsst.sims.utils.HalfSpace.findAllTrixels

    mjd_dict -- a dict keyed on ObsHistID.  Values are the MJD(TAI) of the OpSim pointings.

    filter_dict -- a dict keyed on ObsHistID.  Values are the 'ugrizy' filter of the
    OpSim pointings.
    """

    radius = 2.0  # field of view of a pointing in degrees

    if not os.path.isfile(opsim_db_name):
        raise RuntimeError("%s is not a valid file name" % opsim_db_name)

    if not os.path.isdir(pointing_dir):
        raise RuntimeError("%s is not a valid dir name" % pointing_dir)

    dtype = np.dtype([('obshistid', int), ('mjd', float)])
    obs_data = None
    for file_name in os.listdir(pointing_dir):
        if 'visits' in file_name:
            full_name = os.path.join(pointing_dir, file_name)
            data = np.genfromtxt(full_name, dtype=dtype)
            if obs_data is None:
                obs_data = data['obshistid']
            else:
                obs_data = np.concatenate((obs_data, data['obshistid']), axis=0)

    obs_data = np.sort(obs_data)

    db = DBObject(opsim_db_name, driver='sqlite')
    dtype = np.dtype([('obshistid', int), ('mjd', float), ('band', str, 1),
                      ('ra', float), ('dec', float)])

    htmid_bound_dict = {}
    mjd_dict = {}
    filter_dict = {}

    d_obs = len(obs_data)//5
    for i_start in range(0,len(obs_data), d_obs):
        i_end = i_start + d_obs
        if len(obs_data)-i_start < d_obs:
            i_end = len(obs_data)

        subset = obs_data[i_start:i_end]

        query = 'SELECT obsHistId, expMJD, filter, %s, %s FROM Summary' % (ra_colname, dec_colname)
        query += ' WHERE obsHistID BETWEEN %d and %e' % (subset.min(), subset.max())
        query += ' GROUP BY obsHistID'

        results = db.execute_arbitrary(query, dtype=dtype)

        for ii in range(len(results)):
            obshistid = results['obshistid'][ii]
            if obshistid not in obs_data:
                continue

            hs = halfSpaceFromRaDec(np.degrees(results['ra'][ii]),
                                    np.degrees(results['dec'][ii]),
                                    radius)

            trixel_bounds = hs.findAllTrixels(_truth_trixel_level)
            htmid_bound_dict[obshistid] = trixel_bounds
            mjd_dict[obshistid] = results['mjd'][ii]
            filter_dict[obshistid] = results['band'][ii]

    assert len(obs_data) == len(htmid_bound_dict)

    return htmid_bound_dict, mjd_dict, filter_dict


def write_star_truth_db(obs_md, out_file_name, bp_dict):
    """
    Write out the truth information for static stellar sources.

    Parameters
    ----------
    obs_md is an ObservationMetaData covering the whole survey area

    out_file_name is the name of the sqlite file where the results
    will be written

    bp_dict is a BandpassDict containing lsst.sims.photUtils.Bandpasses
    corresponding to the telescope filters.

    Returns
    -------
    Nothing

    Writes out a table 'star_truth' containing the RA, Dec,
    and magnitudes of all of the stars in the survey area to the
    file specified by out_file_name
    """

    sed_dir = os.environ['SIMS_SED_LIBRARY_DIR']

    db = StarObj(database='LSSTCATSIM', host='fatboy.phys.washington.edu',
                 port=1433, driver='mssql+pymssql')

    colnames = ['ra', 'decl', 'magNorm', 'sedFilename']

    star_iter = db.query_columns(colnames=colnames,
                                 obs_metadata=obs_md,
                                 chunk_size=10000)

    with sqlite3.connect(out_file_name) as conn:
        cursor = conn.cursor()
        cursor.execute('''CREATE TABLE star_truth
                          (ra float, dec float,
                           u float, g float, r float,
                           i float, z float, y float)''')

        conn.commit()

        for star_chunk in star_iter:
            mags = np.NaN*np.ones((len(star_chunk), 6), dtype=float)
            for i_star, star in enumerate(star_chunk):
                full_name = os.path.join(sed_dir,
                                         defaultSpecMap[star['sedFilename']])
                spec = Sed()
                spec.readSED_flambda(full_name)
                fnorm = getImsimFluxNorm(spec, star['magNorm'])
                spec.multiplyFluxNorm(fnorm)
                mag_list = bp_dict.magListForSed(spec)
                mags[i_star] = mag_list

            values = ((star_chunk['ra'][i_obj], star_chunk['dec'][i_obj],
                       mags[i_obj][0], mags[i_obj][1], mags[i_obj][2],
                       mags[i_obj][3], mags[i_obj][4], mags[i_obj][5])
                      for i_obj in range(len(star_chunk)))

            cursor.executemany('''INSERT INTO star_truth
                                  VALUES (?,?,?,?,?,?,?,?)''', values)
            conn.commit()
