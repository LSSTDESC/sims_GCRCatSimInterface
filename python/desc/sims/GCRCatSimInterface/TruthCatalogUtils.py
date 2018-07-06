import numpy as np
import os
import sqlite3
import time

from lsst.sims.catalogs.decorators import cached, compound
from lsst.sims.catalogs.definitions import InstanceCatalog
from lsst.sims.catalogs.db import DBObject
from lsst.sims.utils import findHtmid, halfSpaceFromRaDec
from lsst.sims.photUtils import Sed, getImsimFluxNorm
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
           "get_pointing_htmid"]

_truth_trixel_level = 7

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
                      'umag']#, 'gmag', 'rmag', 'imag',
                      #'zmag', 'ymag']

    _bp_dict = None

    #def get_umag(self):
    #    return np.random.random_sample(len(self.column_by_name('raJ2000')))

    @compound('umag', 'gmag', 'rmag', 'imag', 'zmag', 'ymag')
    def get_SersicMagnitudes(self):
        # list of columns that need to be requested
        # so that the sprinkler will work
        self.column_by_name('majorAxis')
        self.column_by_name('minorAxis')
        self.column_by_name('positionAngle')
        self.column_by_name('is_sprinkled')

        self.column_by_name('internalAv')
        self.column_by_name('internalRv')
        self.column_by_name('galacticAv')
        self.column_by_name('galacticRv')

        sed_name = self.column_by_name('sedFilepath')
        redshift = self.column_by_name('redshift')
        magnorm = self.column_by_name('magNorm')

        sed_dir = os.environ['SIMS_SED_LIBRARY_DIR']

        if not hasattr(self, '_gal_ct'):
            self._gal_ct = 0

        t_start = time.time()
        output = np.zeros((6,len(sed_name)), dtype=float)
        for i_obj in range(len(sed_name)):
            spec = Sed()
            spec.readSED_flambda(os.path.join(sed_dir, sed_name[i_obj]))
            fnorm = getImsimFluxNorm(spec, magnorm[i_obj])
            spec.multiplyFluxNorm(fnorm)
            spec.redshiftSED(redshift[i_obj], dimming=True)
            mag_arr = self._bp_dict.magListForSed(spec)
            output[:,i_obj] = mag_arr
        self._gal_ct += len(sed_name)
        print('did %d upto %d in %e seconds' %
        (len(sed_name), self._gal_ct, time.time()-t_start))

        return output


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
                             out_dir=None, bp_dict=None):
    """
    This method writes out a sqlite database that contains truth information
    on all of the sprinkled sources.  It will return the name of the database
    and a list of the tables in that database.
    ----

    obs is an ObservationMetaData
    """
    if not os.path.isfile(agn_db):
        raise RuntimeError("%s is not a valid file" % agn_db)

    if not os.path.isdir(out_dir):
        raise RuntimeError("%s is not a valid dir" % out_dir)

    twinkles_spec_map.subdir_map['(^specFileGLSN)'] = 'Dynamic'

    BulgeTruth._bp_dict = bp_dict
    DiskTruth._bp_dict = bp_dict

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

    cat.write_catalog(os.path.join(out_dir,'params.txt'), chunk_size=100000)

    txt_name = os.path.join(out_dir, 'params.txt')
    if os.path.exists(txt_name):
        os.unlink(txt_name)

    zpoint_file_name = os.path.join(out_dir, AgnTruth._file_name)
    with sqlite3.connect(zpoint_file_name) as connection:
        cursor = connection.cursor()
        index_cmd = 'CREATE INDEX htmid_indx ON %s (htmid, is_sn, is_agn)' % (AgnTruth._table_name)
        cursor.execute(index_cmd)
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

    radius = 2.0  # field of view in degrees

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
