import numpy as np
import os
import sqlite3

from lsst.sims.catalogs.decorators import cached
from lsst.sims.catalogs.definitions import InstanceCatalog
from lsst.sims.utils import findHtmid
from . import sprinklerCompound_DC2_truth
from . import TwinklesCompoundInstanceCatalog_DC2
from . import SQLSubCatalogMixin
from . import diskDESCQAObject_protoDC2 as diskDESCQAObject
from . import bulgeDESCQAObject_protoDC2 as bulgeDESCQAObject
from . import agnDESCQAObject_protoDC2 as agnDESCQAObject
from . import PhoSimDESCQA
from . import TwinklesCatalogZPoint_DC2
from .TwinklesClasses import twinkles_spec_map

__all__ = ["write_sprinkled_truth"]

class _SprinkledTruth(object):

    def get_htmid(self):
        ra = np.degrees(self.column_by_name('raJ2000'))
        dec = np.degrees(self.column_by_name('decJ2000'))
        return findHtmid(ra, dec, max_level=7)

    def write_header(self, file_handle):
        InstanceCatalog.write_header(self, file_handle)

class _SersicTruth(_SprinkledTruth):
    column_outputs = ['uniqueId', 'galaxy_id', 'htmid',
                      'raJ2000', 'decJ2000', 'redshift',
                      'sedFilepath', 'magNorm',
                      'majorAxis', 'minorAxis',
                      'positionAngle',
                      'internalAv', 'internalRv']


class _ZPointTruth(_SprinkledTruth):
    column_outputs = ['uniqueId', 'galaxy_id', 'htmid',
                      'raJ2000', 'decJ2000', 'redshift',
                      'sedFilepath', 'magNorm',
                      'varParamStr', 'sn_truth_params',
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
    cannot_be_null = ['hasBulge', 'sprinkling_switch', 'sedFilepath']
    _file_name = 'sprinkled_objects.sqlite'
    _table_name = 'bulge'

class DiskTruth(_SersicTruth, SQLSubCatalogMixin, PhoSimDESCQA):
    cannot_be_null = ['hasDisk', 'sprinkling_switch', 'sedFilepath']
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

def write_sprinkled_truth(obs, field_ra=55.064, field_dec=-29.783,
                          agn_db=None, yaml_file='proto-dc2_v4.6.1'):
    """
    obs is an ObservationMetaData
    """
    assert os.path.isfile(agn_db)

    out_dir = 'workspace/catalogs/truth_params/'
    assert os.path.isdir(out_dir)

    twinkles_spec_map.subdir_map['(^specFileGLSN)'] = 'Dynamic'

    cat_class_list = [BulgeTruth, DiskTruth, AgnTruth]
    db_class_list = [bulgeDESCQAObject,
                     diskDESCQAObject,
                     agnDESCQAObject]

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

    zpoint_file_name = os.path.join(out_dir, AgnTruth._file_name)
    with sqlite3.connect(zpoint_file_name) as connection:
        cursor = connection.cursor()
        index_cmd = 'CREATE INDEX htmid_indx ON %s (htmid, is_sn, is_agn)' % (AgnTruth._table_name)
        cursor.execute(index_cmd)
        connection.commit()
