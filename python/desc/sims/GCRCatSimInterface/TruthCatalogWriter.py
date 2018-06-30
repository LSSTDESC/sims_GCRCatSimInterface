import numpy as np
import os

from lsst.sims.catalogs.definitions import InstanceCatalog
from . import sprinklerCompound_DC2_truth
from . import TwinklesCompoundInstanceCatalog_DC2
from . import SubCatalogMixin
from . import diskDESCQAObject_protoDC2 as diskDESCQAObject
from . import bulgeDESCQAObject_protoDC2 as bulgeDESCQAObject
from . import agnDESCQAObject_protoDC2 as agnDESCQAObject
from . import PhoSimDESCQA
from . import TwinklesCatalogZPoint_DC2
from .TwinklesClasses import twinkles_spec_map

__all__ = ["write_sprinkled_truth"]

class _SprinkledTruth(object):

    def write_header(self, file_handle):
        InstanceCatalog.write_header(self, file_handle)

class _SersicTruth(_SprinkledTruth):
    column_outputs = ['uniqueId', 'galaxy_id',
                      'raJ2000', 'decJ2000',
                      'sedFilepath', 'magNorm']


class _ZPointTruth(_SprinkledTruth):
    column_outputs = ['uniqueId', 'galaxy_id',
                      'raJ2000', 'decJ2000',
                      'sedFilepath', 'magNorm',
                      'varParamStr', 'sn_truth_params']

    override_formats = {'varParamStr': '%s', 'sn_truth_params': '%s'}

    def get_sn_truth_params(self):
        n_obj = len(self.column_by_name('raJ2000'))
        return np.array([None]*n_obj).astype(str)


class BulgeTruth(_SersicTruth, SubCatalogMixin, PhoSimDESCQA):
    cannot_be_null = ['hasBulge', 'is_sprinkled', 'sedFilepath']
    subcat_prefix = 'bulge'

class DiskTruth(_SersicTruth, SubCatalogMixin, PhoSimDESCQA):
    cannot_be_null = ['hasDisk', 'is_sprinkled', 'sedFilepath']
    subcat_prefix = 'disk'

class AgnTruth(_ZPointTruth, SubCatalogMixin, TwinklesCatalogZPoint_DC2):
    cannot_be_null = ['is_sprinkled']
    subcat_prefix = 'agn'

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
