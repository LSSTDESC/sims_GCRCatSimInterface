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
                      'raJ2000', 'decJ2000', 'redshift',
                      'sedFilepath', 'magNorm',
                      'majorAxis', 'minorAxis',
                      'positionAngle',
                      'internalAv', 'internalRv']


class _ZPointTruth(_SprinkledTruth):
    column_outputs = ['uniqueId', 'galaxy_id',
                      'raJ2000', 'decJ2000', 'redshift',
                      'sedFilepath', 'magNorm',
                      'varParamStr', 'sn_truth_params',
                      'has_params']

    override_formats = {'varParamStr': '%s', 'sn_truth_params': '%s'}

class BulgeTruth(_SersicTruth, SubCatalogMixin, PhoSimDESCQA):
    cannot_be_null = ['hasBulge', 'sprinkling_switch', 'sedFilepath']
    subcat_prefix = 'bulge'
    _write_subcat_header = True

class DiskTruth(_SersicTruth, SubCatalogMixin, PhoSimDESCQA):
    cannot_be_null = ['hasDisk', 'sprinkling_switch', 'sedFilepath']
    subcat_prefix = 'disk'
    _write_subcat_header = True

class AgnTruth(_ZPointTruth, SubCatalogMixin, TwinklesCatalogZPoint_DC2):
    cannot_be_null = ['sprinkling_switch', 'has_params']
    subcat_prefix = 'agn'
    _write_subcat_header = True

    def get_has_params(self):
        varpar = self.column_by_name('varParamStr').astype(str)
        snpar = self.column_by_name('sn_truth_params').astype(str)
        magnorm = self.column_by_name('magNorm')

        output = []
        for vv, ss, mm in zip(varpar, snpar, magnorm):
            if vv == 'None' and ss == 'None':
                output.append(None)
            elif vv!='None' and not (np.isfinite(mm) or mm < 100.0):
                output.append(None)
            else:
                output.append(True)

        return np.array(output)

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
