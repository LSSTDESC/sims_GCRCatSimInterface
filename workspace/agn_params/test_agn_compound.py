import os
from desc.sims.GCRCatSimInterface import CompoundDESCQAInstanceCatalog
from desc.sims.GCRCatSimInterface import bulgeDESCQAObject_protoDC2
from desc.sims.GCRCatSimInterface import diskDESCQAObject_protoDC2
from desc.sims.GCRCatSimInterface import agnDESCQAObject_protoDC2
from desc.sims.GCRCatSimInterface import GalaxyCompoundDESCQAObject
from lsst.sims.catalogs.definitions import InstanceCatalog
from lsst.sims.utils import arcsecFromRadians
from lsst.sims.utils import ObservationMetaData
import numpy as np

global_agn_params_db = os.path.join(os.environ['SCRATCH'], 'proto_dc2_agn', 'test_agn.db')
assert os.path.exists(global_agn_params_db)

class _testDESCQAObj(object):
    yaml_file_name = 'proto-dc2_v2.1.2'
    field_ra = 23.0
    field_dec = -22.3

class bulgeDESCQAObject_test(_testDESCQAObj, bulgeDESCQAObject_protoDC2):
    objid = 'bulge_descqa'

class diskDESCQAObject_test(_testDESCQAObj, diskDESCQAObject_protoDC2):
    objid = 'disk_descqa'

class agnDESCQAObject_test(_testDESCQAObj, agnDESCQAObject_protoDC2):
    objid = 'agn_descqa'

class CatForBulge(InstanceCatalog):
    column_outputs = ['uniqueId', 'objid', 'galaxy_id', 'raJ2000', 'decJ2000',
                      'minorAxis', 'majorAxis']
    transformations = {'minorAxis': arcsecFromRadians,
                       'majorAxis': arcsecFromRadians,
                       'raJ2000': np.degrees,
                       'decJ2000': np.degrees}

    def get_objid(self):
        return np.array([self.db_obj.objectTypeId]*len(self.column_by_name('raJ2000')))

class CatForDisk(InstanceCatalog):
    column_outputs = ['uniqueId', 'objid', 'galaxy_id', 'raJ2000', 'decJ2000',
                      'minorAxis']
    transformations = {'minorAxis': arcsecFromRadians,
                       'raJ2000': np.degrees,
                       'decJ2000': np.degrees}


    def get_objid(self):
        return np.array([self.db_obj.objectTypeId]*len(self.column_by_name('raJ2000')))

class CatForAgn(InstanceCatalog):
    column_outputs = ['uniqueId', 'objid', 'galaxy_id', 'raJ2000', 'decJ2000', 'varParamStr',
                      'magNorm']

    transformations = {'raJ2000': np.degrees,
                       'decJ2000': np.degrees}


    def get_objid(self):
        return np.array([self.db_obj.objectTypeId]*len(self.column_by_name('raJ2000')))


class LocalGalaxyDESCQAObject(GalaxyCompoundDESCQAObject):
    agn_params_db = global_agn_params_db
    agn_objid = 'agn_descqa'

if __name__ == "__main__":

    obs = ObservationMetaData(pointingRA=22.7, pointingDec=-22.7,
                              boundType='circle', boundLength=0.01)

    cat = CompoundDESCQAInstanceCatalog([CatForBulge, CatForDisk, CatForAgn],
                                        [bulgeDESCQAObject_test,
                                         diskDESCQAObject_test,
                                         agnDESCQAObject_test],
                                        obs_metadata=obs,
                                        compoundDBclass=LocalGalaxyDESCQAObject)

    cat.write_catalog('descqa_compound_cat.txt', chunk_size=10000)
