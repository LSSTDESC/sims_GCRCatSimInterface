from desc.sims.GCRCatSimInterface import CompoundDESCQAInstanceCatalog
from desc.sims.GCRCatSimInterface import bulgeDESCQAObject_protoDC2
from desc.sims.GCRCatSimInterface import diskDESCQAObject_protoDC2
from desc.sims.GCRCatSimInterface import CompoundDESCQAObject
from lsst.sims.catalogs.definitions import InstanceCatalog
from lsst.sims.utils import arcsecFromRadians
from lsst.sims.utils import ObservationMetaData
import numpy as np

class _testDESCQAObj(object):
    yaml_file_name = 'proto-dc2_v2.1.2'
    field_ra = 23.0
    field_dec = -22.3

class bulgeDESCQAObject_test(_testDESCQAObj, bulgeDESCQAObject_protoDC2):
    objid = 'bulge_descqa'

class diskDESCQAObject_test(_testDESCQAObj, diskDESCQAObject_protoDC2):
    objid = 'disk_descqa'

class CatForBulge(InstanceCatalog):
    column_outputs = ['uniqueId', 'raJ2000', 'decJ2000',
                      'minorAxis', 'majorAxis']
    transformations = {'minorAxis': arcsecFromRadians,
                       'majorAxis': arcsecFromRadians,
                       'raJ2000': np.degrees,
                       'decJ2000': np.degrees}

class CatForDisk(InstanceCatalog):
    column_outputs = ['uniqueId', 'raJ2000', 'decJ2000',
                      'minorAxis']
    transformations = {'minorAxis': arcsecFromRadians,
                       'raJ2000': np.degrees,
                       'decJ2000': np.degrees}

if __name__ == "__main__":

    obs = ObservationMetaData(pointingRA=22.7, pointingDec=-22.7,
                              boundType='circle', boundLength=0.01)

    cat = CompoundDESCQAInstanceCatalog([CatForBulge, CatForDisk],
                                        [bulgeDESCQAObject_test,
                                         diskDESCQAObject_test],
                                        obs_metadata=obs,
                                        compoundDBclass=CompoundDESCQAObject)

    cat.write_catalog('descqa_compound_cat.txt', chunk_size=10000)
